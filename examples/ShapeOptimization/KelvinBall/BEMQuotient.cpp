/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "BEM.h"
#include "BEMQuotientKernel.h"
#include "ChamberSymmetry.h"
#include <Rodin/Solver/PartialPivLU.h>
#include <map>
#include <chrono>
#include <iostream>

namespace KelvinBall
{
  BEM::Result BEM::evaluateQuotient(const Mesh& surface,bool derivative,Diagnostics* diagnostics) const
  {
    using Clock=std::chrono::steady_clock;
    if (diagnostics) *diagnostics=Diagnostics{};
    if (surface.getDimension()!=2 || !surface.getCellCount())
      throw std::runtime_error("O24 BEM requires a nonempty chamber interface.");
    Result result; result.panels=surface.getCellCount();
    Real totalArea=0;
    // Audit virtual image edge incidence, without constructing any image mesh.
    using Key=ChamberSymmetry::Key;
    std::map<std::pair<Key,Key>,size_t> edges;
    for (auto cell=surface.getCell();cell;++cell)
    {
      const auto& v=cell->getVertices();
      if (v.size()!=3) throw std::runtime_error("O24 BEM requires affine triangles.");
      const Math::SpatialVector<Real> a=surface.getVertexCoordinates(v(0));
      const Math::SpatialVector<Real> b=surface.getVertexCoordinates(v(1));
      const Math::SpatialVector<Real> c=surface.getVertexCoordinates(v(2));
      const Real area=(b-a).cross(c-a).norm()/2;
      if (!(area>0)) throw std::runtime_error("Degenerate chamber BEM triangle.");
      totalArea+=area;
      for (const auto& g:SewedOutput::getCubeRotations())
        for (size_t i=0;i<3;++i)
        {
          const auto x=ChamberSymmetry::key(Math::SpatialPoint(g*surface.getVertexCoordinates(v(i))));
          const auto y=ChamberSymmetry::key(Math::SpatialPoint(g*surface.getVertexCoordinates(v((i+1)%3))));
          ++edges[std::minmax(x,y)];
        }
    }
    for (const auto& [edge,count]:edges) if (count!=2) ++result.openEdges;
    if (result.openEdges) throw std::runtime_error("O24 chamber interface has unmatched image edges: "+
      std::to_string(result.openEdges));
    result.epsilon=m_epsilonFactor*std::sqrt(totalArea/result.panels);
    const BEMQuotientKernel imageKernel(result.epsilon,m_viscosity);
    const RegularizedStokeslet kernel(result.epsilon,m_viscosity);
    P1 space(surface,3,3);
    TrialFunction trial(space); TestFunction test(space);
    auto form=Integral(Potential(imageKernel,trial),test); form.setOrder(m_order);
    DenseProblem problem(trial,test);
    Solver::PartialPivLU solver(problem);
    const auto basisStart=Clock::now();
    const auto E=ChamberSymmetry(surface).tensorBasis(space);
    if (diagnostics) diagnostics->basisSeconds=std::chrono::duration<Real>(Clock::now()-basisStart).count();
    result.stateDOFs=space.getSize(); result.independentDOFs=E.cols();
    Math::Matrix<Real> rhs(space.getSize(),2),solution(space.getSize(),2),A;
    for (size_t family=0;family<2;++family)
    {
      MatrixFunction rigid(size_t(3),size_t(3),[family](const Point& p) {
        Math::SpatialMatrix<Real> value(3,3); value.setIdentity();
        if (family==1)
        {
          for (size_t i=0;i<3;++i)
          {
            Math::SpatialVector<Real> axis(3);axis.setZero();axis(i)=1;
            const auto column=axis.cross(p.vector());
            for (size_t j=0;j<3;++j) value(j,i)=column(j);
          }
        }
        return value;
      });
      auto linear=Integral(rigid,test); linear.setOrder(m_order);
      problem=form-linear;
      const auto assemblyStart=Clock::now();
      if (family==0) problem.assemble();
      else problem.assemble(AssemblyTarget::RHS);
      if (diagnostics)
      {
        diagnostics->assemblySeconds[family]=std::chrono::duration<Real>(Clock::now()-assemblyStart).count();
        std::cout<<"diagnostic family="<<family<<" assembly_seconds="<<diagnostics->assemblySeconds[family]<<std::endl;
      }
      auto& system=problem.getLinearSystem();
      if (family==0) A=system.getOperator();
      rhs.col(family)=system.getVector();
      const auto reductionStart=Clock::now();
      system.getOperator()=E.transpose()*A*E;
      system.getVector()=E.transpose()*rhs.col(family);
      if (diagnostics) diagnostics->reductionSeconds[family]=std::chrono::duration<Real>(Clock::now()-reductionStart).count();
      const auto& reduced=system.getOperator();
      const Real asymmetry=(reduced-reduced.transpose()).norm()/std::max(reduced.norm(),Real(1e-30));
      if (asymmetry>1e-11) throw std::runtime_error("O24 BEM operator lost reciprocity.");
      if (diagnostics)
      {
        // Match PartialPivLU::solve exactly, exposing its two existing stages.
        Eigen::PartialPivLU<Math::Matrix<Real>> lu;
        const auto factorStart=Clock::now();
        lu.compute(reduced);
        diagnostics->factorSeconds[family]=std::chrono::duration<Real>(Clock::now()-factorStart).count();
        const auto solveStart=Clock::now();
        system.getSolution()=lu.solve(system.getVector());
        diagnostics->solveSeconds[family]=std::chrono::duration<Real>(Clock::now()-solveStart).count();
        std::cout<<"diagnostic family="<<family<<" reduction_seconds="<<diagnostics->reductionSeconds[family]
          <<" factor_seconds="<<diagnostics->factorSeconds[family]<<" solve_seconds="<<diagnostics->solveSeconds[family]<<std::endl;
        if (family==1)
        {
          const auto conditionStart=Clock::now();
          diagnostics->conditionOneEstimate=1/lu.rcond();
          // This extra all-RHS solve is diagnostic overhead, not state back substitution.
          const Math::Matrix<Real> inverse=lu.solve(Math::Matrix<Real>::Identity(reduced.rows(),reduced.cols()));
          diagnostics->operatorNormOne=reduced.cwiseAbs().colwise().sum().maxCoeff();
          diagnostics->inverseNormOne=inverse.cwiseAbs().colwise().sum().maxCoeff();
          diagnostics->conditionOne=diagnostics->operatorNormOne*diagnostics->inverseNormOne;
          const Math::Matrix<Real> defect=reduced*inverse-Math::Matrix<Real>::Identity(reduced.rows(),reduced.cols());
          diagnostics->inverseResidual=defect.cwiseAbs().colwise().sum().maxCoeff()/diagnostics->conditionOne;
          diagnostics->conditionSeconds=std::chrono::duration<Real>(Clock::now()-conditionStart).count();
          if (!inverse.allFinite() || !std::isfinite(diagnostics->conditionOne))
            throw std::runtime_error("Nonfinite BEM condition diagnostic.");
          std::cout<<"diagnostic condition_one="<<diagnostics->conditionOne<<std::endl;
        }
      }
      else solver.solve(system);
      solution.col(family)=E*system.getSolution();
      const Real residual=(E.transpose()*(A*solution.col(family)-rhs.col(family))).norm()/
        std::max((E.transpose()*rhs.col(family)).norm(),Real(1e-30));
      result.residual=std::max(result.residual,residual);
    }
    if (!solution.allFinite() || !std::isfinite(result.residual) || result.residual>1e-9)
      throw std::runtime_error("O24 BEM constrained state residual failed.");
    const Math::Matrix<Real> response=rhs.transpose()*solution;
    result.reciprocity=std::abs(response(0,1)-response(1,0))/std::max(response.norm(),Real(1e-30));
    if (result.reciprocity>1e-9) throw std::runtime_error("O24 BEM response reciprocity failed.");
    // 24 image patches and three rigid-load columns: 24/3 = 8.
    constexpr Real multiplicity=24, columns=3;
    result.k=multiplicity/columns*response(0,0);
    result.q=multiplicity/columns*response(1,1);
    result.c=multiplicity/(2*columns)*(response(0,1)+response(1,0));
    if (!(result.k>0 && result.q>0)) throw std::runtime_error("Invalid O24 BEM resistance.");
    result.rho=std::abs(result.c)/std::sqrt(result.k*result.q);
    result.resistance=Math::Matrix<Real>::Zero(6,6);
    for (size_t i=0;i<3;++i)
    {
      result.resistance(i,i)=result.k; result.resistance(i+3,i+3)=result.q;
      result.resistance(i,i+3)=result.resistance(i+3,i)=result.c;
    }
    if (!derivative) return result;
    if (std::abs(result.c)<1e-13*std::sqrt(result.k*result.q))
      throw std::runtime_error("O24 absolute coupling gradient is undefined at zero.");
    Math::Matrix<Real> W(2,2);
    const Real wc=(result.c>0 ? 1 : -1)/std::sqrt(result.k*result.q);
    W(0,0)=-4*result.rho/result.k; W(1,1)=-4*result.rho/result.q;
    W(0,1)=W(1,0)=4*wc;
    struct Sample
    {
      std::array<Index,3> vertices;
      std::array<Real,3> barycentric;
      Math::SpatialVector<Real> x;
      std::array<Math::SpatialMatrix<Real>,2> density,weighted;
      std::array<Math::SpatialVector<Real>,3> areaGradient;
      Real weight,area;
    };
    std::vector<Sample> samples;
    P1<Real> scalarSpace(surface);
    for (auto cell=surface.getCell();cell;++cell)
    {
      const auto& v=cell->getVertices();
      const Math::SpatialVector<Real> a=surface.getVertexCoordinates(v(0));
      const Math::SpatialVector<Real> b=surface.getVertexCoordinates(v(1));
      const Math::SpatialVector<Real> c=surface.getVertexCoordinates(v(2));
      Math::SpatialVector<Real> n=(b-a).cross(c-a);const Real area=n.norm()/2;n/=2*area;
      const auto& qf=kernel.getQuadrature(m_order,cell->getGeometry());
      const auto& quadrature=cell->getQuadrature(qf);
      const auto& fe=space.getFiniteElement(2,cell->getIndex());
      const auto& sfe=scalarSpace.getFiniteElement(2,cell->getIndex());
      const auto& dofs=space.getDOFs(2,cell->getIndex());
      for (size_t q=0;q<quadrature.getSize();++q)
      {
        Sample s;s.area=area;s.x=quadrature.getPoint(q).vector();
        s.weight=qf.getWeight(q)*quadrature.getPoint(q).getDistortion();
        for (size_t family=0;family<2;++family)
        {
          s.density[family].resize(3,3);s.density[family].setZero();
          for (size_t i=0;i<fe.getCount();++i)
            s.density[family]+=fe.getBasis(i)(qf.getPoint(q))*solution(dofs(i),family);
        }
        s.weighted[0]=W(0,0)*s.density[0]+W(1,0)*s.density[1];
        s.weighted[1]=W(0,1)*s.density[0]+W(1,1)*s.density[1];
        for (size_t i=0;i<3;++i) { s.vertices[i]=v(i);s.barycentric[i]=sfe.getBasis(i)(qf.getPoint(q)); }
        s.areaGradient={Math::SpatialVector<Real>((b-c).cross(n)/2),
          Math::SpatialVector<Real>((c-a).cross(n)/2),Math::SpatialVector<Real>((a-b).cross(n)/2)};
        samples.push_back(std::move(s));
      }
    }
    result.gradient=Math::Matrix<Real>::Zero(surface.getVertexCount(),3);
    Real epsilonGradient=0;
    auto deposit=[&](const Sample& s,const Math::SpatialVector<Real>& dx,Real dArea) {
      for (size_t i=0;i<3;++i) result.gradient.row(s.vertices[i])+=
        (s.barycentric[i]*dx+dArea*s.areaGradient[i]).getData().transpose();
    };
    for (const auto& x:samples)
    {
      Math::SpatialMatrix<Real> rigid(3,3);
      for (size_t i=0;i<3;++i)
      {
        Math::SpatialVector<Real> axis(3);axis.setZero();axis(i)=1;
        const auto column=axis.cross(x.x);
        for (size_t j=0;j<3;++j) rigid(j,i)=column(j);
      }
      const Math::SpatialMatrix<Real> P=2*x.weighted[1];
      Math::SpatialVector<Real> dx(3);
      dx.getData()<<P(1,2)-P(2,1),P(2,0)-P(0,2),P(0,1)-P(1,0);
      const Real value=2*x.weight*(x.weighted[0].trace()+
        x.weighted[1].getData().cwiseProduct(rigid.getData()).sum());
      deposit(x,Math::SpatialVector<Real>(x.weight*dx),value/x.area);
      for (const auto& y:samples)
        for (const auto& g:SewedOutput::getCubeRotations())
        {
          Math::SpatialMatrix<Real> H(3,3);H.setZero();
          for (size_t family=0;family<2;++family)
          {
            const Math::SpatialMatrix<Real> image=g*y.density[family]*g.transpose();
            H=H-x.weight*y.weight*x.weighted[family]*image.transpose();
          }
          const Math::SpatialVector<Real> r=x.x-g*y.x;
          Math::SpatialMatrix<Real> K(3,3);kernel.evaluate(K,r);
          Math::SpatialVector<Real> dr(3);Real de;kernel.differential(r,H,dr,de);
          const Real value=H.getData().cwiseProduct(K.getData()).sum();
          deposit(x,dr,value/x.area);
          deposit(y,Math::SpatialVector<Real>(-g.transpose()*dr),value/y.area);
          epsilonGradient+=de;
        }
    }
    const Real coefficient=epsilonGradient*result.epsilon/(2*totalArea);
    for (auto cell=surface.getCell();cell;++cell)
    {
      const auto& v=cell->getVertices();
      const Math::SpatialVector<Real> a=surface.getVertexCoordinates(v(0));
      const Math::SpatialVector<Real> b=surface.getVertexCoordinates(v(1));
      const Math::SpatialVector<Real> c=surface.getVertexCoordinates(v(2));
      Math::SpatialVector<Real> n=(b-a).cross(c-a);n.normalize();
      result.gradient.row(v(0))+=(coefficient*(b-c).cross(n)/2).getData().transpose();
      result.gradient.row(v(1))+=(coefficient*(c-a).cross(n)/2).getData().transpose();
      result.gradient.row(v(2))+=(coefficient*(a-b).cross(n)/2).getData().transpose();
    }
    if (!result.gradient.allFinite()) throw std::runtime_error("Nonfinite O24 BEM gradient.");
    return result;
  }
}
