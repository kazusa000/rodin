/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "BEM.h"
#include "Sphere.h"
#include "ChamberSymmetry.h"
#include <Rodin/MMG/MeshLoader.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace KelvinBall;
int main(int argc,char** argv)
{
  try
  {
    Configuration configuration;
    configuration.points=7;configuration.hminFactor=.8;configuration.hmaxFactor=1.2;
    configuration.seedAmplitude=.1;configuration.finalize();
    MMG::Mesh mesh;
    const bool reference=argc==2 && std::string(argv[1])=="--reference";
    if (argc==2 && !reference)
    {
      std::ifstream input(argv[1]);if (!input) throw std::runtime_error("Cannot read chamber mesh.");
      MMG::MeshLoader(mesh).load(input);
    }
    else mesh=Sphere(configuration).discretize(true,0,true).mesh;
    prepare(mesh);
    VelocitySpace space(mesh,3);Math::Vector<Real> load;
    const auto result=BEM().evaluateChamber(mesh,space,load,true);
    Real referenceError=0;
    if (reference)
    {
      auto solid=mesh.keep(Obstacle);solid.getConnectivity().compute(2,3);
      const auto skin=solid.skin();const auto patch=skin.keep(Gamma);
      std::map<ChamberSymmetry::Key,Index> indices;
      std::vector<Math::SpatialPoint> coordinates;
      std::vector<std::array<Index,3>> faces;
      for (const auto& g:SewedOutput::getCubeRotations())
      {
        std::vector<Index> map(patch.getVertexCount());
        for (Index v=0;v<patch.getVertexCount();++v)
        {
          const Math::SpatialPoint x=g*patch.getVertexCoordinates(v);
          auto [it,inserted]=indices.emplace(ChamberSymmetry::key(x),coordinates.size());
          if (inserted) coordinates.push_back(x);
          map[v]=it->second;
        }
        for (auto cell=patch.getCell();cell;++cell)
        {
          const auto& v=cell->getVertices();faces.push_back({map[v(0)],map[v(1)],map[v(2)]});
        }
      }
      KelvinBall::Mesh::Builder builder;builder.initialize(3).nodes(coordinates.size());
      for (const auto& x:coordinates) builder.vertex(x);
      for (const auto& f:faces)
      {
        IndexArray v(3);for (size_t j=0;j<3;++j) v(j)=f[j];
        Index index;builder.polytope(Polytope::Type::Triangle,std::move(v),index);
      }
      const auto fullSurface=builder.finalize();
      const auto full=BEM().evaluate(fullSurface,false);
      for (const auto& pair:std::vector<std::pair<Real,Real>>{
        {result.k,full.k},{result.c,full.c},{result.q,full.q},{result.rho,full.rho}})
        referenceError=std::max(referenceError,std::abs(pair.first-pair.second)/std::max(std::abs(pair.second),Real(1e-20)));
      std::cout<<std::setprecision(17)<<"full_surface_reference relative_error="<<referenceError<<std::endl;
    }
    GridFunction direction(space);
    direction=VectorFunction(size_t(3),[](const Point& p) {
      Math::SpatialVector<Real> x=p.vector();
      const Real r=x.norm();if (r==0) return x;
      const Real a=x(0)/r,b=x(1)/r,c=x(2)/r;
      // The even H4 term moves cut vertices as well as patch interiors.
      const Real harmonic=a*a*a*a+b*b*b*b+c*c*c*c-Real(3)/5;
      return Math::SpatialVector<Real>(x*(harmonic+27*a*b*c*(a*a-b*b)*(b*b-c*c)*(c*c-a*a)));
    });
    const Real analytic=load.dot(direction.getData());
    std::cout<<std::setprecision(17)<<"state rho="<<result.rho<<" residual="<<result.residual
      <<" open_edges="<<result.openEdges<<" chamber_panels="<<result.panels
      <<" state_dofs="<<result.stateDOFs<<" independent_dofs="<<result.independentDOFs<<std::endl;
    Real finalError=1;
    for (Real step:{Real(1e-3),Real(5e-4),Real(2.5e-4)})
    {
      Real values[2];size_t side=0;
      for (Real sign:{Real(-1),Real(1)})
      {
        MMG::Mesh perturbed(mesh);
        for (Index v=0;v<mesh.getVertexCount();++v)
        {
          Math::SpatialVector<Real> movement(3);const auto& dofs=space.getDOFs(0,v);
          for (size_t j=0;j<3;++j) movement(j)=direction.getData()(dofs(j));
          perturbed.setVertexCoordinates(v,mesh.getVertexCoordinates(v)+sign*step*movement);
        }
        VelocitySpace movedSpace(perturbed,3);Math::Vector<Real> unused;
        values[side++]=BEM().evaluateChamber(perturbed,movedSpace,unused,false).rho;
      }
      const Real fd=(values[1]-values[0])/(2*step);
      finalError=std::abs(fd-analytic)/std::max(std::abs(analytic),Real(1e-20));
      std::cout<<"finite_difference step="<<step<<" analytic="<<analytic<<" fd="<<fd
        <<" relative_error="<<finalError<<std::endl;
    }
    std::ofstream record("bem-chamber-validation.json");
    const bool passed=result.residual<=1e-9 && result.openEdges==0 && finalError<.002 && referenceError<1e-8;
    record<<std::setprecision(17)<<"{\"passed\":"<<(passed ? "true":"false")
      <<",\"design_domain\":\"O24-chamber\",\"bem_domain\":\"infinite-exterior\""
      <<",\"rho\":"<<result.rho<<",\"residual\":"<<result.residual
      <<",\"chamber_panels\":"<<result.panels<<",\"state_dofs\":"<<result.stateDOFs
      <<",\"independent_dofs\":"<<result.independentDOFs<<",\"open_edges\":"<<result.openEdges
      <<",\"reference_checked\":"<<(reference ? "true":"false")
      <<",\"reference_relative_error\":"<<referenceError
      <<",\"analytic_derivative\":"<<analytic<<",\"fd_relative_error\":"<<finalError<<"}\n";
    if (!record) throw std::runtime_error("Cannot save chamber derivative validation.");
    return passed ? 0 : 2;
  }
  catch (const std::exception& error) { std::cerr<<error.what()<<std::endl;return 1; }
}
