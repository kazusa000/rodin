#include "BEM.h"
#include <Rodin/Solver/PartialPivLU.h>
#include <map>
#include <sstream>
#include <utility>

namespace KelvinBall
{
  BEM::Result BEM::evaluate(const Mesh& surface, bool derivative) const
  {
    if (surface.getDimension() != 2 || surface.getSpaceDimension() != 3 ||
        !surface.getCellCount() || m_order < 1)
      throw std::runtime_error("BEM requires a nonempty triangular surface in three dimensions.");
    Result result;
    result.panels = surface.getCellCount();
    Real totalArea = 0;
    std::map<std::pair<Index, Index>, size_t> edgeCounts;
    for (auto cell = surface.getCell(); cell; ++cell)
    {
      const auto& v = cell->getVertices();
      if (v.size() != 3) throw std::runtime_error("BEM supports affine triangles only.");
      for (size_t i = 0; i < 3; ++i)
        ++edgeCounts[std::minmax(v(i), v((i + 1) % 3))];
      const Math::SpatialVector<Real> a = surface.getVertexCoordinates(v(0));
      const Math::SpatialVector<Real> b = surface.getVertexCoordinates(v(1));
      const Math::SpatialVector<Real> c = surface.getVertexCoordinates(v(2));
      const Real area = (b-a).cross(c-a).norm()/2;
      if (!(area > 0)) throw std::runtime_error("Degenerate BEM panel.");
      totalArea += area;
    }
    for (const auto& [edge, count] : edgeCounts)
    {
      if (count > 2) throw std::runtime_error("Nonmanifold BEM surface.");
      if (count != 2) ++result.openEdges;
    }
    if (result.openEdges) throw std::runtime_error("BEM requires a closed body surface.");
    result.epsilon = m_epsilonFactor * std::sqrt(totalArea / result.panels);
    RegularizedStokeslet kernel(result.epsilon, m_viscosity);
    VelocitySpace space(surface, 3);
    TrialFunction trial(space);
    TestFunction test(space);
    auto form = Integral(Potential(kernel, trial), test);
    form.setOrder(m_order);
    DenseProblem problem(trial, test);
    Solver::PartialPivLU solver(problem);
    Math::Matrix<Real> rhs(space.getSize(), 6), solution(space.getSize(), 6);
    for (size_t load = 0; load < 6; ++load)
    {
      VectorFunction rigid(static_cast<size_t>(3), [load](const Point& p) {
        Math::SpatialVector<Real> velocity(3), axis(3);
        axis.setZero(); axis(load % 3) = 1;
        velocity = load < 3 ? axis : Math::SpatialVector<Real>(axis.cross(p.vector()));
        return velocity;
      });
      auto linear = Integral(rigid, test);
      linear.setOrder(m_order);
      problem = form - linear;
      if (load == 0)
      {
        problem.solve(solver);
        const auto& A = problem.getLinearSystem().getOperator();
        const Real asymmetry = (A-A.transpose()).norm()/std::max(A.norm(), Real(1e-30));
        if (asymmetry > 1e-12)
        {
          std::ostringstream message; message << "BEM operator is not symmetric: " << asymmetry;
          throw std::runtime_error(message.str());
        }
      }
      else
      {
        problem.assemble(AssemblyTarget::RHS);
        auto& system = problem.getLinearSystem();
        system.getSolution() = solver.getSolver().solve(system.getVector());
      }
      const auto& system = problem.getLinearSystem();
      rhs.col(load) = system.getVector();
      solution.col(load) = system.getSolution();
      result.residual = std::max(result.residual,
        (system.getOperator() * solution.col(load) - rhs.col(load)).norm() /
        std::max(rhs.col(load).norm(), Real(1e-30)));
    }
    if (!solution.allFinite() || !std::isfinite(result.residual) || result.residual > 1e-9)
      throw std::runtime_error("Native BEM state failed the original-system residual check.");
    result.resistance = rhs.transpose() * solution;
    result.reciprocity = (result.resistance - result.resistance.transpose()).norm() /
      std::max(result.resistance.norm(), Real(1));
    if (result.reciprocity > 1e-9)
    {
      std::ostringstream message; message << "BEM reciprocity check failed: " << result.reciprocity;
      throw std::runtime_error(message.str());
    }
    result.k = result.resistance.topLeftCorner(3, 3).trace() / 3;
    result.q = result.resistance.bottomRightCorner(3, 3).trace() / 3;
    result.c = (result.resistance.topRightCorner(3, 3).trace() +
                result.resistance.bottomLeftCorner(3, 3).trace()) / 6;
    if (!(result.k > 0 && result.q > 0)) throw std::runtime_error("Invalid BEM resistance.");
    result.rho = std::abs(result.c) / std::sqrt(result.k * result.q);
    if (!derivative) return result;
    if (std::abs(result.c) < 1e-13 * std::sqrt(result.k * result.q))
      throw std::runtime_error("Absolute coupling objective is nondifferentiable at zero coupling.");
    Math::Matrix<Real> W(6, 6);
    W.setZero();
    const Real wk = -result.rho / (2 * result.k), wq = -result.rho / (2 * result.q);
    const Real wc = (result.c > 0 ? 1 : -1) / std::sqrt(result.k * result.q);
    for (size_t i = 0; i < 3; ++i)
    {
      W(i, i) = wk / 3; W(i + 3, i + 3) = wq / 3;
      W(i, i + 3) = W(i + 3, i) = wc / 6;
    }
    struct Sample
    {
      Index cell;
      std::array<Index, 3> vertices;
      std::array<Real, 3> barycentric;
      Math::SpatialVector<Real> x;
      Math::Matrix<Real> density;
      Real weight, area;
      std::array<Math::SpatialVector<Real>, 3> areaGradient;
    };
    std::vector<Sample> samples;
    P1<Real> scalarSpace(surface);
    for (auto cell = surface.getCell(); cell; ++cell)
    {
      const auto& vertices = cell->getVertices();
      const Math::SpatialVector<Real> a = surface.getVertexCoordinates(vertices(0));
      const Math::SpatialVector<Real> b = surface.getVertexCoordinates(vertices(1));
      const Math::SpatialVector<Real> c = surface.getVertexCoordinates(vertices(2));
      Math::SpatialVector<Real> n = (b - a).cross(c - a);
      const Real area = n.norm() / 2;
      if (!(area > 0)) throw std::runtime_error("Degenerate BEM panel.");
      n /= 2 * area;
      const auto& qf = QF::PolytopeQuadratureFormula::get(m_order, cell->getGeometry());
      const auto& quadrature = cell->getQuadrature(qf);
      const auto& fe = space.getFiniteElement(2, cell->getIndex());
      const auto& scalarFE = scalarSpace.getFiniteElement(2, cell->getIndex());
      const auto& dofs = space.getDOFs(2, cell->getIndex());
      for (size_t q = 0; q < quadrature.getSize(); ++q)
      {
        Sample sample;
        sample.cell = cell->getIndex(); sample.area = area;
        sample.x = quadrature.getPoint(q).vector();
        sample.weight = qf.getWeight(q) * quadrature.getPoint(q).getDistortion();
        sample.density.setZero(3, 6);
        for (size_t i = 0; i < fe.getCount(); ++i)
          sample.density += fe.getBasis(i)(qf.getPoint(q)).getData() * solution.row(dofs(i));
        for (size_t i = 0; i < 3; ++i)
        {
          sample.vertices[i] = vertices(i);
          sample.barycentric[i] = scalarFE.getBasis(i)(qf.getPoint(q));
        }
        sample.areaGradient = {Math::SpatialVector<Real>((b-c).cross(n)/2),
          Math::SpatialVector<Real>((c-a).cross(n)/2), Math::SpatialVector<Real>((a-b).cross(n)/2)};
        samples.push_back(std::move(sample));
      }
    }
    result.gradient.setZero(surface.getVertexCount(), 3);
    Real epsilonGradient = 0;
    auto deposit = [&](const Sample& sample, const Math::SpatialVector<Real>& dx, Real dArea) {
      for (size_t i = 0; i < 3; ++i)
        result.gradient.row(sample.vertices[i]) +=
          (sample.barycentric[i] * dx + dArea * sample.areaGradient[i]).getData().transpose();
    };
    for (const auto& x : samples)
    {
      Math::Matrix<Real> rigid(3, 6); rigid.setZero();
      rigid.leftCols(3).setIdentity();
      for (size_t i = 0; i < 3; ++i)
      {
        Math::SpatialVector<Real> axis(3); axis.setZero(); axis(i) = 1;
        rigid.col(i + 3) = axis.cross(x.x).getData();
      }
      const Math::Matrix<Real> P = 2 * x.density * W;
      Math::SpatialVector<Real> dx(3);
      dx.getData() << P(1, 5)-P(2, 4), P(2, 3)-P(0, 5), P(0, 4)-P(1, 3);
      const Real value = (P.cwiseProduct(rigid)).sum() * x.weight;
      deposit(x, Math::SpatialVector<Real>(dx * x.weight), value / x.area);
      for (const auto& y : samples)
      {
        const Math::SpatialMatrix<Real> H = -x.weight * y.weight * x.density * W * y.density.transpose();
        const Math::SpatialVector<Real> r = x.x - y.x;
        Math::SpatialMatrix<Real> K(3, 3); kernel.evaluate(K, r);
        Math::SpatialVector<Real> dr(3); Real de;
        kernel.differential(r, H, dr, de);
        const Real value = H.getData().cwiseProduct(K.getData()).sum();
        deposit(x, dr, value / x.area);
        deposit(y, Math::SpatialVector<Real>(-dr), value / y.area);
        epsilonGradient += de;
      }
    }
    // Epsilon follows the area-based model during deformation, with fixed topology.
    const Real areaCoefficient = epsilonGradient * result.epsilon / (2 * totalArea);
    for (auto cell = surface.getCell(); cell; ++cell)
    {
      const auto& v = cell->getVertices();
      const Math::SpatialVector<Real> a = surface.getVertexCoordinates(v(0));
      const Math::SpatialVector<Real> b = surface.getVertexCoordinates(v(1));
      const Math::SpatialVector<Real> c = surface.getVertexCoordinates(v(2));
      Math::SpatialVector<Real> n = (b-a).cross(c-a); n.normalize();
      result.gradient.row(v(0)) += (areaCoefficient * (b-c).cross(n) / 2).getData().transpose();
      result.gradient.row(v(1)) += (areaCoefficient * (c-a).cross(n) / 2).getData().transpose();
      result.gradient.row(v(2)) += (areaCoefficient * (a-b).cross(n) / 2).getData().transpose();
    }
    if (!result.gradient.allFinite()) throw std::runtime_error("Nonfinite BEM gradient.");
    return result;
  }

  BEM::Result BEM::evaluateChamber(const Mesh& chamber, const VelocitySpace& shapeSpace,
                                  Math::Vector<Real>& load, bool derivative) const
  {
    SewedOutput full(chamber, {Gamma});
    full.getMesh().getConnectivity().compute(2, 3);
    auto solid = full.getMesh().keep(Obstacle);
    solid.getConnectivity().compute(2, 3);
    const auto skin = solid.skin();
    const auto surface = skin.keep(Gamma);
    auto result = evaluate(surface, derivative);
    if (derivative)
    {
      Math::Matrix<Real> fullGradient(full.getMesh().getVertexCount(), 3);
      fullGradient.setZero();
      for (Index i = 0; i < surface.getVertexCount(); ++i)
      {
        const Index sk = surface.getPolytopeMap(0).left[i];
        const Index so = skin.getPolytopeMap(0).left[sk];
        const Index vertex = solid.getPolytopeMap(0).left[so];
        fullGradient.row(vertex) += result.gradient.row(i);
      }
      load = full.pullBackCoordinateGradient(fullGradient, shapeSpace);
    }
    return result;
  }
}
