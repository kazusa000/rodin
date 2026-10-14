#include "Common.h"
#include "SewedOutput.h"
#include "RotatedNitscheIntegrator.h"
#include <iostream>
#include <Rodin/Advection/Lagrangian.h>

using namespace KelvinBall;

int main()
{
  try
  {
    // Two tetrahedra form the canonical O24 chamber, with matching planar traces.
    KelvinBall::Mesh::Builder builder;
    builder.initialize(3).nodes(5);
    for (const Math::SpatialPoint x : std::vector<Math::SpatialPoint>{
      {0,0,0}, {2,0,0}, {2,2,0}, {2,2,2}, {2,2,-2}}) builder.vertex(x);
    for (const auto& tetra : std::vector<std::array<Index,4>>{{0,1,2,3},{0,1,4,2}})
    {
      IndexArray vertices(4);
      for (size_t i=0; i<4; ++i) vertices(i)=tetra[i];
      Index index;
      builder.polytope(Polytope::Type::Tetrahedron, std::move(vertices), index);
      builder.attribute({3,index}, Fluid);
    }
    auto mesh = builder.finalize();
    mesh.getConnectivity().compute(2,3);
    for (auto face=mesh.getBoundary(); face; ++face)
    {
      const auto x=centroid(mesh,*face);
      Attribute attribute;
      if (std::abs(x(0)-2) < 1e-12) attribute=Outer;
      else if (std::abs(x(0)-x(1)) < 1e-12) attribute=x(2)<0 ? SigmaXYMinus : SigmaXYPlus;
      else attribute=x(2)<0 ? SigmaMinus : SigmaPlus;
      mesh.setAttribute({2,face->getIndex()},attribute);
    }
    mesh.getConnectivity().compute(3,2);
    P1 scalarSpace(mesh);
    VelocitySpace vectorSpace(mesh, 3);
    GridFunction scalar(scalarSpace);
    GridFunction vector(vectorSpace);
    for (Eigen::Index i = 0; i < scalar.getData().size(); ++i)
      scalar.getData()(i) = std::sin(Real(i + 1) * .73);
    for (Eigen::Index i = 0; i < vector.getData().size(); ++i)
      vector.getData()(i) = std::sin(Real(i + 1) * .51);
    GridFunction coordinate(scalarSpace);
    coordinate=RealFunction([](const Point& p) { return p.x(); });
    const auto cell=mesh.getCell();
    const Point point(*cell,Math::SpatialPoint{.25,.25,.25});
    const auto backtrace=Flow(Real(-.01),coordinate,VectorFunction{1,0,0},Math::RungeKutta::RK4{});
    if (std::abs(backtrace.getValue(point)-(point.x()-.01))>1e-11)
      throw std::runtime_error("Negative-time pullback must advance the boundary along +velocity.");
    const SewedOutput quotient(mesh);
    VelocitySpace fullSpace(quotient.getMesh(), 3);
    GridFunction fullVector(fullSpace);
    quotient.setVector(fullVector, vector);
    Math::Matrix<Real> fullGradient(quotient.getMesh().getVertexCount(), 3);
    for (Index i = 0; i < quotient.getMesh().getVertexCount(); ++i)
      for (size_t j = 0; j < 3; ++j) fullGradient(i, j) = std::cos(Real(3*i+j+1)*.37);
    const auto load = quotient.pullBackCoordinateGradient(fullGradient, vectorSpace);
    const Real left = load.dot(vector.getData());
    Real right = 0;
    for (Index i = 0; i < quotient.getMesh().getVertexCount(); ++i)
    {
      const auto& dofs = fullSpace.getDOFs(0, i);
      for (size_t j = 0; j < 3; ++j) right += fullGradient(i, j)*fullVector.getData()(dofs(j));
    }
    if (std::abs(left-right) > 1e-11*std::max(Real(1),std::abs(right)))
      throw std::runtime_error("BEM coordinate pullback is not the adjoint of averaged transport.");
    quotient.enforceScalarRepresentatives(scalar);
    quotient.enforceVectorRepresentatives(vector);
    const Math::Vector<Real> firstScalar = scalar.getData(), firstVector = vector.getData();
    quotient.enforceScalarRepresentatives(scalar);
    quotient.enforceVectorRepresentatives(vector);
    if ((firstScalar-scalar.getData()).norm() > 1e-12 ||
        (firstVector-vector.getData()).norm() > 1e-12)
      throw std::runtime_error("Quotient trace projection is not idempotent.");
    RotatedNitscheIntegrator coupling(mesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus}, .01, .25);
    const Real scalarJump = coupling.scalarJump(scalar), vectorJump = coupling.vectorJump(vector);
    if (scalarJump > 1e-11 || vectorJump > 1e-11)
      throw std::runtime_error("Projected canonical-chamber traces do not match.");
    std::cout << "adjoint_error=" << std::abs(left-right) << " scalar_jump=" << scalarJump
      << " vector_jump=" << vectorJump << '\n';
    return 0;
  }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
