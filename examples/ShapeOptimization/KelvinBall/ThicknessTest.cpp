#include "BaselineThickness.h"
#include <iostream>
#include <stdexcept>

using namespace KelvinBall;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }
}

int main()
{
  KelvinBall::Mesh mesh = KelvinBall::Mesh::UniformGrid(Polytope::Type::Tetrahedron, {5, 5, 5});
  mesh.scale(Real(0.5));
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    Math::SpatialPoint x = Math::SpatialPoint::Zero(3);
    for (const auto vertex : cell->getVertices())
      x += mesh.getVertexCoordinates(vertex);
    x /= Real(cell->getVertices().size());
    const bool solid = x(0) > .5 && x(0) < 1.5 && x(1) > .5 && x(1) < 1.5
      && x(2) > .5 && x(2) < 1.5;
    mesh.setAttribute({3, cell->getIndex()}, solid ? Obstacle : Fluid);
  }
  mesh.getConnectivity().compute(2, 3);
  for (auto face = mesh.getPolytope(2); face; ++face)
  {
    const auto& cells = mesh.getConnectivity().getIncidence({2, 3}, face->getIndex());
    if (cells.size() == 2 && mesh.getPolytope(3, cells[0])->getAttribute()
        != mesh.getPolytope(3, cells[1])->getAttribute())
      mesh.setAttribute({2, face->getIndex()}, Gamma);
  }
  VelocitySpace space(mesh, 3);
  Math::Vector<Real> load = Math::Vector<Real>::Zero(space.getSize());
  Math::Vector<Real> nominalLoad = load;
  BaselineThicknessPenalty guarded(mesh, Real(1.5));
  const auto result = guarded.evaluate(mesh, space, load, Real(1.2));
  BaselineThicknessPenalty nominal(mesh, Real(1.2));
  const auto separate = nominal.evaluate(mesh, space, nominalLoad);
  require(result.samples > 0 && result.violating > 0, "no short rays in the unit cube");
  require(result.penalty > 0 && result.nominalPenalty > 0 && load.norm() > 0,
    "missing thickness penalty or paired repair load");
  require(std::abs(result.minimumExit - 1) < 1e-10, "first-exit distance is not cube thickness");
  require(std::abs(result.nominalPenalty - separate.penalty) < 1e-12,
    "nominal penalty depends on guard distance");
  Math::Vector<Real> disabledLoad = Math::Vector<Real>::Zero(space.getSize());
  BaselineThicknessPenalty inactive(mesh, Real(.5));
  const auto inactiveResult = inactive.evaluate(mesh, space, disabledLoad);
  require(inactiveResult.penalty == 0 && disabledLoad.norm() == 0,
    "thick-enough geometry received an artificial repair");
  Math::Vector<Real> expansion = Math::Vector<Real>::Zero(space.getSize());
  for (Index v = 0; v < mesh.getVertexCount(); ++v)
  {
    const auto dofs = space.getDOFs(0, v);
    for (size_t d = 0; d < 3; ++d)
      expansion(dofs(d)) = mesh.getVertexCoordinates(v)(d) - Real(1);
  }
  require(load.dot(expansion) > 0, "paired repair opposes expansion of a thin cube");
  bool rejected = false;
  try { guarded.evaluate(mesh, space, load, Real(2)); }
  catch (const std::runtime_error&) { rejected = true; }
  require(rejected, "nominal distance larger than guard was accepted");
  std::cout << "PASS baseline first-exit, nominal/guard, paired load and inactive case\n";
}
