/* Distributed under the Boost Software License, Version 1.0. */
#include "PressureGauge.h"
#include <iostream>

using namespace KelvinBall;

int main()
{
  KelvinBall::Mesh::Builder builder;
  builder.initialize(3).nodes(12).reserve(0, 12).reserve(3, 3);
  for (Index i = 0; i < 3; ++i)
  {
    const Real shift = 3 * i;
    builder.vertex({shift, 0, 0}).vertex({shift + 1, 0, 0})
      .vertex({shift, 1, 0}).vertex({shift, 0, 1});
    builder.polytope(Polytope::Type::Tetrahedron, {4*i, 4*i+1, 4*i+2, 4*i+3});
  }
  auto mesh = builder.finalize();
  prepare(mesh);
  for (auto face = mesh.getBoundary(); face; ++face)
    mesh.setAttribute({2, face->getIndex()}, Gamma);
  PressureSpace Qh(mesh);
  const auto references = enclosedPressureReferences(Qh);
  if (references.size() != 3)
    throw std::runtime_error("Missing enclosed-component pressure references.");
  TrialFunction p(Qh);
  TestFunction q(Qh);
  Problem problem(p, q);
  problem = Integral(Grad(p), Grad(q));
  problem.assemble();
  auto& system = problem.getLinearSystem();
  const auto original = system.getOperator();
  Math::Vector<Real> exact(Qh.getSize());
  for (Index v = 0; v < mesh.getVertexCount(); ++v)
  {
    const auto& x = mesh.getVertexCoordinates(v);
    exact(Qh.getDOFs(0, v)(0)) = x(0) - 3 * (v / 4) - x(1) + x(2);
  }
  system.getVector() = original * exact;
  const auto rhs = system.getVector();
  IndexMap<Real> gauges;
  for (Index dof : references)
    gauges.emplace(dof, exact(dof));
  system.eliminate(gauges);
  solveDirect(problem);
  if ((original * system.getSolution() - rhs).norm() > 1e-12 ||
      (system.getSolution() - exact).norm() > 1e-12)
    throw std::runtime_error("Pressure gauge changed the manufactured equations.");
  for (auto face = mesh.getBoundary(); face; ++face)
  {
    const Real x = centroid(mesh, *face)(0);
    if (x > 5) mesh.setAttribute({2, face->getIndex()}, Outer);
    else if (x > 2) mesh.setAttribute({2, face->getIndex()}, SigmaPlus);
  }
  if (enclosedPressureReferences(Qh).size() != 1)
    throw std::runtime_error("Outer/cut-connected components must not be pinned.");
  std::cout << "PASS enclosed pressure gauges, original equations and cut exclusion\n";
}
