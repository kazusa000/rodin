/* Distributed under the Boost Software License, Version 1.0. */
#include <fstream>
#include <iomanip>
#include <iostream>
#include <Rodin/MMG.h>
#include "Thickness.h"

using namespace KelvinBall;

// Fixed-connectivity perturbations below are derivative diagnostics ONLY.
int main(int argc, char** argv)
{
  if (argc == 1)
  {
    const auto& rotations = SewedOutput::getCubeRotations();
    for (const Math::SpatialPoint& x : {Math::SpatialPoint{0.83, 0.24, -0.49},
                                      Math::SpatialPoint{-0.17, 0.64, 0.38}})
    {
      size_t count = 0;
      for (const auto& rotation : rotations)
        count += ThicknessPenalty::inChamber(rotation.transpose() * x);
      if (count != 1 || rotations.size() != ChamberMultiplicity)
        throw std::runtime_error("Thickness chamber/group mismatch.");
    }
    if (Tetrahedral && !ThicknessPenalty::inChamber(Math::SpatialPoint{0.3, 0.8, 0.2}))
      throw std::runtime_error("T12 thickness still uses cubic chamber inequalities.");
    std::cout << "PASS thickness chamber group=" << rotations.size() << '\n';
    return 0;
  }
  std::string path;
  Real minimum = 0, epsilon = 0;
  for (int i = 1; i < argc; ++i)
  {
    const std::string arg(argv[i]);
    if (arg.rfind("--mesh=", 0) == 0) path = arg.substr(7);
    else if (arg.rfind("--minimum=", 0) == 0) minimum = std::stod(arg.substr(10));
    else if (arg.rfind("--fd=", 0) == 0) epsilon = std::stod(arg.substr(5));
    else throw std::runtime_error("Unknown thickness diagnostic option.");
  }
  if (path.empty() || !(minimum > 0) || !std::isfinite(minimum) || epsilon < 0)
    throw std::runtime_error("Provide mesh and positive absolute thickness.");
  const auto evaluate = [&](Real displacement) {
    MMG::Mesh mesh;
    mesh.load(path, IO::FileFormat::MEDIT);
    P1 space(mesh, 3);
    Math::Vector<Real> direction = Math::Vector<Real>::Zero(space.getSize());
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const Math::SpatialVector<Real> x = mesh.getVertexCoordinates(v);
      const Math::SpatialVector<Real> value = x * x.squaredNorm();
      const auto dofs = space.getDOFs(0, v);
      for (Index c = 0; c < 3; ++c) direction(dofs(c)) = value(c);
      if (displacement != 0) mesh.setVertexCoordinates(v, x + displacement * value);
    }
    mesh.getConnectivity().compute(2, 3);
    const Location::AABB<MMG::Mesh> locator(mesh);
    Math::Vector<Real> load = Math::Vector<Real>::Zero(space.getSize());
    const ThicknessPenalty penalty(mesh, minimum, true);
    const auto result = penalty.evaluate(mesh, locator, space, 1, load);
    return std::make_pair(result, -load.dot(direction));
  };
  const auto [result, derivative] = evaluate(0);
  std::ofstream out("thickness.json");
  out << std::setprecision(17) << "{\"minimum\":" << minimum
      << ",\"group_order\":" << ChamberMultiplicity
      << ",\"penalty\":" << result.penalty << ",\"deepest\":" << result.deepest
      << ",\"rays\":" << result.rays << ",\"violating\":" << result.violating
      << ",\"cross_cut_samples\":" << result.crossCut
      << ",\"exact_discrete_derivative\":true,\"derivative\":" << derivative;
  if (epsilon > 0)
  {
    const Real fd = (evaluate(epsilon).first.penalty - evaluate(-epsilon).first.penalty) / (2 * epsilon);
    const Real error = std::abs(fd - derivative) / std::max({std::abs(fd), std::abs(derivative), Real(1e-16)});
    out << ",\"fd\":" << fd << ",\"fd_epsilon\":" << epsilon
        << ",\"fd_relative_error\":" << error;
    if (error > 1e-6) { out << ",\"passed\":false}\n"; return 1; }
  }
  out << ",\"passed\":true}\n";
  if (!out) throw std::runtime_error("Failed writing thickness diagnostic.");
  return 0;
}
