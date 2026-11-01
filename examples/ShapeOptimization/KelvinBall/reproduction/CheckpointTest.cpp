/* Distributed under the Boost Software License, Version 1.0. */
#include "Configuration.h"
#include "Checkpoint.h"
#include <iostream>

int main(int argc, char** argv)
{
  using namespace Rodin;
  if (argc == 3 && std::string(argv[1]) == "--load")
  {
    KelvinBall::Checkpoint::loadMesh(argv[2]);
    return 0;
  }
  KelvinBall::Configuration config;
  for (const auto value : {"--graded-chamber", "--h=0.16666666666666666",
       "--outer-radius=15", "--near-field-radius=2", "--chamber-gradation=1.3"})
  {
    if (!config.parse(value))
      throw std::runtime_error("Test option was not parsed.");
  }
  config.finalize();
  const auto axis = config.getGridCoordinates();
  if (config.getGridSpacing() != Real(1) / 6 || axis.back() != 15)
    throw std::runtime_error("L15 changed the reference spacing or wall.");
  for (size_t i = 0; i <= 12; ++i)
  {
    if (axis[i] != static_cast<Real>(i) * (Real(1) / 6))
      throw std::runtime_error("Near-field coordinates changed.");
  }
  for (size_t i = 1; i < axis.size(); ++i)
  {
    if (!(axis[i] > axis[i - 1]) || axis[i] - axis[i - 1] > config.hmax + 1e-14)
      throw std::runtime_error("Invalid graded intervals.");
  }
  auto bad = config;
  bad.nearFieldRadius = 1.01;
  bool rejected = false;
  try { bad.finalize(); }
  catch (const std::runtime_error&) { rejected = true; }
  if (!rejected)
    throw std::runtime_error("Fractional near-field intervals were accepted.");

  // A synthetic topology fixture, not a small-box FEM experiment.
  MMG::Mesh mesh(Geometry::Mesh<Context::Local>::UniformGrid(
    Geometry::Polytope::Type::Tetrahedron, {3, 3, 3}));
  mesh.getConnectivity().compute(2, 3);
  mesh.getConnectivity().compute(1, 3);
  mesh.setCorner(0).setRidge(0).setRequiredVertex(1).setRequiredEdge(1)
    .setRequiredTriangle(0).setRequiredTetrahedron(0);
  mesh.setAttribute({0, 0}, 40);
  mesh.setAttribute({2, 0}, 13);
  KelvinBall::Checkpoint::saveMesh("fixture-original.h5", mesh);
  const auto loaded = KelvinBall::Checkpoint::loadMesh("fixture-original.h5");
  KelvinBall::Checkpoint::saveMesh("fixture-restored.h5", loaded);
  if (loaded.getCorners() != mesh.getCorners() || loaded.getRidges() != mesh.getRidges() ||
      loaded.getRequiredVertices() != mesh.getRequiredVertices() ||
      loaded.getRequiredEdges() != mesh.getRequiredEdges() ||
      loaded.getRequiredTriangles() != mesh.getRequiredTriangles() ||
      loaded.getRequiredTetrahedra() != mesh.getRequiredTetrahedra())
    throw std::runtime_error("Native MMG tag roundtrip changed indices.");
  for (Index i = 0; i < mesh.getVertexCount(); ++i)
  {
    if (loaded.getVertexCoordinates(i) != mesh.getVertexCoordinates(i))
      throw std::runtime_error("Native coordinate roundtrip changed doubles.");
  }
  std::cout << "L15 spacing, grading, invalid configuration and native MMG-tag roundtrip passed.\n";
}
