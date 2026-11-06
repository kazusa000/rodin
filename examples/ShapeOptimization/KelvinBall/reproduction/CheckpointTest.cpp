/* Distributed under the Boost Software License, Version 1.0. */
#include "../Checkpoint.h"
#include <iostream>
#include <Rodin/Variational.h>

int main(int argc, char** argv)
{
  using namespace Rodin;
  if (argc == 3 && std::string(argv[1]) == "--load")
  {
    KelvinBall::Checkpoint::loadMesh(argv[2]);
    return 0;
  }
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
    for (size_t d = 0; d < mesh.getSpaceDimension(); ++d)
    {
      if (loaded.getVertexCoordinates(i)(d) != mesh.getVertexCoordinates(i)(d))
        throw std::runtime_error("Native coordinate roundtrip changed doubles.");
    }
  }
  // An already evaluated direction is checkpoint data, indexed by original vertices.
  Variational::P1 space(mesh, 3);
  Variational::GridFunction direction(space);
  std::vector<IO::HDF5::F64> expected(3 * mesh.getVertexCount());
  for (size_t i = 0; i < expected.size(); ++i)
    expected[i] = Real(i) / 7 - Real(1.25);
  {
    const auto file = IO::HDF5::File(H5Fcreate("direction.h5", H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT));
    const auto group = IO::HDF5::Group(H5Gcreate2(file.get(), "/GridFunction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    const auto values = IO::HDF5::Group(H5Gcreate2(file.get(), "/GridFunction/Values", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    IO::HDF5::writeMatrixDataset(file.get(), "/GridFunction/Values/Data", expected, mesh.getVertexCount(), 3);
  }
  KelvinBall::Checkpoint::restoreVector("direction.h5", direction, mesh.getVertexCount());
  for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
  {
    const auto& dofs = space.getDOFs(0, vertex);
    for (size_t component = 0; component < 3; ++component)
      if (direction.getData()(dofs[component]) != expected[3 * vertex + component])
        throw std::runtime_error("Preserved direction changed a component or vertex index.");
  }
  bool rejected = false;
  try { KelvinBall::Checkpoint::restoreVector("direction.h5", direction, mesh.getVertexCount()+1); }
  catch (const std::runtime_error&) { rejected = true; }
  if (!rejected)
    throw std::runtime_error("Invalid saved direction dimensions were accepted.");
  KelvinBall::Checkpoint::save(24, mesh, loaded, {1, 24, 0.1744360039783153});
  std::cout << "Native checkpoint topology, coordinates and MMG-tag roundtrip passed.\n";
}
