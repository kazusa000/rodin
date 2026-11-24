/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>

#include "Common.h"

namespace KelvinBall
{
#if defined(KELVIN_AXIAL_C3) || defined(KELVIN_AXIAL_C4)
#ifdef KELVIN_AXIAL_C4
  // Quarter-turn about z maps the positive x half-plane to positive y.
  const std::array<RotationPair, CutPairCount> RotationPairs{
    {{SigmaMinus, SigmaPlus, {0, -1, 0, 1, 0, 0, 0, 0, 1}}}};
#else
  // On the (1,1,1) axis, (x,y,z) -> (z,x,y) maps x=z to x=y.
  const std::array<RotationPair, CutPairCount> RotationPairs{
    {{SigmaMinus, SigmaPlus, {0, 0, 1, 1, 0, 0, 0, 1, 0}}}};
#endif
  const FlatSet<Attribute> MasterCuts{SigmaPlus};
  const FlatSet<Attribute> AllCuts{SigmaPlus, SigmaMinus};
  const FlatSet<Attribute> FixedFaces{Outer, SigmaPlus, SigmaMinus};
#else
  const std::array<RotationPair, CutPairCount> RotationPairs{
    {{SigmaXYMinus, SigmaXYPlus, {0, 1, 0, 1, 0, 0, 0, 0, -1}},
      {SigmaMinus, SigmaPlus, {1, 0, 0, 0, 0, -1, 0, 1, 0}}}};

  const FlatSet<Attribute> MasterCuts{SigmaPlus, SigmaXYPlus};
  const FlatSet<Attribute> AllCuts{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
  const FlatSet<Attribute> FixedFaces{Outer, SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus};
#endif

  RotationPair::RotationPair(
    Attribute slave, Attribute master, std::initializer_list<Real> coefficients)
    : slave(slave),
      master(master),
      rotation(3, 3)
  {
    auto coefficient = coefficients.begin();
    for (size_t row = 0; row < 3; ++row)
    {
      for (size_t column = 0; column < 3; ++column)
        rotation(row, column) = *coefficient++;
    }
  }

  Math::SpatialPoint centroid(const Mesh&, const Polytope& face)
  {
    const Geometry::Point point(face, Polytope::Traits(face.getGeometry()).getCentroid());
    return point.getPhysicalCoordinates();
  }

  Real cellSize(const Polytope& cell)
  {
    // A regular tetrahedron of edge h has measure h^3 / (6 sqrt(2)).
    static const Real scale = std::cbrt(6 * std::sqrt(Real(2)));
    return scale * std::cbrt(cell.getMeasure());
  }

  Real meanElementSize(const Mesh& mesh)
  {
    Real sum = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      const auto& vertices = cell->getVertices();
      if (vertices.size() != 4)
        throw std::runtime_error("Mean element size requires tetrahedra.");
      for (size_t i = 0; i < vertices.size(); ++i)
      {
        for (size_t j = i + 1; j < vertices.size(); ++j)
        {
          sum += (mesh.getVertexCoordinates(vertices[i]) -
            mesh.getVertexCoordinates(vertices[j]))
                   .norm();
        }
      }
    }
    if (mesh.getCellCount() == 0)
      throw std::runtime_error("Mean element size requires a nonempty mesh.");
    return sum / (Real(6) * static_cast<Real>(mesh.getCellCount()));
  }

  void splitSelfPairedCut(Mesh& mesh)
  {
    const size_t faceDimension = mesh.getDimension() - 1;
    for (auto face = mesh.getPolytope(faceDimension); face; ++face)
    {
      if (face->getAttribute() == SigmaXYPlus && centroid(mesh, *face).z() < 0)
        mesh.setAttribute({faceDimension, face->getIndex()}, SigmaXYMinus);
    }
  }

  void prepare(Mesh& mesh)
  {
    splitSelfPairedCut(mesh);
    auto& connectivity = mesh.getConnectivity();
    connectivity.discover(3, 2);
    connectivity.discover(3, 1);
    connectivity.restrict(1, 0);
    connectivity.restrict(2, 0);
    connectivity.restrict(2, 3);
  }
}
