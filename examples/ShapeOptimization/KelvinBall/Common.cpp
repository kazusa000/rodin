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
#ifdef KELVIN_TETRAHEDRAL
  // y=z -> x=z, and x=-z -> y=-z. No self-paired x=y cut.
  const std::array<RotationPair, 2> RotationPairs{
    {{SigmaPlus, SigmaMinus, {0, 0, 1, 1, 0, 0, 0, 1, 0}},
      {SigmaXYPlus, SigmaXYMinus, {0, 1, 0, 0, 0, -1, -1, 0, 0}}}};
  const FlatSet<Attribute> MasterCuts{SigmaMinus, SigmaXYMinus};
#else
  const std::array<RotationPair, 2> RotationPairs{
    {{SigmaXYMinus, SigmaXYPlus, {0, 1, 0, 1, 0, 0, 0, 0, -1}},
      {SigmaMinus, SigmaPlus, {1, 0, 0, 0, 0, -1, 0, 1, 0}}}};

  const FlatSet<Attribute> MasterCuts{SigmaPlus, SigmaXYPlus};
#endif

  Real initialRadius(const Math::SpatialPoint& x)
  {
    if (!Tetrahedral || x.norm() == 0)
      return 1;
    const Math::SpatialPoint u = x / x.norm();
    const Real a = u(0), b = u(1), c = u(2);
    const Real h3 = 3 * std::sqrt(Real(3)) * a * b * c;
    const Real h6 = 6 * std::sqrt(Real(3)) *
      (a*a - b*b) * (b*b - c*c) * (c*c - a*a);
    return 1 + Real(0.02) * (h3 + h6) / std::sqrt(Real(2));
  }

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

  ChamberPlane fixedPlane(const Mesh& mesh, const Polytope& face)
  {
    using enum ChamberPlane;
    const auto label = face.getAttribute();
    if (label == Outer)
    {
      if (Tetrahedral)
      {
        const auto c = centroid(mesh, face);
        if (c(1) > c(0))
          return OuterY;
      }
      return OuterX;
    }
    if (label == SigmaPlus) return YMinusZ;
    if (label == SigmaMinus) return Tetrahedral ? XMinusZ : YPlusZ;
    if (label == SigmaXYPlus) return Tetrahedral ? XPlusZ : XY;
    if (label == SigmaXYMinus) return Tetrahedral ? YPlusZ : XY;
    return None;
  }

  Math::SpatialPoint planeNormal(ChamberPlane plane)
  {
    using enum ChamberPlane;
    switch (plane)
    {
      case OuterX: return {1, 0, 0};
      case OuterY: return {0, 1, 0};
      case XY: return {1, -1, 0};
      case YMinusZ: return {0, 1, -1};
      case YPlusZ: return {0, 1, 1};
      case XMinusZ: return {1, 0, -1};
      case XPlusZ: return {1, 0, 1};
      default: throw std::logic_error("Face has no prescribed chamber plane.");
    }
  }

  Real planeOffset(ChamberPlane plane, Real outerRadius)
  {
    return plane == ChamberPlane::OuterX || plane == ChamberPlane::OuterY ? outerRadius : 0;
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
    if (Tetrahedral)
      return;
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
