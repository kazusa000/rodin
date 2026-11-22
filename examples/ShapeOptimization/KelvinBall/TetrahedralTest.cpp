#include "SewedOutput.h"
#include "Sphere.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

using namespace KelvinBall;

namespace
{
  void require(bool condition, const char* message)
  {
    if (!condition)
      throw std::runtime_error(message);
  }

  // Centroids on the uniform grid have coordinates which are thirds of a grid step.
  using Key = std::array<long long, 3>;
  Key key(const Math::SpatialPoint& x, Real spacing)
  {
    return {std::llround(3 * x(0) / spacing), std::llround(3 * x(1) / spacing),
      std::llround(3 * x(2) / spacing)};
  }
}

int main()
{
  const auto& rotations = SewedOutput::getRotations();
  require(ChamberMultiplicity == 12 && rotations.size() == 12, "Wrong group order");
  std::set<Key> orbit;
  for (const auto& rotation : rotations)
  {
    require(rotation.determinant() == 1, "Improper rotation");
    const Math::SpatialPoint v = rotation * Math::SpatialPoint{1, 2, 3};
    orbit.insert(key(v, 3));
    const Math::SpatialPoint corner = rotation * Math::SpatialPoint{1, 1, 1};
    require(corner(0) * corner(1) * corner(2) == 1, "Rotation does not preserve tetrahedron");
  }
  require(orbit.size() == 12, "Duplicate rotations");
  const Math::SpatialPoint diagonal{1, 1, 1};
  const Real amplitude = Real(0.02) / std::sqrt(Real(2));
  require(std::abs(initialRadius(diagonal) - 1 - amplitude) < 1e-13, "Wrong cubic harmonic");
  // At (1,2,0) the cubic harmonic vanishes, leaving only the degree-six term.
  const Math::SpatialPoint planar{1, 2, 0};
  const Real planarH6 = 72 * std::sqrt(Real(3)) / 125;
  require(std::abs(initialRadius(planar) - 1 - amplitude * planarH6) < 1e-13,
    "Wrong degree-six harmonic");
  for (const auto& rotation : rotations)
  {
    for (const Math::SpatialPoint u : {diagonal, planar, Math::SpatialPoint{1, 2, 3}})
    {
      const Math::SpatialPoint v = rotation * u;
      require(std::abs(initialRadius(v) - initialRadius(u)) < 1e-13,
        "Tetrahedral rotation changed the initial radius");
    }
    const Math::SpatialPoint v = rotation * diagonal;
    const Math::SpatialPoint quarterTurn{-v(1), v(0), v(2)};
    require(std::abs(initialRadius(diagonal) - initialRadius(quarterTurn) - 2 * amplitude) < 1e-13,
      "An extra cubic rotation preserved the initial radius");
  }
  std::cout << "PASS: archived 0.02 seed, all 12 tetrahedral rotations preserved, "
               "all 12 extra cubic rotations broken.\n";
  Configuration config;
  config.points = 7;
  config.outerRadius = 10;
  config.finalize();
  const Sphere sphere(config);
  MMG::Mesh chamber(sphere.makeUniformChamber());
  sphere.protectFixedGeometry(chamber, false);
  sphere.projectFixedGeometry(chamber);
  Real volume = 0;
  for (auto cell = chamber.getCell(); cell; ++cell)
  {
    volume += cell->getMeasure();
    const auto x = centroid(chamber, *cell);
    size_t inside = 0;
    for (const auto& rotation : rotations)
    {
      const Math::SpatialPoint y = rotation * x;
      inside += y(0) >= std::abs(y(2)) && y(1) >= std::abs(y(2));
    }
    require(inside == 1, "Rotated chambers overlap or have gaps");
  }
  require(std::abs(12 * volume - std::pow(2 * config.outerRadius, 3)) < 1e-8,
    "Twelve chambers must fill the outer box");

  std::map<Attribute, std::set<Key>> faceCenters;
  std::set<ChamberPlane> planes;
  for (auto face = chamber.getBoundary(); face; ++face)
  {
    const auto plane = fixedPlane(chamber, *face);
    require(plane != ChamberPlane::None, "Missing boundary plane");
    planes.insert(plane);
    for (const auto vertex : face->getVertices())
      require(std::abs(planeNormal(plane).dot(chamber.getVertexCoordinates(vertex)) -
                         planeOffset(plane, config.outerRadius)) < 1e-12,
        "Vertex off its boundary plane");
    faceCenters[*face->getAttribute()].insert(key(centroid(chamber, *face), config.getGridSpacing()));
  }
  require(planes.size() == 6, "Expected two outer faces and four cuts");
  for (const auto& pair : RotationPairs)
  {
    require(MasterCuts.contains(pair.master), "Wrong master cut");
    std::set<Key> mapped;
    for (const auto& coords : faceCenters.at(pair.slave))
    {
      const Math::SpatialPoint x{Real(coords[0]), Real(coords[1]), Real(coords[2])};
      mapped.insert(key(pair.rotation * x, 3));
    }
    require(mapped == faceCenters.at(pair.master), "Paired cut triangulations do not match");
  }
  std::cout << "PASS: 12 proper tetrahedral rotations, box tiling, chamber volume, "
               "six fixed planes, and both rotated cut triangulations.\n";
}
