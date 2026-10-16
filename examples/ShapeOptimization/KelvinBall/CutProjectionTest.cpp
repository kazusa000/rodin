/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "RotatedNitscheIntegrator.h"

/**
 * @brief Checks bounded physical trace projection on the recorded thin cut face.
 * @returns Zero when all exact-hit, physical-bound and trace-value checks pass.
 */
int main()
{
  using namespace KelvinBall;
  const auto require = [](bool value, const char* message) {
    if (!value)
      throw std::runtime_error(message);
  };
  // Unmodified master triangle from the global-8 failure. The expected closest
  // point was independently computed with VTK's physical cell locator.
  KelvinBall::Mesh::Builder builder;
  builder.initialize(3).nodes(4)
    .vertex({0.6926848647402944, 0.6926848647402944, 0.07779284284612764})
    .vertex({0.6958727109776132, 0.6958727109776136, 0.07442924982104593})
    .vertex({0.7001091057425602, 0.7001091057425608, 0.0691800572768777})
    .vertex({0.695, 0.6, 0.07})
    .polytope(Geometry::Polytope::Type::Triangle, {0, 1, 2})
    .polytope(Geometry::Polytope::Type::Tetrahedron, {0, 1, 2, 3});
  auto mesh = builder.finalize();
  mesh.setAttribute({2, 0}, SigmaXYPlus);
  const Math::SpatialPoint query{
    0.69585418458382831, 0.69585418458382808, 0.073360980728514602};
  const Math::SpatialPoint expected{
    0.6961160230351721, 0.6961160230351723, 0.07381239142504428};
  constexpr Real h = Real(2) / 19;
  constexpr Real hausdorff = Real(0.1) * h * h;
  AttributeFaceLocator strict(mesh, FlatSet<Attribute>{SigmaXYPlus});
  require(!strict.locate(SigmaXYPlus, query), "Default-off projection accepted an outside point.");
  AttributeFaceLocator projected(mesh, FlatSet<Attribute>{SigmaXYPlus}, 1e-10, 1e-10, hausdorff);
  const auto mapped = projected.locate(SigmaXYPlus, query);
  require(bool(mapped), "The MMG-tolerance rim projection failed.");
  require((mapped->vector() - expected).norm() < 1e-13, "Physical projection disagrees with VTK.");
  require((query - mapped->vector()).norm() <= hausdorff, "Projection exceeded the absolute geometry budget.");
  Math::SpatialPoint reconstructed;
  mapped->getPolytope().getTransformation().transform(reconstructed, mapped->getReferenceCoordinates());
  require((reconstructed - mapped->vector()).norm() < 1e-14,
    "Projected physical and reference coordinates are inconsistent.");
  const Geometry::Polytope::Traits traits(mapped->getPolytope().getGeometry());
  const auto& halfSpace = traits.getHalfSpace();
  for (Eigen::Index i = 0; i < halfSpace.vector.size(); ++i)
  {
    require(halfSpace.vector[i] - mapped->getReferenceCoordinates().dot(
      halfSpace.matrix.row(i).transpose()) >= -1e-14, "Projected trace extrapolates outside the face.");
  }
  P1<Real, KelvinBall::Mesh> space(mesh);
  GridFunction field(space);
  field = F::x + 2 * F::z;
  require(std::abs(field.getValue(*mapped) - (expected[0] + 2 * expected[2])) < 1e-13,
    "A projected P1 trace was evaluated at the original outside coordinate.");
  require(projected.getProjectionHitCount() == 1, "Projection diagnostics were not recorded.");
  AttributeFaceLocator tooSmall(mesh, FlatSet<Attribute>{SigmaXYPlus}, 1e-10, 1e-10, 0.0005);
  require(!tooSmall.locate(SigmaXYPlus, query), "A gap beyond the geometry budget was accepted.");
  require(!projected.locate(SigmaPlus, query), "Projection crossed a cut attribute.");
  const auto face = mesh.getPolytope(2, 0);
  const Math::SpatialPoint interior = centroid(mesh, *face);
  const auto exact = projected.locate(SigmaXYPlus, interior);
  require(bool(exact) && (exact->vector() - interior).norm() < 1e-15,
    "An exact existing trace point moved.");
  require(projected.getProjectionHitCount() == 1, "An exact hit used the projection fallback.");
  std::cout << "cut_projection_distance=" << projected.getMaximumProjectionDistance()
            << " budget=" << hausdorff << " exact_hit_unchanged=1 outside_rejected=1\n";
  return 0;
}
