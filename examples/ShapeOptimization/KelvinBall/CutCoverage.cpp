/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iomanip>
#include <iostream>
#include <limits>

#include "Common.h"
#include "RotatedNitscheIntegrator.h"

/**
 * @brief Audits every rotated cut quadrature point without solving or moving a mesh.
 * @param argc Number of command-line arguments.
 * @param argv Arguments: the preserved volume mesh to inspect.
 * @returns Zero for complete coverage, one for a miss, two for invalid arguments.
 */
int main(int argc, char** argv)
{
  using namespace KelvinBall;
  if (argc != 2)
    return 2;
  KelvinBall::Mesh volume;
  volume.load(argv[1], IO::FileFormat::MEDIT);
  prepare(volume);
  auto fluid = volume.trim(Obstacle);
  splitSelfPairedCut(fluid);
  // These are the historical optimizer's unchanged trace tolerances.
  constexpr Real physicalTolerance = 0.01;
  constexpr Real referenceTolerance = 0.25;
  const RotatedNitscheIntegrator coupling(
    fluid, MasterCuts, physicalTolerance, referenceTolerance);
  size_t totalMissed = 0;
  std::cout << std::setprecision(17) << "{\"pairs\":[";
  bool firstPair = true;
  for (const auto& pair : RotationPairs)
  {
    if (!firstPair)
      std::cout << ',';
    firstPair = false;
    typename SubMesh<Context::Local>::Builder builder;
    builder.initialize(fluid);
    for (auto face = fluid.getPolytope(2); face; ++face)
    {
      if (face->getAttribute() == pair.master)
        builder.include(2, face->getIndex());
    }
    auto master = builder.finalize();
    Location::AABB<decltype(master)> exhaustive(master);
    exhaustive.setTolerance(physicalTolerance).setReferenceTolerance(referenceTolerance)
      .setExhaustiveFallback(true);
    size_t found = 0;
    size_t missed = 0;
    std::cout << "{\"slave\":" << pair.slave << ",\"master\":" << pair.master
              << ",\"misses\":[";
    for (auto face = fluid.getPolytope(2); face; ++face)
    {
      if (face->getAttribute() != pair.slave)
        continue;
      const auto& qf = QF::PolytopeQuadratureFormula::get(4, face->getGeometry());
      const auto& quadrature = face->getQuadrature(qf);
      for (size_t qp = 0; qp < quadrature.getSize(); ++qp)
      {
        const Math::SpatialPoint point = pair.rotation * quadrature.getPoint(qp).vector();
        if (coupling.getLocator().locate(pair.master, point))
        {
          ++found;
          continue;
        }
        if (missed++)
          std::cout << ',';
        Real bestMargin = -std::numeric_limits<Real>::infinity();
        Index bestFace = 0;
        Math::SpatialPoint bestReference;
        Real planeResidual = 0;
        for (auto candidate = master.getPolytope(2); candidate; ++candidate)
        {
          Math::SpatialPoint reference, projected;
          candidate->getTransformation().inverse(reference, point);
          const Geometry::Polytope::Traits traits(candidate->getGeometry());
          const auto& hs = traits.getHalfSpace();
          Real margin = std::numeric_limits<Real>::infinity();
          for (Eigen::Index i = 0; i < hs.vector.size(); ++i)
            margin = std::min(margin, hs.vector[i] - reference.dot(hs.matrix.row(i).transpose()));
          if (margin > bestMargin)
          {
            bestMargin = margin;
            bestFace = candidate->getIndex();
            bestReference = reference;
            candidate->getTransformation().transform(projected, reference);
            planeResidual = (point - projected).norm();
          }
        }
        std::cout << "{\"slave_face\":" << face->getIndex() << ",\"quadrature_index\":" << qp
                  << ",\"point\":[" << point[0] << ',' << point[1] << ',' << point[2]
                  << "],\"exhaustive_found\":" << (exhaustive.locate(2, point) ? "true" : "false")
                  << ",\"best_reference_margin\":" << bestMargin
                  << ",\"best_master_face\":" << bestFace
                  << ",\"best_master_reference\":[" << bestReference[0] << ',' << bestReference[1]
                  << "],\"best_plane_residual\":" << planeResidual << ",\"best_master_vertices\":[";
        const auto best = master.getPolytope(2, bestFace);
        bool firstVertex = true;
        for (const Index vertex : best->getVertices())
        {
          if (!firstVertex)
            std::cout << ',';
          firstVertex = false;
          const auto& coordinates = master.getVertexCoordinates(vertex);
          std::cout << '[' << coordinates[0] << ',' << coordinates[1] << ',' << coordinates[2] << ']';
        }
        std::cout << "]}";
      }
    }
    totalMissed += missed;
    std::cout << "],\"found\":" << found << ",\"missed\":" << missed << '}';
  }
  std::cout << "],\"total_missed\":" << totalMissed << "}\n";
  return totalMissed ? 1 : 0;
}
