/*
 * 3Dv8_test: 3Dv6 null-space gradient flow with
 *   - volume equality constraint
 *   - surface-area upper-bound inequality constraint
 *   - centroid-based rigid-body rotations
 *   - optional whole-mesh translation via SHIFT_{X,Y,Z}
 *   - same-phase internal face cleanup before solves
 *   - no MEDIT round-trip in the main iteration state
 *   - checked PETSc KSP solves for fail-fast diagnostics
 *
 * This is intentionally a small, test-only stability variant of 3Dv6. MEDIT
 * meshes are still written as artifacts, but the next iteration continues from
 * the in-memory mesh to avoid text I/O precision loss.
 */

#include <Rodin/MMG.h>
#include <Rodin/PETSc.h>
#include <Rodin/IO/XDMF.h>

#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Advection/Lagrangian.h>

#include <petscksp.h>
#include <petscsys.h>

#include "RuntimeConfig.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef LEVELSET_STOKES_MESH_FILE
#define LEVELSET_STOKES_MESH_FILE "3dshapes/sphere_init_aligned.mesh"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using MeshType = Rodin::MMG::Mesh;

// Phase labels
static constexpr Attribute Obstacle   = 2;
static constexpr Attribute Fluid      = 3;
static constexpr Attribute GammaShape = 13;

// Stokes viscosity
static constexpr Real mu = 1.0;

// Optimization parameters
static constexpr size_t defaultMaxIt   = 5;
static constexpr Real   defaultHMax    = 0.2;
static constexpr Real   defaultAlpha   = 0.2;
static constexpr Real   defaultDt      = 0.5 * (defaultHMax - 0.1 * defaultHMax);
static constexpr Real   regularization = 1e-12;
static constexpr size_t defaultConvergenceWindow = 5;
static constexpr Real   defaultConvergenceRtolJraw = 5e-2;
static constexpr Real   defaultNsAlphaJ = 0.5;
static constexpr Real   defaultNsAlphaC = 0.5;
static constexpr Real   defaultSurfaceAreaFactor = 1.05;
static constexpr Real   areaActiveToleranceRatio = 1e-4;
static constexpr Real   gramDetToleranceRatio = 1e-12;
static constexpr Real   defaultAreaCorrectionGain = 0.1;
static constexpr Real   areaMultiplierTolerance = 1e-10;
static constexpr Real defaultPcFactorShiftAmount = 1e-10;
static constexpr const char* defaultKspType = "gmres";
static constexpr const char* defaultKspPcType = "bjacobi";
static constexpr const char* defaultKspPcFactorSolverType = "";
static constexpr size_t defaultKspMaxIt = 500;
static constexpr Real defaultKspRtol = 1e-8;
static constexpr const char* defaultStokesKspType = "fgmres";
static constexpr const char* defaultStokesPcType = "fieldsplit";
static constexpr const char* defaultStokesFieldSplitType = "schur";
static constexpr const char* defaultStokesSchurFactType = "lower";
static constexpr const char* defaultStokesSchurPrecondition = "selfp";
static constexpr const char* defaultStokesVelocityPcType = "lu";
static constexpr const char* defaultStokesPressurePcType = "jacobi";
static constexpr Real defaultStokesRtol = 1e-5;
static constexpr size_t defaultStokesMaxIt = 2000;
static constexpr Real defaultQualityMinTetVolume = 1e-15;
static constexpr Real defaultQualityMinAltitudeRatio = 1e-10;
static constexpr Real defaultQualityDuplicateTol = 1e-13;
static constexpr size_t defaultRemeshQualityMaxRetries = 6;
static constexpr Real defaultJLineSearchRtol = 1e-2;
static constexpr Real defaultAdaptGradation = 1.3;
static constexpr const char* algorithmName = "3Dv12_test";

struct MeshQualityConfig
{
  Real minTetVolume = defaultQualityMinTetVolume;
  Real minAltitudeRatio = defaultQualityMinAltitudeRatio;
  Real duplicateTol = defaultQualityDuplicateTol;
};

struct MeshQualityStats
{
  size_t vertexCount = 0;
  size_t cellCount = 0;
  size_t duplicateGroups = 0;
  size_t duplicateVertices = 0;
  size_t badVolumeCells = 0;
  size_t badAltitudeCells = 0;
  Real minTetVolume = std::numeric_limits<Real>::infinity();
  Real minAltitudeRatio = std::numeric_limits<Real>::infinity();

  bool accepted() const
  {
    return duplicateGroups == 0 && badVolumeCells == 0 && badAltitudeCells == 0;
  }
};

struct InterfaceFaceStats
{
  size_t gammaFaces = 0;
  size_t gammaBoundaryFaces = 0;
  size_t gammaBadAdjacencyFaces = 0;
  size_t gammaSamePhaseFaces = 0;
  size_t gammaFluidObstacleFaces = 0;
  Real minGammaArea = std::numeric_limits<Real>::infinity();
  size_t tinyGammaAreaFaces = 0;
};

struct Vec3
{
  Real x = 0.0;
  Real y = 0.0;
  Real z = 0.0;
};

struct MeshCheckpointStats
{
  size_t vertexCount = 0;
  size_t cellCount = 0;
  size_t faceCount = 0;
  double obstacleVolume = std::numeric_limits<double>::quiet_NaN();
  double interfaceArea = std::numeric_limits<double>::quiet_NaN();
  Vec3 centroid;
  Vec3 bboxMin;
  Vec3 bboxMax;
  MeshQualityStats quality;
};

static Vec3 toVec3(const Math::SpatialPoint& p)
{
  Vec3 out;
  if (p.size() > 0)
    out.x = p(0);
  if (p.size() > 1)
    out.y = p(1);
  if (p.size() > 2)
    out.z = p(2);
  return out;
}

static Vec3 operator-(const Vec3& a, const Vec3& b)
{
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static Vec3 cross(const Vec3& a, const Vec3& b)
{
  return {
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  };
}

static Real dot(const Vec3& a, const Vec3& b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

static Real norm(const Vec3& a)
{
  return std::sqrt(dot(a, a));
}

static Real tetVolume(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d)
{
  return std::abs(dot(b - a, cross(c - a, d - a))) / 6.0;
}

static Real triangleArea(const Vec3& a, const Vec3& b, const Vec3& c)
{
  return 0.5 * norm(cross(b - a, c - a));
}

static InterfaceFaceStats analyzeGammaShapeFaces(MeshType& mesh)
{
  auto& conn = mesh.getConnectivity();
  conn.compute(2, 3);
  conn.compute(3, 2);

  InterfaceFaceStats stats;
  for (auto fit = mesh.getFace(); fit; ++fit)
  {
    const auto faceAttr = fit->getAttribute();
    if (!faceAttr || *faceAttr != GammaShape)
      continue;

    ++stats.gammaFaces;
    if (fit->isBoundary())
      ++stats.gammaBoundaryFaces;

    const auto& vertices = fit->getVertices();
    if (vertices.size() == 3)
    {
      const Real area =
        triangleArea(
          toVec3(mesh.getVertexCoordinates(vertices[0])),
          toVec3(mesh.getVertexCoordinates(vertices[1])),
          toVec3(mesh.getVertexCoordinates(vertices[2])));
      stats.minGammaArea = std::min(stats.minGammaArea, area);
      if (!std::isfinite(area) || area <= Real(1e-14))
        ++stats.tinyGammaAreaFaces;
    }

    const auto& cells = conn.getIncidence({2, 3}, fit->getIndex());
    if (cells.size() != 2)
    {
      ++stats.gammaBadAdjacencyFaces;
      continue;
    }

    const auto c0 = mesh.getAttribute(3, cells[0]);
    const auto c1 = mesh.getAttribute(3, cells[1]);
    if (!c0 || !c1)
    {
      ++stats.gammaBadAdjacencyFaces;
      continue;
    }

    if (*c0 == *c1)
    {
      ++stats.gammaSamePhaseFaces;
      continue;
    }

    if ((*c0 == Fluid && *c1 == Obstacle) || (*c0 == Obstacle && *c1 == Fluid))
      ++stats.gammaFluidObstacleFaces;
    else
      ++stats.gammaBadAdjacencyFaces;
  }

  if (stats.gammaFaces == 0)
    stats.minGammaArea = 0.0;
  return stats;
}

static void logGammaShapeStats(MeshType& mesh, const std::string& label)
{
  const auto stats = analyzeGammaShapeFaces(mesh);
  Alert::Info()
    << "   | GammaShape stats " << label
    << ": faces=" << stats.gammaFaces
    << ", fluid_obstacle=" << stats.gammaFluidObstacleFaces
    << ", boundary=" << stats.gammaBoundaryFaces
    << ", bad_adjacency=" << stats.gammaBadAdjacencyFaces
    << ", same_phase=" << stats.gammaSamePhaseFaces
    << ", min_area=" << stats.minGammaArea
    << ", tiny_area_faces=" << stats.tinyGammaAreaFaces
    << Alert::Raise;
}

static Real tetMinAltitudeRatio(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d)
{
  const Real volume = tetVolume(a, b, c, d);
  const std::array<Real, 6> edges{
    norm(a - b), norm(a - c), norm(a - d),
    norm(b - c), norm(b - d), norm(c - d)
  };
  const Real maxEdge = *std::max_element(edges.begin(), edges.end());
  if (maxEdge <= 0.0)
    return 0.0;

  const std::array<Real, 4> faceAreas{
    triangleArea(b, c, d),
    triangleArea(a, c, d),
    triangleArea(a, b, d),
    triangleArea(a, b, c)
  };

  Real minAltitude = std::numeric_limits<Real>::infinity();
  for (const Real area : faceAreas)
  {
    if (area <= 0.0)
      return 0.0;
    minAltitude = std::min(minAltitude, 3.0 * volume / area);
  }

  return minAltitude / maxEdge;
}

static MeshQualityStats evaluateMeshQuality(
  const MeshType& mesh,
  const MeshQualityConfig& config)
{
  MeshQualityStats stats;

  std::vector<std::array<long long, 3>> duplicateKeys;
  for (auto vit = mesh.getVertex(); vit; ++vit)
  {
    ++stats.vertexCount;
    const auto p = mesh.getVertexCoordinates(vit->getIndex());
    for (Eigen::Index i = 0; i < p.size(); ++i)
    {
      if (!std::isfinite(p(i)))
      {
        ++stats.duplicateGroups;
        ++stats.duplicateVertices;
        break;
      }
    }

    const Real tol = std::max(config.duplicateTol, Real(1e-300));
    duplicateKeys.push_back({
      static_cast<long long>(std::llround((p.size() > 0 ? p(0) : 0.0) / tol)),
      static_cast<long long>(std::llround((p.size() > 1 ? p(1) : 0.0) / tol)),
      static_cast<long long>(std::llround((p.size() > 2 ? p(2) : 0.0) / tol))
    });
  }

  std::sort(duplicateKeys.begin(), duplicateKeys.end());
  for (size_t i = 0; i < duplicateKeys.size();)
  {
    size_t j = i + 1;
    while (j < duplicateKeys.size() && duplicateKeys[j] == duplicateKeys[i])
      ++j;
    if (j - i > 1)
    {
      ++stats.duplicateGroups;
      stats.duplicateVertices += j - i;
    }
    i = j;
  }

  for (auto cit = mesh.getCell(); cit; ++cit)
  {
    ++stats.cellCount;
    const auto& vertices = cit->getVertices();
    if (vertices.size() != 4)
    {
      ++stats.badVolumeCells;
      ++stats.badAltitudeCells;
      continue;
    }

    const Vec3 a = toVec3(mesh.getVertexCoordinates(vertices[0]));
    const Vec3 b = toVec3(mesh.getVertexCoordinates(vertices[1]));
    const Vec3 c = toVec3(mesh.getVertexCoordinates(vertices[2]));
    const Vec3 d = toVec3(mesh.getVertexCoordinates(vertices[3]));

    const Real volume = tetVolume(a, b, c, d);
    const Real altitudeRatio = tetMinAltitudeRatio(a, b, c, d);
    stats.minTetVolume = std::min(stats.minTetVolume, volume);
    stats.minAltitudeRatio = std::min(stats.minAltitudeRatio, altitudeRatio);

    if (!std::isfinite(volume) || volume <= config.minTetVolume)
      ++stats.badVolumeCells;
    if (!std::isfinite(altitudeRatio) || altitudeRatio <= config.minAltitudeRatio)
      ++stats.badAltitudeCells;
  }

  if (stats.cellCount == 0)
  {
    stats.minTetVolume = 0.0;
    stats.minAltitudeRatio = 0.0;
    ++stats.badVolumeCells;
    ++stats.badAltitudeCells;
  }

  return stats;
}

static void logMeshQuality(
  const MeshQualityStats& stats,
  const std::string& label,
  bool accepted)
{
  Alert::Info()
    << "   | Mesh quality " << label
    << ": accepted=" << (accepted ? 1 : 0)
    << ", vertices=" << stats.vertexCount
    << ", cells=" << stats.cellCount
    << ", min_tet_volume=" << stats.minTetVolume
    << ", min_altitude_ratio=" << stats.minAltitudeRatio
    << ", bad_volume_cells=" << stats.badVolumeCells
    << ", bad_altitude_cells=" << stats.badAltitudeCells
    << ", duplicate_groups=" << stats.duplicateGroups
    << ", duplicate_vertices=" << stats.duplicateVertices
    << Alert::Raise;
}

static double computeInterfaceArea(MeshType& mesh)
{
  P0g p0(mesh);
  TestFunction z0(p0);
  LinearForm lfArea(z0);
  lfArea = FaceIntegral(z0).over({GammaShape});
  lfArea.assemble();

  GridFunction one(p0);
  one = 1.0;
  return lfArea(one);
}

static size_t clearSamePhaseInternalFaces(
  MeshType& mesh,
  const std::string& label = "",
  bool alwaysLog = false)
{
  auto& conn = mesh.getConnectivity();
  conn.compute(2, 3);
  conn.compute(3, 2);

  size_t cleared = 0;
  for (auto fit = mesh.getFace(); fit; ++fit)
  {
    const auto faceAttr = fit->getAttribute();
    if (!faceAttr || fit->isBoundary())
      continue;

    const auto& cells = conn.getIncidence({2, 3}, fit->getIndex());
    if (cells.size() != 2)
      continue;

    const auto c0 = mesh.getAttribute(3, cells[0]);
    const auto c1 = mesh.getAttribute(3, cells[1]);
    if (c0 && c1 && *c0 == *c1)
    {
      mesh.setAttribute(fit.key(), {});
      ++cleared;
    }
  }

  if (cleared > 0 || alwaysLog)
  {
    Alert::Info()
      << "   | Same-phase internal face cleanup"
      << (label.empty() ? "" : " ")
      << label
      << ": cleared=" << cleared
      << Alert::Raise;
  }

  return cleared;
}

static void saveArtifactMesh(
  MeshType& mesh,
  const std::filesystem::path& path,
  const std::string& label)
{
  clearSamePhaseInternalFaces(mesh, label + " before save", true);
  mesh.save(path.string(), IO::FileFormat::MEDIT);

  Alert::Info()
    << "   | " << label
    << " saved as MEDIT artifact; continuing with in-memory mesh"
    << Alert::Raise;
}

template <class ProblemType>
static void solveKSPChecked(
  ProblemType& problem,
  const std::string& label,
  PetscReal rtol,
  PetscInt maxIt,
  const Optional<std::string>& prefix = {})
{
  Solver::KSP solver(problem);
  solver.setPrefix(prefix);
  solver.setTolerances(rtol, PETSC_DEFAULT, PETSC_DEFAULT, maxIt);
  solver.solve();

  ::KSPConvergedReason reason = KSP_CONVERGED_ITERATING;
  PetscInt iterations = 0;
  PetscReal residual = 0.0;

  PetscErrorCode ierr = KSPGetConvergedReason(solver.getHandle(), &reason);
  if (ierr != PETSC_SUCCESS)
    throw std::runtime_error(label + " KSPGetConvergedReason failed.");

  ierr = KSPGetIterationNumber(solver.getHandle(), &iterations);
  if (ierr != PETSC_SUCCESS)
    throw std::runtime_error(label + " KSPGetIterationNumber failed.");

  ierr = KSPGetResidualNorm(solver.getHandle(), &residual);
  if (ierr != PETSC_SUCCESS)
    throw std::runtime_error(label + " KSPGetResidualNorm failed.");

  ::PC pc = PETSC_NULLPTR;
  ::PCFailedReason pcReason = PC_NOERROR;
  ierr = KSPGetPC(solver.getHandle(), &pc);
  if (ierr == PETSC_SUCCESS && pc)
    ierr = PCGetFailedReason(pc, &pcReason);

  Alert::Info()
    << "   | " << label
    << " KSP reason=" << static_cast<int>(reason)
    << ", its=" << iterations
    << ", residual=" << residual
    << ", pcReason=" << static_cast<int>(pcReason)
    << Alert::Raise;

  if (reason <= 0)
  {
    throw std::runtime_error(
      label + " KSP did not converge: reason=" + std::to_string(static_cast<int>(reason)) +
      ", its=" + std::to_string(static_cast<long long>(iterations)) +
      ", residual=" + std::to_string(static_cast<double>(residual)) +
      ", pcReason=" + std::to_string(static_cast<int>(pcReason)));
  }
}

template <class GF>
static bool isFiniteGridFunction(const GF& gf)
{
  const auto& data = gf.getData();
  for (auto i = decltype(data.size())(0); i < data.size(); ++i)
  {
    if (!std::isfinite(static_cast<double>(data[i])))
      return false;
  }
  return true;
}

static Math::SpatialPoint zeroPoint(size_t dim)
{
  Math::SpatialPoint p(dim);
  p.setZero();
  return p;
}

static void translateMesh(MeshType& mesh, const std::array<Real, 3>& shift)
{
  const Real shiftNorm =
    std::sqrt(shift[0] * shift[0] + shift[1] * shift[1] + shift[2] * shift[2]);
  if (shiftNorm == 0.0)
    return;

  for (auto it = mesh.getVertex(); it; ++it)
  {
    auto coords = mesh.getVertexCoordinates(it->getIndex());
    coords(0) += shift[0];
    if (coords.size() > 1)
      coords(1) += shift[1];
    if (coords.size() > 2)
      coords(2) += shift[2];
    mesh.setVertexCoordinates(it->getIndex(), coords);
  }
}

static Math::SpatialPoint computeRegionCentroid(const MeshType& mesh, Attribute region)
{
  Math::SpatialPoint centroid = zeroPoint(mesh.getSpaceDimension());
  Real totalMeasure = 0.0;

  for (auto it = mesh.getCell(); it; ++it)
  {
    if (it->getAttribute() != region)
      continue;

    const Real measure = it->getMeasure();
    if (!(measure > 0.0))
      continue;

    const auto rc = Geometry::Polytope::Traits(it->getGeometry()).getCentroid();
    const Geometry::Point p(*it, rc);
    centroid(0) += measure * p.x();
    if (centroid.size() > 1)
      centroid(1) += measure * p.y();
    if (centroid.size() > 2)
      centroid(2) += measure * p.z();
    totalMeasure += measure;
  }

  if (!(totalMeasure > 0.0))
    throw std::runtime_error("Failed to compute obstacle centroid: obstacle measure is zero.");

  centroid /= totalMeasure;
  return centroid;
}

static MeshCheckpointStats computeMeshCheckpointStats(
  MeshType& mesh,
  const MeshQualityConfig& qualityConfig)
{
  MeshCheckpointStats stats;
  stats.quality = evaluateMeshQuality(mesh, qualityConfig);
  stats.vertexCount = stats.quality.vertexCount;
  stats.cellCount = stats.quality.cellCount;
  stats.obstacleVolume = mesh.getVolume(Obstacle);
  stats.interfaceArea = computeInterfaceArea(mesh);
  stats.centroid = toVec3(computeRegionCentroid(mesh, Obstacle));

  stats.bboxMin = {
    std::numeric_limits<Real>::infinity(),
    std::numeric_limits<Real>::infinity(),
    std::numeric_limits<Real>::infinity()
  };
  stats.bboxMax = {
    -std::numeric_limits<Real>::infinity(),
    -std::numeric_limits<Real>::infinity(),
    -std::numeric_limits<Real>::infinity()
  };

  bool hasVertex = false;
  for (auto vit = mesh.getVertex(); vit; ++vit)
  {
    hasVertex = true;
    const auto p = toVec3(mesh.getVertexCoordinates(vit->getIndex()));
    stats.bboxMin.x = std::min(stats.bboxMin.x, p.x);
    stats.bboxMin.y = std::min(stats.bboxMin.y, p.y);
    stats.bboxMin.z = std::min(stats.bboxMin.z, p.z);
    stats.bboxMax.x = std::max(stats.bboxMax.x, p.x);
    stats.bboxMax.y = std::max(stats.bboxMax.y, p.y);
    stats.bboxMax.z = std::max(stats.bboxMax.z, p.z);
  }
  if (!hasVertex)
  {
    stats.bboxMin = {};
    stats.bboxMax = {};
  }

  for (auto fit = mesh.getFace(); fit; ++fit)
    ++stats.faceCount;

  return stats;
}

static void logMeshCheckpointStats(
  const MeshCheckpointStats& stats,
  const std::string& label)
{
  Alert::Info()
    << "   | Mesh checkpoint " << label
    << ": vertices=" << stats.vertexCount
    << ", cells=" << stats.cellCount
    << ", faces=" << stats.faceCount
    << ", Vobs=" << stats.obstacleVolume
    << ", A=" << stats.interfaceArea
    << ", centroid=(" << stats.centroid.x << ", " << stats.centroid.y << ", " << stats.centroid.z << ")"
    << ", bbox_min=(" << stats.bboxMin.x << ", " << stats.bboxMin.y << ", " << stats.bboxMin.z << ")"
    << ", bbox_max=(" << stats.bboxMax.x << ", " << stats.bboxMax.y << ", " << stats.bboxMax.z << ")"
    << ", min_tet_volume=" << stats.quality.minTetVolume
    << ", min_altitude_ratio=" << stats.quality.minAltitudeRatio
    << Alert::Raise;
}

static void logMeshCheckpointCompare(
  const MeshCheckpointStats& previousEval,
  const MeshCheckpointStats& currentInput,
  size_t previousIteration,
  size_t currentIteration)
{
  const Vec3 centroidDiff{
    currentInput.centroid.x - previousEval.centroid.x,
    currentInput.centroid.y - previousEval.centroid.y,
    currentInput.centroid.z - previousEval.centroid.z
  };
  const Vec3 bboxMinDiff{
    currentInput.bboxMin.x - previousEval.bboxMin.x,
    currentInput.bboxMin.y - previousEval.bboxMin.y,
    currentInput.bboxMin.z - previousEval.bboxMin.z
  };
  const Vec3 bboxMaxDiff{
    currentInput.bboxMax.x - previousEval.bboxMax.x,
    currentInput.bboxMax.y - previousEval.bboxMax.y,
    currentInput.bboxMax.z - previousEval.bboxMax.z
  };

  Alert::Info()
    << "   | Checkpoint compare: candidate_eval_optimized iteration " << previousIteration
    << " vs iteration_input_optimized_from_raw iteration " << currentIteration
    << ": d_vertices=" << static_cast<long long>(currentInput.vertexCount) - static_cast<long long>(previousEval.vertexCount)
    << ", d_cells=" << static_cast<long long>(currentInput.cellCount) - static_cast<long long>(previousEval.cellCount)
    << ", d_faces=" << static_cast<long long>(currentInput.faceCount) - static_cast<long long>(previousEval.faceCount)
    << ", dV=" << currentInput.obstacleVolume - previousEval.obstacleVolume
    << ", dA=" << currentInput.interfaceArea - previousEval.interfaceArea
    << ", centroid_diff_norm=" << norm(centroidDiff)
    << ", bbox_min_diff_norm=" << norm(bboxMinDiff)
    << ", bbox_max_diff_norm=" << norm(bboxMaxDiff)
    << ", d_min_tet_volume=" << currentInput.quality.minTetVolume - previousEval.quality.minTetVolume
    << ", d_min_altitude_ratio=" << currentInput.quality.minAltitudeRatio - previousEval.quality.minAltitudeRatio
    << Alert::Raise;
}

static auto makeRigidVelocity(
  const std::array<Real, 3>& translation,
  const std::array<Real, 3>& omega,
  const Math::SpatialPoint& center)
{
  const Real cx = center.size() > 0 ? center(0) : 0.0;
  const Real cy = center.size() > 1 ? center(1) : 0.0;
  const Real cz = center.size() > 2 ? center(2) : 0.0;

  return VectorFunction{
    [=](const Geometry::Point& p)
    {
      const Real x = p.x() - cx;
      const Real y = p.y() - cy;
      const Real z = p.z() - cz;
      return translation[0] + omega[1] * z - omega[2] * y;
    },
    [=](const Geometry::Point& p)
    {
      const Real x = p.x() - cx;
      const Real y = p.y() - cy;
      const Real z = p.z() - cz;
      return translation[1] + omega[2] * x - omega[0] * z;
    },
    [=](const Geometry::Point& p)
    {
      const Real x = p.x() - cx;
      const Real y = p.y() - cy;
      const Real z = p.z() - cz;
      return translation[2] + omega[0] * y - omega[1] * x;
    }
  };
}

static std::array<Real, 3> basis(int axis)
{
  if (axis < 0 || axis > 2)
    throw std::runtime_error("Axis must be 0, 1 or 2.");

  std::array<Real, 3> e{0.0, 0.0, 0.0};
  e[static_cast<size_t>(axis)] = 1.0;
  return e;
}

enum class ObjectiveMode
{
  K,
  C,
  Q
};

enum class ObjectiveSense
{
  Min,
  Max
};

enum class AreaDebugMode
{
  Off,
  Descent,
  Ascent
};

static ObjectiveMode parseObjectiveMode(const char* value)
{
  const std::string mode = value ? value : "K";
  if (mode == "K" || mode == "k")
    return ObjectiveMode::K;
  if (mode == "C" || mode == "c")
    return ObjectiveMode::C;
  if (mode == "Q" || mode == "q")
    return ObjectiveMode::Q;

  throw std::runtime_error(
    std::string("Unsupported ") + algorithmName + " objective mode. Use OBJECTIVE_MODE=K, OBJECTIVE_MODE=C or OBJECTIVE_MODE=Q.");
}

static ObjectiveSense parseObjectiveSense(const char* value)
{
  const std::string sense = value ? value : "min";
  if (sense == "min" || sense == "MIN" || sense == "Min")
    return ObjectiveSense::Min;
  if (sense == "max" || sense == "MAX" || sense == "Max")
    return ObjectiveSense::Max;

  throw std::runtime_error(
    "Unsupported objective sense. Use OBJECTIVE_SENSE=min or OBJECTIVE_SENSE=max.");
}

static AreaDebugMode parseAreaDebugMode(const char* value)
{
  const std::string mode = value ? value : "off";
  if (mode == "off" || mode == "OFF" || mode == "Off")
    return AreaDebugMode::Off;
  if (mode == "descent" || mode == "DESCENT" || mode == "Descent")
    return AreaDebugMode::Descent;
  if (mode == "ascent" || mode == "ASCENT" || mode == "Ascent")
    return AreaDebugMode::Ascent;

  throw std::runtime_error(
    "Unsupported AREA_DEBUG_MODE. Use off, descent or ascent.");
}

static void setPetscOption(const std::string& name, const std::string& value)
{
  PetscOptionsSetValue(nullptr, name.c_str(), value.c_str());
}

static void configureStokesFieldSplitOptions(
  const std::string& prefix,
  const std::string& velocitySplitName,
  const std::string& pressureSplitName,
  const std::string& kspType,
  const std::string& pcType,
  const std::string& fieldSplitType,
  const std::string& schurFactType,
  const std::string& schurPrecondition,
  const std::string& velocityPcType,
  const std::string& pressurePcType,
  Real pcFactorShiftAmount)
{
  const std::string p = "-" + prefix;
  setPetscOption(p + "ksp_type", kspType);
  setPetscOption(p + "pc_type", pcType);

  if (pcType == "fieldsplit")
  {
    setPetscOption(p + "pc_fieldsplit_type", fieldSplitType);
    if (fieldSplitType == "schur")
    {
      setPetscOption(p + "pc_fieldsplit_schur_fact_type", schurFactType);
      setPetscOption(p + "pc_fieldsplit_schur_precondition", schurPrecondition);
    }

    setPetscOption(p + "fieldsplit_" + velocitySplitName + "_ksp_type", "preonly");
    setPetscOption(p + "fieldsplit_" + velocitySplitName + "_pc_type", velocityPcType);
    setPetscOption(p + "fieldsplit_" + pressureSplitName + "_ksp_type", "preonly");
    setPetscOption(p + "fieldsplit_" + pressureSplitName + "_pc_type", pressurePcType);
  }

  if (pcType != "bjacobi" && pcType != "jacobi" && pcType != "none" && pcType != "fieldsplit")
  {
    const std::string shift = std::to_string(pcFactorShiftAmount);
    setPetscOption(p + "pc_factor_shift_type", "NONZERO");
    setPetscOption(p + "pc_factor_shift_amount", shift);
  }
}

template <class ProblemType>
static void setStokesVelocityPressureSplits(
  ProblemType& problem,
  const std::string& velocityName,
  const std::string& pressureName)
{
  const auto& offsets = problem.getTrialOffsets();
  if (offsets.size() < 2)
    throw std::runtime_error("Stokes field split requires at least velocity and pressure fields.");

  auto& axb = problem.getLinearSystem();
  using Split = Rodin::PETSc::Math::LinearSystem::FieldSplits::Split;

  const PetscInt velocityStart = static_cast<PetscInt>(offsets[0]);
  const PetscInt velocitySize = static_cast<PetscInt>(offsets[1] - offsets[0]);
  const PetscInt pressureStart = static_cast<PetscInt>(offsets[1]);
  const PetscInt pressureSize =
    static_cast<PetscInt>(problem.getTotalTrialSize() - offsets[1]);

  ::IS velocityIS = PETSC_NULLPTR;
  ::IS pressureIS = PETSC_NULLPTR;
  PetscErrorCode ierr = ISCreateStride(
    axb.getCommunicator(), velocitySize, velocityStart, 1, &velocityIS);
  if (ierr != PETSC_SUCCESS)
    throw std::runtime_error("Failed to create velocity field split.");

  ierr = ISCreateStride(
    axb.getCommunicator(), pressureSize, pressureStart, 1, &pressureIS);
  if (ierr != PETSC_SUCCESS)
  {
    ISDestroy(&velocityIS);
    throw std::runtime_error("Failed to create pressure field split.");
  }

  std::vector<Split> splits;
  splits.push_back(Split{velocityName, velocityIS});
  // Pressure and the global mean-pressure Lagrange multiplier are contiguous.
  // Grouping them gives PETSc the two-field saddle-point structure expected
  // by Schur field-split preconditioning.
  splits.push_back(Split{pressureName, pressureIS});
  axb.setFieldSplits(Rodin::PETSc::Math::LinearSystem::FieldSplits{std::move(splits)});
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Force line-buffered stdout so output appears promptly even under ld-linux loader
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  std::cout << std::unitbuf;

  const Real defaultHMaxRuntime = LevelSetStokes::Runtime::envDouble("HMAX", defaultHMax);
  const Real hminRatio = LevelSetStokes::Runtime::envDouble("HMIN_RATIO", 0.1);
  const Real hausdRatio = LevelSetStokes::Runtime::envDouble("HAUSD_RATIO", 0.1);
  const Real retryHMinRatio = hminRatio > 0.33 ? hminRatio : 0.33;

  const Real hmax0 = defaultHMaxRuntime;

  Real hmax  = hmax0;

  try
  {
    const char* meshFile = LevelSetStokes::Runtime::envCString("MESH", LEVELSET_STOKES_MESH_FILE);
    const size_t maxIt   = LevelSetStokes::Runtime::envSizeT("MAX_ITERS", defaultMaxIt);
    const size_t convergenceWindow =
      LevelSetStokes::Runtime::envSizeT("CONVERGENCE_WINDOW", defaultConvergenceWindow);
    const Real convergenceRtolJraw =
      LevelSetStokes::Runtime::envDouble("CONVERGENCE_RTOL_JRAW", defaultConvergenceRtolJraw);
    const Real nsAlphaJParam =
      LevelSetStokes::Runtime::envDouble("NS_ALPHA_J", defaultNsAlphaJ);
    const Real nsAlphaCParam =
      LevelSetStokes::Runtime::envDouble("NS_ALPHA_C", defaultNsAlphaC);
    const Real surfaceAreaFactor =
      LevelSetStokes::Runtime::envDouble("SURFACE_AREA_FACTOR", defaultSurfaceAreaFactor);
    const AreaDebugMode areaDebugMode =
      parseAreaDebugMode(LevelSetStokes::Runtime::envCString("AREA_DEBUG_MODE", "off"));
    const Real areaCorrectionGain =
      LevelSetStokes::Runtime::envDouble("AREA_CORRECTION_GAIN", defaultAreaCorrectionGain);
    const Real pcFactorShiftAmount =
      LevelSetStokes::Runtime::envDouble("KSP_PC_SHIFT_AMOUNT", defaultPcFactorShiftAmount);
    const std::string kspType =
      LevelSetStokes::Runtime::envCString("KSP_TYPE", defaultKspType);
    const std::string kspPcType =
      LevelSetStokes::Runtime::envCString("KSP_PC_TYPE", defaultKspPcType);
    const std::string kspPcFactorSolverType =
      LevelSetStokes::Runtime::envCString("KSP_PC_FACTOR_SOLVER_TYPE", defaultKspPcFactorSolverType);
    const size_t kspMaxIt =
      LevelSetStokes::Runtime::envSizeT("KSP_MAX_ITERS", defaultKspMaxIt);
    const Real kspRtol =
      LevelSetStokes::Runtime::envDouble("KSP_RTOL", defaultKspRtol);
    const Real stokesRtol =
      LevelSetStokes::Runtime::envDouble("STOKES_RTOL", defaultStokesRtol);
    const size_t stokesMaxIt =
      LevelSetStokes::Runtime::envSizeT("STOKES_MAX_ITERS", defaultStokesMaxIt);
    const std::string stokesKspType =
      LevelSetStokes::Runtime::envCString("STOKES_KSP_TYPE", defaultStokesKspType);
    const std::string stokesPcType =
      LevelSetStokes::Runtime::envCString("STOKES_PC_TYPE", defaultStokesPcType);
    const std::string stokesFieldSplitType =
      LevelSetStokes::Runtime::envCString("STOKES_FIELDSPLIT_TYPE", defaultStokesFieldSplitType);
    const std::string stokesSchurFactType =
      LevelSetStokes::Runtime::envCString("STOKES_SCHUR_FACT_TYPE", defaultStokesSchurFactType);
    const std::string stokesSchurPrecondition =
      LevelSetStokes::Runtime::envCString("STOKES_SCHUR_PRECONDITION", defaultStokesSchurPrecondition);
    const std::string stokesVelocityPcType =
      LevelSetStokes::Runtime::envCString("STOKES_VELOCITY_PC_TYPE", defaultStokesVelocityPcType);
    const std::string stokesPressurePcType =
      LevelSetStokes::Runtime::envCString("STOKES_PRESSURE_PC_TYPE", defaultStokesPressurePcType);
    const size_t remeshQualityMaxRetries =
      LevelSetStokes::Runtime::envSizeT("REMESH_QUALITY_MAX_RETRIES", defaultRemeshQualityMaxRetries);
    const Real jLineSearchRtol =
      LevelSetStokes::Runtime::envDouble("J_LINESEARCH_RTOL", defaultJLineSearchRtol);
    const bool adaptSizeMapEnabled =
      LevelSetStokes::Runtime::envDouble("ADAPT_SIZE_MAP", 1.0) != 0.0;
    const Real adaptHNearParam =
      LevelSetStokes::Runtime::envDouble("ADAPT_H_NEAR", -1.0);
    const Real adaptHFarParam =
      LevelSetStokes::Runtime::envDouble("ADAPT_H_FAR", -1.0);
    const Real adaptRNearParam =
      LevelSetStokes::Runtime::envDouble("ADAPT_R_NEAR", -1.0);
    const Real adaptRFarParam =
      LevelSetStokes::Runtime::envDouble("ADAPT_R_FAR", -1.0);
    const Real adaptGradation =
      LevelSetStokes::Runtime::envDouble("ADAPT_GRADATION", defaultAdaptGradation);
    const MeshQualityConfig meshQualityConfig{
      LevelSetStokes::Runtime::envDouble("QUALITY_MIN_TET_VOLUME", defaultQualityMinTetVolume),
      LevelSetStokes::Runtime::envDouble("QUALITY_MIN_ALTITUDE_RATIO", defaultQualityMinAltitudeRatio),
      LevelSetStokes::Runtime::envDouble("QUALITY_DUPLICATE_TOL", defaultQualityDuplicateTol)
    };
    PetscOptionsSetValue(nullptr, "-ksp_type", kspType.c_str());
    PetscOptionsSetValue(nullptr, "-pc_type", kspPcType.c_str());
    if (!kspPcFactorSolverType.empty())
      PetscOptionsSetValue(nullptr, "-pc_factor_mat_solver_type", kspPcFactorSolverType.c_str());
    if (kspPcType != "bjacobi" && kspPcType != "jacobi" && kspPcType != "none")
    {
      const std::string pcFactorShiftAmountValue = std::to_string(pcFactorShiftAmount);
      PetscOptionsSetValue(nullptr, "-pc_factor_shift_type", "NONZERO");
      PetscOptionsSetValue(nullptr, "-pc_factor_shift_amount", pcFactorShiftAmountValue.c_str());
    }
    configureStokesFieldSplitOptions(
      "state_",
      "u",
      "p",
      stokesKspType,
      stokesPcType,
      stokesFieldSplitType,
      stokesSchurFactType,
      stokesSchurPrecondition,
      stokesVelocityPcType,
      stokesPressurePcType,
      pcFactorShiftAmount);
    configureStokesFieldSplitOptions(
      "adjoint_",
      "v",
      "q",
      stokesKspType,
      stokesPcType,
      stokesFieldSplitType,
      stokesSchurFactType,
      stokesSchurPrecondition,
      stokesVelocityPcType,
      stokesPressurePcType,
      pcFactorShiftAmount);
    const ObjectiveMode objectiveMode =
      parseObjectiveMode(LevelSetStokes::Runtime::envCString("OBJECTIVE_MODE", "K"));
    const ObjectiveSense objectiveSense =
      parseObjectiveSense(LevelSetStokes::Runtime::envCString("OBJECTIVE_SENSE", "min"));
    const std::array<Real, 3> shift{
      LevelSetStokes::Runtime::envDouble("SHIFT_X", 0.0),
      LevelSetStokes::Runtime::envDouble("SHIFT_Y", 0.0),
      LevelSetStokes::Runtime::envDouble("SHIFT_Z", 0.0)
    };

    // Objective component K_ij, default K_11 in 0-based storage -> axis 1,1 as in original file.
    const int iAxis = LevelSetStokes::Runtime::envInt("IAXIS", 1);
    const int jAxis = LevelSetStokes::Runtime::envInt("JAXIS", 1);

    Real hmin  = hminRatio * hmax;
    Real hausd = hausdRatio * hmin;

    const Real stepK = LevelSetStokes::Runtime::envDouble("STEP_K", 0.1);
    const Real alpha = 16 * (hmax0 - hmin);

    // Load and pre-optimize initial mesh
    MeshType th;
    th.load(meshFile, IO::FileFormat::MEDIT);
    clearSamePhaseInternalFaces(th, "after initial load", true);

    MMG::Optimizer()
      .setHMax(hmax)
      .setHMin(hmin)
      .setHausdorff(hausd)
      .optimize(th);
    clearSamePhaseInternalFaces(th, "after initial optimizer", true);

    translateMesh(th, shift);
    saveArtifactMesh(th, "Omega0.mesh", "Omega0.mesh");

    // XDMF output
    std::filesystem::create_directories("out");
    const std::filesystem::path checkpointDir = std::filesystem::path("out") / "checkpoints";
    std::filesystem::create_directories(checkpointDir);
    IO::XDMF xdmf(std::string("out/") + algorithmName);

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fObjRaw("obj_raw.txt");
    std::ofstream fObjBefore("obj_before.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fNS("ns.txt");
    std::ofstream fCentroid("centroid.txt");
    std::ofstream fLineSearch("line_search.txt");

    if (!fObj || !fObjRaw || !fObjBefore || !fVol || !fNS || !fCentroid || !fLineSearch)
      throw std::runtime_error("Failed to open output history files.");

    // Constraint targets measured on the pre-optimized initial mesh.
    const double targetObstacleVolume = th.getVolume(Obstacle);
    const double initialShapeArea = computeInterfaceArea(th);
    const double targetSurfaceArea = surfaceAreaFactor * initialShapeArea;
    std::vector<double> jrawHistory;
    jrawHistory.reserve(maxIt);
    bool hasPreviousEvalCheckpoint = false;
    size_t previousEvalCheckpointIteration = 0;
    MeshCheckpointStats previousEvalCheckpointStats;
    bool hasLastAcceptedObjective = false;
    double lastAcceptedJraw = 0.0;
    double lastAcceptedJ = 0.0;
    const auto ei = basis(iAxis);
    const auto ej = basis(jAxis);
    const bool rotationalState = objectiveMode == ObjectiveMode::C || objectiveMode == ObjectiveMode::Q;
    const bool rotationalAdjoint = objectiveMode == ObjectiveMode::Q;
    const std::array<Real, 3> stateTranslation =
      objectiveMode == ObjectiveMode::K ? ej : std::array<Real, 3>{0.0, 0.0, 0.0};
    const std::array<Real, 3> stateOmega =
      rotationalState ? ej : std::array<Real, 3>{0.0, 0.0, 0.0};
    const std::array<Real, 3> adjointTranslation =
      rotationalAdjoint ? std::array<Real, 3>{0.0, 0.0, 0.0} : ei;
    const std::array<Real, 3> adjointOmega =
      rotationalAdjoint ? ei : std::array<Real, 3>{0.0, 0.0, 0.0};

    // Initial connectivity
    {
      auto& conn = th.getConnectivity();
      conn.compute(2, 3);
      conn.compute(3, 2);
      conn.compute(2, 1);
      conn.compute(1, 0);
      conn.compute(0, 0);
    }

    // Detect the fixed exterior boundary attributes once from the initial mesh
    FlatSet<Attribute> baseOuterBdr;
    for (auto it = th.getFace(); !it.end(); ++it)
    {
      const auto& face = *it;
      if (!face.isBoundary())
        continue;

      const auto a = face.getAttribute();
      if (!a || *a == GammaShape)
        continue;

      baseOuterBdr.insert(*a);
    }

    struct ObjectiveSample
    {
      double Jraw = std::numeric_limits<double>::quiet_NaN();
      double J = std::numeric_limits<double>::quiet_NaN();
      double obstacleVolume = std::numeric_limits<double>::quiet_NaN();
      double shapeArea = std::numeric_limits<double>::quiet_NaN();
      Math::SpatialPoint centroid;
    };

    auto evaluateObjectiveOnMesh =
      [&](MeshType& mesh, const std::string& label) -> ObjectiveSample
      {
        {
          auto& conn = mesh.getConnectivity();
          conn.compute(2, 3);
          conn.compute(3, 2);
          conn.compute(2, 1);
          conn.compute(1, 0);
          conn.compute(0, 0);
        }

        const auto sampleCentroid = computeRegionCentroid(mesh, Obstacle);
        const auto sampleZeroCenter = zeroPoint(mesh.getSpaceDimension());
        const auto sampleStateVelocity = makeRigidVelocity(
          stateTranslation,
          stateOmega,
          rotationalState ? sampleCentroid : sampleZeroCenter);
        const auto sampleAdjointVelocity = makeRigidVelocity(
          adjointTranslation,
          adjointOmega,
          rotationalAdjoint ? sampleCentroid : sampleZeroCenter);

        auto sampleFluidMesh = mesh.trim(Obstacle);
        const size_t sampleD = mesh.getSpaceDimension();

        auto sampleVelocitySpace = H1(std::integral_constant<size_t, 2>{}, sampleFluidMesh, sampleD);
        auto samplePressureSpace = H1(std::integral_constant<size_t, 1>{}, sampleFluidMesh);
        P0g sampleGlobalP0(sampleFluidMesh);

        FlatSet<Attribute> sampleFluidShapeBdr{GammaShape};
        FlatSet<Attribute> sampleShapeInterface{GammaShape};

        PETSc::Variational::TrialFunction sampleUp(sampleVelocitySpace); sampleUp.setName("u");
        PETSc::Variational::TrialFunction samplePp(samplePressureSpace); samplePp.setName("p");
        PETSc::Variational::TrialFunction sampleLmb(sampleGlobalP0);     sampleLmb.setName("lambda");

        PETSc::Variational::TestFunction sampleVp(sampleVelocitySpace);
        PETSc::Variational::TestFunction sampleQp(samplePressureSpace);
        PETSc::Variational::TestFunction sampleMu0(sampleGlobalP0);

        Problem sampleState(sampleUp, samplePp, sampleLmb, sampleVp, sampleQp, sampleMu0);
        sampleState =
            Integral(mu * Jacobian(sampleUp), Jacobian(sampleVp))
          - Integral(samplePp, Div(sampleVp))
          + Integral(Div(sampleUp), sampleQp)
          + Integral(sampleLmb, sampleQp)
          + Integral(samplePp, sampleMu0)
          + regularization * Integral(samplePp, sampleQp)
          + regularization * Integral(sampleLmb, sampleMu0)
          + DirichletBC(sampleUp, sampleStateVelocity).on(sampleFluidShapeBdr)
          + DirichletBC(sampleUp, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

        auto sampleU = GridFunction(sampleVelocitySpace);
        auto sampleP = GridFunction(samplePressureSpace);
        sampleState.assemble();
        bool sampleStateFallbackNeeded = false;
        try
        {
          Solver::KSP(sampleState).solve();
        }
        catch (const std::exception& e)
        {
          sampleStateFallbackNeeded = true;
          Alert::Warning()
            << "   | Fast line-search state " << label
            << " failed: " << e.what()
            << ". Retrying with Schur+LU fallback."
            << Alert::Raise;
        }

        if (!sampleStateFallbackNeeded)
        {
          sampleU = sampleUp.getSolution();
          sampleP = samplePp.getSolution();
          sampleStateFallbackNeeded =
            !isFiniteGridFunction(sampleU) || !isFiniteGridFunction(sampleP);
        }

        if (sampleStateFallbackNeeded)
        {
          Alert::Warning()
            << "   | Fast line-search state " << label
            << " produced non-finite values; retrying with Schur+LU fallback."
            << Alert::Raise;

          PETSc::Variational::TrialFunction sampleUpFallback(sampleVelocitySpace);
          sampleUpFallback.setName("u");
          PETSc::Variational::TrialFunction samplePpFallback(samplePressureSpace);
          samplePpFallback.setName("p");
          PETSc::Variational::TrialFunction sampleLmbFallback(sampleGlobalP0);
          sampleLmbFallback.setName("lambda");

          PETSc::Variational::TestFunction sampleVpFallback(sampleVelocitySpace);
          PETSc::Variational::TestFunction sampleQpFallback(samplePressureSpace);
          PETSc::Variational::TestFunction sampleMu0Fallback(sampleGlobalP0);

          Problem sampleStateFallback(
            sampleUpFallback,
            samplePpFallback,
            sampleLmbFallback,
            sampleVpFallback,
            sampleQpFallback,
            sampleMu0Fallback);
          sampleStateFallback =
              Integral(mu * Jacobian(sampleUpFallback), Jacobian(sampleVpFallback))
            - Integral(samplePpFallback, Div(sampleVpFallback))
            + Integral(Div(sampleUpFallback), sampleQpFallback)
            + Integral(sampleLmbFallback, sampleQpFallback)
            + Integral(samplePpFallback, sampleMu0Fallback)
            + regularization * Integral(samplePpFallback, sampleQpFallback)
            + regularization * Integral(sampleLmbFallback, sampleMu0Fallback)
            + DirichletBC(sampleUpFallback, sampleStateVelocity).on(sampleFluidShapeBdr)
            + DirichletBC(sampleUpFallback, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

          sampleStateFallback.assemble();
          setStokesVelocityPressureSplits(sampleStateFallback, "u", "p");
          solveKSPChecked(
            sampleStateFallback,
            "line-search state fallback " + label,
            stokesRtol,
            static_cast<PetscInt>(stokesMaxIt),
            std::string("state_"));

          sampleU = sampleUpFallback.getSolution();
          sampleP = samplePpFallback.getSolution();
          if (!isFiniteGridFunction(sampleU) || !isFiniteGridFunction(sampleP))
            throw std::runtime_error("line-search state fallback " + label + " produced non-finite values.");
        }

        auto sampleJu = Jacobian(sampleU);
        auto sampleNObj = -FaceNormal(sampleFluidMesh);
        auto sampleSigma = mu * (sampleJu + sampleJu.T()) - sampleP * IdentityMatrix(sampleD);
        auto sampleTraction = sampleSigma * sampleNObj;

        P0g sampleP0Obj(sampleFluidMesh);
        TestFunction sampleZ0Obj(sampleP0Obj);
        LinearForm sampleLfObj(sampleZ0Obj);
        sampleLfObj =
          FaceIntegral(Dot(sampleTraction, sampleAdjointVelocity), sampleZ0Obj).over(sampleShapeInterface);
        sampleLfObj.assemble();

        GridFunction sampleOneObj(sampleP0Obj);
        sampleOneObj = 1.0;

        const double sampleJraw = -sampleLfObj(sampleOneObj);
        const double sampleSenseSign = objectiveSense == ObjectiveSense::Max ? -1.0 : 1.0;
        const double sampleJ = sampleSenseSign * sampleJraw;
        if (!std::isfinite(sampleJraw) || !std::isfinite(sampleJ))
          throw std::runtime_error("line-search objective " + label + " is non-finite.");

        return ObjectiveSample{
          sampleJraw,
          sampleJ,
          mesh.getVolume(Obstacle),
          computeInterfaceArea(mesh),
          sampleCentroid
        };
      };

    Alert::Info()
      << algorithmName << " configuration:"
      << " mode=" << (objectiveMode == ObjectiveMode::K ? "K" : (objectiveMode == ObjectiveMode::C ? "C" : "Q"))
      << ", sense=" << (objectiveSense == ObjectiveSense::Min ? "min" : "max")
      << ", iAxis=" << iAxis
      << ", jAxis=" << jAxis
      << ", shift=(" << shift[0] << ", " << shift[1] << ", " << shift[2] << ")"
      << ", surfaceAreaFactor=" << surfaceAreaFactor
      << ", areaCorrectionGain=" << areaCorrectionGain
      << ", kspType=" << kspType
      << ", kspPcType=" << kspPcType
      << ", kspPcFactorSolverType=" << kspPcFactorSolverType
      << ", stokesKspType=" << stokesKspType
      << ", stokesPcType=" << stokesPcType
      << ", stokesFieldSplitType=" << stokesFieldSplitType
      << ", stokesSchurFactType=" << stokesSchurFactType
      << ", stokesSchurPrecondition=" << stokesSchurPrecondition
      << ", stokesVelocityPcType=" << stokesVelocityPcType
      << ", stokesPressurePcType=" << stokesPressurePcType
      << ", stokesRtol=" << stokesRtol
      << ", stokesMaxIt=" << stokesMaxIt
      << ", pcFactorShiftAmount=" << pcFactorShiftAmount
      << ", kspMaxIt=" << kspMaxIt
      << ", kspRtol=" << kspRtol
      << ", remeshQualityMaxRetries=" << remeshQualityMaxRetries
      << ", jLineSearchRtol=" << jLineSearchRtol
      << ", adaptSizeMap=" << (adaptSizeMapEnabled ? 1 : 0)
      << ", adaptHNear=" << adaptHNearParam
      << ", adaptHFar=" << adaptHFarParam
      << ", adaptRNear=" << adaptRNearParam
      << ", adaptRFar=" << adaptRFarParam
      << ", adaptGradation=" << adaptGradation
      << ", qualityMinTetVolume=" << meshQualityConfig.minTetVolume
      << ", qualityMinAltitudeRatio=" << meshQualityConfig.minAltitudeRatio
      << ", qualityDuplicateTol=" << meshQualityConfig.duplicateTol
      << Alert::Raise;

    for (size_t it = 0; it < maxIt; ++it)
    {
      Alert::Info() << "----- Iteration: " << it << Alert::Raise;

      // ----------------------------------------------------------------------
      // Optimize and validate current mesh
      // ----------------------------------------------------------------------
      Alert::Info()
        << "   | Optimizing iteration input mesh: hmax=" << hmax
        << ", hmin=" << hmin
        << ", hausd=" << hausd
        << Alert::Raise;
      clearSamePhaseInternalFaces(th, "iteration input", true);
      try
      {
        MMG::Optimizer()
          .setHMax(hmax)
          .setHMin(hmin)
          .setHausdorff(hausd)
          .setAngleDetection(false)
          .optimize(th);
        clearSamePhaseInternalFaces(th, "after iteration input optimizer", true);
      }
      catch (const Alert::Exception&)
      {
        hmax /= 2;
        hmin = retryHMinRatio * hmax;
        hausd = hausdRatio * hmin;

        Alert::Warning()
          << "Iteration input optimizer failed at iteration " << it
          << ". Reducing hmax to " << hmax
          << " and retrying."
          << Alert::Raise;
        continue;
      }

      const auto currentQuality = evaluateMeshQuality(th, meshQualityConfig);
      logMeshQuality(currentQuality, "after iteration input optimizer", currentQuality.accepted());

      if (hasPreviousEvalCheckpoint)
      {
        const auto checkpointPath =
          checkpointDir / ("iteration_input_optimized_from_raw." + std::to_string(it) + ".mesh");
        saveArtifactMesh(
          th,
          checkpointPath,
          checkpointPath.string());
        const auto currentInputCheckpointStats =
          computeMeshCheckpointStats(th, meshQualityConfig);
        logMeshCheckpointStats(
          currentInputCheckpointStats,
          "iteration_input_optimized_from_raw." + std::to_string(it));
        logMeshCheckpointCompare(
          previousEvalCheckpointStats,
          currentInputCheckpointStats,
          previousEvalCheckpointIteration,
          it);
        hasPreviousEvalCheckpoint = false;
      }

      // ----------------------------------------------------------------------
      // Connectivity
      // ----------------------------------------------------------------------
      {
        auto& conn = th.getConnectivity();
        conn.compute(2, 3);
        conn.compute(3, 2);
        conn.compute(2, 1);
        conn.compute(1, 0);
        conn.compute(0, 0);
      }

      const auto obstacleCentroid = computeRegionCentroid(th, Obstacle);
      const auto zeroCenter = zeroPoint(th.getSpaceDimension());
      const auto stateVelocity = makeRigidVelocity(
        stateTranslation,
        stateOmega,
        rotationalState ? obstacleCentroid : zeroCenter);
      const auto adjointVelocity = makeRigidVelocity(
        adjointTranslation,
        adjointOmega,
        rotationalAdjoint ? obstacleCentroid : zeroCenter);

      Alert::Info()
        << "   | Obstacle centroid = ("
        << obstacleCentroid(0) << ", "
        << (obstacleCentroid.size() > 1 ? obstacleCentroid(1) : 0.0) << ", "
        << (obstacleCentroid.size() > 2 ? obstacleCentroid(2) : 0.0) << ")"
        << Alert::Raise;

      // ----------------------------------------------------------------------
      // Trim obstacle to get fluid mesh
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Trimming fluid mesh." << Alert::Raise;
      auto fluidMesh = th.trim(Obstacle);
      fluidMesh.save("Omega.mesh", IO::FileFormat::MEDIT);

      const size_t d = th.getSpaceDimension();

      // FE spaces
      auto velocitySpace = H1(std::integral_constant<size_t, 2>{}, fluidMesh, d);
      auto pressureSpace = H1(std::integral_constant<size_t, 1>{}, fluidMesh);
      P0g globalP0(fluidMesh);

      FlatSet<Attribute> fluidShapeBdr{GammaShape};
      FlatSet<Attribute> shapeInterface{GammaShape};

      // ----------------------------------------------------------------------
      // State Stokes solve
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Solving state Stokes." << Alert::Raise;

      PETSc::Variational::TrialFunction up(velocitySpace); up.setName("u");
      PETSc::Variational::TrialFunction pp(pressureSpace); pp.setName("p");
      PETSc::Variational::TrialFunction lmb(globalP0);     lmb.setName("lambda");

      PETSc::Variational::TestFunction vp(velocitySpace);
      PETSc::Variational::TestFunction qp(pressureSpace);
      PETSc::Variational::TestFunction mu0(globalP0);

      Problem state(up, pp, lmb, vp, qp, mu0);
      state =
          Integral(mu * Jacobian(up), Jacobian(vp))
        - Integral(pp, Div(vp))
        + Integral(Div(up), qp)
        + Integral(lmb, qp)
        + Integral(pp, mu0)
        + regularization * Integral(pp, qp)
        + regularization * Integral(lmb, mu0)
        + DirichletBC(up, stateVelocity).on(fluidShapeBdr)
        + DirichletBC(up, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Alert::Info() << "   | Assembling Stokes." << Alert::Raise;
      state.assemble();

      Alert::Info() << "   | Solving Stokes." << Alert::Raise;

      auto u = GridFunction(velocitySpace);
      auto p = GridFunction(pressureSpace);
      bool stateFallbackNeeded = false;
      try
      {
        Solver::KSP(state).solve();
      }
      catch (const std::exception& e)
      {
        stateFallbackNeeded = true;
        Alert::Warning()
          << "   | Fast state Stokes failed: " << e.what()
          << ". Retrying with Schur+LU fallback."
          << Alert::Raise;
      }

      if (!stateFallbackNeeded)
      {
        u = up.getSolution();
        p = pp.getSolution();
        stateFallbackNeeded = !isFiniteGridFunction(u) || !isFiniteGridFunction(p);
      }

      if (stateFallbackNeeded)
      {
        Alert::Warning()
          << "   | Fast state Stokes produced non-finite values; retrying with Schur+LU fallback."
          << Alert::Raise;

        PETSc::Variational::TrialFunction upFallback(velocitySpace); upFallback.setName("u");
        PETSc::Variational::TrialFunction ppFallback(pressureSpace); ppFallback.setName("p");
        PETSc::Variational::TrialFunction lmbFallback(globalP0);     lmbFallback.setName("lambda");

        PETSc::Variational::TestFunction vpFallback(velocitySpace);
        PETSc::Variational::TestFunction qpFallback(pressureSpace);
        PETSc::Variational::TestFunction mu0Fallback(globalP0);

        Problem stateFallback(upFallback, ppFallback, lmbFallback, vpFallback, qpFallback, mu0Fallback);
        stateFallback =
            Integral(mu * Jacobian(upFallback), Jacobian(vpFallback))
          - Integral(ppFallback, Div(vpFallback))
          + Integral(Div(upFallback), qpFallback)
          + Integral(lmbFallback, qpFallback)
          + Integral(ppFallback, mu0Fallback)
          + regularization * Integral(ppFallback, qpFallback)
          + regularization * Integral(lmbFallback, mu0Fallback)
          + DirichletBC(upFallback, stateVelocity).on(fluidShapeBdr)
          + DirichletBC(upFallback, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

        stateFallback.assemble();
        setStokesVelocityPressureSplits(stateFallback, "u", "p");
        solveKSPChecked(
          stateFallback,
          "state Stokes fallback",
          stokesRtol,
          static_cast<PetscInt>(stokesMaxIt),
          std::string("state_"));

        u = upFallback.getSolution();
        p = ppFallback.getSolution();
        if (!isFiniteGridFunction(u) || !isFiniteGridFunction(p))
          throw std::runtime_error("state Stokes fallback produced non-finite values.");
      }

      // ----------------------------------------------------------------------
      // Adjoint Stokes solve
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Solving adjoint Stokes." << Alert::Raise;

      PETSc::Variational::TrialFunction ua(velocitySpace); ua.setName("v");
      PETSc::Variational::TrialFunction pa(pressureSpace); pa.setName("q");
      PETSc::Variational::TrialFunction la(globalP0);      la.setName("lambda_adj");

      PETSc::Variational::TestFunction va(velocitySpace);
      PETSc::Variational::TestFunction qa(pressureSpace);
      PETSc::Variational::TestFunction ma(globalP0);

      Problem adj(ua, pa, la, va, qa, ma);
      adj =
          Integral(mu * Jacobian(ua), Jacobian(va))
        - Integral(pa, Div(va))
        + Integral(Div(ua), qa)
        + Integral(la, qa)
        + Integral(pa, ma)
        + regularization * Integral(pa, qa)
        + regularization * Integral(la, ma)
        + DirichletBC(ua, adjointVelocity).on(fluidShapeBdr)
        + DirichletBC(ua, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      adj.assemble();
      auto v = GridFunction(velocitySpace);
      auto q = GridFunction(pressureSpace);
      bool adjointFallbackNeeded = false;
      try
      {
        Solver::KSP(adj).solve();
      }
      catch (const std::exception& e)
      {
        adjointFallbackNeeded = true;
        Alert::Warning()
          << "   | Fast adjoint Stokes failed: " << e.what()
          << ". Retrying with Schur+LU fallback."
          << Alert::Raise;
      }

      if (!adjointFallbackNeeded)
      {
        v = ua.getSolution();
        q = pa.getSolution();
        adjointFallbackNeeded = !isFiniteGridFunction(v) || !isFiniteGridFunction(q);
      }

      if (adjointFallbackNeeded)
      {
        Alert::Warning()
          << "   | Fast adjoint Stokes produced non-finite values; retrying with Schur+LU fallback."
          << Alert::Raise;

        PETSc::Variational::TrialFunction uaFallback(velocitySpace); uaFallback.setName("v");
        PETSc::Variational::TrialFunction paFallback(pressureSpace); paFallback.setName("q");
        PETSc::Variational::TrialFunction laFallback(globalP0);      laFallback.setName("lambda_adj");

        PETSc::Variational::TestFunction vaFallback(velocitySpace);
        PETSc::Variational::TestFunction qaFallback(pressureSpace);
        PETSc::Variational::TestFunction maFallback(globalP0);

        Problem adjFallback(uaFallback, paFallback, laFallback, vaFallback, qaFallback, maFallback);
        adjFallback =
            Integral(mu * Jacobian(uaFallback), Jacobian(vaFallback))
          - Integral(paFallback, Div(vaFallback))
          + Integral(Div(uaFallback), qaFallback)
          + Integral(laFallback, qaFallback)
          + Integral(paFallback, maFallback)
          + regularization * Integral(paFallback, qaFallback)
          + regularization * Integral(laFallback, maFallback)
          + DirichletBC(uaFallback, adjointVelocity).on(fluidShapeBdr)
          + DirichletBC(uaFallback, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

        adjFallback.assemble();
        setStokesVelocityPressureSplits(adjFallback, "v", "q");
        solveKSPChecked(
          adjFallback,
          "adjoint Stokes fallback",
          stokesRtol,
          static_cast<PetscInt>(stokesMaxIt),
          std::string("adjoint_"));

        v = uaFallback.getSolution();
        q = paFallback.getSolution();
        if (!isFiniteGridFunction(v) || !isFiniteGridFunction(q))
          throw std::runtime_error("adjoint Stokes fallback produced non-finite values.");
      }

      // ----------------------------------------------------------------------
      // Objective evaluation
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Evaluating objective." << Alert::Raise;

      auto ju = Jacobian(u);

      auto nObj = -FaceNormal(fluidMesh);
      auto sigma = mu * (ju + ju.T()) - p * IdentityMatrix(d);
      auto traction = sigma * nObj;

      P0g p0Obj(fluidMesh);
      TestFunction z0Obj(p0Obj);
      LinearForm lfObj(z0Obj);
      lfObj = FaceIntegral(Dot(traction, adjointVelocity), z0Obj).over(shapeInterface);
      lfObj.assemble();

      GridFunction oneObj(p0Obj);
      oneObj = 1.0;

      const double Jraw = -lfObj(oneObj);
      const double senseSign = objectiveSense == ObjectiveSense::Max ? -1.0 : 1.0;
      const double J = senseSign * Jraw;

      Alert::Info()
        << "   | Objective raw: " << Jraw
        << ", objective opt: " << J
        << Alert::Raise;

      // ----------------------------------------------------------------------
      // Shape gradient
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Computing shape gradient." << Alert::Raise;

      auto juShape = Jacobian(u);
      auto jvShape = Jacobian(v);

      auto eu = 0.5 * (juShape + juShape.T());
      auto ev = 0.5 * (jvShape + jvShape.T());

      auto Graw = 2.0 * mu * Dot(eu, ev);
      auto G = senseSign * Graw;

      const double obstacleVolume = th.getVolume(Obstacle);
      const double volumeViolation = obstacleVolume - targetObstacleVolume;

      // ----------------------------------------------------------------------
      // Signed-distance reconstruction and curvature projection
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Distancing obstacle phase." << Alert::Raise;

      P1 sh(th);
      GridFunction dist(sh);
      Distance::Eikonal(dist)
        .setInterior(Fluid)
        .setInterface(shapeInterface)
        .solve()
        .sign();

      const Real distMax = dist.max();
      GridFunction negDist(sh);
      negDist = -dist;
      const Real distMin = -negDist.max();
      const Real maxAbsDist = std::max(std::abs(distMin), std::abs(distMax));

      Alert::Info() << "   | Computing curvature from signed distance." << Alert::Raise;

      P1 kappaSpace(th);
      PETSc::Variational::TrialFunction kappaT(kappaSpace); kappaT.setName("kappa");
      PETSc::Variational::TestFunction kappaW(kappaSpace);
      const Real ellKappa = 2.0 * hmin;

      Problem curvProj(kappaT, kappaW);
      curvProj =
          Integral(kappaT, kappaW)
        + ellKappa * ellKappa * Integral(Grad(kappaT), Grad(kappaW))
        + Integral(Grad(dist), Grad(kappaW));
      curvProj.assemble();
      solveKSPChecked(
        curvProj,
        "curvature projection",
        kspRtol,
        static_cast<PetscInt>(kspMaxIt));

      Alert::Info() << "   | Curvature checkpoint: extracting kappa solution." << Alert::Raise;
      GridFunction kappa(kappaSpace);
      kappa = kappaT.getSolution();

      Alert::Info() << "   | Curvature checkpoint: computing kappa extrema." << Alert::Raise;
      const Real kappaMax = kappa.max();
      GridFunction negKappa(kappaSpace);
      negKappa = -kappa;
      const Real kappaMin = -negKappa.max();
      const Real maxAbsKappa = std::max(std::abs(kappaMin), std::abs(kappaMax));
      if (!isFiniteGridFunction(kappa))
        throw std::runtime_error("curvature projection produced non-finite kappa values.");

      Alert::Info() << "   | Curvature checkpoint: assembling objective-gradient mean." << Alert::Raise;
      P0g p0ObjectiveStats(fluidMesh);
      TestFunction z0ObjectiveStats(p0ObjectiveStats);

      LinearForm lfG(z0ObjectiveStats);
      lfG = FaceIntegral(G, z0ObjectiveStats).over(shapeInterface);
      lfG.assemble();

      GridFunction oneObjectiveStats(p0ObjectiveStats);
      oneObjectiveStats = 1.0;

      Alert::Info() << "   | Curvature checkpoint: computing interface area." << Alert::Raise;
      const double shapeArea = computeInterfaceArea(th);
      const double areaViolation = shapeArea - targetSurfaceArea;
      bool areaActive = areaViolation >= -areaActiveToleranceRatio * targetSurfaceArea;
      const double meanG = (shapeArea > 0.0) ? lfG(oneObjectiveStats) / shapeArea : 0.0;

      Alert::Info() << "   | Curvature checkpoint: tracing obstacle face normal." << Alert::Raise;
      auto n = FaceNormal(th);
      n.traceOf(Obstacle);

      Alert::Info() << "   | Curvature checkpoint: complete." << Alert::Raise;
      Alert::Info() << "   | Solving Hilbert identification problems for DJ, DV and DA." << Alert::Raise;

      P1 vh(th, d);

      // Discrete V-inner product:
      //   aV(u, v) = ∫_D alpha^2 ∇u:∇v + ∫_D u·v
      TrialFunction innerU(vh);
      TestFunction innerV(vh);
      BilinearForm innerVForm(innerU, innerV);
      innerVForm =
          Integral(alpha * alpha * Jacobian(innerU), Jacobian(innerV))
        + Integral(innerU, innerV);
      innerVForm.assemble();

      // gradJ in V: <gradJ, psi>_V = DJ(psi) = ∫_Γ G (psi·n)
      PETSc::Variational::TrialFunction gradJTrial(vh); gradJTrial.setName("gradJ");
      PETSc::Variational::TestFunction gradJTest(vh);

      Problem identifyJ(gradJTrial, gradJTest);
      identifyJ =
          Integral(alpha * alpha * Jacobian(gradJTrial), Jacobian(gradJTest))
        + Integral(gradJTrial, gradJTest)
        - FaceIntegral(G * Dot(n, gradJTest)).over(shapeInterface)
        + DirichletBC(gradJTrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      solveKSPChecked(
        identifyJ,
        "identify gradJ",
        kspRtol,
        static_cast<PetscInt>(kspMaxIt));

      // gradV in V: <gradV, psi>_V = DV(psi) = ∫_Γ (psi·n)
      PETSc::Variational::TrialFunction gradVTrial(vh); gradVTrial.setName("gradV");
      PETSc::Variational::TestFunction gradVTest(vh);

      Problem identifyV(gradVTrial, gradVTest);
      identifyV =
          Integral(alpha * alpha * Jacobian(gradVTrial), Jacobian(gradVTest))
        + Integral(gradVTrial, gradVTest)
        - FaceIntegral(Dot(n, gradVTest)).over(shapeInterface)
        + DirichletBC(gradVTrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      solveKSPChecked(
        identifyV,
        "identify gradV",
        kspRtol,
        static_cast<PetscInt>(kspMaxIt));

      // The weak kappa projection above has the opposite sign from the
      // obstacle outward-normal area derivative.
      // gradA in V: <gradA, psi>_V = DA(psi) = -∫_Γ kappa * (psi·n)
      PETSc::Variational::TrialFunction gradATrial(vh); gradATrial.setName("gradA");
      PETSc::Variational::TestFunction gradATest(vh);

      Problem identifyA(gradATrial, gradATest);
      identifyA =
          Integral(alpha * alpha * Jacobian(gradATrial), Jacobian(gradATest))
        + Integral(gradATrial, gradATest)
        + FaceIntegral(kappa * Dot(n, gradATest)).over(shapeInterface)
        + DirichletBC(gradATrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      solveKSPChecked(
        identifyA,
        "identify gradA",
        kspRtol,
        static_cast<PetscInt>(kspMaxIt));

      GridFunction gradJ(vh);
      gradJ = gradJTrial.getSolution();
      GridFunction gradV(vh);
      gradV = gradVTrial.getSolution();
      GridFunction gradA(vh);
      gradA = gradATrial.getSolution();

      const Real innerJV = innerVForm(gradJ, gradV);
      const Real innerJA = innerVForm(gradJ, gradA);
      const Real innerVV = innerVForm(gradV, gradV);
      const Real innerVA = innerVForm(gradV, gradA);
      const Real innerAA = innerVForm(gradA, gradA);
      const Real determinant = innerVV * innerAA - innerVA * innerVA;
      const Real determinantThreshold =
        gramDetToleranceRatio * std::max<Real>(1.0, innerVV * innerAA);
      const Real determinantRelative =
        determinant / std::max<Real>(innerVV * innerAA, Real(1e-30));
      const bool gramDegenerate = determinant <= determinantThreshold;
      Real areaCorrectionViolation = areaViolation;
      bool areaRangeActive = areaActive;
      bool areaNullActive = areaActive && !gramDegenerate;

      Real solveVV = innerVV;
      Real solveVA = innerVA;
      Real solveAA = innerAA;
      Real solveDeterminant = determinant;
      const bool useAreaSystem = areaRangeActive && !gramDegenerate;

      if (!useAreaSystem)
        areaNullActive = false;

      GridFunction xiJ(vh);
      xiJ = gradJ;
      GridFunction xiC(vh);
      xiC = VectorFunction{0.0, 0.0, 0.0};

      Real projCoeffV = 0.0;
      Real projCoeffA = 0.0;
      Real rangeCoeffV = 0.0;
      Real rangeCoeffA = 0.0;
      Real trialProjCoeffV = 0.0;
      Real trialProjCoeffA = 0.0;
      Real areaDualMu = 0.0;

      if (useAreaSystem)
      {
        trialProjCoeffV = (solveAA * innerJV - solveVA * innerJA) / solveDeterminant;
        trialProjCoeffA = (solveVV * innerJA - solveVA * innerJV) / solveDeterminant;
        areaDualMu = -trialProjCoeffA;

        if (areaNullActive)
        {
          projCoeffV = trialProjCoeffV;
          projCoeffA = trialProjCoeffA;
        }
        else
        {
          projCoeffV = innerJV / std::max(innerVV, Real(1e-30));
        }

        rangeCoeffV = (solveAA * volumeViolation - solveVA * areaCorrectionViolation) / solveDeterminant;
        rangeCoeffA = (solveVV * areaCorrectionViolation - solveVA * volumeViolation) / solveDeterminant;
      }
      else
      {
        projCoeffV = innerJV / std::max(innerVV, Real(1e-30));
        rangeCoeffV = volumeViolation / std::max(innerVV, Real(1e-30));
      }

      {
        GridFunction tmp(vh);
        tmp = gradV;
        tmp *= projCoeffV;
        xiJ -= tmp;
      }
      if (areaNullActive && useAreaSystem)
      {
        GridFunction tmp(vh);
        tmp = gradA;
        tmp *= projCoeffA;
        xiJ -= tmp;
      }

      xiC = gradV;
      xiC *= rangeCoeffV;
      if (areaRangeActive && useAreaSystem)
      {
        GridFunction tmp(vh);
        tmp = gradA;
        tmp *= rangeCoeffA;
        xiC += tmp;
      }

      P1 shNormSpace(th);
      GridFunction normGradJ(shNormSpace);
      normGradJ = Frobenius(gradJ);
      GridFunction normGradV(shNormSpace);
      normGradV = Frobenius(gradV);
      GridFunction normGradA(shNormSpace);
      normGradA = Frobenius(gradA);
      GridFunction normXiJ(shNormSpace);
      normXiJ = Frobenius(xiJ);
      GridFunction normXiC(shNormSpace);
      normXiC = Frobenius(xiC);

      const Real maxGradJ = normGradJ.max();
      const Real maxGradV = normGradV.max();
      const Real maxGradA = normGradA.max();
      const Real maxXiJ = normXiJ.max();
      const Real maxXiC = normXiC.max();

      const Real volumeOrthResidual = innerVForm(xiJ, gradV);
      const Real areaOrthResidual = innerVForm(xiJ, gradA);

      Alert::Info()
        << "   | Identification diag: <gradJ,gradV>=" << innerJV
        << ", <gradJ,gradA>=" << innerJA
        << ", <gradV,gradV>=" << innerVV
        << ", <gradV,gradA>=" << innerVA
        << ", <gradA,gradA>=" << innerAA
        << ", distMin=" << distMin
        << ", distMax=" << distMax
        << ", max|dist|=" << maxAbsDist
        << ", ellKappa=" << ellKappa
        << ", kappaMin=" << kappaMin
        << ", kappaMax=" << kappaMax
        << ", max|kappa|=" << maxAbsKappa
        << ", det=" << determinant
        << ", detRel=" << determinantRelative
        << ", areaActive=" << areaActive
        << ", areaRangeActive=" << areaRangeActive
        << ", areaNullActive=" << areaNullActive
        << ", areaDualMu=" << areaDualMu
        << ", areaCorrection=" << areaCorrectionViolation
        << ", max|gradJ|=" << maxGradJ
        << ", max|gradV|=" << maxGradV
        << ", max|gradA|=" << maxGradA
        << ", max|xiJ|=" << maxXiJ
        << ", <xiJ,gradV>_V=" << volumeOrthResidual
        << ", <xiJ,gradA>_V=" << areaOrthResidual
        << ", mean(G)=" << meanG
        << ", Vobs=" << obstacleVolume
        << ", Vtarget=" << targetObstacleVolume
        << ", cV=" << volumeViolation
        << ", A=" << shapeArea
        << ", Amax=" << targetSurfaceArea
        << ", cA=" << areaViolation
        << Alert::Raise;

      if (areaRangeActive && !useAreaSystem)
      {
        Alert::Warning()
          << "   | area constraint fallback to volume-only: det(G)=" << determinant
          << ", detSolve=" << solveDeterminant
          << ", threshold=" << determinantThreshold
          << ", detRel=" << determinantRelative
          << Alert::Raise;
      }

      Alert::Info()
        << "   | Range space diag: projCoeffV=" << projCoeffV
        << ", projCoeffA=" << projCoeffA
        << ", rangeCoeffV=" << rangeCoeffV
        << ", rangeCoeffA=" << rangeCoeffA
        << ", max|xiC|=" << maxXiC
        << Alert::Raise;
      //   αJ = nsAlphaJParam * hmin / ||ξJ||∞
      //   αC = min(0.9, nsAlphaCParam * hmin / max(1e-9, ||ξC||∞))
      const Real alphaJ = nsAlphaJParam * hmin / std::max(maxXiJ, Real(1e-30));
      const Real alphaC = std::min(Real(0.9), nsAlphaCParam * hmin / std::max(maxXiC, Real(1e-9)));

      //   θ = -(αJ · ξJ + αC · ξC)
      GridFunction thetaField(vh);
      Real areaDebugAlpha = 0.0;
      if (areaDebugMode == AreaDebugMode::Off)
      {
        thetaField = xiJ;
        thetaField *= -alphaJ;
        {
          GridFunction tmp(vh);
          tmp = xiC;
          tmp *= -alphaC;
          thetaField += tmp;
        }

        Alert::Info()
          << "  alphaJ=" << alphaJ
          << ", alphaC=" << alphaC
          << ", stepK=" << stepK
          << ", dispJ=" << stepK * alphaJ * maxXiJ
          << ", dispC=" << stepK * alphaC * maxXiC
          << Alert::Raise;
      }
      else
      {
        areaDebugAlpha = nsAlphaCParam * hmin / std::max(maxGradA, Real(1e-9));
        thetaField = gradA;
        thetaField *=
          (areaDebugMode == AreaDebugMode::Descent ? -areaDebugAlpha : areaDebugAlpha);

        Alert::Warning()
          << "   | AREA_DEBUG_MODE active: mode="
          << (areaDebugMode == AreaDebugMode::Descent ? "descent" : "ascent")
          << ", stepK=" << stepK
          << ", alphaA=" << areaDebugAlpha
          << ", dispA=" << stepK * areaDebugAlpha * maxGradA
          << Alert::Raise;
      }

      // ----------------------------------------------------------------------
      // Advection
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Advecting level set." << Alert::Raise;

      // αJ and αC encode the spatially adaptive direction scaling.
      // STEP_K remains the global advection step size multiplier.
      P1 shSpeed(th);
      GridFunction speed(shSpeed);
      speed = Frobenius(thetaField);
      const Real vmax = speed.max();

      Alert::Info()
        << "   | Advection: vmax=" << vmax
        << ", ||theta||∞=" << vmax
        << Alert::Raise;

      TrialFunction advect(sh);
      TestFunction test(sh);
      Advection::Lagrangian(advect, test, dist, thetaField).step(stepK);

      // ----------------------------------------------------------------------
      // XDMF snapshot before remeshing
      // ----------------------------------------------------------------------
      domainGrid.clear();
      domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);
      domainGrid.add("dist",   dist,   IO::XDMF::Center::Node);
      domainGrid.add("theta",  thetaField,  IO::XDMF::Center::Node);
      domainGrid.add("advect", advect.getSolution(), IO::XDMF::Center::Node);

      stateGrid.clear();
      stateGrid.setMesh(fluidMesh, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", u, IO::XDMF::Center::Node);
      stateGrid.add("p", p, IO::XDMF::Center::Node);
      stateGrid.add("v", v, IO::XDMF::Center::Node);
      stateGrid.add("q", q, IO::XDMF::Center::Node);

      xdmf.write(it).flush();

      // ----------------------------------------------------------------------
      // Remeshing
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Remeshing." << Alert::Raise;

      bool remeshAccepted = false;
      MeshType acceptedMesh;
      ObjectiveSample acceptedObjective;
      Real acceptedStepK = 0.0;
      Real acceptedHMax = 0.0;
      Real acceptedHMin = 0.0;
      Real acceptedHausd = 0.0;
      double acceptedReferenceJraw = lastAcceptedJraw;
      double acceptedReferenceJ = lastAcceptedJ;
      size_t acceptedAttempt = 0;
      if (!hasLastAcceptedObjective)
      {
        lastAcceptedJraw = Jraw;
        lastAcceptedJ = J;
        hasLastAcceptedObjective = true;
        Alert::Info()
          << "   | Initializing last accepted objective: Jraw=" << lastAcceptedJraw
          << ", Jopt=" << lastAcceptedJ
          << Alert::Raise;
      }

      const Real lineSearchReferenceJ = static_cast<Real>(lastAcceptedJ);
      const Real jLineSearchTolerance =
        jLineSearchRtol * std::max<Real>(1.0, std::abs(lineSearchReferenceJ));
      const size_t maxRemeshAttempts = std::max<size_t>(1, remeshQualityMaxRetries);
      const Real evalOptHMax = hmax0;
      const Real evalOptHMin = hminRatio * evalOptHMax;
      const Real evalOptHausd = hausdRatio * evalOptHMin;
      for (size_t attempt = 0; attempt < maxRemeshAttempts; ++attempt)
      {
        Real attemptHMax = hmax;
        Real attemptHMin = hmin;
        Real attemptHausd = hausd;
        Real attemptStepK = stepK;
        if (attempt == 1)
        {
          attemptStepK = Real(0.5) * stepK;
        }
        else if (attempt == 2)
        {
          attemptStepK = Real(0.25) * stepK;
        }
        else if (attempt == 3)
        {
          attemptStepK = Real(0.5) * stepK;
          attemptHMin = std::max(attemptHMin, Real(0.2) * attemptHMax);
          attemptHausd = std::max(attemptHausd, Real(0.2) * attemptHMin);
        }
        else if (attempt == 4)
        {
          attemptStepK = Real(0.25) * stepK;
          attemptHMin = std::max(attemptHMin, Real(0.3) * attemptHMax);
          attemptHausd = std::max(attemptHausd, Real(0.25) * attemptHMin);
        }
        else if (attempt >= 5)
        {
          attemptStepK = std::max(Real(0.05) * stepK, stepK / std::pow(Real(2.0), static_cast<Real>(attempt)));
          attemptHMax *= Real(0.75);
          attemptHMin = std::max(attemptHMin, Real(0.3) * attemptHMax);
          attemptHausd = std::max(attemptHausd, Real(0.25) * attemptHMin);
        }

        Alert::Info()
          << "   | Remesh attempt " << attempt
          << ": stepK=" << attemptStepK
          << ", hmax=" << attemptHMax
          << ", hmin=" << attemptHMin
          << ", hausd=" << attemptHausd
          << Alert::Raise;

        try
        {
          TrialFunction attemptAdvect(sh);
          TestFunction attemptTest(sh);
          Advection::Lagrangian(attemptAdvect, attemptTest, dist, thetaField).step(attemptStepK);

          MeshType candidateRaw =
            MMG::LevelSetDiscretizer()
              .split(Fluid,    {Fluid, Obstacle})
              .split(Obstacle, {Fluid, Obstacle})
              .setRMC(1e-5)
              .setHMax(attemptHMax)
              .setHMin(attemptHMin)
              .setHausdorff(attemptHausd)
              .setAngleDetection(false)
              .setBoundaryReference(GammaShape)
              .setBaseReferences(baseOuterBdr)
              .discretize(attemptAdvect.getSolution());

          // MMG may leave a layer of stale interface faces tagged as obstacle.
          // Clear them so the next optimization step sees the intended interface.
          for (auto fit = candidateRaw.getFace(); fit; ++fit)
          {
            if (fit->getAttribute() == Geometry::Attribute(Obstacle))
              candidateRaw.setAttribute(fit.key(), {});
          }
          clearSamePhaseInternalFaces(candidateRaw, "after remesh raw", true);

          if (adaptSizeMapEnabled)
          {
            {
              auto& candidateConn = candidateRaw.getConnectivity();
              candidateConn.compute(2, 3);
              candidateConn.compute(3, 2);
              candidateConn.compute(2, 1);
              candidateConn.compute(1, 0);
              candidateConn.compute(0, 0);
            }

            const Real adaptHNear =
              adaptHNearParam > 0.0 ? adaptHNearParam : attemptHMin;
            const Real adaptHFar =
              adaptHFarParam > 0.0 ? adaptHFarParam : attemptHMax;
            const Real adaptRNear =
              adaptRNearParam > 0.0 ? adaptRNearParam : Real(2.0) * adaptHNear;
            const Real adaptRFar =
              adaptRFarParam > 0.0
                ? adaptRFarParam
                : std::max(attemptHMax, adaptRNear + adaptHNear);
            const Real adaptWidth =
              std::max(adaptRFar - adaptRNear, Real(1e-12));

            P1 adaptSpace(candidateRaw);
            GridFunction adaptDist(adaptSpace);
            Distance::Eikonal(adaptDist)
              .setInterior(Fluid)
              .setInterface(shapeInterface)
              .solve()
              .sign();

            MMG::RealGridFunction adaptSizeMap(adaptSpace);
            auto& sizeData = adaptSizeMap.getData();
            const auto& distData = adaptDist.getData();
            if (sizeData.size() != distData.size())
              throw std::runtime_error("Adapt size map and distance field have incompatible sizes.");

            for (auto i = decltype(sizeData.size())(0); i < sizeData.size(); ++i)
            {
              const Real d = std::abs(static_cast<Real>(distData[i]));
              const Real t = std::clamp((d - adaptRNear) / adaptWidth, Real(0.0), Real(1.0));
              sizeData[i] = adaptHNear + t * (adaptHFar - adaptHNear);
            }

            Alert::Info()
              << "   | Adapting candidate size map: hNear=" << adaptHNear
              << ", hFar=" << adaptHFar
              << ", rNear=" << adaptRNear
              << ", rFar=" << adaptRFar
              << ", gradation=" << adaptGradation
              << Alert::Raise;

            MMG::Adapt()
              .setHMax(adaptHFar)
              .setHMin(adaptHNear)
              .setHausdorff(attemptHausd)
              .setGradation(std::max(Real(1.0), adaptGradation))
              .setAngleDetection(false)
              .adapt(candidateRaw, adaptSizeMap);

            for (auto fit = candidateRaw.getFace(); fit; ++fit)
            {
              if (fit->getAttribute() == Geometry::Attribute(Obstacle))
                candidateRaw.setAttribute(fit.key(), {});
            }
            clearSamePhaseInternalFaces(candidateRaw, "after remesh adapt", true);
          }

          const auto remeshQuality = evaluateMeshQuality(candidateRaw, meshQualityConfig);
          logMeshQuality(
            remeshQuality,
            "raw candidate after remesh attempt " + std::to_string(attempt),
            remeshQuality.accepted());

          Alert::Info()
            << "   | Optimizing eval candidate copy before objective: hmax=" << evalOptHMax
            << ", hmin=" << evalOptHMin
            << ", hausd=" << evalOptHausd
            << Alert::Raise;

          MeshType candidateEval = candidateRaw;
          try
          {
            MMG::Optimizer()
              .setHMax(evalOptHMax)
              .setHMin(evalOptHMin)
              .setHausdorff(evalOptHausd)
              .setAngleDetection(false)
              .optimize(candidateEval);
          }
          catch (const Alert::Exception&)
          {
            Alert::Warning()
              << "Eval candidate optimizer failed at iteration " << it
              << ", attempt=" << attempt
              << Alert::Raise;
            fLineSearch
              << it << " " << attempt << " " << attemptStepK << " "
              << attemptHMax << " " << attemptHMin << " " << attemptHausd << " "
              << "eval_optimizer_reject nan nan " << Jraw << " " << J << " "
              << lastAcceptedJraw << " " << lastAcceptedJ << " "
              << jLineSearchTolerance << " "
              << remeshQuality.minTetVolume << " " << remeshQuality.minAltitudeRatio << " "
              << remeshQuality.badVolumeCells << " " << remeshQuality.badAltitudeCells << "\n";
            fLineSearch.flush();
            continue;
          }
          clearSamePhaseInternalFaces(candidateEval, "after eval candidate optimizer", true);

          const auto evalQuality = evaluateMeshQuality(candidateEval, meshQualityConfig);
          logMeshQuality(
            evalQuality,
            "eval candidate after optimizer attempt " + std::to_string(attempt),
            evalQuality.accepted());

          ObjectiveSample candidateObjective;
          try
          {
            candidateObjective =
              evaluateObjectiveOnMesh(candidateEval, "iteration " + std::to_string(it) +
                " attempt " + std::to_string(attempt) + " eval-optimized");
          }
          catch (const std::exception& e)
          {
            Alert::Warning()
              << "Line-search candidate rejected by objective evaluation at iteration " << it
              << ", attempt=" << attempt
              << ": " << e.what()
              << Alert::Raise;
            fLineSearch
              << it << " " << attempt << " " << attemptStepK << " "
              << attemptHMax << " " << attemptHMin << " " << attemptHausd << " "
              << "objective_eval_reject nan nan " << Jraw << " " << J << " "
              << lastAcceptedJraw << " " << lastAcceptedJ << " "
              << jLineSearchTolerance << " "
              << remeshQuality.minTetVolume << " " << remeshQuality.minAltitudeRatio << " "
              << remeshQuality.badVolumeCells << " " << remeshQuality.badAltitudeCells << "\n";
            fLineSearch.flush();
            continue;
          }

          const bool objectiveAccepted =
            candidateObjective.J <= lineSearchReferenceJ + jLineSearchTolerance;
          Alert::Info()
            << "   | J line-search attempt " << attempt
            << ": Jraw_before=" << Jraw
            << ", Jraw_after=" << candidateObjective.Jraw
            << ", J_before=" << J
            << ", J_after=" << candidateObjective.J
            << ", Jraw_last_accepted=" << lastAcceptedJraw
            << ", J_last_accepted=" << lastAcceptedJ
            << ", tol=" << jLineSearchTolerance
            << ", accepted=" << (objectiveAccepted ? 1 : 0)
            << Alert::Raise;

          fLineSearch
            << it << " " << attempt << " " << attemptStepK << " "
            << attemptHMax << " " << attemptHMin << " " << attemptHausd << " "
            << (objectiveAccepted ? "accepted" : "j_reject") << " "
            << candidateObjective.Jraw << " " << candidateObjective.J << " "
            << Jraw << " " << J << " " << jLineSearchTolerance << " "
            << lastAcceptedJraw << " " << lastAcceptedJ << " "
            << remeshQuality.minTetVolume << " " << remeshQuality.minAltitudeRatio << " "
            << remeshQuality.badVolumeCells << " " << remeshQuality.badAltitudeCells << "\n";
          fLineSearch.flush();

          if (!objectiveAccepted)
            continue;

          const auto rawCheckpointPath =
            checkpointDir / ("candidate_raw_accepted." + std::to_string(it) + ".mesh");
          const auto evalCheckpointPath =
            checkpointDir / ("candidate_eval_optimized." + std::to_string(it) + ".mesh");

          saveArtifactMesh(
            candidateRaw,
            rawCheckpointPath,
            rawCheckpointPath.string());
          saveArtifactMesh(
            candidateEval,
            evalCheckpointPath,
            evalCheckpointPath.string());

          const auto rawCheckpointStats =
            computeMeshCheckpointStats(candidateRaw, meshQualityConfig);
          const auto evalCheckpointStats =
            computeMeshCheckpointStats(candidateEval, meshQualityConfig);
          logMeshCheckpointStats(
            rawCheckpointStats,
            "candidate_raw_accepted." + std::to_string(it));
          logMeshCheckpointStats(
            evalCheckpointStats,
            "candidate_eval_optimized." + std::to_string(it));

          previousEvalCheckpointStats = evalCheckpointStats;
          previousEvalCheckpointIteration = it;
          hasPreviousEvalCheckpoint = true;

          acceptedMesh = std::move(candidateRaw);
          acceptedObjective = candidateObjective;
          acceptedStepK = attemptStepK;
          acceptedHMax = attemptHMax;
          acceptedHMin = attemptHMin;
          acceptedHausd = attemptHausd;
          acceptedReferenceJraw = lastAcceptedJraw;
          acceptedReferenceJ = lastAcceptedJ;
          acceptedAttempt = attempt;
          remeshAccepted = true;
          break;
        }
        catch (const Alert::Exception&)
        {
          Alert::Warning()
            << "Remeshing attempt failed at iteration " << it
            << ", attempt=" << attempt
            << Alert::Raise;
        }
      }

      if (!remeshAccepted)
      {
        hmax  /= 2;
        hmin   = 0.33 * hmax;
        hausd  = 0.1 * hmin;

        Alert::Warning()
          << "All remesh candidates failed or were rejected. Keeping previous mesh and retrying with hmax=" << hmax
          << ", hmin=" << hmin
          << ", hausd=" << hausd
          << Alert::Raise;
        continue;
      }

      th = std::move(acceptedMesh);
      saveArtifactMesh(
        th,
        std::filesystem::path("out") / ("Omega." + std::to_string(it) + ".mesh"),
        "out/Omega." + std::to_string(it) + ".mesh");

      const double acceptedVolumeViolation = acceptedObjective.obstacleVolume - targetObstacleVolume;
      const double acceptedAreaViolation = acceptedObjective.shapeArea - targetSurfaceArea;
      const bool acceptedAreaActive =
        acceptedAreaViolation >= -areaActiveToleranceRatio * targetSurfaceArea;

      lastAcceptedJraw = acceptedObjective.Jraw;
      lastAcceptedJ = acceptedObjective.J;

      hmax  = hmax0;
      hmin  = hminRatio * hmax;
      hausd = hausdRatio * hmin;

      // ----------------------------------------------------------------------
      // Logging
      // ----------------------------------------------------------------------
      fObj << acceptedObjective.J << "\n";
      fObjRaw << acceptedObjective.Jraw << "\n";
      fObjBefore << Jraw << " " << J << "\n";
      fVol << acceptedObjective.obstacleVolume << " " << acceptedVolumeViolation << " "
           << acceptedObjective.shapeArea << " " << acceptedAreaViolation << " "
           << (acceptedAreaActive ? 1 : 0) << "\n";
      fNS  << alphaJ << " " << alphaC << " " << (areaActive ? 1 : 0) << " " << projCoeffV << " " << projCoeffA << " " << rangeCoeffV << " " << rangeCoeffA << " " << maxXiJ << " " << maxXiC << "\n";
      fCentroid
        << acceptedObjective.centroid(0) << " "
        << (acceptedObjective.centroid.size() > 1 ? acceptedObjective.centroid(1) : 0.0) << " "
        << (acceptedObjective.centroid.size() > 2 ? acceptedObjective.centroid(2) : 0.0) << "\n";

      fObj.flush();
      fObjRaw.flush();
      fObjBefore.flush();
      fVol.flush();
      fNS.flush();
      fCentroid.flush();

      Alert::Info()
        << "   | Accepted line-search attempt=" << acceptedAttempt
        << ", acceptedStepK=" << acceptedStepK
        << ", acceptedHMax=" << acceptedHMax
        << ", acceptedHMin=" << acceptedHMin
        << ", acceptedHausd=" << acceptedHausd
        << ", Jraw_before=" << Jraw
        << ", Jraw_after=" << acceptedObjective.Jraw
        << ", Jopt_before=" << J
        << ", Jopt_after=" << acceptedObjective.J
        << ", Jraw_reference=" << acceptedReferenceJraw
        << ", J_reference=" << acceptedReferenceJ
        << ", Vobs=" << acceptedObjective.obstacleVolume
        << ", cV=" << acceptedVolumeViolation
        << ", A=" << acceptedObjective.shapeArea
        << ", Amax=" << targetSurfaceArea
        << ", cA=" << acceptedAreaViolation
        << ", areaActive=" << acceptedAreaActive
        << ", areaNullActive=" << areaNullActive
        << ", mean(G)=" << meanG
        << ", alphaJ=" << alphaJ
        << ", alphaC=" << alphaC
        << ", projCoeffV=" << projCoeffV
        << ", projCoeffA=" << projCoeffA
        << ", rangeCoeffV=" << rangeCoeffV
        << ", rangeCoeffA=" << rangeCoeffA
        << ", max|xiJ|=" << maxXiJ
        << ", max|xiC|=" << maxXiC
        << Alert::Raise;

      jrawHistory.push_back(acceptedObjective.Jraw);
      if (convergenceWindow >= 2 && convergenceRtolJraw > 0.0 && jrawHistory.size() >= convergenceWindow)
      {
        const double referenceJraw = jrawHistory[jrawHistory.size() - convergenceWindow];
        const double scale =
          std::max(1.0, std::max(std::abs(acceptedObjective.Jraw), std::abs(referenceJraw)));
        const double relativeChange = std::abs(acceptedObjective.Jraw - referenceJraw) / scale;
        if (relativeChange < convergenceRtolJraw)
        {
          Alert::Info()
            << "   | Converged by Jraw criterion: window=" << convergenceWindow
            << ", rel_change=" << relativeChange
            << ", tol=" << convergenceRtolJraw
            << Alert::Raise;
          break;
        }
      }
    }

    xdmf.close();

    Alert::Success()
      << "Saved final mesh sequence to out/Omega.*.mesh"
      << Alert::Raise;

    PetscFinalize();
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << algorithmName << " failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
