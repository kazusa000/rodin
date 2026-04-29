/*
 * Standalone post-processing smooth pass for final 3D obstacle meshes.
 *
 * This executable reuses the final level-set smoothing idea from 3Dv3+,
 * but runs independently from the optimization loop. It reads an existing
 * obstacle mesh, rebuilds a signed distance field, applies scalar Helmholtz
 * smoothing, optionally shifts the zero iso-surface to recover volume, and
 * discretizes a new smoothed mesh.
 */

#include <Rodin/MMG.h>
#include <Rodin/PETSc.h>

#include <Rodin/Assembly.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/Geometry.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

#include <petscsys.h>

#include "RuntimeConfig.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef LEVELSET_STOKES_MESH_FILE
#define LEVELSET_STOKES_MESH_FILE "3dshapes/sphere_init_aligned.mesh"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

using MeshType = Rodin::MMG::Mesh;

static constexpr Attribute Obstacle = 2;
static constexpr Attribute Fluid = 3;
static constexpr Attribute GammaShape = 13;

static constexpr Real defaultHMax = 0.4;
static constexpr Real defaultFinalHmaxFactor = 0.1;
static constexpr Real defaultFinalHminRatio = 0.1;
static constexpr Real defaultFinalHausdRatio = 3.0;
static constexpr Real defaultFinalRmc = 1e-4;
static constexpr size_t defaultSmoothSteps = 1;
static constexpr Real defaultSmoothEpsFactor = 1.0;
static constexpr Real defaultSmoothIsoShift = 0.0;
static constexpr Real defaultFeatureSmoothKappaFactor = 1.0;
static constexpr Real defaultFeatureSmoothMinWeight = 0.08;

struct ShapeStats
{
  double volume = 0.0;
  double area = 0.0;
  double kappaMean = 0.0;
  double kappaStd = 0.0;
  double kappaMaxAbs = 0.0;
};

struct MeshPhaseStats
{
  size_t obstacleCells = 0;
  size_t fluidCells = 0;
  size_t interfaceFaces = 0;
};

static void computeConnectivity(MeshType& mesh)
{
  auto& conn = mesh.getConnectivity();
  conn.compute(2, 3);
  conn.compute(3, 2);
  conn.compute(2, 1);
  conn.compute(1, 0);
  conn.compute(0, 0);
}

static void rebuildInterfaceAttribute(MeshType& mesh)
{
  for (auto it = mesh.getFace(); !it.end(); ++it)
  {
    if (it->getAttribute() == Geometry::Attribute(GammaShape))
      mesh.setAttribute(it.key(), {});
  }

  mesh.trace(
    Map<Pair<Attribute, Attribute>, Attribute>{
      {{Obstacle, Fluid}, GammaShape},
    }
  );
}

static FlatSet<Attribute> detectBaseOuterBoundary(MeshType& mesh)
{
  FlatSet<Attribute> baseOuterBdr;
  for (auto it = mesh.getFace(); !it.end(); ++it)
  {
    const auto& face = *it;
    if (!face.isBoundary())
      continue;

    const auto a = face.getAttribute();
    if (!a || *a == GammaShape)
      continue;

    baseOuterBdr.insert(*a);
  }

  return baseOuterBdr;
}

static ShapeStats computeShapeStats(MeshType& mesh, const auto& distField, Real ellKappa)
{
  P1 kappaSpace(mesh);
  PETSc::Variational::TrialFunction kappaT(kappaSpace); kappaT.setName("kappa_postsmooth_diag");
  PETSc::Variational::TestFunction kappaW(kappaSpace);

  Problem curvProj(kappaT, kappaW);
  curvProj =
      Integral(kappaT, kappaW)
    + ellKappa * ellKappa * Integral(Grad(kappaT), Grad(kappaW))
    + Integral(Grad(distField), Grad(kappaW));
  curvProj.assemble();
  Solver::KSP(curvProj).solve();

  GridFunction kappa(kappaSpace);
  kappa = kappaT.getSolution();

  GridFunction absKappa(kappaSpace);
  absKappa = Abs(kappa);

  P0g p0Stats(mesh);
  TestFunction z0Stats(p0Stats);

  LinearForm lfArea(z0Stats);
  lfArea = FaceIntegral(z0Stats).over({GammaShape});
  lfArea.assemble();

  LinearForm lfKappa(z0Stats);
  lfKappa = FaceIntegral(kappa, z0Stats).over({GammaShape});
  lfKappa.assemble();

  LinearForm lfKappaSq(z0Stats);
  lfKappaSq = FaceIntegral(kappa * kappa, z0Stats).over({GammaShape});
  lfKappaSq.assemble();

  GridFunction oneStats(p0Stats);
  oneStats = 1.0;

  ShapeStats stats;
  stats.volume = mesh.getVolume(Obstacle);
  stats.area = lfArea(oneStats);
  const double mean = stats.area > 0.0 ? lfKappa(oneStats) / stats.area : 0.0;
  const double meanSq = stats.area > 0.0 ? lfKappaSq(oneStats) / stats.area : 0.0;
  stats.kappaMean = mean;
  stats.kappaStd = std::sqrt(std::max(0.0, meanSq - mean * mean));
  stats.kappaMaxAbs = absKappa.max();
  return stats;
}

static MeshPhaseStats computePhaseStats(MeshType& mesh)
{
  MeshPhaseStats stats;
  for (auto it = mesh.getCell(); !it.end(); ++it)
  {
    const auto attr = it->getAttribute();
    if (!attr)
      continue;
    if (*attr == Obstacle)
      stats.obstacleCells++;
    else if (*attr == Fluid)
      stats.fluidCells++;
  }

  for (auto it = mesh.getFace(); !it.end(); ++it)
  {
    const auto attr = it->getAttribute();
    if (attr && *attr == GammaShape)
      stats.interfaceFaces++;
  }

  return stats;
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  std::cout << std::unitbuf;

  try
  {
    const char* meshFile =
      LevelSetStokes::Runtime::envCString("MESH", LEVELSET_STOKES_MESH_FILE);
    const char* outputMeshFile =
      LevelSetStokes::Runtime::envCString("OUTPUT_MESH", "out/Omega.postsmooth.mesh");
    const Real hmax =
      LevelSetStokes::Runtime::envDouble("HMAX", defaultHMax);
    const Real finalHmaxFactor =
      LevelSetStokes::Runtime::envDouble("FINAL_HMAX_FACTOR", defaultFinalHmaxFactor);
    const Real finalHminRatio =
      LevelSetStokes::Runtime::envDouble("FINAL_HMIN_RATIO", defaultFinalHminRatio);
    const Real finalHausdRatio =
      LevelSetStokes::Runtime::envDouble("FINAL_HAUSD_RATIO", defaultFinalHausdRatio);
    const Real finalRmc =
      LevelSetStokes::Runtime::envDouble("FINAL_RMC", defaultFinalRmc);
    const size_t smoothSteps =
      LevelSetStokes::Runtime::envSizeT("SMOOTH_STEPS", defaultSmoothSteps);
    const Real smoothEpsFactor =
      LevelSetStokes::Runtime::envDouble("SMOOTH_EPS_FACTOR", defaultSmoothEpsFactor);
    const Real smoothIsoShift =
      LevelSetStokes::Runtime::envDouble("SMOOTH_ISO_SHIFT", defaultSmoothIsoShift);
    const std::string smoothMode =
      LevelSetStokes::Runtime::envCString("SMOOTH_MODE", "global");
    const bool featureSmooth =
      smoothMode == "feature" || smoothMode == "FEATURE" || smoothMode == "Feature";
    if (!featureSmooth && smoothMode != "global" && smoothMode != "GLOBAL" && smoothMode != "Global")
      throw std::runtime_error("Unsupported SMOOTH_MODE. Use global or feature.");
    const Real featureSmoothKappaFactor =
      LevelSetStokes::Runtime::envDouble("FEATURE_SMOOTH_KAPPA_FACTOR", defaultFeatureSmoothKappaFactor);
    const Real featureSmoothMinWeight = std::clamp(
      LevelSetStokes::Runtime::envDouble("FEATURE_SMOOTH_MIN_WEIGHT", defaultFeatureSmoothMinWeight),
      Real(0.0),
      Real(1.0));

    const Real finalDiscHmax = std::max(Real(1e-6), finalHmaxFactor * hmax);
    const Real finalDiscHmin = finalHminRatio * finalDiscHmax;
    const Real finalDiscHausd = finalHausdRatio * finalDiscHmin;
    const Real smoothEps = std::pow(smoothEpsFactor * finalDiscHmax, 2);
    const Real diagEllKappa = std::max(Real(1e-8), 2.0 * finalDiscHmax);

    std::filesystem::create_directories("out");

    std::ofstream fSmooth("postsmooth.txt");
    if (!fSmooth)
      throw std::runtime_error("Failed to open postsmooth.txt.");

    fSmooth
      << "# smooth_steps smooth_eps smooth_iso_shift_manual smooth_iso_shift_auto smooth_iso_shift_correct smooth_iso_shift_total "
      << "V_target V_smooth V_before V_after dV_rel "
      << "A_smooth A_before A_after dA_rel "
      << "kappa_mean_before kappa_mean_smooth kappa_mean_after "
      << "kappa_std_before kappa_std_smooth kappa_std_after "
      << "kappa_max_before kappa_max_smooth kappa_max_after\n";

    MeshType th;
    th.load(meshFile, IO::FileFormat::MEDIT);
    computeConnectivity(th);
    rebuildInterfaceAttribute(th);
    const FlatSet<Attribute> baseOuterBdr = detectBaseOuterBoundary(th);
    const FlatSet<Attribute> finalShapeInterface{GammaShape};

    auto discretizePhi = [&](const auto& phiField)
    {
      MeshType mesh =
        MMG::LevelSetDiscretizer()
          .split(Fluid, {Fluid, Obstacle})
          .split(Obstacle, {Fluid, Obstacle})
          .setRMC(finalRmc)
          .setHMax(finalDiscHmax)
          .setHMin(finalDiscHmin)
          .setHausdorff(finalDiscHausd)
          .setAngleDetection(false)
          .setBoundaryReference(GammaShape)
          .setBaseReferences(baseOuterBdr)
          .discretize(phiField);

      for (auto fit = mesh.getFace(); fit; ++fit)
      {
        if (fit->getAttribute() == Geometry::Attribute(Obstacle))
          mesh.setAttribute(fit.key(), {});
      }

      computeConnectivity(mesh);
      rebuildInterfaceAttribute(mesh);
      return mesh;
    };

    P1 finalSh(th);
    GridFunction finalDist(finalSh);
    Distance::Eikonal(finalDist)
      .setInterior(Fluid)
      .setInterface(finalShapeInterface)
      .solve()
      .sign();

    finalDist.save("out/Omega.postsmooth.dist.gf", IO::FileFormat::MFEM);
    const ShapeStats beforeStats = computeShapeStats(th, finalDist, diagEllKappa);

    GridFunction currentPhi(finalSh);
    currentPhi = finalDist;

    GridFunction smoothWeight(finalSh);
    smoothWeight = 1.0;

    if (featureSmooth)
    {
      PETSc::Variational::TrialFunction kappaT(finalSh); kappaT.setName("kappa_postsmooth_feature");
      PETSc::Variational::TestFunction kappaW(finalSh);

      Problem curvProj(kappaT, kappaW);
      curvProj =
          Integral(kappaT, kappaW)
        + diagEllKappa * diagEllKappa * Integral(Grad(kappaT), Grad(kappaW))
        + Integral(Grad(finalDist), Grad(kappaW));
      curvProj.assemble();
      Solver::KSP(curvProj).solve();

      GridFunction kappa(finalSh);
      kappa = kappaT.getSolution();

      const Real kappaScale =
        std::max(Real(1e-12), featureSmoothKappaFactor / std::max(finalDiscHmax, Real(1e-12)));
      smoothWeight =
        featureSmoothMinWeight
        + (1.0 - featureSmoothMinWeight) / (1.0 + Pow(Abs(kappa) / kappaScale, 2));

      GridFunction negWeight(finalSh);
      negWeight = -smoothWeight;
      smoothWeight.save("out/Omega.postsmooth.smooth_weight.gf", IO::FileFormat::MFEM);

      Alert::Info()
        << "Post smooth feature weight: kappa_scale=" << kappaScale
        << ", min_weight=" << featureSmoothMinWeight
        << ", weight_min=" << -negWeight.max()
        << ", weight_max=" << smoothWeight.max()
        << Alert::Raise;
    }

    for (size_t smoothIt = 0; smoothIt < smoothSteps; ++smoothIt)
    {
      PETSc::Variational::TrialFunction smoothPhi(finalSh); smoothPhi.setName("phi_postsmooth");
      PETSc::Variational::TestFunction smoothW(finalSh);

      Problem smoothProblem(smoothPhi, smoothW);
      smoothProblem =
          smoothEps * Integral(smoothWeight * Grad(smoothPhi), Grad(smoothW))
        + Integral(smoothPhi, smoothW)
        - Integral(currentPhi, smoothW);

      Solver::KSP(smoothProblem).solve();
      currentPhi = smoothPhi.getSolution();
    }

    currentPhi.save("out/Omega.postsmooth.smooth.gf", IO::FileFormat::MFEM);

    MeshType smoothTh = discretizePhi(currentPhi);
    const MeshPhaseStats smoothPhaseStats = computePhaseStats(smoothTh);

    Alert::Info()
      << "Post smooth zero-iso mesh: obstacle_cells=" << smoothPhaseStats.obstacleCells
      << ", fluid_cells=" << smoothPhaseStats.fluidCells
      << ", interface_faces=" << smoothPhaseStats.interfaceFaces
      << Alert::Raise;

    P1 smoothSh(smoothTh);
    GridFunction smoothDist(smoothSh);
    Distance::Eikonal(smoothDist)
      .setInterior(Fluid)
      .setInterface(finalShapeInterface)
      .solve()
      .sign();

    const ShapeStats smoothStats = computeShapeStats(smoothTh, smoothDist, diagEllKappa);
    const double targetVolume = beforeStats.volume;
    const double autoIsoShift =
      smoothStats.area > 0.0 ? (targetVolume - smoothStats.volume) / smoothStats.area : 0.0;

    double isoShiftCorrection = 0.0;
    double totalIsoShift = autoIsoShift + smoothIsoShift;

    auto buildShiftedResult = [&](double isoShift)
    {
      GridFunction phi(smoothSh);
      phi = smoothDist + isoShift;

      MeshType mesh = discretizePhi(phi);
      P1 meshSh(mesh);
      GridFunction meshDist(meshSh);
      Distance::Eikonal(meshDist)
        .setInterior(Fluid)
        .setInterface(finalShapeInterface)
        .solve()
        .sign();
      return std::make_pair(std::move(mesh), computeShapeStats(mesh, meshDist, diagEllKappa));
    };

    auto optimizeFinalMesh = [&](MeshType& mesh)
    {
      Alert::Info()
        << "Post smooth final optimizer: hmax=" << finalDiscHmax
        << ", hmin=" << finalDiscHmin
        << ", hausd=" << finalDiscHausd
        << Alert::Raise;

      MMG::Optimizer()
        .setHMax(finalDiscHmax)
        .setHMin(finalDiscHmin)
        .setHausdorff(finalDiscHausd)
        .setAngleDetection(false)
        .optimize(mesh);

      computeConnectivity(mesh);
      rebuildInterfaceAttribute(mesh);
    };

    auto computeOptimizedStats = [&](MeshType& mesh)
    {
      P1 meshSh(mesh);
      GridFunction meshDist(meshSh);
      Distance::Eikonal(meshDist)
        .setInterior(Fluid)
        .setInterface(finalShapeInterface)
        .solve()
        .sign();
      return computeShapeStats(mesh, meshDist, diagEllKappa);
    };

    auto shiftedResult = buildShiftedResult(totalIsoShift);
    MeshType finalTh = std::move(shiftedResult.first);
    ShapeStats afterStats = shiftedResult.second;

    if (afterStats.area > 0.0)
    {
      isoShiftCorrection = (targetVolume - afterStats.volume) / afterStats.area;
      totalIsoShift += isoShiftCorrection;
      shiftedResult = buildShiftedResult(totalIsoShift);
      finalTh = std::move(shiftedResult.first);
      afterStats = shiftedResult.second;
    }

    const double volumeBeforeFinalOptimizer = afterStats.volume;
    const double areaBeforeFinalOptimizer = afterStats.area;
    optimizeFinalMesh(finalTh);
    afterStats = computeOptimizedStats(finalTh);

    GridFunction shiftedPhi(smoothSh);
    shiftedPhi = smoothDist + totalIsoShift;
    shiftedPhi.save("out/Omega.postsmooth.levelset.gf", IO::FileFormat::MFEM);

    const double volumeRelChange =
      beforeStats.volume > 0.0
        ? std::abs(afterStats.volume - beforeStats.volume) / beforeStats.volume
        : 0.0;
    const double areaRelChange =
      beforeStats.area > 0.0
        ? std::abs(afterStats.area - beforeStats.area) / beforeStats.area
        : 0.0;

    finalTh.save(outputMeshFile, IO::FileFormat::MEDIT);

    fSmooth
      << smoothSteps << ' '
      << smoothEps << ' '
      << smoothIsoShift << ' '
      << autoIsoShift << ' '
      << isoShiftCorrection << ' '
      << totalIsoShift << ' '
      << targetVolume << ' '
      << smoothStats.volume << ' '
      << beforeStats.volume << ' '
      << afterStats.volume << ' '
      << volumeRelChange << ' '
      << smoothStats.area << ' '
      << beforeStats.area << ' '
      << afterStats.area << ' '
      << areaRelChange << ' '
      << beforeStats.kappaMean << ' '
      << smoothStats.kappaMean << ' '
      << afterStats.kappaMean << ' '
      << beforeStats.kappaStd << ' '
      << smoothStats.kappaStd << ' '
      << afterStats.kappaStd << ' '
      << beforeStats.kappaMaxAbs << ' '
      << smoothStats.kappaMaxAbs << ' '
      << afterStats.kappaMaxAbs << '\n';
    fSmooth.flush();

    Alert::Info()
      << "Post smooth configuration: hmax=" << finalDiscHmax
      << ", hmin=" << finalDiscHmin
      << ", hausd=" << finalDiscHausd
      << ", rmc=" << finalRmc
      << ", smooth_steps=" << smoothSteps
      << ", smooth_eps=" << smoothEps
      << ", manual_iso_shift=" << smoothIsoShift
      << ", smooth_mode=" << (featureSmooth ? "feature" : "global")
      << ", output_mesh=" << outputMeshFile
      << ", feature_kappa_factor=" << featureSmoothKappaFactor
      << ", feature_min_weight=" << featureSmoothMinWeight
      << Alert::Raise;

    Alert::Info()
      << "Post smooth iso shift: auto=" << autoIsoShift
      << ", correct=" << isoShiftCorrection
      << ", total=" << totalIsoShift
      << Alert::Raise;

    Alert::Info()
      << "Post smooth diagnostics: V_before=" << beforeStats.volume
      << ", V_after=" << afterStats.volume
      << ", dV_rel=" << volumeRelChange
      << ", A_before=" << beforeStats.area
      << ", A_after=" << afterStats.area
      << ", dA_rel=" << areaRelChange
      << Alert::Raise;

    Alert::Info()
      << "Post smooth final optimizer diagnostics: V_before_opt=" << volumeBeforeFinalOptimizer
      << ", V_after_opt=" << afterStats.volume
      << ", A_before_opt=" << areaBeforeFinalOptimizer
      << ", A_after_opt=" << afterStats.area
      << Alert::Raise;

    Alert::Info()
      << "Post smooth curvature: "
      << "mean_before=" << beforeStats.kappaMean
      << ", mean_smooth=" << smoothStats.kappaMean
      << ", mean_after=" << afterStats.kappaMean
      << ", std_before=" << beforeStats.kappaStd
      << ", std_smooth=" << smoothStats.kappaStd
      << ", std_after=" << afterStats.kappaStd
      << ", max_before=" << beforeStats.kappaMaxAbs
      << ", max_smooth=" << smoothStats.kappaMaxAbs
      << ", max_after=" << afterStats.kappaMaxAbs
      << Alert::Raise;

    Alert::Success()
      << "Saved smoothed mesh to " << outputMeshFile
      << Alert::Raise;

    PetscFinalize();
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "3DPostSmooth failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
