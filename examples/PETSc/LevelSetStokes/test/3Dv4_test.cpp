/*
 * 3Dv4: null-space gradient flow with
 *   - volume equality constraint
 *   - surface-area upper-bound inequality constraint
 *
 * The active-set treatment is intentionally simplified for a single
 * inequality: when the area bound is active, it is treated as a second
 * tangent constraint; otherwise the algorithm falls back to the 3Dv2
 * volume-only null-space flow.
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
static constexpr Real   defaultConvergenceRtolJraw = 5e-3;
static constexpr Real   defaultNsAlphaJ = 0.5;
static constexpr Real   defaultNsAlphaC = 0.5;
static constexpr Real   defaultSurfaceAreaFactor = 1.05;
static constexpr int    defaultFinalRefine = 1;
static constexpr Real   defaultFinalHmaxFactor = 0.1;
static constexpr Real   defaultFinalHminRatio = 0.1;
static constexpr Real   defaultFinalHausdRatio = 3.0;
static constexpr Real   defaultFinalRmc = 1e-4;
static constexpr size_t defaultSmoothSteps = 1;
static constexpr Real   defaultSmoothEpsFactor = 1.0;
static constexpr Real   defaultSmoothIsoShift = 0.0;
static constexpr Real   areaActiveToleranceRatio = 1e-4;
static constexpr Real   gramDetToleranceRatio = 1e-12;

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
    "Unsupported 3Dv4 objective mode. Use OBJECTIVE_MODE=K, OBJECTIVE_MODE=C or OBJECTIVE_MODE=Q.");
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
    const bool finalRefineEnabled =
      LevelSetStokes::Runtime::envInt("FINAL_REFINE", defaultFinalRefine) != 0;
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
    const AreaDebugMode areaDebugMode =
      parseAreaDebugMode(LevelSetStokes::Runtime::envCString("AREA_DEBUG_MODE", "off"));
    const ObjectiveMode objectiveMode =
      parseObjectiveMode(LevelSetStokes::Runtime::envCString("OBJECTIVE_MODE", "K"));
    const ObjectiveSense objectiveSense =
      parseObjectiveSense(LevelSetStokes::Runtime::envCString("OBJECTIVE_SENSE", "min"));

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

    MMG::Optimizer()
      .setHMax(hmax)
      .setHMin(hmin)
      .setHausdorff(hausd)
      .optimize(th);

    th.save("Omega0.mesh", IO::FileFormat::MEDIT);

    // XDMF output
    std::filesystem::create_directories("out");
    IO::XDMF xdmf("out/3Dv4");

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fObjRaw("obj_raw.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fNS("ns.txt");
    std::ofstream fSmooth("phi_smooth.txt");

    if (!fObj || !fObjRaw || !fVol || !fNS || !fSmooth)
      throw std::runtime_error("Failed to open output history files.");

    fSmooth
      << "# smooth_steps smooth_eps smooth_iso_shift_manual smooth_iso_shift_auto smooth_iso_shift_correct smooth_iso_shift_total "
      << "V_target V_smooth "
      << "V_before V_after dV_rel "
      << "A_smooth "
      << "A_before A_after dA_rel "
      << "kappa_mean_before kappa_mean_smooth kappa_mean_after "
      << "kappa_std_before kappa_std_smooth kappa_std_after "
      << "kappa_max_before kappa_max_smooth kappa_max_after\n";

    // Constraint targets measured on the pre-optimized initial mesh.
    const double targetObstacleVolume = th.getVolume(Obstacle);
    const double initialShapeArea = computeInterfaceArea(th);
    const double targetSurfaceArea = surfaceAreaFactor * initialShapeArea;
    struct ShapeDiagnostics
    {
      double volume = 0.0;
      double area = 0.0;
      double kappaMean = 0.0;
      double kappaStd = 0.0;
      double kappaMaxAbs = 0.0;
    };

    auto computeShapeDiagnostics = [&](MeshType& mesh, const auto& distField) -> ShapeDiagnostics
    {
      P1 kappaSpace(mesh);
      PETSc::Variational::TrialFunction kappaT(kappaSpace); kappaT.setName("kappa");
      PETSc::Variational::TestFunction kappaW(kappaSpace);

      Problem curvProj(kappaT, kappaW);
      curvProj =
          Integral(kappaT, kappaW)
        + Integral(Grad(distField), Grad(kappaW));
      curvProj.assemble();
      Solver::KSP(curvProj).solve();

      auto kappa = GridFunction(kappaSpace);
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

      const double area = lfArea(oneStats);
      const double mean = area > 0.0 ? lfKappa(oneStats) / area : 0.0;
      const double meanSq = area > 0.0 ? lfKappaSq(oneStats) / area : 0.0;

      ShapeDiagnostics stats;
      stats.volume = mesh.getVolume(Obstacle);
      stats.area = area;
      stats.kappaMean = mean;
      stats.kappaStd = std::sqrt(std::max(0.0, meanSq - mean * mean));
      stats.kappaMaxAbs = absKappa.max();
      return stats;
    };
    std::vector<double> jrawHistory;
    jrawHistory.reserve(maxIt);
    // Objective basis vectors
    const auto ei = basis(iAxis);
    const auto ej = basis(jAxis);
    const bool rotationalState = objectiveMode == ObjectiveMode::C || objectiveMode == ObjectiveMode::Q;
    const bool rotationalAdjoint = objectiveMode == ObjectiveMode::Q;
    const Real stateOmegaX = rotationalState ? ej[0] : 0.0;
    const Real stateOmegaY = rotationalState ? ej[1] : 0.0;
    const Real stateOmegaZ = rotationalState ? ej[2] : 0.0;
    const Real stateTx = objectiveMode == ObjectiveMode::K ? ej[0] : 0.0;
    const Real stateTy = objectiveMode == ObjectiveMode::K ? ej[1] : 0.0;
    const Real stateTz = objectiveMode == ObjectiveMode::K ? ej[2] : 0.0;

    const auto stateVelocity = VectorFunction{
      [=](const Geometry::Point& p)
      {
        return stateTx + stateOmegaY * p.z() - stateOmegaZ * p.y();
      },
      [=](const Geometry::Point& p)
      {
        return stateTy + stateOmegaZ * p.x() - stateOmegaX * p.z();
      },
      [=](const Geometry::Point& p)
      {
        return stateTz + stateOmegaX * p.y() - stateOmegaY * p.x();
      }
    };

    const auto adjointVelocity = VectorFunction{
      [=](const Geometry::Point& p)
      {
        return rotationalAdjoint ? ei[1] * p.z() - ei[2] * p.y() : ei[0];
      },
      [=](const Geometry::Point& p)
      {
        return rotationalAdjoint ? ei[2] * p.x() - ei[0] * p.z() : ei[1];
      },
      [=](const Geometry::Point& p)
      {
        return rotationalAdjoint ? ei[0] * p.y() - ei[1] * p.x() : ei[2];
      }
    };

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

    const FlatSet<Attribute> finalShapeInterface{GammaShape};
    bool hasCompletedIteration = false;
    size_t lastCompletedIteration = 0;
    auto runFinalSmooth = [&](size_t iterationIndex)
    {
      if (!finalRefineEnabled)
        return;

      const Real finalDiscHmax = std::max(Real(1e-6), finalHmaxFactor * hmax);
      const Real finalDiscHmin = finalHminRatio * finalDiscHmax;
      const Real finalDiscHausd = finalHausdRatio * finalDiscHmin;
      const Real smoothEps = std::pow(smoothEpsFactor * finalDiscHmax, 2);

      Alert::Info()
        << "   | Final phi smooth discretize: hmax=" << finalDiscHmax
        << ", hmin=" << finalDiscHmin
        << ", hausd=" << finalDiscHausd
        << ", rmc=" << finalRmc
        << Alert::Raise;

      Alert::Info()
        << "   | Final phi smooth Helmholtz: steps=" << smoothSteps
        << ", eps=" << smoothEps
        << ", iso_shift=" << smoothIsoShift
        << Alert::Raise;

      try
      {
        auto discretizeFinalPhi = [&](const auto& phiField)
        {
          MeshType mesh =
            MMG::LevelSetDiscretizer()
              .split(Fluid,    {Fluid, Obstacle})
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

          auto& conn = mesh.getConnectivity();
          conn.compute(2, 3);
          conn.compute(3, 2);
          conn.compute(2, 1);
          conn.compute(1, 0);
          conn.compute(0, 0);
          return mesh;
        };

        {
          auto& conn = th.getConnectivity();
          conn.compute(2, 3);
          conn.compute(3, 2);
          conn.compute(2, 1);
          conn.compute(1, 0);
          conn.compute(0, 0);
        }

        P1 finalSh(th);
        GridFunction finalDist(finalSh);
        Distance::Eikonal(finalDist)
          .setInterior(Fluid)
          .setInterface(finalShapeInterface)
          .solve()
          .sign();

        finalDist.save("out/Omega.final.dist.gf", IO::FileFormat::MFEM);
        const ShapeDiagnostics beforeStats = computeShapeDiagnostics(th, finalDist);

        GridFunction currentPhi(finalSh);
        currentPhi = finalDist;

        for (size_t smoothIt = 0; smoothIt < smoothSteps; ++smoothIt)
        {
          PETSc::Variational::TrialFunction smoothPhi(finalSh); smoothPhi.setName("phi_smooth");
          PETSc::Variational::TestFunction smoothW(finalSh);

          Problem smoothProblem(smoothPhi, smoothW);
          smoothProblem =
              smoothEps * Integral(Grad(smoothPhi), Grad(smoothW))
            + Integral(smoothPhi, smoothW)
            - Integral(currentPhi, smoothW);

          Solver::KSP(smoothProblem).solve();
          currentPhi = smoothPhi.getSolution();
        }

        currentPhi.save("out/Omega.final.smooth.gf", IO::FileFormat::MFEM);

        MeshType smoothTh = discretizeFinalPhi(currentPhi);

        P1 smoothSh(smoothTh);
        GridFunction smoothDist(smoothSh);
        Distance::Eikonal(smoothDist)
          .setInterior(Fluid)
          .setInterface(finalShapeInterface)
          .solve()
          .sign();

        const ShapeDiagnostics smoothStats = computeShapeDiagnostics(smoothTh, smoothDist);
        const double targetVolume = beforeStats.volume;
        const double autoIsoShift =
          smoothStats.area > 0.0
            ? (targetVolume - smoothStats.volume) / smoothStats.area
            : 0.0;
        double isoShiftCorrection = 0.0;
        double totalIsoShift = autoIsoShift + smoothIsoShift;

        auto buildShiftedResult = [&](double isoShift)
        {
          GridFunction phi(smoothSh);
          phi = smoothDist + isoShift;

          MeshType mesh = discretizeFinalPhi(phi);

          P1 meshSh(mesh);
          GridFunction meshDist(meshSh);
          Distance::Eikonal(meshDist)
            .setInterior(Fluid)
            .setInterface(finalShapeInterface)
            .solve()
            .sign();

          return std::make_pair(std::move(mesh), computeShapeDiagnostics(mesh, meshDist));
        };

        Alert::Info()
          << "   | Final phi smooth discretize."
          << Alert::Raise;

        auto shiftedResult = buildShiftedResult(totalIsoShift);
        MeshType finalTh = std::move(shiftedResult.first);
        ShapeDiagnostics afterStats = shiftedResult.second;

        if (afterStats.area > 0.0)
        {
          isoShiftCorrection = (targetVolume - afterStats.volume) / afterStats.area;
          totalIsoShift += isoShiftCorrection;
          shiftedResult = buildShiftedResult(totalIsoShift);
          finalTh = std::move(shiftedResult.first);
          afterStats = shiftedResult.second;
        }

        GridFunction shiftedPhi(smoothSh);
        shiftedPhi = smoothDist + totalIsoShift;
        shiftedPhi.save("out/Omega.final.levelset.gf", IO::FileFormat::MFEM);

        const double volumeRelChange =
          beforeStats.volume > 0.0
            ? std::abs(afterStats.volume - beforeStats.volume) / beforeStats.volume
            : 0.0;
        const double areaRelChange =
          beforeStats.area > 0.0
            ? std::abs(afterStats.area - beforeStats.area) / beforeStats.area
            : 0.0;

        finalTh.save("out/Omega.final.mesh", IO::FileFormat::MEDIT);
        finalTh.save("out/Omega." + std::to_string(iterationIndex + 1) + ".mesh", IO::FileFormat::MEDIT);

        const double finalVolume = finalTh.getVolume(Obstacle);

        auto finalFluidMesh = finalTh.trim(Obstacle);
        domainGrid.clear();
        domainGrid.setMesh(finalTh, IO::XDMF::MeshPolicy::Transient);
        stateGrid.clear();
        stateGrid.setMesh(finalFluidMesh, IO::XDMF::MeshPolicy::Transient);
        xdmf.write(iterationIndex + 1).flush();

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
          << "   | Final phi smooth iso shift: "
          << "manual=" << smoothIsoShift
          << ", auto=" << autoIsoShift
          << ", correct=" << isoShiftCorrection
          << ", total=" << totalIsoShift
          << ", V_target=" << targetVolume
          << ", V_smooth=" << smoothStats.volume
          << ", A_smooth=" << smoothStats.area
          << Alert::Raise;

        Alert::Info()
          << "   | Final phi smooth curvature(smooth zero-iso): "
          << "mean_smooth=" << smoothStats.kappaMean
          << ", std_smooth=" << smoothStats.kappaStd
          << ", max_smooth=" << smoothStats.kappaMaxAbs
          << Alert::Raise;

        Alert::Info()
          << "   | Final phi smooth diagnostics: "
          << "V_before=" << beforeStats.volume
          << ", V_after=" << afterStats.volume
          << ", dV_rel=" << volumeRelChange
          << ", A_before=" << beforeStats.area
          << ", A_after=" << afterStats.area
          << ", dA_rel=" << areaRelChange
          << Alert::Raise;

        Alert::Info()
          << "   | Final phi smooth curvature: "
          << "mean_before=" << beforeStats.kappaMean
          << ", mean_after=" << afterStats.kappaMean
          << ", std_before=" << beforeStats.kappaStd
          << ", std_after=" << afterStats.kappaStd
          << ", max_before=" << beforeStats.kappaMaxAbs
          << ", max_after=" << afterStats.kappaMaxAbs
          << Alert::Raise;

        Alert::Info()
          << "   | Final phi smooth saved: out/Omega.final.mesh and out/Omega." << (iterationIndex + 1)
          << ".mesh"
          << ", Vobs_final=" << finalVolume
          << Alert::Raise;
      }
      catch (const std::exception& e)
      {
        Alert::Warning()
          << "Final phi smooth failed: " << e.what()
          << ". Keeping coarse final mesh."
          << Alert::Raise;
      }
    };

    for (size_t it = 0; it < maxIt; ++it)
    {
      Alert::Info() << "----- Iteration: " << it << Alert::Raise;

      // ----------------------------------------------------------------------
      // Optimize current mesh
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Optimizing the domain." << Alert::Raise;
      try
      {
        MMG::Optimizer()
          .setHMax(hmax)
          .setHMin(hmin)
          .setHausdorff(hausd)
          .setAngleDetection(false)
          .optimize(th);

        hmax  = hmax0;
        hmin  = hminRatio * hmax;
        hausd = hausdRatio * hmin;
      }
      catch (const Alert::Exception&)
      {
        hmax  /= 2;
        hmin   = retryHMinRatio * hmax;
        hausd  = hausdRatio * hmin;

        Alert::Warning()
          << "Mesh optimization failed at iteration " << it
          << ". Reducing hmax to " << hmax
          << " and retrying."
          << Alert::Raise;
        continue;
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

      Solver::KSP(state).solve();

      auto u = GridFunction(velocitySpace);
      auto p = GridFunction(pressureSpace);
      u = up.getSolution();
      p = pp.getSolution();

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

      Solver::KSP(adj).solve();

      auto v = GridFunction(velocitySpace);
      auto q = GridFunction(pressureSpace);
      v = ua.getSolution();
      q = pa.getSolution();

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
      Solver::KSP(curvProj).solve();

      GridFunction kappa(kappaSpace);
      kappa = kappaT.getSolution();

      const Real kappaMax = kappa.max();
      GridFunction negKappa(kappaSpace);
      negKappa = -kappa;
      const Real kappaMin = -negKappa.max();
      const Real maxAbsKappa = std::max(std::abs(kappaMin), std::abs(kappaMax));

      P0g p0SurfaceStats(th);
      TestFunction z0SurfaceStats(p0SurfaceStats);
      LinearForm lfKappaSq(z0SurfaceStats);
      lfKappaSq = FaceIntegral(kappa * kappa, z0SurfaceStats).over(shapeInterface);
      lfKappaSq.assemble();

      GridFunction oneSurfaceStats(p0SurfaceStats);
      oneSurfaceStats = 1.0;
      const Real interfaceKappaSq = lfKappaSq(oneSurfaceStats);

      P0g p0ObjectiveStats(fluidMesh);
      TestFunction z0ObjectiveStats(p0ObjectiveStats);

      LinearForm lfG(z0ObjectiveStats);
      lfG = FaceIntegral(G, z0ObjectiveStats).over(shapeInterface);
      lfG.assemble();

      GridFunction oneObjectiveStats(p0ObjectiveStats);
      oneObjectiveStats = 1.0;

      const double shapeArea = computeInterfaceArea(th);
      const double areaViolation = shapeArea - targetSurfaceArea;
      const bool areaActive = areaViolation >= -areaActiveToleranceRatio * targetSurfaceArea;
      const double meanG = (shapeArea > 0.0) ? lfG(oneObjectiveStats) / shapeArea : 0.0;

      auto n = FaceNormal(th);
      n.traceOf(Obstacle);


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

      Solver::KSP(identifyJ).solve();

      // gradV in V: <gradV, psi>_V = DV(psi) = ∫_Γ (psi·n)
      PETSc::Variational::TrialFunction gradVTrial(vh); gradVTrial.setName("gradV");
      PETSc::Variational::TestFunction gradVTest(vh);

      Problem identifyV(gradVTrial, gradVTest);
      identifyV =
          Integral(alpha * alpha * Jacobian(gradVTrial), Jacobian(gradVTest))
        + Integral(gradVTrial, gradVTest)
        - FaceIntegral(Dot(n, gradVTest)).over(shapeInterface)
        + DirichletBC(gradVTrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP(identifyV).solve();

      // gradA in V: <gradA, psi>_V = DA(psi) = ∫_Γ kappa * (psi·n)
      PETSc::Variational::TrialFunction gradATrial(vh); gradATrial.setName("gradA");
      PETSc::Variational::TestFunction gradATest(vh);

      Problem identifyA(gradATrial, gradATest);
      identifyA =
          Integral(alpha * alpha * Jacobian(gradATrial), Jacobian(gradATest))
        + Integral(gradATrial, gradATest)
        + FaceIntegral(kappa * Dot(n, gradATest)).over(shapeInterface)
        + DirichletBC(gradATrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP(identifyA).solve();

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
      const bool gramDegenerate = determinant <= determinantThreshold;

      GridFunction xiJ(vh);
      xiJ = gradJ;
      GridFunction xiC(vh);
      xiC = VectorFunction{0.0, 0.0, 0.0};

      Real projCoeffV = 0.0;
      Real projCoeffA = 0.0;
      Real rangeCoeffV = 0.0;
      Real rangeCoeffA = 0.0;

      if (areaActive && !gramDegenerate)
      {
        projCoeffV = (innerAA * innerJV - innerVA * innerJA) / determinant;
        projCoeffA = (innerVV * innerJA - innerVA * innerJV) / determinant;
        rangeCoeffV = (innerAA * volumeViolation - innerVA * areaViolation) / determinant;
        rangeCoeffA = (innerVV * areaViolation - innerVA * volumeViolation) / determinant;
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
      if (areaActive && !gramDegenerate)
      {
        GridFunction tmp(vh);
        tmp = gradA;
        tmp *= projCoeffA;
        xiJ -= tmp;
      }

      xiC = gradV;
      xiC *= rangeCoeffV;
      if (areaActive && !gramDegenerate)
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
        << ", int_Gamma(kappa^2)=" << interfaceKappaSq
        << ", det=" << determinant
        << ", areaActive=" << areaActive
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

      if (areaActive && gramDegenerate)
      {
        Alert::Warning()
          << "   | area constraint fallback to volume-only: det(G)=" << determinant
          << ", threshold=" << determinantThreshold
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

      try
      {
        th =
          MMG::LevelSetDiscretizer()
            .split(Fluid,    {Fluid, Obstacle})
            .split(Obstacle, {Fluid, Obstacle})
            .setRMC(1e-5)
            .setHMax(hmax)
            .setHMin(hmin)
            .setHausdorff(hausd)
            .setAngleDetection(false)
            .setBoundaryReference(GammaShape)
            .setBaseReferences(baseOuterBdr)
            .discretize(advect.getSolution());

        // MMG may leave a layer of stale interface faces tagged as obstacle.
        // Clear them so the next optimization step sees the intended interface.
        for (auto fit = th.getFace(); fit; ++fit)
        {
          if (fit->getAttribute() == Geometry::Attribute(Obstacle))
            th.setAttribute(fit.key(), {});
        }
      }
      catch (const Alert::Exception&)
      {
        hmax  /= 2;
        hmin   = 0.33 * hmax;
        hausd  = 0.1 * hmin;

        Alert::Warning()
          << "Remeshing failed. Retrying with hmax=" << hmax
          << Alert::Raise;
        continue;
      }

      // ----------------------------------------------------------------------
      // Logging
      // ----------------------------------------------------------------------
      fObj << J << "\n";
      fObjRaw << Jraw << "\n";
      fVol << obstacleVolume << " " << volumeViolation << " " << shapeArea << " " << areaViolation << " " << (areaActive ? 1 : 0) << "\n";
      fNS  << alphaJ << " " << alphaC << " " << (areaActive ? 1 : 0) << " " << projCoeffV << " " << projCoeffA << " " << rangeCoeffV << " " << rangeCoeffA << " " << maxXiJ << " " << maxXiC << "\n";

      fObj.flush();
      fObjRaw.flush();
      fVol.flush();
      fNS.flush();

      Alert::Info()
        << "   | Jraw=" << Jraw
        << ", Jopt=" << J
        << ", Vobs=" << obstacleVolume
        << ", cV=" << volumeViolation
        << ", A=" << shapeArea
        << ", Amax=" << targetSurfaceArea
        << ", cA=" << areaViolation
        << ", areaActive=" << areaActive
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

      th.save("out/Omega." + std::to_string(it) + ".mesh", IO::FileFormat::MEDIT);
      hasCompletedIteration = true;
      lastCompletedIteration = it;

      jrawHistory.push_back(Jraw);
      if (convergenceWindow >= 2 && convergenceRtolJraw > 0.0 && jrawHistory.size() >= convergenceWindow)
      {
        const double referenceJraw = jrawHistory[jrawHistory.size() - convergenceWindow];
        const double scale =
          std::max(1.0, std::max(std::abs(Jraw), std::abs(referenceJraw)));
        const double relativeChange = std::abs(Jraw - referenceJraw) / scale;
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

    if (hasCompletedIteration)
      runFinalSmooth(lastCompletedIteration);

    xdmf.close();

    Alert::Success()
      << "Saved final mesh sequence to out/Omega.*.mesh"
      << Alert::Raise;

    PetscFinalize();
    return 0;
  }
  catch (const std::exception& e)
  {
    std::cerr << "3Dv4 failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
