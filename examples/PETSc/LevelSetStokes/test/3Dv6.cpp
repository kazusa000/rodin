/*
 * 3Dv6: null-space gradient flow with
 *   - volume equality constraint
 *   - surface-area upper-bound inequality constraint
 *   - centroid-based rigid-body rotations
 *   - optional whole-mesh translation via SHIFT_{X,Y,Z}
 *
 * This combines the area-constrained flow of 3Dv4 with the centroid-based
 * rotational semantics of 3Dv5. The active-set treatment remains simplified
 * for a single inequality: when the area bound is active, it is treated as a
 * second tangent constraint; otherwise the algorithm falls back to the 3Dv2
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
static constexpr Real   defaultHilbertAlphaFactor = 16.0;
static constexpr Real   areaActiveToleranceRatio = 1e-4;
static constexpr Real   gramDetToleranceRatio = 1e-12;
static constexpr Real   defaultAreaCorrectionGain = 0.1;
static constexpr Real   areaMultiplierTolerance = 1e-10;
static constexpr const char* algorithmName = "3Dv6";

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
    "Unsupported 3Dv6 objective mode. Use OBJECTIVE_MODE=K, OBJECTIVE_MODE=C or OBJECTIVE_MODE=Q.");
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
    const Real hilbertAlphaFactor =
      LevelSetStokes::Runtime::envDouble("HILBERT_ALPHA_FACTOR", defaultHilbertAlphaFactor);
    if (!(hilbertAlphaFactor > 0.0))
      throw std::runtime_error("HILBERT_ALPHA_FACTOR must be positive.");
    const AreaDebugMode areaDebugMode =
      parseAreaDebugMode(LevelSetStokes::Runtime::envCString("AREA_DEBUG_MODE", "off"));
    const Real areaCorrectionGain =
      LevelSetStokes::Runtime::envDouble("AREA_CORRECTION_GAIN", defaultAreaCorrectionGain);
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
    const Real alpha = hilbertAlphaFactor * (hmax0 - hmin);

    // Load and pre-optimize initial mesh
    MeshType th;
    th.load(meshFile, IO::FileFormat::MEDIT);

    MMG::Optimizer()
      .setHMax(hmax)
      .setHMin(hmin)
      .setHausdorff(hausd)
      .optimize(th);

    translateMesh(th, shift);
    th.save("Omega0.mesh", IO::FileFormat::MEDIT);

    // XDMF output
    std::filesystem::create_directories("out");
    IO::XDMF xdmf(std::string("out/") + algorithmName);

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fObjRaw("obj_raw.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fNS("ns.txt");
    std::ofstream fCentroid("centroid.txt");

    if (!fObj || !fObjRaw || !fVol || !fNS || !fCentroid)
      throw std::runtime_error("Failed to open output history files.");

    // Constraint targets measured on the pre-optimized initial mesh.
    const double targetObstacleVolume = th.getVolume(Obstacle);
    const double initialShapeArea = computeInterfaceArea(th);
    const double targetSurfaceArea = surfaceAreaFactor * initialShapeArea;
    std::vector<double> jrawHistory;
    jrawHistory.reserve(maxIt);
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

    Alert::Info()
      << algorithmName << " configuration:"
      << " mode=" << (objectiveMode == ObjectiveMode::K ? "K" : (objectiveMode == ObjectiveMode::C ? "C" : "Q"))
      << ", sense=" << (objectiveSense == ObjectiveSense::Min ? "min" : "max")
      << ", iAxis=" << iAxis
      << ", jAxis=" << jAxis
      << ", shift=(" << shift[0] << ", " << shift[1] << ", " << shift[2] << ")"
      << ", surfaceAreaFactor=" << surfaceAreaFactor
      << ", hilbertAlphaFactor=" << hilbertAlphaFactor
      << ", hilbertAlpha=" << alpha
      << ", areaCorrectionGain=" << areaCorrectionGain
      << Alert::Raise;

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
      bool areaActive = areaViolation >= -areaActiveToleranceRatio * targetSurfaceArea;
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
        << ", int_Gamma(kappa^2)=" << interfaceKappaSq
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
      fCentroid
        << obstacleCentroid(0) << " "
        << (obstacleCentroid.size() > 1 ? obstacleCentroid(1) : 0.0) << " "
        << (obstacleCentroid.size() > 2 ? obstacleCentroid(2) : 0.0) << "\n";

      fObj.flush();
      fObjRaw.flush();
      fVol.flush();
      fNS.flush();
      fCentroid.flush();

      Alert::Info()
        << "   | Jraw=" << Jraw
        << ", Jopt=" << J
        << ", Vobs=" << obstacleVolume
        << ", cV=" << volumeViolation
        << ", A=" << shapeArea
        << ", Amax=" << targetSurfaceArea
        << ", cA=" << areaViolation
        << ", areaActive=" << areaActive
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

      th.save("out/Omega." + std::to_string(it) + ".mesh", IO::FileFormat::MEDIT);

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
