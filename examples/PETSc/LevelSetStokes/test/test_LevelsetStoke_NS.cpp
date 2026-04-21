/*
 * Null-space gradient flow variant of LevelSetStokes3DObstacle.
 *
 * Replaces the Augmented Lagrangian volume constraint handling with the
 * null-space gradient flow method (Feppon, Allaire, Dapogny 2020).
 * The descent direction is decomposed into:
 *   ξ_J: null-space component (objective descent within feasible directions)
 *   ξ_C: range-space component (Gauss-Newton constraint correction)
 * with adaptive normalization — no penalty parameters to tune.
 *
 * Semantic choice (unchanged):
 *   - obstacle cells: attribute 2
 *   - fluid cells:    attribute 3
 *   - shape interface: attribute 13
 * The optimization target is the obstacle shape.
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

// Null-space gradient flow step-size constants (O(1), insensitive to tuning)
static constexpr Real A_J = 0.5;
static constexpr Real A_C = 0.5;

struct ArrayDiagnostics
{
  bool finite = true;
  long long badIndex = -1;
  double badValue = 0.0;
  double norm2 = 0.0;
  double normInf = 0.0;
  double min = 0.0;
  double max = 0.0;
};

struct MeshQualityDiagnostics
{
  size_t cellCount = 0;
  size_t badAspectCount = 0;
  double minMeasure = std::numeric_limits<double>::infinity();
  double maxMeasure = 0.0;
  double minEdge = std::numeric_limits<double>::infinity();
  double maxEdge = 0.0;
  double maxAspect = 0.0;
};

template <class Derived>
static ArrayDiagnostics inspectArray(const Eigen::MatrixBase<Derived>& vec)
{
  ArrayDiagnostics out;
  if (vec.size() == 0)
    return out;

  out.norm2 = vec.norm();
  out.normInf = vec.template lpNorm<Eigen::Infinity>();
  out.min = vec.minCoeff();
  out.max = vec.maxCoeff();

  for (Eigen::Index i = 0; i < vec.size(); ++i)
  {
    const double value = vec.derived().coeff(i);
    if (!std::isfinite(value))
    {
      out.finite = false;
      out.badIndex = static_cast<long long>(i);
      out.badValue = value;
      break;
    }
  }
  return out;
}

template <class Derived>
static void printArrayDiagnostics(const std::string& name, const Eigen::MatrixBase<Derived>& vec)
{
  const auto d = inspectArray(vec);
  std::cout
    << "Info:    | " << name
    << ": finite=" << (d.finite ? "true" : "false")
    << ", norm2=" << d.norm2
    << ", normInf=" << d.normInf
    << ", min=" << d.min
    << ", max=" << d.max;
  if (!d.finite)
    std::cout << ", badIndex=" << d.badIndex << ", badValue=" << d.badValue;
  std::cout << '\n';
}

static void printKSPSummary(const std::string& name, ::KSP ksp)
{
  KSPConvergedReason reason;
  PetscInt its = 0;
  PetscReal residual = 0.0;
  KSPGetConvergedReason(ksp, &reason);
  KSPGetIterationNumber(ksp, &its);
  KSPGetResidualNorm(ksp, &residual);

  const char* reasonLabel = "unknown";
  switch (reason)
  {
    case KSP_CONVERGED_RTOL_NORMAL: reasonLabel = "CONVERGED_RTOL_NORMAL"; break;
    case KSP_CONVERGED_RTOL: reasonLabel = "CONVERGED_RTOL"; break;
    case KSP_CONVERGED_ATOL: reasonLabel = "CONVERGED_ATOL"; break;
    case KSP_CONVERGED_ITS: reasonLabel = "CONVERGED_ITS"; break;
    case KSP_DIVERGED_ITS: reasonLabel = "DIVERGED_ITS"; break;
    case KSP_DIVERGED_DTOL: reasonLabel = "DIVERGED_DTOL"; break;
    case KSP_DIVERGED_BREAKDOWN: reasonLabel = "DIVERGED_BREAKDOWN"; break;
    case KSP_DIVERGED_BREAKDOWN_BICG: reasonLabel = "DIVERGED_BREAKDOWN_BICG"; break;
    case KSP_DIVERGED_NONSYMMETRIC: reasonLabel = "DIVERGED_NONSYMMETRIC"; break;
    case KSP_DIVERGED_INDEFINITE_PC: reasonLabel = "DIVERGED_INDEFINITE_PC"; break;
    case KSP_DIVERGED_NANORINF: reasonLabel = "DIVERGED_NANORINF"; break;
    case KSP_DIVERGED_INDEFINITE_MAT: reasonLabel = "DIVERGED_INDEFINITE_MAT"; break;
    case KSP_DIVERGED_PC_FAILED: reasonLabel = "DIVERGED_PC_FAILED"; break;
    default: break;
  }

  Alert::Info()
    << "   | " << name
    << ": reason=" << static_cast<int>(reason)
    << " (" << reasonLabel << ")"
    << ", iterations=" << its
    << ", residual=" << residual
    << Alert::Raise;
}

template <class MeshLike>
static MeshQualityDiagnostics inspectMeshQuality(const MeshLike& mesh)
{
  MeshQualityDiagnostics out;

  for (auto it = mesh.getCell(); it; ++it)
  {
    const auto& cell = *it;
    if (cell.getGeometry() != Geometry::Polytope::Type::Tetrahedron)
      continue;

    out.cellCount++;
    const double measure = cell.getMeasure();
    out.minMeasure = std::min(out.minMeasure, measure);
    out.maxMeasure = std::max(out.maxMeasure, measure);

    std::vector<Geometry::Vertex> vertices;
    vertices.reserve(4);
    for (auto vit = cell.getVertex(); vit; ++vit)
      vertices.push_back(*vit);

    if (vertices.size() < 4)
      continue;

    double localMinEdge = std::numeric_limits<double>::infinity();
    double localMaxEdge = 0.0;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      for (size_t j = i + 1; j < vertices.size(); ++j)
      {
        const double edge = (vertices[i].getCoordinates() - vertices[j].getCoordinates()).norm();
        localMinEdge = std::min(localMinEdge, edge);
        localMaxEdge = std::max(localMaxEdge, edge);
      }
    }

    out.minEdge = std::min(out.minEdge, localMinEdge);
    out.maxEdge = std::max(out.maxEdge, localMaxEdge);

    const double aspect =
      localMinEdge > 0.0 ? localMaxEdge / localMinEdge : std::numeric_limits<double>::infinity();
    out.maxAspect = std::max(out.maxAspect, aspect);
    if (!std::isfinite(aspect) || aspect > 10.0)
      out.badAspectCount++;
  }

  if (out.cellCount == 0)
  {
    out.minMeasure = 0.0;
    out.minEdge = 0.0;
  }

  return out;
}

template <class MeshLike>
static void printMeshQuality(const std::string& name, const MeshLike& mesh)
{
  const auto q = inspectMeshQuality(mesh);
  Alert::Info()
    << "   | " << name
    << ": cells=" << q.cellCount
    << ", minMeasure=" << q.minMeasure
    << ", maxMeasure=" << q.maxMeasure
    << ", minEdge=" << q.minEdge
    << ", maxEdge=" << q.maxEdge
    << ", maxAspect=" << q.maxAspect
    << ", badAspectCount=" << q.badAspectCount
    << Alert::Raise;
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
    "Unsupported 3D NS objective mode. Use OBJECTIVE_MODE=K, OBJECTIVE_MODE=C or OBJECTIVE_MODE=Q.");
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

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Force line-buffered stdout so output appears promptly even under ld-linux loader
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  std::cout << std::unitbuf;

    const Real defaultHMaxRuntime = LevelSetStokes::Runtime::envDouble("HMAX", defaultHMax);
    const Real hminRatio = LevelSetStokes::Runtime::envDouble("HMIN_RATIO", 0.1);
    const Real hausdRatio = LevelSetStokes::Runtime::envDouble("HAUSD_RATIO", 0.1);
    const Real rmc = LevelSetStokes::Runtime::envDouble("RMC", 1e-5);
    const Real retryHMinRatio = hminRatio > 0.33 ? hminRatio : 0.33;

  const Real hmax0 = defaultHMaxRuntime;

  Real hmax  = hmax0;

  try
  {
    const char* meshFile = LevelSetStokes::Runtime::envCString("MESH", LEVELSET_STOKES_MESH_FILE);
    const size_t maxIt   = LevelSetStokes::Runtime::envSizeT("MAX_ITERS", defaultMaxIt);
    const ObjectiveMode objectiveMode =
      parseObjectiveMode(LevelSetStokes::Runtime::envCString("OBJECTIVE_MODE", "K"));
    const ObjectiveSense objectiveSense =
      parseObjectiveSense(LevelSetStokes::Runtime::envCString("OBJECTIVE_SENSE", "min"));

    // Objective component K_ij, default K_11 in 0-based storage -> axis 1,1 as in original file.
    const int iAxis = LevelSetStokes::Runtime::envInt("IAXIS", 1);
    const int jAxis = LevelSetStokes::Runtime::envInt("JAXIS", 1);

    Real hmin  = hminRatio * hmax;
    Real hausd = hausdRatio * hmin;

    const Real k = LevelSetStokes::Runtime::envDouble("STEP_K", 0.1);
    const Real dt    = k * (hmax0 - hmin);
    const Real alpha = 16 * (hmax0 - hmin);

    Alert::Info()
      << "Configuration: HMAX=" << hmax
      << ", HMIN_RATIO=" << hminRatio
      << ", HAUSD_RATIO=" << hausdRatio
      << ", RMC=" << rmc
      << Alert::Raise;

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
    IO::XDMF xdmf("out/LevelSetStokes3DObstacle_test");

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fObjRaw("obj_raw.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fNS("ns.txt");

    if (!fObj || !fObjRaw || !fVol || !fNS)
      throw std::runtime_error("Failed to open output history files.");

    // Volume constraint target
    const double targetObstacleVolume = th.getVolume(Obstacle);

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
      printMeshQuality("Domain mesh quality", th);
      printMeshQuality("Fluid mesh quality", fluidMesh);

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

      Solver::KSP stateSolver(state);
      stateSolver.solve();
      printKSPSummary("State KSP", stateSolver.getHandle());

      auto u = GridFunction(velocitySpace);
      auto p = GridFunction(pressureSpace);
      u = up.getSolution();
      p = pp.getSolution();
      printArrayDiagnostics("u", u.getData());
      printArrayDiagnostics("p", p.getData());

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

      Solver::KSP adjSolver(adj);
      adjSolver.solve();
      printKSPSummary("Adjoint KSP", adjSolver.getHandle());

      auto v = GridFunction(velocitySpace);
      auto q = GridFunction(pressureSpace);
      v = ua.getSolution();
      q = pa.getSolution();
      printArrayDiagnostics("v", v.getData());
      printArrayDiagnostics("q", q.getData());

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
      printArrayDiagnostics("lfObj", lfObj.getVector());

      GridFunction oneObj(p0Obj);
      oneObj = 1.0;
      printArrayDiagnostics("oneObj", oneObj.getData());

      const double Jraw = -lfObj(oneObj);
      const double senseSign = objectiveSense == ObjectiveSense::Max ? -1.0 : 1.0;
      const double J = senseSign * Jraw;

      if (!std::isfinite(Jraw) || !std::isfinite(J))
      {
        Alert::Warning()
          << "Objective became non-finite at iteration " << it
          << ". Saving debug mesh snapshot."
          << Alert::Raise;
        th.save("debug_nan_th.mesh", IO::FileFormat::MEDIT);
        fluidMesh.save("debug_nan_fluid.mesh", IO::FileFormat::MEDIT);
      }

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
      const double violation = obstacleVolume - targetObstacleVolume;

      P0g p0Stats(fluidMesh);
      TestFunction z0Stats(p0Stats);

      LinearForm lfArea(z0Stats);
      lfArea = FaceIntegral(z0Stats).over(shapeInterface);
      lfArea.assemble();

      LinearForm lfG(z0Stats);
      lfG = FaceIntegral(G, z0Stats).over(shapeInterface);
      lfG.assemble();

      GridFunction oneStats(p0Stats);
      oneStats = 1.0;

      const double shapeArea = lfArea(oneStats);
      const double meanG = (shapeArea > 0.0) ? lfG(oneStats) / shapeArea : 0.0;

      // Use the fluid-domain boundary normal and flip it so it matches the
      // obstacle-outward convention used in the objective/shape gradient.
      auto n = -FaceNormal(fluidMesh);


      Alert::Info() << "   | Solving Hilbert identification problems for DJ and Dg." << Alert::Raise;

      P1 vh(fluidMesh, d);

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

      Solver::KSP identifyJSolver(identifyJ);
      identifyJSolver.solve();
      printKSPSummary("IdentifyJ KSP", identifyJSolver.getHandle());

      // gradg in V: <gradg, psi>_V = Dg(psi) = ∫_Γ (psi·n)
      PETSc::Variational::TrialFunction gradgTrial(vh); gradgTrial.setName("gradg");
      PETSc::Variational::TestFunction gradgTest(vh);

      Problem identifyg(gradgTrial, gradgTest);
      identifyg =
          Integral(alpha * alpha * Jacobian(gradgTrial), Jacobian(gradgTest))
        + Integral(gradgTrial, gradgTest)
        - FaceIntegral(Dot(n, gradgTest)).over(shapeInterface)
        + DirichletBC(gradgTrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP identifygSolver(identifyg);
      identifygSolver.solve();
      printKSPSummary("Identifyg KSP", identifygSolver.getHandle());

      GridFunction gradJ(vh);
      gradJ = gradJTrial.getSolution();
      printArrayDiagnostics("gradJ", gradJ.getData());
      GridFunction gradg(vh);
      gradg = gradgTrial.getSolution();
      printArrayDiagnostics("gradg", gradg.getData());

      // Single equality constraint:
      //   xiJ = gradJ - <gradJ, gradg>_V / <gradg, gradg>_V * gradg
      const Real innerJG = innerVForm(gradJ, gradg);
      const Real innerGG = innerVForm(gradg, gradg);
      const Real projCoeff = innerJG / std::max(innerGG, Real(1e-30));

      GridFunction xiJ(vh);
      xiJ = gradJ;
      {
        GridFunction tmp(vh);
        tmp = gradg;
        tmp *= projCoeff;
        xiJ -= tmp;
      }

      P1 shNormSpace(fluidMesh);
      GridFunction normGradJ(shNormSpace);
      normGradJ = Frobenius(gradJ);
      GridFunction normGradg(shNormSpace);
      normGradg = Frobenius(gradg);
      GridFunction normXiJ(shNormSpace);
      normXiJ = Frobenius(xiJ);

      const Real maxGradJ = normGradJ.max();
      const Real maxGradg = normGradg.max();
      const Real maxXiJ = normXiJ.max();
      const Real orthResidual = innerVForm(xiJ, gradg);

      LinearForm lfNullCheck(z0Stats);
      lfNullCheck = FaceIntegral(Dot(n, xiJ), z0Stats).over(shapeInterface);
      lfNullCheck.assemble();
      const Real nullResidual = lfNullCheck(oneStats);

      Alert::Info()
        << "   | Identification diag: <gradJ,gradg>=" << innerJG
        << ", <gradg,gradg>=" << innerGG
        << ", projCoeff=" << projCoeff
        << ", max|gradJ|=" << maxGradJ
        << ", max|gradg|=" << maxGradg
        << ", max|xiJ|=" << maxXiJ
        << ", <xiJ,gradg>_V=" << orthResidual
        << ", Dg(xiJ)=" << nullResidual
        << ", mean(G)=" << meanG
        << ", violation=" << violation
        << Alert::Raise;

      // Range space step (Gauss-Newton constraint correction), formula (2.8):
      //   xiC = (g(x) / <gradg, gradg>_V) * gradg
      const Real rangeCoeff = violation / std::max(innerGG, Real(1e-30));

      GridFunction xiC(vh);
      xiC = gradg;
      xiC *= rangeCoeff;

      GridFunction normXiC(shNormSpace);
      normXiC = Frobenius(xiC);
      const Real maxXiC = normXiC.max();

      Alert::Info()
        << "   | Range space diag: rangeCoeff=" << rangeCoeff
        << ", max|xiC|=" << maxXiC
        << Alert::Raise;
      //   αJ = A_J * hmin / ||ξJ||∞
      //   αC = min(0.9, A_C * hmin / max(1e-9, ||ξC||∞))
      const Real alphaJ = A_J * hmin / std::max(maxXiJ, Real(1e-30));
      const Real alphaC = std::min(Real(0.9), A_C * hmin / std::max(maxXiC, Real(1e-9)));

      //   θ = -(αJ · ξJ + αC · ξC)
      GridFunction thetaFluid(vh);
      thetaFluid = xiJ;
      thetaFluid *= -alphaJ;
      {
        GridFunction tmp(vh);
        tmp = xiC;
        tmp *= -alphaC;
        thetaFluid += tmp;
      }

      // Advection still happens on the full level-set mesh th, so extend the
      // fluid-domain descent field by zero inside the obstacle.
      P1 thetaSpace(th, d);
      GridFunction thetaField(thetaSpace);
      thetaField.project(Region::Cells, VectorFunction{0.0, 0.0, 0.0});
      thetaField.project(
        Region::Cells,
        VectorFunction{
          [&](const Geometry::Point& p) { return thetaFluid.getValue(p)(0); },
          [&](const Geometry::Point& p) { return thetaFluid.getValue(p)(1); },
          [&](const Geometry::Point& p) { return thetaFluid.getValue(p)(2); }
        },
        Fluid);

      Alert::Info()
        << "  alphaJ=" << alphaJ
        << ", alphaC=" << alphaC
        << ", dispJ=" << alphaJ * maxXiJ
        << ", dispC=" << alphaC * maxXiC
        << Alert::Raise;

      // ----------------------------------------------------------------------
      // Signed-distance reconstruction
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Distancing obstacle phase." << Alert::Raise;

      P1 sh(th);
      GridFunction dist(sh);
      Distance::Eikonal(dist)
        .setInterior(Fluid)
        .setInterface(shapeInterface)
        .solve()
        .sign();

      // ----------------------------------------------------------------------
      // Advection
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Advecting level set." << Alert::Raise;

      // Δt = 1 per paper §5.3.2: step size is encoded in αJ, αC.
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
      Advection::Lagrangian(advect, test, dist, thetaField).step(1.0);

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
            .setRMC(rmc)
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

        printMeshQuality("Remeshed domain quality", th);
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
      fVol << obstacleVolume << " " << violation << "\n";
      fNS  << alphaJ << " " << alphaC << " " << projCoeff << " " << rangeCoeff << " " << maxXiJ << " " << maxXiC << "\n";

      fObj.flush();
      fObjRaw.flush();
      fVol.flush();
      fNS.flush();

      Alert::Info()
        << "   | Jraw=" << Jraw
        << ", Jopt=" << J
        << ", Vobs=" << obstacleVolume
        << ", c=" << violation
        << ", mean(G)=" << meanG
        << ", alphaJ=" << alphaJ
        << ", alphaC=" << alphaC
        << ", projCoeff=" << projCoeff
        << ", max|xiJ|=" << maxXiJ
        << ", max|xiC|=" << maxXiC
        << Alert::Raise;

      th.save("out/Omega." + std::to_string(it) + ".mesh", IO::FileFormat::MEDIT);
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
    std::cerr << "LevelSetStokes3DObstacle_NS failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
