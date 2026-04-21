/*
 * Paper-style direct-mesh deformation test for LevelSetStokes3DObstacle.
 *
 * This target keeps the FEM state/adjoint solves from the Rodin examples but
 * replaces the level-set update with a direct vertex displacement:
 *   1. compute the shape gradient G on the current interface,
 *   2. build the augmented-Lagrangian scalar phi = G - ell + b * c,
 *   3. solve a vector Laplace problem for theta in the fluid domain,
 *   4. move mesh vertices by x <- x - tau * theta(x),
 *   5. backtrack tau until the augmented Lagrangian decreases.
 *
 * The goal is to approximate the update mechanism described in Moreau,
 * Ishimoto, Privat, while still using the existing Rodin FEM pipeline.
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
static constexpr size_t defaultMaxIt        = 20;
static constexpr Real   defaultHMax         = 0.2;
static constexpr Real   regularization      = 1e-12;
static constexpr Real   defaultTau0         = 1e-3;
static constexpr Real   defaultTauMin       = 1e-6;
static constexpr size_t defaultMaxBacktrack = 8;
static constexpr Real   defaultEll0         = 0.0;
static constexpr Real   defaultB0           = 1.0;
static constexpr Real   defaultBTarget      = 50.0;
static constexpr Real   defaultALGrowth     = 1.03;

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
        "Unsupported paper-style objective mode. Use OBJECTIVE_MODE=K, OBJECTIVE_MODE=C or OBJECTIVE_MODE=Q.");
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
      parseObjectiveMode(LevelSetStokes::Runtime::envCString("OBJECTIVE_MODE", "C"));
    const ObjectiveSense objectiveSense =
      parseObjectiveSense(LevelSetStokes::Runtime::envCString("OBJECTIVE_SENSE", "min"));

    // Objective component K_ij, default K_11 in 0-based storage -> axis 1,1 as in original file.
    const int iAxis = LevelSetStokes::Runtime::envInt("IAXIS", 1);
    const int jAxis = LevelSetStokes::Runtime::envInt("JAXIS", 1);

    Real hmin  = hminRatio * hmax;
    Real hausd = hausdRatio * hmin;

    const Real tau0 = LevelSetStokes::Runtime::envDouble("TAU0", defaultTau0);
    const Real tauMin = LevelSetStokes::Runtime::envDouble("TAU_MIN", defaultTauMin);
    const size_t maxBacktrack =
      LevelSetStokes::Runtime::envSizeT("MAX_BACKTRACK", defaultMaxBacktrack);
    Real ell = LevelSetStokes::Runtime::envDouble("ELL0", defaultEll0);
    Real b = LevelSetStokes::Runtime::envDouble("B0", defaultB0);
    const Real bTarget = LevelSetStokes::Runtime::envDouble("B_TARGET", defaultBTarget);
    const Real alGrowth = LevelSetStokes::Runtime::envDouble("AL_GROWTH", defaultALGrowth);

    Alert::Info()
      << "Configuration: HMAX=" << hmax
      << ", HMIN_RATIO=" << hminRatio
      << ", HAUSD_RATIO=" << hausdRatio
      << ", RMC=" << rmc
      << ", TAU0=" << tau0
      << ", ELL0=" << ell
      << ", B0=" << b
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
    IO::XDMF xdmf("out/LevelSetStokes3DObstacle_paperC11");

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fObjRaw("obj_raw.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fPaper("paper.txt");

    if (!fObj || !fObjRaw || !fVol || !fPaper)
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

    const size_t d = th.getSpaceDimension();

    struct ObjectiveEval
    {
      double Jraw = 0.0;
      double Jopt = 0.0;
      double volume = 0.0;
      double violation = 0.0;
      double augmented = 0.0;
      bool finite = true;
    };

    const auto evaluateCandidate =
      [&](const MeshType& mesh, double ellValue, double bValue, bool verbose = false) -> ObjectiveEval
      {
        ObjectiveEval out;
        auto fluid = mesh.trim(Obstacle);
        {
          auto& conn = fluid.getConnectivity();
          conn.compute(2, 3);
          conn.compute(3, 2);
          conn.compute(2, 1);
          conn.compute(1, 0);
          conn.compute(0, 0);
        }

        auto velocitySpaceEval = H1(std::integral_constant<size_t, 2>{}, fluid, d);
        auto pressureSpaceEval = H1(std::integral_constant<size_t, 1>{}, fluid);
        P0g globalP0Eval(fluid);

        FlatSet<Attribute> fluidShapeBdrEval{GammaShape};
        FlatSet<Attribute> shapeInterfaceEval{GammaShape};

        PETSc::Variational::TrialFunction upEval(velocitySpaceEval);
        PETSc::Variational::TrialFunction ppEval(pressureSpaceEval);
        PETSc::Variational::TrialFunction lmbEval(globalP0Eval);
        PETSc::Variational::TestFunction vpEval(velocitySpaceEval);
        PETSc::Variational::TestFunction qpEval(pressureSpaceEval);
        PETSc::Variational::TestFunction mu0Eval(globalP0Eval);

        Problem stateEval(upEval, ppEval, lmbEval, vpEval, qpEval, mu0Eval);
        stateEval =
            Integral(mu * Jacobian(upEval), Jacobian(vpEval))
          - Integral(ppEval, Div(vpEval))
          + Integral(Div(upEval), qpEval)
          + Integral(lmbEval, qpEval)
          + Integral(ppEval, mu0Eval)
          + regularization * Integral(ppEval, qpEval)
          + regularization * Integral(lmbEval, mu0Eval)
          + DirichletBC(upEval, stateVelocity).on(fluidShapeBdrEval)
          + DirichletBC(upEval, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

        Solver::KSP stateEvalSolver(stateEval);
        stateEvalSolver.solve();
        if (verbose)
          printKSPSummary("Candidate State KSP", stateEvalSolver.getHandle());

        auto uEval = GridFunction(velocitySpaceEval);
        auto pEval = GridFunction(pressureSpaceEval);
        uEval = upEval.getSolution();
        pEval = ppEval.getSolution();

        PETSc::Variational::TrialFunction uaEval(velocitySpaceEval);
        PETSc::Variational::TrialFunction paEval(pressureSpaceEval);
        PETSc::Variational::TrialFunction laEval(globalP0Eval);
        PETSc::Variational::TestFunction vaEval(velocitySpaceEval);
        PETSc::Variational::TestFunction qaEval(pressureSpaceEval);
        PETSc::Variational::TestFunction maEval(globalP0Eval);

        Problem adjEval(uaEval, paEval, laEval, vaEval, qaEval, maEval);
        adjEval =
            Integral(mu * Jacobian(uaEval), Jacobian(vaEval))
          - Integral(paEval, Div(vaEval))
          + Integral(Div(uaEval), qaEval)
          + Integral(laEval, qaEval)
          + Integral(paEval, maEval)
          + regularization * Integral(paEval, qaEval)
          + regularization * Integral(laEval, maEval)
          + DirichletBC(uaEval, adjointVelocity).on(fluidShapeBdrEval)
          + DirichletBC(uaEval, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

        Solver::KSP adjEvalSolver(adjEval);
        adjEvalSolver.solve();
        if (verbose)
          printKSPSummary("Candidate Adjoint KSP", adjEvalSolver.getHandle());

        auto vEval = GridFunction(velocitySpaceEval);
        vEval = uaEval.getSolution();

        auto juEval = Jacobian(uEval);
        auto nObjEval = -FaceNormal(fluid);
        auto sigmaEval = mu * (juEval + juEval.T()) - pEval * IdentityMatrix(d);
        auto tractionEval = sigmaEval * nObjEval;

        P0g p0ObjEval(fluid);
        TestFunction z0ObjEval(p0ObjEval);
        LinearForm lfObjEval(z0ObjEval);
        lfObjEval = FaceIntegral(Dot(tractionEval, adjointVelocity), z0ObjEval).over(shapeInterfaceEval);
        lfObjEval.assemble();

        GridFunction oneObjEval(p0ObjEval);
        oneObjEval = 1.0;

        out.Jraw = -lfObjEval(oneObjEval);
        const double senseSignEval =
          objectiveSense == ObjectiveSense::Max ? -1.0 : 1.0;
        out.Jopt = senseSignEval * out.Jraw;
        out.volume = mesh.getVolume(Obstacle);
        out.violation = out.volume - targetObstacleVolume;
        out.augmented = out.Jopt - ellValue * out.violation
                      + 0.5 * bValue * out.violation * out.violation;
        out.finite = std::isfinite(out.Jraw) && std::isfinite(out.Jopt)
                  && std::isfinite(out.augmented);
        return out;
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
      printMeshQuality("Domain mesh quality", th);
      printMeshQuality("Fluid mesh quality", fluidMesh);

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
      // Shape gradient and paper-style descent direction
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
      const double currentAugmented =
        J - ell * violation + 0.5 * b * violation * violation;

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

      auto n = -FaceNormal(fluidMesh);
      auto phi = G - ell + b * violation;

      Alert::Info() << "   | Solving paper-style Laplace descent field." << Alert::Raise;

      P1 vh(fluidMesh, d);
      PETSc::Variational::TrialFunction thetaTrial(vh); thetaTrial.setName("theta");
      PETSc::Variational::TestFunction thetaTest(vh);
      Problem descent(thetaTrial, thetaTest);
      descent =
          Integral(Jacobian(thetaTrial), Jacobian(thetaTest))
        + regularization * Integral(thetaTrial, thetaTest)
        - FaceIntegral(phi * Dot(n, thetaTest)).over(shapeInterface)
        + DirichletBC(thetaTrial, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP descentSolver(descent);
      descentSolver.solve();
      printKSPSummary("Descent KSP", descentSolver.getHandle());

      GridFunction thetaFluid(vh);
      thetaFluid = thetaTrial.getSolution();
      printArrayDiagnostics("thetaFluid", thetaFluid.getData());

      P1 shNormSpace(fluidMesh);
      GridFunction normTheta(shNormSpace);
      normTheta = Frobenius(thetaFluid);
      const Real maxTheta = normTheta.max();

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
        << "   | mean(G)=" << meanG
        << ", phiShift=" << (-ell + b * violation)
        << ", ||theta||∞=" << maxTheta
        << ", Laug=" << currentAugmented
        << Alert::Raise;

      // ----------------------------------------------------------------------
      // Snapshot of current iterate and descent field
      // ----------------------------------------------------------------------
      domainGrid.clear();
      domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);
      domainGrid.add("theta", thetaField, IO::XDMF::Center::Node);

      stateGrid.clear();
      stateGrid.setMesh(fluidMesh, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", u, IO::XDMF::Center::Node);
      stateGrid.add("p", p, IO::XDMF::Center::Node);
      stateGrid.add("v", v, IO::XDMF::Center::Node);
      stateGrid.add("q", q, IO::XDMF::Center::Node);

      xdmf.write(it).flush();

      // ----------------------------------------------------------------------
      // Backtracking line search with direct mesh displacement
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Backtracking direct mesh displacement." << Alert::Raise;

      bool accepted = false;
      Real tau = tau0;
      MeshType acceptedMesh;
      ObjectiveEval acceptedEval;
      const auto& thetaData = thetaField.getData();
      const size_t vertexCount = th.getVertexCount();
      const size_t scalarSize = vertexCount;
      if (thetaData.size() != static_cast<Eigen::Index>(scalarSize * d))
      {
        throw std::runtime_error("Unexpected thetaField data size for P1 vector field.");
      }

      for (size_t bt = 0; bt < maxBacktrack && tau >= tauMin; ++bt)
      {
        MeshType candidate = th;
        for (Index vid = 0; vid < static_cast<Index>(vertexCount); ++vid)
        {
          auto coords = th.getVertexCoordinates(vid);
          for (size_t c = 0; c < d; ++c)
            coords(c) -= tau * thetaData[vid + static_cast<Index>(c * scalarSize)];
          candidate.setVertexCoordinates(vid, coords);
        }
        candidate.flush();

        try
        {
          MMG::Optimizer()
            .setHMax(hmax)
            .setHMin(hmin)
            .setHausdorff(hausd)
            .setAngleDetection(false)
            .optimize(candidate);
          auto& candidateConn = candidate.getConnectivity();
          candidateConn.compute(2, 3);
          candidateConn.compute(3, 2);
          candidateConn.compute(2, 1);
          candidateConn.compute(1, 0);
          candidateConn.compute(0, 0);
        }
        catch (const Alert::Exception&)
        {
          Alert::Warning()
            << "   | Candidate mesh invalid for tau=" << tau
            << ", retrying with smaller step."
            << Alert::Raise;
          tau *= 0.5;
          continue;
        }

        const auto cand = evaluateCandidate(candidate, ell, b, false);
        Alert::Info()
          << "   | Candidate tau=" << tau
          << ", Jraw=" << cand.Jraw
          << ", Jopt=" << cand.Jopt
          << ", Laug=" << cand.augmented
          << Alert::Raise;

        if (cand.finite && cand.augmented < currentAugmented)
        {
          accepted = true;
          acceptedMesh = std::move(candidate);
          acceptedEval = cand;
          break;
        }

        tau *= 0.5;
      }

      if (!accepted)
      {
        Alert::Warning()
          << "Failed to find a decreasing direct-displacement step at iteration " << it
          << ". Stopping."
          << Alert::Raise;
        break;
      }

      th = std::move(acceptedMesh);
      ell = ell - b * acceptedEval.violation;
      if (b < bTarget)
        b = std::min(bTarget, alGrowth * b);

      fObj << acceptedEval.Jopt << "\n";
      fObjRaw << acceptedEval.Jraw << "\n";
      fVol << acceptedEval.volume << " " << acceptedEval.violation << "\n";
      fPaper << tau << " " << ell << " " << b << " "
             << maxTheta << " " << currentAugmented << " "
             << acceptedEval.augmented << "\n";

      fObj.flush();
      fObjRaw.flush();
      fVol.flush();
      fPaper.flush();

      Alert::Info()
        << "   | Accepted tau=" << tau
        << ", next Jraw=" << acceptedEval.Jraw
        << ", next Jopt=" << acceptedEval.Jopt
        << ", next Vobs=" << acceptedEval.volume
        << ", next c=" << acceptedEval.violation
        << ", ell=" << ell
        << ", b=" << b
        << Alert::Raise;

      printMeshQuality("Accepted domain quality", th);
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
    std::cerr << "test_PaperStyleC11 failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
