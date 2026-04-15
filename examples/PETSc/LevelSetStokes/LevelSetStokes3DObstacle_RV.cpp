/*
 * Stokes shape-optimization with REDUCED VOLUME constraint.
 *
 * Based on LevelSetStokes3DObstacle.cpp, with the key improvement:
 *   volume-only constraint  -->  reduced volume  nu = 6*sqrt(pi)*V / A^(3/2)
 *
 * This prevents degenerate (infinitely flat / infinitely slender) optimal
 * shapes that arise when only volume is constrained (see Paper 1, Sec 5.1).
 *
 * Semantic choice:
 *   - obstacle cells:  attribute 2
 *   - fluid cells:     attribute 3
 *   - shape interface: attribute 13
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
static constexpr size_t defaultMaxIt   = 60;
static constexpr Real   defaultHMax    = 0.2;
static constexpr Real   defaultAlpha   = 0.2;
static constexpr Real   defaultDt      = 0.5 * (defaultHMax - 0.1 * defaultHMax);
static constexpr Real   regularization = 1e-12;

static std::array<Real, 3> basis(int axis)
{
  if (axis < 0 || axis > 2)
    throw std::runtime_error("Axis must be 0, 1 or 2.");

  std::array<Real, 3> e{0.0, 0.0, 0.0};
  e[static_cast<size_t>(axis)] = 1.0;
  return e;
}

/// Compute the reduced volume: nu = 6*sqrt(pi)*V / A^(3/2)
static double reducedVolume(double V, double A)
{
  return 6.0 * std::sqrt(M_PI) * V / std::pow(A, 1.5);
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

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

    // Objective component K_ij, default K_11
    const int iAxis = LevelSetStokes::Runtime::envInt("IAXIS", 1);
    const int jAxis = LevelSetStokes::Runtime::envInt("JAXIS", 1);

    Real hmin  = hminRatio * hmax;
    Real hausd = hausdRatio * hmin;

    const Real k = LevelSetStokes::Runtime::envDouble("STEP_K", 0.1);
    const Real dt    = k * (hmax0 - hmin);
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
    IO::XDMF xdmf("out/LevelSetStokes3DObstacle_RV");

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj_rv.txt");
    std::ofstream fVol("vol_rv.txt");
    std::ofstream fAL("al_rv.txt");

    if (!fObj || !fVol || !fAL)
      throw std::runtime_error("Failed to open output history files.");

    // ------------------------------------------------------------------
    // Compute initial surface area and reduced volume target
    // ------------------------------------------------------------------
    const double initObstacleVolume = th.getVolume(Obstacle);

    // Compute initial GammaShape surface area via face integral
    {
      auto& conn = th.getConnectivity();
      conn.compute(2, 3);
      conn.compute(3, 2);
      conn.compute(2, 1);
      conn.compute(1, 0);
      conn.compute(0, 0);
    }

    P0g p0Init(th);
    TestFunction z0Init(p0Init);
    LinearForm lfAreaInit(z0Init);
    lfAreaInit = FaceIntegral(z0Init).over({GammaShape});
    lfAreaInit.assemble();

    GridFunction oneInit(p0Init);
    oneInit = 1.0;
    const double initShapeArea = lfAreaInit(oneInit);

    const double targetReducedVol = reducedVolume(initObstacleVolume, initShapeArea);

    Alert::Info()
      << "Initial: V=" << initObstacleVolume
      << ", A=" << initShapeArea
      << ", nu=" << targetReducedVol
      << Alert::Raise;

    // Augmented Lagrangian parameters — TWO constraints:
    //   Constraint 1 (volume):         C_V  = V - V0 = 0
    //   Constraint 2 (reduced volume): C_nu = nu - nu0 = 0
    // Volume prevents uniform shrinkage; reduced volume prevents degenerate shapes.
    double lambdaV  = 0.0;      // Lagrange multiplier for volume
    double penaltyV = LevelSetStokes::Runtime::envDouble("PENALTY_V", 5.0e3);   // penalty for volume
    double lambdaNu  = 0.0;    // Lagrange multiplier for reduced volume
    double penaltyNu = LevelSetStokes::Runtime::envDouble("PENALTY_NU", 5.0e3);  // penalty for reduced volume (nu violation ~0.02 => shift~100)
    double prevVolViolation = 0.0;
    double prevNuViolation  = 0.0;
    const double targetObstacleVolume = initObstacleVolume;

    // Objective basis vectors
    const auto ei = basis(iAxis);
    const auto ej = basis(jAxis);

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

      // ------------------------------------------------------------------
      // Optimize current mesh
      // ------------------------------------------------------------------
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

      // ------------------------------------------------------------------
      // Connectivity
      // ------------------------------------------------------------------
      {
        auto& conn = th.getConnectivity();
        conn.compute(2, 3);
        conn.compute(3, 2);
        conn.compute(2, 1);
        conn.compute(1, 0);
        conn.compute(0, 0);
      }

      // ------------------------------------------------------------------
      // Trim obstacle to get fluid mesh
      // ------------------------------------------------------------------
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

      // ------------------------------------------------------------------
      // State Stokes solve
      // ------------------------------------------------------------------
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
        + DirichletBC(up, VectorFunction{ej[0], ej[1], ej[2]}).on(fluidShapeBdr)
        + DirichletBC(up, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Alert::Info() << "   | Assembling Stokes." << Alert::Raise;
      state.assemble();

      Alert::Info() << "   | Solving Stokes." << Alert::Raise;

      Solver::KSP(state).solve();

      auto u = GridFunction(velocitySpace);
      auto p = GridFunction(pressureSpace);
      u = up.getSolution();
      p = pp.getSolution();

      // ------------------------------------------------------------------
      // Adjoint Stokes solve
      // ------------------------------------------------------------------
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
        + DirichletBC(ua, VectorFunction{ei[0], ei[1], ei[2]}).on(fluidShapeBdr)
        + DirichletBC(ua, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP(adj).solve();

      auto v = GridFunction(velocitySpace);
      auto q = GridFunction(pressureSpace);
      v = ua.getSolution();
      q = pa.getSolution();

      // ------------------------------------------------------------------
      // Objective evaluation
      // ------------------------------------------------------------------
      Alert::Info() << "   | Evaluating objective." << Alert::Raise;

      auto ju = Jacobian(u);

      auto nObj = FaceNormal(th);
      nObj.traceOf(Obstacle);

      const auto eiVec = VectorFunction{ei[0], ei[1], ei[2]};
      auto sigma = mu * (ju + ju.T()) - p * IdentityMatrix(d);
      auto traction = sigma * nObj;

      P0g p0Obj(th);
      TestFunction z0Obj(p0Obj);
      LinearForm lfObj(z0Obj);
      lfObj = FaceIntegral(Dot(traction, eiVec), z0Obj).over(shapeInterface);
      lfObj.assemble();

      GridFunction oneObj(p0Obj);
      oneObj = 1.0;

      const double J = -lfObj(oneObj);

      Alert::Info() << "   | Objective: " << J << Alert::Raise;

      // ------------------------------------------------------------------
      // Signed-distance reconstruction
      // ------------------------------------------------------------------
      Alert::Info() << "   | Distancing obstacle phase." << Alert::Raise;

      P1 sh(th);
      GridFunction dist(sh);
      Distance::Eikonal(dist)
        .setInterior(Fluid)
        .setInterface(shapeInterface)
        .solve()
        .sign();

      // ------------------------------------------------------------------
      // Curvature computation via weak Laplacian of signed distance
      //
      //   For a signed distance function (|grad(phi)|=1):
      //     kappa = div(grad(phi)) = Laplacian(phi)
      //
      //   Weak projection:  find kappa in P1 s.t.
      //     int kappa * w dx = - int grad(dist) . grad(w) dx
      //
      //   i.e. kappa = -Laplacian(dist) (with sign convention:
      //   positive kappa = convex obstacle surface seen from fluid)
      // ------------------------------------------------------------------
      Alert::Info() << "   | Computing curvature from signed distance." << Alert::Raise;

      P1 kappaSpace(th);
      TrialFunction kappaT(kappaSpace);
      TestFunction  kappaW(kappaSpace);

      Problem curvProj(kappaT, kappaW);
      curvProj =
          Integral(kappaT, kappaW)
        + Integral(Grad(dist), Grad(kappaW));

      curvProj.assemble();
      Solver::CG(curvProj).solve();

      auto kappa = GridFunction(kappaSpace);
      kappa = kappaT.getSolution();

      // ------------------------------------------------------------------
      // Shape gradient + reduced volume constraint
      // ------------------------------------------------------------------
      Alert::Info() << "   | Computing shape gradient." << Alert::Raise;

      auto juShape = Jacobian(u);
      auto jvShape = Jacobian(v);

      auto eu = 0.5 * (juShape + juShape.T());
      auto ev = 0.5 * (jvShape + jvShape.T());

      auto G = 2.0 * mu * Dot(eu, ev);

      // Compute current obstacle volume and surface area
      const double obstacleVolume = th.getVolume(Obstacle);

      P0g p0Stats(th);
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

      // --- Constraint 1: volume ---
      const double volViolation = obstacleVolume - targetObstacleVolume;
      // Shape derivative of volume: dC_V/d(theta_n) = 1  (constant on interface)
      const double volShift = -lambdaV + penaltyV * volViolation;

      // --- Constraint 2: reduced volume ---
      const double nu = reducedVolume(obstacleVolume, shapeArea);
      const double nuViolation = nu - targetReducedVol;
      // Shape derivative of reduced volume:
      //   dC_nu/d(theta_n) = nu * (1/V  -  3*kappa / (2*A))
      // Split into uniform part (from volume derivative) and curvature part (from area derivative)
      const double nuCoeff = -lambdaNu + penaltyNu * nuViolation;
      const double nuUniformPart   = nuCoeff * nu / obstacleVolume;
      const double nuCurvaturePart = -1.5 * nuCoeff * nu / shapeArea;

      // Combined shift:  volShift * 1  +  nuUniformPart  +  nuCurvaturePart * kappa(x)
      const double uniformPart  = volShift + nuUniformPart;
      const double curvaturePart = nuCurvaturePart;

      auto n = FaceNormal(th);
      n.traceOf(Obstacle);

      // ------------------------------------------------------------------
      // Hilbert extension
      //
      // The RHS on GammaShape is:
      //   (G  +  uniformPart  +  curvaturePart * kappa) * (n . psi)
      //
      // where uniformPart and curvaturePart combine TWO constraints:
      //   uniformPart  = volShift                (from dC_V/d(theta_n) = 1)
      //                + nuCoeff * nu / V        (from volume part of dC_nu)
      //   curvaturePart * kappa                  (from area part of dC_nu)
      // ------------------------------------------------------------------
      Alert::Info() << "   | Solving vector Hilbert extension." << Alert::Raise;

      P1 vh(th, d);
      PETSc::Variational::TrialFunction theta(vh); theta.setName("theta");
      PETSc::Variational::TestFunction psi(vh);

      Problem hilbert(theta, psi);
      hilbert =
          Integral(alpha * alpha * Jacobian(theta), Jacobian(psi))
        + Integral(theta, psi)
        + FaceIntegral((G + uniformPart + curvaturePart * kappa) * Dot(n, psi))
            .over(shapeInterface)
        + DirichletBC(theta, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP(hilbert).solve();

      // ------------------------------------------------------------------
      // Advection
      // ------------------------------------------------------------------
      Alert::Info() << "   | Advecting level set." << Alert::Raise;

      auto& V = theta.getSolution();

      P1 shNorm(th);
      GridFunction speed(shNorm);
      speed = Frobenius(V);

      // Adaptive step: bound max displacement to dt * 1.0 (i.e. normalize so max speed <= 1)
      // and further clamp so displacement <= hmin to avoid mesh inversion
      const Real vmax = speed.max();
      Real effectiveDt = dt;
      if (vmax > 0.0)
      {
        const Real maxDisplacement = vmax * dt;
        if (maxDisplacement > hmin)
          effectiveDt = hmin / vmax;
      }
      if (vmax > 1.0)
        V /= vmax;

      Alert::Info() << "   | vmax=" << vmax << ", effectiveDt=" << effectiveDt << Alert::Raise;

      TrialFunction advect(sh);
      TestFunction test(sh);
      Advection::Lagrangian(advect, test, dist, V).step(effectiveDt);

      // ------------------------------------------------------------------
      // XDMF snapshot before remeshing
      // ------------------------------------------------------------------
      domainGrid.clear();
      domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);
      domainGrid.add("dist",     dist,                 IO::XDMF::Center::Node);
      domainGrid.add("theta",    theta.getSolution(),  IO::XDMF::Center::Node);
      domainGrid.add("speed",    speed,                IO::XDMF::Center::Node);
      domainGrid.add("advect",   advect.getSolution(), IO::XDMF::Center::Node);
      domainGrid.add("curvature", kappa,               IO::XDMF::Center::Node);

      stateGrid.clear();
      stateGrid.setMesh(fluidMesh, IO::XDMF::MeshPolicy::Transient);
      stateGrid.add("u", u, IO::XDMF::Center::Node);
      stateGrid.add("p", p, IO::XDMF::Center::Node);
      stateGrid.add("v", v, IO::XDMF::Center::Node);
      stateGrid.add("q", q, IO::XDMF::Center::Node);

      xdmf.write(it).flush();

      // ------------------------------------------------------------------
      // Remeshing
      // ------------------------------------------------------------------
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

      // ------------------------------------------------------------------
      // Augmented Lagrangian update (dual constraints)
      // ------------------------------------------------------------------
      const double Jaug = J
        - lambdaV  * volViolation + 0.5 * penaltyV  * volViolation * volViolation
        - lambdaNu * nuViolation  + 0.5 * penaltyNu * nuViolation  * nuViolation;

      fObj << J << " " << Jaug << "\n";
      fVol << obstacleVolume << " " << shapeArea << " " << nu
           << " " << volViolation << " " << nuViolation << "\n";
      fAL  << lambdaV << " " << penaltyV << " " << lambdaNu << " " << penaltyNu << "\n";

      fObj.flush();
      fVol.flush();
      fAL.flush();

      // Update volume constraint multipliers
      lambdaV += -penaltyV * volViolation;
      // Only increase penalty if violation did not improve sufficiently
      if (std::abs(volViolation) > 0.01 * targetObstacleVolume &&
          std::abs(volViolation) > 0.5 * std::abs(prevVolViolation))
        penaltyV = std::min(1.0e6, 1.5 * penaltyV);
      prevVolViolation = volViolation;

      // Update reduced volume constraint multipliers
      lambdaNu += -penaltyNu * nuViolation;
      if (std::abs(nuViolation) > 0.01 * targetReducedVol &&
          std::abs(nuViolation) > 0.5 * std::abs(prevNuViolation))
        penaltyNu = std::min(1.0e6, 1.5 * penaltyNu);
      prevNuViolation = nuViolation;

      Alert::Info()
        << "   | J=" << J
        << ", Jaug=" << Jaug
        << ", Vobs=" << obstacleVolume
        << ", V0=" << targetObstacleVolume
        << ", cV=" << volViolation
        << ", A=" << shapeArea
        << ", nu=" << nu
        << ", nu0=" << targetReducedVol
        << ", cNu=" << nuViolation
        << ", mean(G)=" << meanG
        << ", lamV=" << lambdaV
        << ", bV=" << penaltyV
        << ", lamNu=" << lambdaNu
        << ", bNu=" << penaltyNu
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
    std::cerr << "LevelSetStokes3DObstacle_RV failed: " << e.what() << '\n';
    PetscFinalize();
    return 1;
  }
}
