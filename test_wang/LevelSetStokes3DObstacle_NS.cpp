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

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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
static constexpr Real   defaultHMax    = 0.08;
static constexpr Real   defaultAlpha   = 0.2;
static constexpr Real   defaultDt      = 0.5 * (defaultHMax - 0.33 * defaultHMax);
static constexpr Real   regularization = 1e-12;

// Null-space gradient flow step-size constants (O(1), insensitive to tuning)
static constexpr Real A_J = 1.0;
static constexpr Real A_C = 1.0;

static std::array<Real, 3> basis(int axis)
{
  if (axis < 0 || axis > 2)
    throw std::runtime_error("Axis must be 0, 1 or 2.");

  std::array<Real, 3> e{0.0, 0.0, 0.0};
  e[static_cast<size_t>(axis)] = 1.0;
  return e;
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const Real hmax0 = defaultHMax;

  Real hmax  = hmax0;

  try
  {
    const char* meshFile = "sphere_init_aligned.mesh";
    const size_t maxIt   = defaultMaxIt;

    // Objective component K_ij, default K_11 in 0-based storage -> axis 1,1 as in original file.
    const int iAxis = 0;
    const int jAxis = 0;

    Real hmin  = 0.33 * hmax;
    Real hausd = 0.1 * hmin;

    const Real k = 0.5;
    const Real dt    = k * (defaultHMax - hmin);
    const Real alpha = 4 * dt;

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
    IO::XDMF xdmf("out/LevelSetStokes3DObstacle_NS");

    auto domainGrid = xdmf.grid("domain");
    domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);

    auto stateGrid = xdmf.grid("state");

    // Histories
    std::ofstream fObj("obj.txt");
    std::ofstream fVol("vol.txt");
    std::ofstream fNS("ns.txt");

    if (!fObj || !fVol || !fNS)
      throw std::runtime_error("Failed to open output history files.");

    // Volume constraint target
    const double targetObstacleVolume = th.getVolume(Obstacle);

    // Objective basis vectors
    const auto ei = basis(iAxis);
    const auto ej = basis(jAxis);

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
        hmin  = 0.1 * hmax;
        hausd = 0.1 * hmin;
      }
      catch (const Alert::Exception&)
      {
        hmax  /= 2;
        hmin   = 0.33 * hmax;
        hausd  = 0.1 * hmin;

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
        + DirichletBC(ua, VectorFunction{ei[0], ei[1], ei[2]}).on(fluidShapeBdr)
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

      // ----------------------------------------------------------------------
      // Shape gradient
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Computing shape gradient." << Alert::Raise;

      auto juShape = Jacobian(u);
      auto jvShape = Jacobian(v);

      auto eu = 0.5 * (juShape + juShape.T());
      auto ev = 0.5 * (jvShape + jvShape.T());

      auto G = 2.0 * mu * Dot(eu, ev);

      const double obstacleVolume = th.getVolume(Obstacle);
      const double violation = obstacleVolume - targetObstacleVolume;

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

      auto n = FaceNormal(th);
      n.traceOf(Obstacle);

      // ----------------------------------------------------------------------
      // Null-space gradient flow: two Hilbert extensions + projection
      // ----------------------------------------------------------------------
      Alert::Info() << "   | Solving Hilbert extensions (objective + volume)." << Alert::Raise;

      P1 vh(th, d);

      // Extension 1: objective gradient θ_J
      PETSc::Variational::TrialFunction thetaJ(vh); thetaJ.setName("thetaJ");
      PETSc::Variational::TestFunction psiJ(vh);

      Problem hilbertJ(thetaJ, psiJ);
      hilbertJ =
          Integral(alpha * alpha * Jacobian(thetaJ), Jacobian(psiJ))
        + Integral(thetaJ, psiJ)
        + FaceIntegral(G * Dot(n, psiJ)).over(shapeInterface)
        + DirichletBC(thetaJ, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP(hilbertJ).solve();

      // Extension 2: volume gradient θ_g (constant 1.0 on shape boundary)
      PETSc::Variational::TrialFunction thetaG(vh); thetaG.setName("thetaG");
      PETSc::Variational::TestFunction psiG(vh);

      Problem hilbertG(thetaG, psiG);
      hilbertG =
          Integral(alpha * alpha * Jacobian(thetaG), Jacobian(psiG))
        + Integral(thetaG, psiG)
        + FaceIntegral(1.0 * Dot(n, psiG)).over(shapeInterface)
        + DirichletBC(thetaG, VectorFunction{0.0, 0.0, 0.0}).on(baseOuterBdr);

      Solver::KSP(hilbertG).solve();

      // Copy PETSc solutions into plain GridFunctions for arithmetic
      GridFunction gfJ(vh);
      gfJ = thetaJ.getSolution();
      GridFunction gfG(vh);
      gfG = thetaG.getSolution();

      // Compute inner products via Eigen data vectors
      // (nodal L² inner product as proxy for H¹ metric encoded by the Hilbert extension)
      const Real innerJG = gfJ.getData().dot(gfG.getData());
      const Real innerGG = gfG.getData().dot(gfG.getData());

      const Real projCoeff  = innerJG / std::max(innerGG, 1e-30);
      const Real rangeCoeff = violation / std::max(innerGG, 1e-30);

      // ξ_J = θ_J - projCoeff · θ_g  (null-space: objective descent without affecting volume)
      GridFunction xiJ(vh);
      xiJ = gfJ;
      {
        GridFunction tmp(vh);
        tmp = gfG;
        tmp *= projCoeff;
        xiJ -= tmp;
      }

      // ξ_C = rangeCoeff · θ_g  (range-space: Gauss-Newton constraint correction)
      GridFunction xiC(vh);
      xiC = gfG;
      xiC *= rangeCoeff;

      // Adaptive normalization (step size proportional to hmin)
      P1 shNormSpace(th);
      GridFunction normXiJ(shNormSpace);
      normXiJ = Frobenius(xiJ);
      GridFunction normXiC(shNormSpace);
      normXiC = Frobenius(xiC);

      const Real maxXiJ = normXiJ.max();
      const Real maxXiC = normXiC.max();

      // Paper (5.13-5.14): α normalizes so displacement = α·||ξ||∞ ≈ hmin.
      // Cap αC at 0.9 (paper default) so constraint step ≤ 0.9·||ξC||∞.
      const Real alphaJ = A_J * hmin / std::max(maxXiJ, 1e-30);
      const Real alphaC = std::min(Real(0.9), A_C * hmin / std::max(maxXiC, 1e-30));

      Alert::Info()
        << "   | NS diag: ||xiJ||=" << maxXiJ
        << ", ||xiC||=" << maxXiC
        << ", alphaJ=" << alphaJ
        << ", alphaC=" << alphaC
        << ", dispJ=" << alphaJ * maxXiJ
        << ", dispC=" << alphaC * maxXiC
        << Alert::Raise;

      // Combined descent direction: θ = -α_J · ξ_J - α_C · ξ_C
      GridFunction theta(vh);
      theta = xiJ;
      theta *= -alphaJ;
      {
        GridFunction tmp(vh);
        tmp = xiC;
        tmp *= -alphaC;
        theta += tmp;
      }

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

      // The adaptive normalization (alphaJ, alphaC) already targets
      // ||theta||∞ ≈ hmin, so we advect for unit time (displacement = 1 * theta).
      // No additional vmax normalization needed.
      TrialFunction advect(sh);
      TestFunction test(sh);
      Advection::Lagrangian(advect, test, dist, theta).step(1.0);

      // ----------------------------------------------------------------------
      // XDMF snapshot before remeshing
      // ----------------------------------------------------------------------
      domainGrid.clear();
      domainGrid.setMesh(th, IO::XDMF::MeshPolicy::Transient);
      domainGrid.add("dist",   dist,   IO::XDMF::Center::Node);
      domainGrid.add("theta",  theta,  IO::XDMF::Center::Node);
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
            // .setHausdorff(hausd)
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

      // ----------------------------------------------------------------------
      // Logging
      // ----------------------------------------------------------------------
      fObj << J << "\n";
      fVol << obstacleVolume << " " << violation << "\n";
      fNS  << alphaJ << " " << alphaC << " " << projCoeff << " " << rangeCoeff << "\n";

      fObj.flush();
      fVol.flush();
      fNS.flush();

      Alert::Info()
        << "   | J=" << J
        << ", Vobs=" << obstacleVolume
        << ", c=" << violation
        << ", mean(G)=" << meanG
        << ", alphaJ=" << alphaJ
        << ", alphaC=" << alphaC
        << ", proj=" << projCoeff
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
