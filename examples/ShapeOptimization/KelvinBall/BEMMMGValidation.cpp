#include "BEM.h"
#include "Sphere.h"
#include "NodalTransport.h"
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/MMG.h>
#include <Rodin/MMG/MeshPrinter.h>
#include <Rodin/Solver/SparseLU.h>
#ifdef RODIN_USE_MUMPS
#include <Rodin/Solver/MUMPS.h>
#endif
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace KelvinBall;

namespace
{
  void save(const MMG::Mesh& mesh, const char* path)
  {
    std::ofstream out(path); out << std::setprecision(17);
    MMG::MeshPrinter(mesh).print(out);
    if (!out) throw std::runtime_error("Cannot save full-domain MMG evidence.");
  }

  void prepare(MMG::Mesh& mesh)
  {
    mesh.getConnectivity().compute(2,3);
    for (auto cell=mesh.getCell(); cell; ++cell) mesh.setAttribute({3,cell->getIndex()},Fluid);
    for (auto face=mesh.getFace(); face; ++face)
      if (face->getAttribute()!=Outer) mesh.setAttribute({2,face->getIndex()},{});
    mesh.getRequiredTriangles().clear();
    for (auto face=mesh.getBoundary(); face; ++face)
    {
      mesh.setAttribute({2,face->getIndex()},Outer);
      mesh.setRequiredTriangle(face->getIndex());
      for (Index vertex:face->getVertices()) mesh.setRequiredVertex(vertex);
    }
  }

  BEM::Result evaluate(MMG::Mesh& mesh, const VelocitySpace& space, Math::Vector<Real>& load,
                       bool derivative)
  {
    mesh.getConnectivity().compute(2,3);
    auto solid=mesh.keep(Obstacle); solid.getConnectivity().compute(2,3);
    const auto skin=solid.skin();
    auto result=BEM().evaluate(skin,derivative);
    load.setZero(space.getSize());
    if (derivative)
      for (Index vertex=0; vertex<skin.getVertexCount(); ++vertex)
      {
        const Index parent=solid.getPolytopeMap(0).left[skin.getPolytopeMap(0).left[vertex]];
        const auto& dofs=space.getDOFs(0,parent);
        for (size_t j=0; j<3; ++j) load(dofs(j))+=result.gradient(vertex,j);
      }
    return result;
  }
}

int main()
{
  try
  {
    std::cout << std::setprecision(17);
    constexpr size_t points=13;
    constexpr Real h=Real(4)/(points-1), dt=.01*h;
    auto cube=KelvinBall::Mesh::UniformGrid(Polytope::Type::Tetrahedron,{points,points,points});
    cube.scale(h);
    for (Index i=0; i<cube.getVertexCount(); ++i)
      cube.setVertexCoordinates(i,cube.getVertexCoordinates(i)-Math::SpatialPoint{2,2,2});
    MMG::Mesh background(std::move(cube)); prepare(background);
    Configuration configuration; configuration.seedAmplitude=.02;
    const Sphere seed(configuration);
    P1 backgroundSpace(background);
    GridFunction initialPhi(backgroundSpace);
    initialPhi=RealFunction([&seed](const Point& p) { return seed.levelSet(p.vector()); });
    MMG::LevelSetDiscretizer cutter;
    cutter.split(Fluid,{Obstacle,Fluid}).setBoundaryReference(Gamma)
      .setHMin(.8*h).setHMax(1.2*h).setHausdorff(.1*h*h)
      .setGradation(2).setAngleDetection(false);
    MMG::Mesh mesh=cutter.discretize(initialPhi);
    auto& connectivity=mesh.getConnectivity();
    connectivity.compute(3,2); connectivity.compute(3,1);
    connectivity.compute(1,0); connectivity.compute(2,0);
    connectivity.compute(2,3); connectivity.compute(0,0);
    save(mesh,"full-initial.mesh");
    VelocitySpace space(mesh,3);
    Math::Vector<Real> physicalLoad;
    const auto initial=evaluate(mesh,space,physicalLoad,true);
    const Real initialVolume=mesh.getVolume(Obstacle);
    std::cout << "initial rho=" << initial.rho << " panels=" << initial.panels << std::endl;

    GridFunction rhoGradient(space), volumeGradient(space);
    TrialFunction trial(space); TestFunction test(space);
    auto normal=FaceNormal(mesh); normal.traceOf(Fluid);
    LinearForm volumeForm(test);
    volumeForm=FaceIntegral(RealFunction{-1},Dot(normal,test)).over(Gamma);
    volumeForm.assemble();
    {
      auto metric=Integral(16*h*h*Jacobian(trial),Jacobian(test))+Integral(trial,test);
      Problem problem(trial,test);
      problem=metric+DirichletBC(trial,VectorFunction{0,0,0}).on(Outer);
      problem.assemble();
      using System=std::remove_reference_t<decltype(problem.getLinearSystem())>;
#ifdef RODIN_USE_MUMPS
      Solver::MUMPS<System> solver(problem);
#else
      Solver::SparseLU<System> solver(problem);
#endif
      for (size_t state=0; state<2; ++state)
      {
        problem=metric-Integral(VectorFunction{0,0,0},test)+
          DirichletBC(trial,VectorFunction{0,0,0}).on(Outer);
        problem.assemble(AssemblyTarget::RHS);
        auto& system=problem.getLinearSystem();
        if (state==0) system.getVector()-=physicalLoad;
        else system.getVector()+=volumeForm.getVector();
        problem.solve(solver);
        const Real residual=(system.getOperator()*system.getSolution()-system.getVector()).norm()/
          std::max(system.getVector().norm(),Real(1));
        if (!solver.success() || !std::isfinite(residual) || residual>1e-9)
          throw std::runtime_error("Full-domain native H1 identification failed.");
        if (state==0) rhoGradient=trial.getSolution();
        else volumeGradient=trial.getSolution();
      }
    }
    const Real metric=volumeForm(volumeGradient);
    if (std::abs(metric)<1e-20) throw std::runtime_error("Zero native volume metric.");
    GridFunction theta(space);
    theta=rhoGradient-volumeForm(rhoGradient)/metric*volumeGradient;
    P1 scalarSpace(mesh);
    GridFunction magnitude(scalarSpace);
    magnitude=Frobenius(theta);
    theta/=std::max(magnitude.max(),Real(1e-30));
    const Real slope=-physicalLoad.dot(theta.getData());
    if (!(slope>0)) throw std::runtime_error("Native BEM transport direction is not ascent.");
    GridFunction distance(scalarSpace);
    Distance::Eikonal(distance).setInterior(Obstacle).setInterface(Gamma).solve().sign();
    GridFunction advected(scalarSpace);
    const auto flow=Flow(-dt,distance,theta,Math::RungeKutta::RK4{});
    const NodalTransport transport(mesh,distance,flow,advected);
    const Real increment=(advected.getData()-distance.getData()).lpNorm<Eigen::Infinity>();
    if (!(increment>1e-12)) throw std::runtime_error("No actual level-set transport occurred.");
    prepare(mesh);
    MMG::Mesh candidate=cutter.discretize(advected);
    save(candidate,"full-updated.mesh");
    VelocitySpace candidateSpace(candidate,3);
    Math::Vector<Real> unused;
    const auto updated=evaluate(candidate,candidateSpace,unused,false);
    if (std::abs(updated.rho-initial.rho)<1e-12)
      throw std::runtime_error("Reconstructed BEM design did not change.");
    std::ofstream record("bem-mmg-validation.json"); record << std::setprecision(17);
    record << "{\"passed\":true,\"route\":\"full-volume-levelset-mmg\",\"symmetry_enforced\":false,\"points\":" << points
      << ",\"seed_amplitude\":0.02,\"dt\":" << dt << ",\"initial_rho\":" << initial.rho
      << ",\"updated_rho\":" << updated.rho << ",\"improved\":" << (updated.rho>initial.rho ? "true":"false")
      << ",\"predicted_delta_rho\":" << dt*slope << ",\"volume_slope\":" << volumeForm(theta)
      << ",\"relative_volume_change\":" << candidate.getVolume(Obstacle)/initialVolume-1
      << ",\"transport_increment\":" << increment << ",\"transport_attempted\":" << transport.getAttemptedCount()
      << ",\"transport_fallback\":" << transport.getFallbackCount()
      << ",\"initial_panels\":" << initial.panels << ",\"updated_panels\":" << updated.panels
      << ",\"initial_epsilon\":" << initial.epsilon << ",\"updated_epsilon\":" << updated.epsilon
      << ",\"updated_residual\":" << updated.residual << ",\"updated_open_edges\":" << updated.openEdges << "}\n";
    if (!record) throw std::runtime_error("Cannot save real MMG update evidence.");
    std::cout << "real_mmg_update rho_before=" << initial.rho << " rho_after=" << updated.rho << std::endl;
    return 0;
  }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
