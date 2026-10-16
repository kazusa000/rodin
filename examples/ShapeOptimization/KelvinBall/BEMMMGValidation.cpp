#include "BEM.h"
#include "Sphere.h"
#include "NodalTransport.h"
#include <Rodin/Advection/Lagrangian.h>
#include <Rodin/Distance/Eikonal.h>
#include <Rodin/MMG.h>
#include <Rodin/MMG/MeshPrinter.h>
#include <Rodin/MMG/MeshLoader.h>
#include <Rodin/Solver/SparseLU.h>
#ifdef RODIN_USE_MUMPS
#include <Rodin/Solver/MUMPS.h>
#endif
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>

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

int main(int argc, char** argv)
{
  try
  {
    std::cout << std::setprecision(17);
    constexpr size_t points=13;
    constexpr Real h=Real(4)/(points-1);
    Real seedAmplitude=.02, volumeTolerance=.005;
    Real requestedTargetVolume=0, expectedRho=0;
    bool enforceVolumeTolerance=true;
    std::string initialMesh;
    size_t updates=0, maxBacktracks=6;
    // Operational bound for finite background campaigns, not a convergence criterion.
    constexpr size_t maximumUpdates=1000;
    // Weak Armijo gain and halving follow the existing Rodin globalization policy.
    constexpr Real armijo=1e-4, reduction=.5;
    std::vector<Real> stepFactors{.01};
    for (int argument=1; argument<argc; ++argument)
    {
      const std::string option(argv[argument]);
      if (option.starts_with("--seed-amplitude="))
        seedAmplitude=std::stod(option.substr(17));
      else if (option.starts_with("--updates=")) updates=std::stoul(option.substr(10));
      else if (option.starts_with("--volume-tolerance=")) volumeTolerance=std::stod(option.substr(19));
      else if (option=="--no-volume-rejection") enforceVolumeTolerance=false;
      else if (option.starts_with("--max-backtracks=")) maxBacktracks=std::stoul(option.substr(17));
      else if (option.starts_with("--initial-mesh=")) initialMesh=option.substr(15);
      else if (option.starts_with("--target-volume=")) requestedTargetVolume=std::stod(option.substr(16));
      else if (option.starts_with("--expected-rho=")) expectedRho=std::stod(option.substr(15));
      else if (option.starts_with("--step-factors="))
      {
        stepFactors.clear();
        std::istringstream values(option.substr(15));
        std::string value;
        while (std::getline(values,value,',')) stepFactors.push_back(std::stod(value));
      }
      else throw std::runtime_error("Unknown full-domain BEM validation option: "+option);
    }
    if (!std::isfinite(seedAmplitude) || seedAmplitude<0 || seedAmplitude>1 || stepFactors.empty())
      throw std::runtime_error("Use a seed amplitude in [0,1] and at least one step factor.");
    for (const Real factor:stepFactors)
      if (!std::isfinite(factor) || factor<0 || factor>1)
        throw std::runtime_error("Step factors must be in [0,1]; zero is a reconstruction control.");
    if (updates>maximumUpdates || maxBacktracks>10 || !std::isfinite(volumeTolerance) ||
        volumeTolerance<=0 || volumeTolerance>.05)
      throw std::runtime_error("Use at most 1000 updates, ten backtracks and a volume tolerance in (0,0.05].");
    if (!std::isfinite(requestedTargetVolume) || requestedTargetVolume<0 ||
        !std::isfinite(expectedRho) || expectedRho<0)
      throw std::runtime_error("Target volume and expected rho must be finite and nonnegative.");
    if (!initialMesh.empty() && (requestedTargetVolume<=0 || expectedRho<=0))
      throw std::runtime_error("Restart requires the original target volume and accepted rho.");
    if (updates>0 && (stepFactors.size()!=1 || stepFactors.front()<=0))
      throw std::runtime_error("Consecutive updates need exactly one positive initial step factor.");
    if (updates>0 && (std::filesystem::exists("bem-mmg-history.json") ||
                     std::filesystem::exists("accepted-0.mesh")))
      throw std::runtime_error("Refusing to overwrite an existing accepted history.");
    MMG::LevelSetDiscretizer cutter;
    cutter.split(Fluid,{Obstacle,Fluid}).setBoundaryReference(Gamma)
      .setHMin(.8*h).setHMax(1.2*h).setHausdorff(.1*h*h)
      .setGradation(2).setAngleDetection(false);
    MMG::Mesh mesh;
    if (initialMesh.empty())
    {
      auto cube=KelvinBall::Mesh::UniformGrid(Polytope::Type::Tetrahedron,{points,points,points});
      cube.scale(h);
      for (Index i=0; i<cube.getVertexCount(); ++i)
        cube.setVertexCoordinates(i,cube.getVertexCoordinates(i)-Math::SpatialPoint{2,2,2});
      MMG::Mesh background(std::move(cube)); prepare(background);
      Configuration configuration; configuration.seedAmplitude=seedAmplitude;
      const Sphere seed(configuration);
      P1 backgroundSpace(background);
      GridFunction initialPhi(backgroundSpace);
      initialPhi=RealFunction([&seed](const Point& p) { return seed.levelSet(p.vector()); });
      mesh=cutter.discretize(initialPhi);
    }
    else
    {
      std::ifstream input(initialMesh);
      if (!input) throw std::runtime_error("Cannot read accepted restart mesh.");
      MMG::MeshLoader(mesh).load(input);
    }
    save(mesh,"full-initial.mesh");
    const Real startingVolume=mesh.getVolume(Obstacle);
    const Real targetVolume=requestedTargetVolume>0 ? requestedTargetVolume : startingVolume;
    if (!std::isfinite(startingVolume) || startingVolume<=0)
      throw std::runtime_error("Initial obstacle volume must be finite and positive.");
    if (enforceVolumeTolerance && updates>0 && std::abs(startingVolume/targetVolume-1)>volumeTolerance)
      throw std::runtime_error("Starting mesh violates the original volume tolerance.");
    std::cout << "initial_geometry_ready restart=" << (!initialMesh.empty())
      << " volume=" << startingVolume << " target_volume=" << targetVolume
      << " expected_rho=" << expectedRho << " requested_updates=" << updates
      << " volume_rejection_enabled=" << enforceVolumeTolerance << std::endl;
    Real seedRho=0, lastRho=0, lastVolume=startingVolume;
    size_t acceptedUpdates=0;
    const auto writeSummary=[&](bool complete, const char* reason)
    {
      std::ofstream record("bem-mmg-history.json"); record << std::setprecision(17);
      record << "{\"complete\":" << (complete ? "true":"false")
        << ",\"route\":\"full-volume-native-bem-mmg\",\"symmetry_enforced\":false"
        << ",\"seed_amplitude\":" << seedAmplitude << ",\"initial_step_factor\":" << stepFactors.front()
        << ",\"requested_updates\":" << updates << ",\"accepted_updates\":" << acceptedUpdates
        << ",\"seed_rho\":" << seedRho << ",\"last_rho\":" << lastRho
        << ",\"target_volume\":" << targetVolume << ",\"last_volume\":" << lastVolume
        << ",\"restart\":" << (initialMesh.empty() ? "false":"true")
        << ",\"starting_volume\":" << startingVolume << ",\"expected_start_rho\":" << expectedRho
        << ",\"volume_projection_enabled\":true,\"volume_rejection_enabled\":" << (enforceVolumeTolerance ? "true":"false")
        << ",\"volume_tolerance\":";
      if (enforceVolumeTolerance) record << volumeTolerance;
      else record << "null";
      record << ",\"armijo\":" << armijo
        << ",\"step_reduction\":" << reduction << ",\"max_backtracks\":" << maxBacktracks
        << ",\"stop_reason\":\"" << reason << "\"}\n";
      if (!record) throw std::runtime_error("Cannot save native optimization summary.");
    };
    if (updates>0) save(mesh,"accepted-0.mesh");
    for (size_t iteration=0; iteration<std::max(updates,size_t(1)); ++iteration)
    {
      // All trial fields borrow this working copy; the accepted mesh remains intact.
      MMG::Mesh domain(mesh);
      auto& connectivity=domain.getConnectivity();
      connectivity.compute(3,2); connectivity.compute(3,1);
      connectivity.compute(1,0); connectivity.compute(2,0);
      connectivity.compute(2,3); connectivity.compute(0,0);
      VelocitySpace space(domain,3);
      Math::Vector<Real> physicalLoad;
      const auto initial=evaluate(domain,space,physicalLoad,true);
      const Real initialVolume=domain.getVolume(Obstacle);
      std::cout << "initial rho=" << initial.rho << " panels=" << initial.panels << std::endl;
      if (iteration==0)
      {
        if (expectedRho>0 && std::abs(initial.rho-expectedRho)>1e-10*std::max(expectedRho,Real(1e-30)))
          throw std::runtime_error("Restart BEM state disagrees with the saved accepted rho.");
        seedRho=initial.rho; lastRho=initial.rho;
      }
      else if (std::abs(initial.rho-lastRho)>1e-10*std::max(initial.rho,Real(1e-30)))
        throw std::runtime_error("Fresh BEM state disagrees with the accepted checkpoint.");
      if (updates>0) writeSummary(false,"running");

      GridFunction rhoGradient(space), volumeGradient(space);
      TrialFunction trial(space); TestFunction test(space);
      auto normal=FaceNormal(domain); normal.traceOf(Fluid);
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
          if (state==0) system.getVector()+=physicalLoad;
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
      P1 scalarSpace(domain);
      GridFunction magnitude(scalarSpace);
      magnitude=Frobenius(theta);
      theta/=std::max(magnitude.max(),Real(1e-30));
      const Real slope=physicalLoad.dot(theta.getData());
      if (!(slope>0)) throw std::runtime_error("Native BEM transport direction is not ascent.");
      GridFunction distance(scalarSpace);
      Distance::Eikonal(distance).setInterior(Obstacle).setInterface(Gamma).solve().sign();
      const Real volumeSlope=volumeForm(theta);
      prepare(domain);
      bool accepted=false;
      const std::filesystem::path iterationDirectory="iteration-"+std::to_string(iteration);
      if (updates>0 && !std::filesystem::create_directory(iterationDirectory))
        throw std::runtime_error("Refusing to overwrite an existing optimization iteration.");
      const size_t trials=updates>0 ? maxBacktracks+1 : stepFactors.size();
      for (size_t step=0; step<trials; ++step)
      {
        const Real factor=updates>0 ? stepFactors.front()*std::pow(reduction,step) : stepFactors[step];
        const Real dt=factor*h;
        const std::filesystem::path directory=updates>0 ? iterationDirectory/("trial-"+std::to_string(step)) :
          (stepFactors.size()==1 ? std::filesystem::path(".") : std::filesystem::path("step-"+std::to_string(step)));
        if ((updates>0 || stepFactors.size()>1) && !std::filesystem::create_directory(directory))
          throw std::runtime_error("Refusing to overwrite existing step evidence.");
        try
        {
          GridFunction advected(scalarSpace);
          size_t attempted=0, fallback=0;
          if (dt==0) advected=distance;
          else
          {
            const auto flow=Flow(-dt,distance,theta,Math::RungeKutta::RK4{});
            const NodalTransport transport(domain,distance,flow,advected);
            attempted=transport.getAttemptedCount(); fallback=transport.getFallbackCount();
          }
          const Real increment=(advected.getData()-distance.getData()).lpNorm<Eigen::Infinity>();
          if (dt>0 && !(increment>1e-12)) throw std::runtime_error("No actual level-set transport occurred.");
          MMG::Mesh candidate=cutter.discretize(advected);
          save(candidate,(directory/"full-updated.mesh").c_str());
          const Real candidateVolume=candidate.getVolume(Obstacle);
          if (!std::isfinite(candidateVolume) || candidateVolume<=0)
            throw std::runtime_error("Candidate obstacle volume must be finite and positive.");
          const Real volumeDrift=candidateVolume/targetVolume-1;
          if (enforceVolumeTolerance && updates>0 && std::abs(volumeDrift)>volumeTolerance)
          {
            std::ofstream rejected(directory/"bem-mmg-validation.json"); rejected << std::setprecision(17);
            rejected << "{\"passed\":false,\"state_solved\":false,\"accepted\":false,"
              << "\"rejection_reason\":\"volume_tolerance\",\"iteration\":" << iteration
              << ",\"attempt\":" << step << ",\"step_factor\":" << factor
              << ",\"dt\":" << dt << ",\"target_volume\":" << targetVolume
              << ",\"candidate_volume\":" << candidateVolume
              << ",\"cumulative_volume_drift\":" << volumeDrift << "}\n";
            if (!rejected) throw std::runtime_error("Cannot save rejected volume evidence.");
            std::cout << "rejected_volume iteration=" << iteration << " factor=" << factor
              << " drift=" << volumeDrift << std::endl;
            continue;
          }
          BEM::Result updated;
          {
            VelocitySpace candidateSpace(candidate,3);
            Math::Vector<Real> unused;
            updated=evaluate(candidate,candidateSpace,unused,false);
          }
          const bool sufficient=updated.rho>initial.rho && updated.rho>=initial.rho+armijo*dt*slope;
          const bool take=updates>0 && sufficient;
          std::ofstream record(directory/"bem-mmg-validation.json"); record << std::setprecision(17);
          record << "{\"passed\":true,\"route\":\"full-volume-levelset-mmg\",\"symmetry_enforced\":false,\"points\":" << points
            << ",\"seed_amplitude\":" << seedAmplitude << ",\"step_factor\":" << factor
            << ",\"zero_step_control\":" << (dt==0 ? "true":"false")
            << ",\"dt\":" << dt << ",\"initial_rho\":" << initial.rho
            << ",\"initial_c\":" << initial.c << ",\"updated_c\":" << updated.c
            << ",\"updated_rho\":" << updated.rho << ",\"improved\":" << (updated.rho>initial.rho ? "true":"false")
            << ",\"predicted_delta_rho\":" << dt*slope << ",\"volume_slope\":" << volumeSlope
            << ",\"relative_volume_change\":" << candidate.getVolume(Obstacle)/initialVolume-1
            << ",\"transport_increment\":" << increment << ",\"transport_attempted\":" << attempted
            << ",\"transport_fallback\":" << fallback
            << ",\"initial_panels\":" << initial.panels << ",\"updated_panels\":" << updated.panels
            << ",\"initial_epsilon\":" << initial.epsilon << ",\"updated_epsilon\":" << updated.epsilon
            << ",\"updated_residual\":" << updated.residual << ",\"updated_open_edges\":" << updated.openEdges
            << ",\"state_solved\":true,\"accepted\":" << (take ? "true":"false")
            << ",\"iteration\":" << iteration << ",\"attempt\":" << step
            << ",\"target_volume\":" << targetVolume << ",\"candidate_volume\":" << candidateVolume
            << ",\"cumulative_volume_drift\":" << volumeDrift
            << ",\"volume_projection_enabled\":true,\"volume_rejection_enabled\":" << (enforceVolumeTolerance ? "true":"false")
            << ",\"armijo_required_rho\":" << initial.rho+armijo*dt*slope
            << ",\"rejection_reason\":\"" << (take ? "none" : (updates>0 ? "insufficient_rho_gain":"diagnostic_only"))
            << "\"}\n";
          if (!record) throw std::runtime_error("Cannot save real MMG update evidence.");
          std::cout << "real_mmg_update seed=" << seedAmplitude << " dt=" << dt
            << " rho_before=" << initial.rho << " rho_after=" << updated.rho << std::endl;
          if (take)
          {
            save(candidate,("accepted-"+std::to_string(iteration+1)+".mesh").c_str());
            mesh=std::move(candidate);
            accepted=true; ++acceptedUpdates; lastRho=updated.rho; lastVolume=candidateVolume;
            writeSummary(acceptedUpdates==updates,acceptedUpdates==updates ? "requested_updates_completed":"running");
            std::cout << "accepted_update index=" << acceptedUpdates << " rho=" << lastRho
              << " volume_drift=" << volumeDrift << " factor=" << factor << std::endl;
            break;
          }
        }
        catch (const std::exception& error)
        {
          if (updates==0 || accepted) throw;
          std::ofstream detail(directory/"error.txt"); detail << error.what() << '\n';
          std::ofstream rejected(directory/"bem-mmg-validation.json");
          rejected << "{\"passed\":false,\"accepted\":false,\"rejection_reason\":\"geometry_or_state_failure\","
            << "\"iteration\":" << iteration << ",\"attempt\":" << step << ",\"step_factor\":" << factor << "}\n";
          if (!detail || !rejected) throw std::runtime_error("Cannot save failed native trial.");
          std::cerr << "rejected_trial iteration=" << iteration << " attempt=" << step << " " << error.what() << std::endl;
        }
      }
      if (updates>0 && !accepted)
      {
        writeSummary(false,"line_search_exhausted");
        return 2;
      }
    }
    return 0;
  }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
