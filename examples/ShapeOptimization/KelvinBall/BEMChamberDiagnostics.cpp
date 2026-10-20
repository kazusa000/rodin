/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "BEM.h"
#include "Sphere.h"
#include <Rodin/MMG/MeshLoader.h>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace KelvinBall;

int main(int argc,char** argv)
{
  try
  {
    constexpr Real stateTolerance=1e-9;
    const bool selfCheck=argc==2 && std::string(argv[1])=="--self-check";
    if (!selfCheck && argc!=3)
      throw std::runtime_error("Usage: KelvinBallBEMChamberDiagnostics CHAMBER.mesh EXPECTED_RHO | --self-check");
    MMG::Mesh mesh;
    if (selfCheck)
    {
      Configuration configuration;
      configuration.points=7;configuration.hminFactor=.8;configuration.hmaxFactor=1.2;
      configuration.seedAmplitude=.1;configuration.finalize();
      mesh=Sphere(configuration).discretize(true,0,true).mesh;
    }
    else
    {
      std::ifstream input(argv[1]);
      if (!input) throw std::runtime_error("Cannot read accepted chamber mesh.");
      MMG::MeshLoader(mesh).load(input);
    }
    prepare(mesh);
    auto solid=mesh.keep(Obstacle);
    solid.getConnectivity().compute(2,3);
    const auto skin=solid.skin();
    const auto surface=skin.keep(Gamma);
    BEM bem; // Production defaults: order=3, epsilonFactor=0.6, viscosity=Mu.
    Real expected=0;
    if (selfCheck) expected=bem.evaluateQuotient(surface,false).rho;
    else expected=std::stod(argv[2]);
    BEM::Diagnostics diagnostics;
    const auto start=std::chrono::steady_clock::now();
    const auto result=bem.evaluateQuotient(surface,false,&diagnostics);
    const Real totalSeconds=std::chrono::duration<Real>(std::chrono::steady_clock::now()-start).count();
    const Real rhoError=std::abs(result.rho-expected)/std::max(std::abs(expected),Real(1e-30));
    const bool passed=result.openEdges==0 && result.residual<=stateTolerance &&
      rhoError<=stateTolerance && diagnostics.inverseResidual<=stateTolerance &&
      diagnostics.conditionOne>=1;
    std::ofstream record("bem-chamber-diagnostics.json");
    record<<std::setprecision(17)<<"{\"passed\":"<<(passed ? "true":"false")
      <<",\"self_check\":"<<(selfCheck ? "true":"false")
      <<",\"gradient_computed\":false,\"shape_updates\":0"
      <<",\"design_domain\":\"O24-chamber\",\"bem_domain\":\"infinite-exterior\""
      <<",\"order\":3,\"epsilon_factor\":0.6,\"viscosity\":"<<Mu
      <<",\"rho\":"<<result.rho<<",\"expected_rho\":"<<expected<<",\"rho_relative_error\":"<<rhoError
      <<",\"k\":"<<result.k<<",\"c\":"<<result.c<<",\"q\":"<<result.q
      <<",\"residual\":"<<result.residual<<",\"reciprocity\":"<<result.reciprocity
      <<",\"epsilon\":"<<result.epsilon<<",\"panels\":"<<result.panels
      <<",\"state_dofs\":"<<result.stateDOFs<<",\"independent_dofs\":"<<result.independentDOFs
      <<",\"open_edges\":"<<result.openEdges<<",\"basis_seconds\":"<<diagnostics.basisSeconds
      <<",\"families\":[";
    for (size_t family=0;family<2;++family)
    {
      if (family) record<<",";
      record<<"{\"family\":\""<<(family==0 ? "translation":"rotation")<<"\""
        <<",\"assembly_seconds\":"<<diagnostics.assemblySeconds[family]
        <<",\"reduction_seconds\":"<<diagnostics.reductionSeconds[family]
        <<",\"lu_seconds\":"<<diagnostics.factorSeconds[family]
        <<",\"back_substitution_seconds\":"<<diagnostics.solveSeconds[family]<<"}";
    }
    record<<"],\"condition_norm\":\"induced matrix 1-norm\""
      <<",\"condition_matrix\":\"E^T A E, original production basis, no scaling or preconditioning\""
      <<",\"condition_method\":\"norm(A,1)*norm(LU.solve(I),1)\""
      <<",\"condition_one\":"<<diagnostics.conditionOne
      <<",\"condition_one_lu_estimate\":"<<diagnostics.conditionOneEstimate
      <<",\"operator_norm_one\":"<<diagnostics.operatorNormOne
      <<",\"inverse_norm_one\":"<<diagnostics.inverseNormOne
      <<",\"normalized_inverse_residual_one\":"<<diagnostics.inverseResidual
      <<",\"condition_diagnostic_seconds\":"<<diagnostics.conditionSeconds
      <<",\"total_seconds\":"<<totalSeconds<<"}\n";
    if (!record) throw std::runtime_error("Cannot save BEM diagnostic record.");
    std::cout<<std::setprecision(17)<<"rho="<<result.rho<<" condition_one="<<diagnostics.conditionOne
      <<" passed="<<passed<<std::endl;
    return passed ? 0 : 2;
  }
  catch (const std::exception& error)
  {
    std::cerr<<error.what()<<std::endl;
    return 1;
  }
}
