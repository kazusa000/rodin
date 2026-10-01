/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <fstream>
#include <iomanip>
#include "Metrics.h"
#include "Sphere.h"

using namespace KelvinBall;

// Fixed-geometry boundary-condition experiment. Never advances the shape or
// treats truncated-domain dissipation with nonzero outer trace as body drag.
int main(int argc, char** argv)
{
  Configuration configuration;
  std::string meshFile, outerFile;
  bool prepareOnly = false;
  for (int i = 1; i < argc; ++i)
  {
    const std::string option(argv[i]);
    if (configuration.parse(option)) continue;
    if (option == "--prepare-only") prepareOnly = true;
    else if (option.rfind("--mesh=", 0) == 0) meshFile = option.substr(7);
    else if (option.rfind("--outer-values=", 0) == 0) outerFile = option.substr(15);
    else throw std::runtime_error("Unknown boundary experiment option: " + option);
  }
  configuration.finalize();
  KelvinBall::Mesh chamber;
  if (meshFile.empty()) chamber = Sphere(configuration).discretize().mesh;
  else chamber.load(meshFile, IO::FileFormat::MEDIT);
  prepare(chamber);
  chamber.save("chamber.mesh", IO::FileFormat::MEDIT);
  SubMesh fluid = chamber.trim(Obstacle);
  splitSelfPairedCut(fluid);
  prepare(fluid);
  fluid.save("fluid.mesh", IO::FileFormat::MEDIT);
  if (prepareOnly) return 0;
  const OuterVelocity outer(outerFile);
  VelocitySpace Vh(fluid, 3);
  PressureSpace Qh(fluid);
  RotatedNitscheIntegrator coupling(fluid, MasterCuts, 0.01, 0.25);
  GridFunction uT0(Vh), uT1(Vh), uT2(Vh), uR0(Vh), uR1(Vh), uR2(Vh);
  GridFunction pT0(Qh), pT1(Qh), pT2(Qh), pR0(Qh), pR1(Qh), pR2(Qh);
  const Metrics metrics({configuration.getH(), configuration.nitschePenalty,
    configuration.stabilizationFactor}, &outer);
  const auto energy = metrics.evaluateChamber(Vh, Qh, coupling,
    uT0,uT1,uT2,uR0,uR1,uR2,pT0,pT1,pT2,pR0,pR1,pR2);
  const std::array<const decltype(uT0)*,6> velocities{&uT0,&uT1,&uT2,&uR0,&uR1,&uR2};
  const std::array<const decltype(pT0)*,6> pressures{&pT0,&pT1,&pT2,&pR0,&pR1,&pR2};
  std::ofstream fields("states.txt");
  fields << std::setprecision(17);
  for (Index vertex = 0; vertex < fluid.getVertexCount(); ++vertex)
  {
    const auto& x = fluid.getVertexCoordinates(vertex);
    const auto vd = Vh.getDOFs(0, vertex);
    const auto pd = Qh.getDOFs(0, vertex);
    fields << x(0) << ' ' << x(1) << ' ' << x(2);
    for (const auto* u : velocities)
      for (size_t j=0; j<3; ++j) fields << ' ' << u->getData()(vd(j));
    for (const auto* p : pressures) fields << ' ' << p->getData()(pd(0));
    fields << '\n';
  }
  std::ofstream diagnostics("energy-diagnostic.txt");
  diagnostics << std::setprecision(17) << energy.k << ' ' << energy.c << ' '
    << energy.q << ' ' << energy.rho << ' ' << energy.nitscheJump << '\n';
  if (!fields || !diagnostics) throw std::runtime_error("Failed writing state diagnostics.");
  return 0;
}
