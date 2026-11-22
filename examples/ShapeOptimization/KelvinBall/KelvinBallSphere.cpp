/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/* Validate the spherical Kelvin ball on its chamber and sewn domains. */
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include <Rodin/Alert/Raise.h>
#include <Rodin/Alert/Success.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>

#include "Metrics.h"
#include "SewedOutput.h"
#include "Sphere.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace KelvinBall
{
#ifdef KELVIN_AXIAL_C3
  // Twist a three-lobed radial displacement along the body diagonal. Build
  // both hands from identical nodes and cells so remeshing cannot mask a sign error.
  static Mesh handedTestMesh(const Mesh& source, Real outerRadius, bool mirror)
  {
    Mesh::Builder builder;
    builder.initialize(3).nodes(source.getVertexCount());
    for (Index i = 0; i < source.getVertexCount(); ++i)
    {
      Math::SpatialPoint x = source.getVertexCoordinates(i);
      const Real r = x.norm();
      if (r > 0 && r < std::min(Real(3), outerRadius))
      {
        const Real height = x.dot(axialDirection()) / r;
        const Real u = (x(0) - x(1)) / (r * std::sqrt(Real(2)));
        const Real v = (x(0) + x(1) - 2 * x(2)) / (r * std::sqrt(Real(6)));
        const Real phase = 1.8 * height;
        const Real lobes = u * (u*u - 3*v*v) * std::cos(phase) +
          v * (3*u*u - v*v) * std::sin(phase);
        const Real cutoff = std::max(Real(0), Real(1) - std::pow(r / std::min(Real(3), outerRadius), 2));
        x *= 1 + 0.2 * lobes * cutoff * cutoff;
      }
      if (mirror)
        std::swap(x(1), x(2));
      builder.vertex(x);
    }
    for (auto cell = source.getCell(); cell; ++cell)
    {
      IndexArray vertices = cell->getVertices();
      if (mirror)
        std::swap(vertices(0), vertices(1));
      Index index;
      builder.polytope(cell->getGeometry(), std::move(vertices), index);
      builder.attribute({3, index}, cell->getAttribute());
    }
    for (auto face = source.getBoundary(); face; ++face)
    {
      IndexArray vertices = face->getVertices();
      if (mirror)
        std::swap(vertices(0), vertices(1));
      Index index;
      builder.polytope(face->getGeometry(), std::move(vertices), index);
      builder.attribute({2, index}, face->getAttribute());
    }
    return builder.finalize();
  }
#endif
  class SphereValidation
  {
    public:
      SphereValidation(int argc, char** argv)
        : m_argc(argc),
          m_argv(argv)
      {}

      int run();

    private:
      void print(const char* name, const Values& values) const
      {
        std::cout << name << ": k=" << values.k << ", c=" << values.c
                  << ", q=" << values.q << ", rho=" << values.rho;
        if (values.nitscheJump > 0)
          std::cout << ", nitsche_jump=" << values.nitscheJump;
        std::cout << '\n';
      }

      Real relativeDifference(Real chamber, Real full) const
      {
        return std::abs(chamber - full) /
          std::max({std::abs(chamber), std::abs(full), Real(1e-30)});
      }

      void checkReconstruction(Mesh& mesh) const
      {
        mesh.getConnectivity().compute(2, 3);
        size_t openCuts = 0;
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          const auto attribute = face->getAttribute();
          if (!attribute || (*attribute != Gamma && *attribute != Outer))
            ++openCuts;
        }
        if (openCuts != 0)
          throw std::runtime_error(
            "The reconstructed mesh has nonconforming internal chamber cuts.");
      }

      int m_argc;
      char** m_argv;
  };
}

int KelvinBall::SphereValidation::run()
{
  const int argc = m_argc;
  char** argv = m_argv;
  Configuration configuration;
  bool saveMesh = false;
#ifdef KELVIN_AXIAL_C3
  bool checkHandedness = false;
#endif
  for (int argument = 1; argument < argc; ++argument)
  {
    const StringView option(argv[argument]);
    if (configuration.parse(option))
      continue;
#ifdef KELVIN_AXIAL_C3
    if (option == "--check-handedness")
    {
      checkHandedness = true;
      continue;
    }
#endif
    if (option == "--save-mesh")
      saveMesh = true;
    else if (option == "--help")
    {
#ifdef KELVIN_AXIAL_C3
      std::cout << "  --check-handedness          Compare a twisted C3 body with its mirror.\n";
#endif
      std::cout
        << "Usage: " << argv[0] << " [options]\n"
        << "  --n=<points>              Background points per edge (default: 13).\n"
        << "  --h=<size>                Requested initial grid spacing; alternative to "
           "--n.\n"
        << "  --hmin-factor=<value>     Minimum MMG size / reference spacing (default: "
           "0.1).\n"
        << "  --hmax-factor=<value>     Maximum MMG size / reference spacing (default: "
           "10).\n"
        << "  --outer-radius=<value>    Chamber outer radius (default: 2).\n"
        << "  --penalty=<value>         Rotational Nitsche penalty (default: 320).\n"
        << "  --stabilization=<value>   P1--P1 stabilization factor (default: 0.05).\n"
        << "  --save-mesh               Write KelvinBallSphere.mesh.\n";
      return 0;
    }
    else
      throw std::runtime_error("Unknown KelvinBallSphere option: " + std::string(option));
  }
  configuration.finalize();

  auto sphere = Sphere(configuration).discretize(true);
  auto& chamber = sphere.mesh;
  const Parameters parameters{meanElementSize(chamber), configuration.nitschePenalty,
    configuration.stabilizationFactor};
  if (saveMesh)
    chamber.save("KelvinBallSphere.mesh", IO::FileFormat::MEDIT);
  prepare(chamber);
  SubMesh chamberFluid = chamber.trim(Obstacle);
  Metrics metrics(parameters);
  const Values chamberValues = metrics.evaluateChamber(chamberFluid);

  SewedOutput reconstructed(chamberFluid, FlatSet<Attribute>{Gamma, Outer});
  checkReconstruction(reconstructed.getMesh());
  const Values fullValues = metrics.evaluateSewed(reconstructed.getMesh());

  std::cout << std::setprecision(15);
  print("chamber", chamberValues);
  print("full", fullValues);
  std::cout << "relative differences: k="
            << relativeDifference(chamberValues.k, fullValues.k)
            << ", c=" << relativeDifference(chamberValues.c, fullValues.c)
            << ", q=" << relativeDifference(chamberValues.q, fullValues.q)
            << ", rho=" << relativeDifference(chamberValues.rho, fullValues.rho) << '\n';

#ifdef KELVIN_AXIAL_C3
  const Real resistanceScale = std::sqrt(fullValues.k * fullValues.q);
  const Real couplingDifference = std::abs(chamberValues.c - fullValues.c) / resistanceScale;
  std::cout << "scaled chamber/full C difference: " << couplingDifference << '\n';
  if (relativeDifference(chamberValues.k, fullValues.k) > 0.02 ||
      relativeDifference(chamberValues.q, fullValues.q) > 0.02 ||
      couplingDifference > 0.005 || std::max(chamberValues.rho, fullValues.rho) > 0.01)
    throw std::runtime_error("The axial sphere chamber/full comparison failed.");
  if (checkHandedness)
  {
    auto rightMesh = handedTestMesh(reconstructed.getMesh(), configuration.outerRadius, false);
    const auto right = metrics.evaluateSewed(rightMesh);
    auto leftMesh = handedTestMesh(reconstructed.getMesh(), configuration.outerRadius, true);
    const auto left = metrics.evaluateSewed(leftMesh);
    print("hand A", right);
    print("mirror", left);
    if (!(right.c * left.c < 0 && std::min(right.rho, left.rho) > 1e-8) ||
        std::abs(right.c + left.c) > 1e-8 * std::sqrt(right.k * right.q) ||
        relativeDifference(right.k, left.k) > 1e-8 || relativeDifference(right.q, left.q) > 1e-8)
      throw std::runtime_error("Mirroring the axial body did not reverse its signed coupling.");
  }
#endif

  IO::XDMF output("KelvinBallSphere");
  output.setMesh(reconstructed.getMesh());
  output.write(0).flush();
  output.close();

  Alert::Success() << "Wrote KelvinBallSphere.xdmf" << Alert::Raise;
  return 0;
}

int main(int argc, char** argv)
{
  return KelvinBall::SphereValidation(argc, argv).run();
}
