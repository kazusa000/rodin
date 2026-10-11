/* Distributed under the Boost Software License, Version 1.0. */
#include "Sphere.h"
#include "SewedOutput.h"
#include <iostream>
using namespace KelvinBall;
int main()
{
  Configuration zero, perturbed;
  perturbed.seedAmplitude = 0.02;
  Sphere original(zero), seed(perturbed);
  const Math::SpatialPoint x{0.81, 0.47, -0.29};
  if (original.levelSet(x) != x.norm() - 1)
    throw std::runtime_error("Default sphere changed.");
  for (const auto& rotation : SewedOutput::getCubeRotations())
    if (std::abs(seed.levelSet(rotation*x)-seed.levelSet(x)) > 1e-14)
      throw std::runtime_error("Seed is not O24 invariant.");
  Math::SpatialPoint mirror = x;
  mirror(2) *= -1;
  if (std::abs(seed.levelSet(x)-seed.levelSet(mirror)) < 1e-5)
    throw std::runtime_error("Seed lost chirality.");
  std::cout << "PASS fixed chiral O24 seed and default sphere\n";
}
