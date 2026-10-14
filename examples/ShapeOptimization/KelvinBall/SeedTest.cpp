/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Sphere.h"
#include "SewedOutput.h"
#include <iostream>
#include <stdexcept>

using namespace KelvinBall;

int main()
{
  Configuration zero, perturbed;
  perturbed.parse("--seed-amplitude=0.02");
  perturbed.finalize();
  Sphere original(zero), seed(perturbed);
  const Math::SpatialPoint point{0.81, 0.47, -0.29};
  if (original.levelSet(point) != point.norm() - 1)
    throw std::runtime_error("Default sphere changed.");
  for (const auto& rotation : SewedOutput::getCubeRotations())
  {
    if (std::abs(seed.levelSet(rotation * point) - seed.levelSet(point)) > 1e-14)
      throw std::runtime_error("Seed is not O24 invariant.");
  }
  Math::SpatialPoint mirror = point;
  mirror(2) *= -1;
  if (std::abs(seed.levelSet(point) - seed.levelSet(mirror)) < 1e-5)
    throw std::runtime_error("Seed lost chirality.");
  if (std::abs(seed.levelSet(point) + seed.levelSet(mirror) -
               2 * (point.norm() - 1)) > 1e-14)
    throw std::runtime_error("Seed perturbation is not reflection odd.");
  Math::SpatialPoint origin(3);
  origin.setZero();
  if (seed.levelSet(origin) != -1)
    throw std::runtime_error("Seed is undefined at its center.");
  std::cout << "PASS fixed 0.02 H3 H6 chiral O24 seed and unchanged default sphere\n";
}
