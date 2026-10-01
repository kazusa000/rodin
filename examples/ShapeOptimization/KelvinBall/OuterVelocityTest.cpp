/* Distributed under the Boost Software License, Version 1.0. */
#include <cstdio>
#include <fstream>
#include "OuterVelocity.h"

using namespace KelvinBall;

int main()
{
  auto require = [](bool value) {
    if (!value) throw std::runtime_error("Outer velocity regression failed.");
  };
  const Math::SpatialPoint surface{1, 0, 0}, far{2, 0, 0};
  const OuterVelocity sphere("analytic-sphere"), zero;
  require(zero.value(far, 0).norm() == 0);
  require(std::abs(sphere.value(surface, 0)(0) - 1) < 1e-14);
  require(std::abs(sphere.value(surface, 1)(1) - 1) < 1e-14);
  require(std::abs(sphere.value(far, 0)(0) - 0.6875) < 1e-14);
  require(std::abs(sphere.value(far, 1)(1) - 0.40625) < 1e-14);
  require(std::abs(sphere.value(far, 5)(1) - 0.25) < 1e-14);
  const char* path = "outer-velocity-test.txt";
  {
    std::ofstream table(path);
    table << "1\n2 0 0";
    for (size_t j=0; j<18; ++j) table << ' ' << j;
    table << '\n';
  }
  const OuterVelocity table(path);
  require(table.value(far, 5)(2) == 17);
  bool missing = false;
  try { table.value(surface, 0); }
  catch (const std::runtime_error&) { missing = true; }
  require(missing);
  std::remove(path);
  return 0;
}
