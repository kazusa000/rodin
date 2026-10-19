/*
 * Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Configuration.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
  using namespace KelvinBall;
  const auto require = [](bool value, const char* message) {
    if (!value)
      throw std::runtime_error(message);
  };
  Configuration uniform;
  uniform.finalize();
  const auto core = uniform.getGridCoordinates();
  Configuration graded;
  graded.parse("--graded-chamber");
  graded.parse("--h=0.16666666666666666");
  graded.parse("--outer-radius=15");
  graded.finalize();
  const auto axis = graded.getGridCoordinates();
  require(graded.getGridSpacing() == uniform.getGridSpacing(), "Reference scale changed");
  require(axis.back() == 15 && axis.size() < 91, "Exterior endpoint or grading failed");
  for (size_t i = 0; i < core.size(); ++i)
    require(axis[i] == core[i], "Near-field coordinates changed");
  constexpr Real roundoff = 1e-12;
  for (size_t i = 1; i < axis.size(); ++i)
    require(axis[i] > axis[i - 1] &&
      axis[i] - axis[i - 1] <= graded.hmax + roundoff, "Invalid exterior interval");
  require(graded.hmin == uniform.hmin && graded.hmax == uniform.hmax,
    "MMG size bounds changed with box size");
  Configuration limiting = graded;
  limiting.outerRadius = 2;
  limiting.finalize();
  require(limiting.getGridCoordinates() == core, "Limiting small box changed");
  for (const char* option : {"--near-field-radius=2.01", "--chamber-gradation=1", "--hmax-factor=0.5"})
  {
    Configuration invalid = graded;
    invalid.parse(option);
    bool rejected = false;
    try { invalid.finalize(); }
    catch (const std::runtime_error&) { rejected = true; }
    require(rejected, "Invalid graded configuration accepted");
  }
  Configuration missing;
  missing.parse("--graded-chamber");
  bool rejected = false;
  try { missing.finalize(); }
  catch (const std::runtime_error&) { rejected = true; }
  require(rejected, "Graded mode accepted without explicit reference h");
  std::cout << "h0=" << graded.getGridSpacing() << " axis_points=" << axis.size()
            << " core_points=" << core.size() << " outer=" << axis.back() << '\n';
}
