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
  const auto require = [](bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
  };
  Configuration old;
  old.parse("--n=20");
  old.finalize();
  const auto core = old.getGridCoordinates();
  require(old.getH() == Real(2) / 19 && core.size() == 20 && core.back() == 2,
          "Historical uniform grid changed.");
  Configuration graded;
  graded.parse("--graded-chamber");
  graded.parse("--h=0.10526315789473684");
  graded.parse("--outer-radius=15");
  graded.finalize();
  const auto coordinates = graded.getGridCoordinates();
  require(graded.getH() == old.getH(), "Reference h changed with outer box.");
  require(coordinates.back() == 15 && coordinates.size() < 60, "Invalid far field.");
  for (size_t i = 0; i < core.size(); ++i)
    require(coordinates[i] == core[i], "Uniform core coordinates are not exactly preserved.");
  for (size_t i = 1; i < coordinates.size(); ++i)
    require(coordinates[i] > coordinates[i - 1] &&
            coordinates[i] - coordinates[i - 1] <= Real(10) * old.getH() + 1e-12,
            "Invalid graded interval.");
  require(100 * graded.getH() == 100 * old.getH() &&
          Real(.5) * graded.getH() == Real(.5) * old.getH() &&
          2 * graded.getH() == 2 * old.getH(), "Physical optimization parameters changed.");
  const Math::SpatialPoint near{2, 1, 0};
  const Math::SpatialPoint far{15, 1, 0};
  require(graded.getGradedSize(near, old.getH()) == old.getH(), "Core target coarsened.");
  require(graded.getGradedSize(far, old.getH()) == 10 * old.getH(), "Far target cap changed.");
  require(graded.getGradedSize(far, old.getH() / 2) == 5 * old.getH(), "MMG retry scaling changed.");
  Configuration missing;
  missing.parse("--graded-chamber");
  bool rejected = false;
  try { missing.finalize(); } catch (const std::runtime_error&) { rejected = true; }
  require(rejected, "Missing reference spacing was accepted.");
  Configuration invalid = graded;
  invalid.nearFieldRadius = Real(2.01);
  rejected = false;
  try { invalid.finalize(); } catch (const std::runtime_error&) { rejected = true; }
  require(rejected, "A nonmatching uniform core was accepted.");
  std::cout << "h=" << graded.getH() << " core_points=" << core.size()
            << " total_axis_points=" << coordinates.size() << " outer=" << coordinates.back() << '\n';
}
