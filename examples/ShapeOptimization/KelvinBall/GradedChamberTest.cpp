#include "Configuration.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
  using namespace KelvinBall;
  const auto require = [](bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
  };
  Configuration baseline;
  baseline.finalize();
  const auto original = baseline.getGridCoordinates();
  require(original.size() == 13 && original.back() == 2, "Original grid changed.");
  Configuration graded;
  graded.parse("--graded-chamber");
  graded.parse("--outer-radius=15");
  graded.parse("--h=0.16666666666666666");
  graded.parse("--near-field-radius=2");
  graded.parse("--chamber-gradation=1.3");
  graded.finalize();
  const auto coordinates = graded.getGridCoordinates();
  require(graded.getGridSpacing() == baseline.getGridSpacing(), "Reference spacing changed.");
  require(graded.hmin == baseline.hmin && graded.hmax == baseline.hmax,
          "Near/far MMG size changed with box.");
  require(coordinates.back() == 15 && coordinates.size() < 91, "Graded extent or size invalid.");
  for (size_t i = 0; i < original.size(); ++i)
    require(std::abs(coordinates[i] - original[i]) < 1e-12, "Near-field coordinates changed.");
  for (size_t i = 1; i < coordinates.size(); ++i)
  {
    const auto width = coordinates[i] - coordinates[i - 1];
    require(width > 0 && width <= graded.hmax + 1e-12, "Invalid graded interval.");
  }
  Configuration missing;
  missing.parse("--graded-chamber");
  bool refused = false;
  try { missing.finalize(); } catch (const std::runtime_error&) { refused = true; }
  require(refused, "Missing explicit near-field spacing accepted.");
  Configuration historical;
  historical.parse("--graded-chamber");
  historical.parse("--outer-radius=15");
  historical.parse("--h=0.10526315789473684");
  historical.finalize();
  const Math::SpatialPoint interface{1, Real(0.5), Real(0.2)};
  const Math::SpatialPoint core{2, 2, 2};
  const Math::SpatialPoint far{15, 15, 15};
  require(historical.getGradedSize(interface, historical.getH()) == historical.getH(),
          "Historical interface target size was not restored.");
  require(historical.getGradedSize(core, historical.getH()) == historical.getH(),
          "The old L2 region was coarsened.");
  require(historical.getGradedSize(far, historical.getH()) == historical.hmax,
          "Far-field target cap changed.");
  require(std::abs(historical.getGradedSize(far, historical.getH() / 2) -
                   historical.hmax / 2) < 1e-14, "MMG retry scale changed.");
  std::cout << "graded_axis_points=" << coordinates.size()
            << " h0=" << graded.getGridSpacing()
            << " hmin=" << graded.hmin << " hmax=" << graded.hmax << '\n';
}
