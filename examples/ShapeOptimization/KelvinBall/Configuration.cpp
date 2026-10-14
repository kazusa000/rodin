/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Configuration.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace KelvinBall
{
  bool Configuration::parse(std::string_view option)
  {
    if (option.rfind("--n=", 0) == 0)
    {
      points = std::stoul(std::string(option.substr(4)));
      pointsSpecified = true;
    }
    else if (option.rfind("--h=", 0) == 0)
    {
      requestedH = std::stod(std::string(option.substr(4)));
      hSpecified = true;
    }
    else if (option.rfind("--outer-radius=", 0) == 0)
      outerRadius = std::stod(std::string(option.substr(15)));
    else if (option == "--graded-chamber")
      gradedChamber = true;
    else if (option.rfind("--near-field-radius=", 0) == 0)
      nearFieldRadius = std::stod(std::string(option.substr(20)));
    else if (option.rfind("--chamber-gradation=", 0) == 0)
      chamberGradation = std::stod(std::string(option.substr(20)));
    else if (option.rfind("--seed-amplitude=", 0) == 0)
      seedAmplitude = std::stod(std::string(option.substr(17)));
    else if (option.rfind("--penalty=", 0) == 0)
      nitschePenalty = std::stod(std::string(option.substr(10)));
    else if (option.rfind("--stabilization=", 0) == 0)
      stabilizationFactor = std::stod(std::string(option.substr(16)));
    else if (option.rfind("--background-hmin=", 0) == 0)
      backgroundHMin = std::stod(std::string(option.substr(18)));
    else if (option.rfind("--background-hmax=", 0) == 0)
      backgroundHMax = std::stod(std::string(option.substr(18)));
    else if (option.rfind("--background-hausdorff=", 0) == 0)
      backgroundHausdorff = std::stod(std::string(option.substr(23)));
    else if (option.rfind("--background-gradation=", 0) == 0)
      backgroundGradation = std::stod(std::string(option.substr(23)));
    else if (option == "--mmg-adapt")
      adapt = true;
    else if (option.rfind("--mmg-adapt-interface-size=", 0) == 0)
      adaptInterfaceSize = std::stod(std::string(option.substr(27)));
    else if (option.rfind("--mmg-adapt-far-size=", 0) == 0)
      adaptFarSize = std::stod(std::string(option.substr(21)));
    else if (option.rfind("--mmg-adapt-width=", 0) == 0)
      adaptWidth = std::stod(std::string(option.substr(18)));
    else if (option.rfind("--mmg-adapt-gradation=", 0) == 0)
      adaptGradation = std::stod(std::string(option.substr(22)));
    else if (option.rfind("--mmg-snap=", 0) == 0)
      mmgSnap = std::stod(std::string(option.substr(11)));
    else if (option.rfind("--mmg-retries=", 0) == 0)
      mmgRetries = std::stoul(std::string(option.substr(14)));
    else
      return false;
    return true;
  }

  void Configuration::finalize()
  {
    if (pointsSpecified && hSpecified)
      throw std::runtime_error("Use either --n or --h, not both.");
    if (!(outerRadius > 1))
      throw std::runtime_error(
        "The chamber outer radius must exceed the unit sphere radius.");
    if (hSpecified)
    {
      if (!(requestedH > 0))
        throw std::runtime_error("The requested mesh size must be positive.");
      points = static_cast<size_t>(std::ceil(outerRadius / requestedH)) + 1;
    }
    if (points < 5)
      throw std::runtime_error("The chamber uniform grid requires --n >= 5.");
    if (!std::isfinite(seedAmplitude) || seedAmplitude < 0 || seedAmplitude > Real(0.1))
      throw std::runtime_error("The seed amplitude must be in [0,0.1].");
    if (gradedChamber)
    {
      if (!hSpecified || !(nearFieldRadius > 1) || nearFieldRadius > outerRadius ||
          !(chamberGradation > 1))
        throw std::runtime_error("A graded chamber requires --h, 1 < core <= box, and growth > 1.");
      points = getGridCoordinates().size();
    }
    hmin = Real(0.1) * getH();
    hmax = Real(10) * getH();
    if (!(nitschePenalty > 0))
      throw std::runtime_error("The Nitsche penalty must be positive.");
    if (stabilizationFactor < 0)
      throw std::runtime_error("The pressure stabilization must be nonnegative.");
    if (!(backgroundHMin > 0) || !(backgroundHMax >= backgroundHMin))
      throw std::runtime_error(
        "The background sizes must satisfy 0 < --background-hmin <= --background-hmax.");
    if (!(backgroundHausdorff > 0))
      throw std::runtime_error("The background Hausdorff tolerance must be positive.");
    if (!(backgroundGradation > 1))
      throw std::runtime_error("The background gradation must exceed one.");
    if (!(adaptInterfaceSize > 0) || !(adaptFarSize > 0))
      throw std::runtime_error("The adaptation sizes must be positive.");
    if (!(adaptWidth > 0))
      throw std::runtime_error("The adaptation width must be positive.");
    if (!(adaptGradation > 1))
      throw std::runtime_error("The adaptation gradation must exceed one.");
    if (!(mmgSnap >= 0) || !(mmgSnap < 0.5))
      throw std::runtime_error(
        "The snapping fraction must satisfy 0 <= --mmg-snap < 0.5.");
  }

  Real Configuration::getH() const
  {
    if (gradedChamber)
      return requestedH;
    return outerRadius / static_cast<Real>(points - 1);
  }

  std::vector<Real> Configuration::getGridCoordinates() const
  {
    std::vector<Real> coordinates{Real(0)};
    if (!gradedChamber)
    {
      for (size_t i = 1; i < points; ++i)
        coordinates.push_back(i * getH());
      return coordinates;
    }
    Real step = getH();
    constexpr Real closureTolerance = 1e-12;
    for (const Real end : {nearFieldRadius, outerRadius})
    {
      while (coordinates.back() < end)
      {
        if (end == outerRadius && end != nearFieldRadius)
          step = std::min(step * chamberGradation, Real(10) * getH());
        const Real remaining = end - coordinates.back();
        const Real next = remaining <= step * (Real(1) + closureTolerance)
          ? end : coordinates.back() + step;
        if (!(next > coordinates.back()))
          throw std::runtime_error("The graded coordinates cannot advance.");
        coordinates.push_back(next);
      }
    }
    return coordinates;
  }

  Real Configuration::getGradedSize(const Math::SpatialPoint& point, Real h) const
  {
    const Real extent = point.cwiseAbs().maxCoeff();
    // A geometric progression of edge lengths yields a size increasing linearly
    // with distance. The entire old box retains its interface target h.
    const Real distance = std::max(Real(0), extent - nearFieldRadius);
    const Real scale = h / getH();
    return std::min(Real(10) * h,
      adaptInterfaceSize * h + scale * (chamberGradation - Real(1)) * distance);
  }
}
