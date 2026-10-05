#include "SewedOutput.h"
#include "Common.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <limits>
#include <Rodin/IO/MEDIT.h>

using namespace KelvinBall;

void require(bool condition, const char* message)
{
  if (!condition)
    throw std::runtime_error(message);
}

int main()
{
  // Saved transport replays must not perturb coordinates during MEDIT input.
  const Math::SpatialPoint saved{0.9999999999999998, 1.0554945123781563,
    4.440892098500626e-16};
  std::ostringstream row;
  row << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << saved(0) << ' ' << saved(1) << ' ' << saved(2) << " 13";
  const std::string text = row.str();
  const auto parsed = IO::MEDIT::ParseVertex(3)(text.begin(), text.end());
  require(bool(parsed), "MEDIT vertex parse failed");
  require((parsed->vertex.array() == saved.array()).all(),
    "MEDIT max_digits10 coordinate round trip changed saved geometry");
  require(parsed->attribute == 13, "MEDIT vertex attribute changed");
  const auto& rotations = SewedOutput::getCubeRotations();
  require(rotations.size() == ChamberMultiplicity, "incorrect group order");
  const Math::SpatialPoint x{0.21, 0.53, 0.82};
  size_t chambers = 0;
  for (const auto& rotation : rotations)
  {
    require(std::abs(rotation.determinant() - 1) < 1e-14, "improper rotation");
    require((rotation.transpose()*rotation - Math::SpatialMatrix<Real>::Identity(3,3)).norm()
      < 1e-14, "non-orthogonal rotation");
    const Math::SpatialPoint y = rotation * x;
    require(std::abs(initialRadius(x) - initialRadius(y)) < 1e-14,
      "seed is not invariant under its group");
    const bool inside = Tetrahedral ? y(0) >= std::abs(y(2)) && y(1) >= std::abs(y(2))
      : y(0) >= y(1) && y(1) >= std::abs(y(2));
    chambers += inside;
  }
  require(chambers == 1, "chamber does not tile a generic orbit once");
  if (Tetrahedral)
  {
    Math::SpatialPoint quarter(3), mirror(3), boundary(3);
    quarter = Math::SpatialPoint{-x(1), x(0), x(2)};
    mirror = Math::SpatialPoint{-x(0), x(1), x(2)};
    require(std::abs(initialRadius(x)-initialRadius(quarter)) > 1e-4,
      "seed accidentally retains cubic symmetry");
    require(std::abs(initialRadius(x)-initialRadius(mirror)) > 1e-4,
      "seed accidentally retains reflection symmetry");
    boundary = Math::SpatialPoint{0.8, 0.3, 0.3};
    const Math::SpatialPoint first = RotationPairs[0].rotation * boundary;
    require(std::abs(first(0)-first(2)) < 1e-14 && first(1)>=first(0),
      "incorrect positive cut pairing");
    boundary = Math::SpatialPoint{0.3, 0.8, -0.3};
    const Math::SpatialPoint second = RotationPairs[1].rotation * boundary;
    require(std::abs(second(1)+second(2)) < 1e-14 && second(0)>=second(1),
      "incorrect negative cut pairing");
  }
  else
    require(initialRadius(x) == 1, "cubic seed changed");
  std::cout << "PASS group=" << ChamberMultiplicity << " seed and chamber\n";
}
