#include "SewedOutput.h"
#include "Common.h"
#include "CrossedEdgeSnapGuard.h"
#include "RotationalQuadDiagonal.h"
#include "Sphere.h"
#include "PeriodicCuts.h"
#include <Rodin/Location.h>
#include <Rodin/Variational.h>
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
  // Manufactured transfer: on a Cartesian simplex grid, interpolation of
  // |x|^2 is the sum of its independent 1D linear interpolants. This field
  // also has exact T/O traces, so periodic projection cannot alter the result.
  Configuration coarseConfiguration,fineConfiguration;
  coarseConfiguration.points=5; fineConfiguration.points=9;
  MMG::Mesh coarse(Sphere(coarseConfiguration).makeUniformChamber());
  MMG::Mesh fine(Sphere(fineConfiguration).makeUniformChamber());
  Variational::P1 coarseSpace(coarse),fineSpace(fine);
  Variational::GridFunction quadratic(coarseSpace),transferred(fineSpace);
  quadratic=Variational::RealFunction([](const Geometry::Point& point) {
    return point.getPhysicalCoordinates().squaredNorm();
  });
  const Location::AABB<MMG::Mesh> sourceLocator(coarse);
  transferred=Variational::RealFunction([&](const Geometry::Point& point) {
    const auto located=sourceLocator.locate(3,point.getPhysicalCoordinates());
    require(bool(located),"manufactured background transfer missed an interior node");
    return quadratic.getValue(*located);
  });
  PeriodicCuts(fine).project(transferred);
  const Real coarseH=coarseConfiguration.getH();
  for (Index vertex=0;vertex<fine.getVertexCount();++vertex) {
    const auto x=fine.getVertexCoordinates(vertex);
    Real exact=0;
    for (size_t j=0;j<3;++j) {
      const Real left=std::floor(x(j)/coarseH)*coarseH;
      const Real fraction=(x(j)-left)/coarseH;
      exact+=(1-fraction)*left*left+fraction*(left+coarseH)*(left+coarseH);
    }
    require(std::abs(transferred[vertex]-exact)<1e-12,
      "background field transfer does not match the manufactured P1 interpolant");
  }
  require(!sourceLocator.locate(3,Math::SpatialPoint{3,3,0}),
    "background transfer accepted a point outside the source domain");
  const std::array<Math::SpatialPoint,4> quad{Math::SpatialPoint{.7,.2,.3},
    Math::SpatialPoint{1.1,.2,.3},Math::SpatialPoint{1.05,.43,.3},
    Math::SpatialPoint{.71,.4,.3}};
  const int chosen=rotationalQuadDiagonal(quad,1e-11);
  require(chosen>=0,"generic quad unnecessarily retained its centre fan");
  const auto endpoints=[](const auto& q,int diagonal) {
    return diagonal==0 ? std::array<Math::SpatialPoint,2>{q[0],q[2]}
      : std::array<Math::SpatialPoint,2>{q[1],q[3]};
  };
  const auto expected=endpoints(quad,chosen);
  for (const auto& rotation:SewedOutput::getCubeRotations())
    for (int offset=0;offset<4;++offset)
      for (const int orientation:{-1,1}) {
        std::array<Math::SpatialPoint,4> moved;
        for (int i=0;i<4;++i) moved[i]=rotation*quad[(offset+4+orientation*i)%4];
        const int selected=rotationalQuadDiagonal(moved,1e-11);
        require(selected>=0,"rotated generic quad lost its diagonal");
        const auto actual=endpoints(moved,selected);
        const Math::SpatialPoint a=rotation*expected[0],b=rotation*expected[1];
        require(((actual[0]-a).norm()<1e-14 && (actual[1]-b).norm()<1e-14) ||
          ((actual[0]-b).norm()<1e-14 && (actual[1]-a).norm()<1e-14),
          "quad diagonal is not rotation/renumbering covariant");
      }
  const std::array<Math::SpatialPoint,4> square{Math::SpatialPoint{1,1,0},
    Math::SpatialPoint{-1,1,0},Math::SpatialPoint{-1,-1,0},Math::SpatialPoint{1,-1,0}};
  require(rotationalQuadDiagonal(square,1e-11)==-1,
    "symmetric quad broke diagonal-exchange symmetry");
  bool invalidQuad=false;
  try { rotationalQuadDiagonal(quad,0); }
  catch(const std::runtime_error&) { invalidQuad=true; }
  require(invalidQuad,"quad diagonal accepted zero tolerance");
  invalidQuad=false;
  try { rotationalQuadDiagonal({quad[0],quad[0],quad[0],quad[0]},1e-11); }
  catch(const std::runtime_error&) { invalidQuad=true; }
  require(invalidQuad,"quad diagonal accepted degenerate geometry");
  // Regression: independent near-edge snaps must not zero both ends of an
  // originally crossed edge, even when its periodic partners are elsewhere.
  const std::vector<Real> original{-2e-4,6e-5,-2e-4,6e-5,1};
  std::vector<Real> snapped{0,0,0,0,1};
  const std::vector<Index> orbit{0,1,0,1,4};
  const std::set<std::array<Index,2>> edges{{0,1},{2,3}};
  require(guardCrossedEdgeSnaps(original,snapped,orbit,edges)==2,
    "conflicting crossed-edge snaps were not restored");
  require(snapped[0]==original[0] && snapped[2]==original[2] &&
    snapped[1]==0 && snapped[3]==0 && snapped[4]==1,
    "snap guard changed unrelated data or broke periodic orbit restoration");
  require(guardCrossedEdgeSnaps(original,snapped,orbit,edges)==0,
    "snap guard is not idempotent");
  std::vector<Real> sameSign{0,0};
  require(guardCrossedEdgeSnaps(std::vector<Real>{1,2},sameSign,
    std::vector<Index>{0,1},std::set<std::array<Index,2>>{{0,1}})==0,
    "snap guard restored an edge which was not originally crossed");
  bool rejected=false;
  try { guardCrossedEdgeSnaps(original,snapped,std::vector<Index>{0},edges); }
  catch (const std::runtime_error&) { rejected=true; }
  require(rejected,"snap guard accepted inconsistent data sizes");
  // Saved transport replays must not perturb coordinates during MEDIT input.
  const Math::SpatialPoint saved{0.9999999999999998, 1.0554945123781563,
    4.440892098500626e-16};
  std::ostringstream row;
  row << std::setprecision(std::numeric_limits<Real>::max_digits10)
      << saved(0) << ' ' << saved(1) << ' ' << saved(2) << " 13";
  const std::string text = row.str();
  const auto parsed = IO::MEDIT::ParseVertex(3)(text.begin(), text.end());
  require(bool(parsed), "MEDIT vertex parse failed");
  for (size_t component = 0; component < 3; ++component)
    require(parsed->vertex(component) == saved(component),
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
