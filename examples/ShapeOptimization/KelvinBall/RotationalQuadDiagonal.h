/* Distributed under the Boost Software License, Version 1.0. */
#ifndef KELVIN_ROTATIONAL_QUAD_DIAGONAL_H
#define KELVIN_ROTATIONAL_QUAD_DIAGONAL_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include "Common.h"
#include "SewedOutput.h"

namespace KelvinBall
{
  /** Select a diagonal of a cyclic convex planar quadrilateral without moving
   * its vertices. Return 0 for (0,2), 1 for (1,3), or -1 for the old centre fan
   * when a group symmetry exchanges the diagonals. Maximise the smaller scaled
   * triangle quality 4 sqrt(3) A / sum(edge^2). Decision-only coordinate keys
   * use the existing periodic matching precision, canonicalised over the same
   * proper rotation group, so rotations and vertex renumbering select the same
   * geometric subdivision. The real, unrounded vertices define every cell.
   */
  inline int rotationalQuadDiagonal(const std::array<Math::SpatialPoint,4>& points,
      Real tolerance)
  {
    if (!(tolerance>0) || !std::isfinite(tolerance))
      throw std::runtime_error("Quad diagonal requires a finite positive matching tolerance.");
    using Key=std::array<long long,3>;
    const auto key=[&](const Math::SpatialPoint& x) {
      Key result;
      for (size_t j=0;j<3;++j) {
        const Real scaled=x(j)/tolerance;
        if (!std::isfinite(scaled) || std::abs(scaled)>Real(std::numeric_limits<long long>::max()/4))
          throw std::runtime_error("Quad diagonal has invalid coordinates.");
        result[j]=std::llround(scaled);
      }
      return result;
    };
    const auto triangleQuality=[&](std::array<size_t,3> ids) {
      std::array<Key,3> canonical;
      bool first=true;
      for (const auto& rotation:SewedOutput::getCubeRotations()) {
        std::array<Key,3> candidate;
        for (size_t i=0;i<3;++i) candidate[i]=key(rotation*points[ids[i]]);
        std::sort(candidate.begin(),candidate.end());
        if (first || candidate<canonical) {canonical=candidate;first=false;}
      }
      Math::SpatialPoint u(3),v(3);
      for (size_t j=0;j<3;++j) {
        u(j)=Real(canonical[1][j]-canonical[0][j]);
        v(j)=Real(canonical[2][j]-canonical[0][j]);
      }
      const Real length=u.squaredNorm()+v.squaredNorm()+(u-v).squaredNorm();
      return length>0 ? 2*std::sqrt(Real(3))*u.cross(v).norm()/length : Real(0);
    };
    const auto diagonalKey=[&](size_t a,size_t b) {
      std::array<Key,2> canonical;
      bool first=true;
      for (const auto& rotation:SewedOutput::getCubeRotations()) {
        std::array<Key,2> candidate{key(rotation*points[a]),key(rotation*points[b])};
        std::sort(candidate.begin(),candidate.end());
        if (first || candidate<canonical) {canonical=candidate;first=false;}
      }
      return canonical;
    };
    const auto a=diagonalKey(0,2),b=diagonalKey(1,3);
    const Real q0=std::min(triangleQuality({0,1,2}),triangleQuality({0,2,3}));
    const Real q1=std::min(triangleQuality({0,1,3}),triangleQuality({1,2,3}));
    if (!(std::max(q0,q1)>0))
      throw std::runtime_error("Quad diagonal has a degenerate polygon.");
    if (a==b) return -1;
    return q0==q1 ? (a<b?0:1) : (q0>q1?0:1);
  }
}
#endif
