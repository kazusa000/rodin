#ifndef KELVIN_CROSSED_EDGE_SNAP_GUARD_H
#define KELVIN_CROSSED_EDGE_SNAP_GUARD_H

#include <array>
#include <cmath>
#include <set>
#include <stdexcept>
#include <vector>

namespace KelvinBall
{
  /** Undo conflicting zero snaps, not a repair of the resulting interface.
   * An originally crossed edge must not have both endpoints snapped to zero:
   * this can collapse distinct zero contours into a nonmanifold zero edge.
   * Restore the endpoint further from zero, together with its periodic orbit.
   * No coordinate, cell or original scalar value is changed. The full closure
   * check remains mandatory; this guard does not claim to cover all topology.
   */
  template<class Scalar, class Index>
  size_t guardCrossedEdgeSnaps(const std::vector<Scalar>& original,
      std::vector<Scalar>& snapped, const std::vector<Index>& representatives,
      const std::set<std::array<Index, 2>>& edges)
  {
    if (original.size()!=snapped.size() || original.size()!=representatives.size())
      throw std::runtime_error("Crossed-edge snap guard has mismatched input sizes.");
    std::set<Index> restore;
    for (const auto& edge:edges)
    {
      const Index a=edge[0], b=edge[1];
      if (a>=original.size() || b>=original.size())
        throw std::runtime_error("Crossed-edge snap guard has invalid connectivity.");
      const bool crossed=(original[a]<0 && original[b]>0) ||
                         (original[a]>0 && original[b]<0);
      if (crossed && snapped[a]==0 && snapped[b]==0)
        restore.insert(representatives[std::abs(original[a])>=std::abs(original[b])?a:b]);
    }
    size_t count=0;
    for (Index i=0;i<original.size();++i)
      if (restore.count(representatives[i]))
      {
        count+=snapped[i]!=original[i];
        snapped[i]=original[i];
      }
    return count;
  }
}
#endif
