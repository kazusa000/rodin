/* Distributed under the Boost Software License, Version 1.0. */
#ifndef KELVIN_BALL_PRESSURE_GAUGE_H
#define KELVIN_BALL_PRESSURE_GAUGE_H

#include "Common.h"

namespace KelvinBall
{
  /**
   * @brief One pressure reference per enclosed P1 fluid component.
   *
   * A component bounded only by the rigid Gamma surface has a constant
   * pressure nullspace. Components touching the outer box or rotated cuts
   * are not pinned: their pressure families may be coupled by the cuts.
   * Vertex CCL follows P1 connectivity, including cells sharing a vertex.
   * Requires vertex adjacency prepared by prepare(). No geometry is changed.
   */
  inline IndexSet enclosedPressureReferences(const PressureSpace& space)
  {
    const auto& mesh = space.getMesh();
    std::vector<Boolean> excluded(mesh.getVertexCount(), false);
    const auto& connectivity = mesh.getConnectivity();
    for (auto face = mesh.getBoundary(); face; ++face)
    {
      if (face->getAttribute() == Gamma)
        continue;
      for (Index vertex : connectivity.getIncidence({2, 0}, face->getIndex()))
        excluded[vertex] = true;
    }
    IndexSet references;
    const auto components = mesh.ccl(0,
      [](const Polytope&, const Polytope&) { return true; });
    for (const auto& component : components.getComponents())
    {
      bool enclosed = !component.empty();
      for (Index vertex : component)
        enclosed = enclosed && !excluded[vertex];
      if (enclosed)
        references.insert(space.getDOFs(0, *component.begin())(0));
    }
    return references;
  }
}

#endif
