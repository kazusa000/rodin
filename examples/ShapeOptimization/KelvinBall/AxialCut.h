/*
 * Distributed under the Boost Software License, Version 1.0.
 * See the accompanying LICENSE file.
 */
#ifndef RODIN_KELVINBALL_AXIALCUT_H
#define RODIN_KELVINBALL_AXIALCUT_H

#include "Common.h"
#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <stdexcept>

namespace KelvinBall
{
  /// Match both faces once per mesh and enforce the same displacement on their nodes.
  class AxialCut
  {
    public:
      explicit AxialCut(Mesh& mesh)
      {
        using Key = std::array<long long, 3>;
        using Triangle = std::array<Index, 3>;
        const auto key = [](const Math::SpatialPoint& x) -> Key {
          return {std::llround(x(0) * 1e10), std::llround(x(1) * 1e10),
                  std::llround(x(2) * 1e10)};
        };
        std::map<Key, Index> master;
        std::map<Triangle, Attribute> triangles;
        std::vector<std::pair<Triangle, Attribute>> slaveFaces;
        mesh.getConnectivity().compute(2, 3);
        for (auto face = mesh.getBoundary(); face; ++face)
        {
          const auto& nodes = face->getVertices();
          const auto label = face->getAttribute();
          if (label == Outer)
            m_wall.insert(nodes.begin(), nodes.end());
          if (label != SigmaPlus && label != SigmaMinus)
            continue;
          Triangle triangle{nodes(0), nodes(1), nodes(2)};
          const auto& incident = mesh.getConnectivity().getIncidence({2, 3}, face->getIndex());
          const Attribute material = mesh.getCellAttribute(incident[0]).value();
          if (label == SigmaPlus)
          {
            for (Index i : triangle)
              master.emplace(key(mesh.getVertexCoordinates(i)), i);
            std::sort(triangle.begin(), triangle.end());
            triangles.emplace(triangle, material);
          }
          else
            slaveFaces.emplace_back(triangle, material);
        }
        if (triangles.empty() || triangles.size() != slaveFaces.size())
          throw std::runtime_error("The C3 cuts have different triangulations.");
        for (auto [triangle, material] : slaveFaces)
        {
          for (Index& i : triangle)
          {
            const auto position = master.find(key(m_rotation * mesh.getVertexCoordinates(i)));
            if (position == master.end())
              throw std::runtime_error("A C3 cut node has no rotated partner.");
            m_pairs.emplace(i, position->second);
            i = position->second;
          }
          std::sort(triangle.begin(), triangle.end());
          auto partner = triangles.find(triangle);
          if (partner == triangles.end() || partner->second != material)
            throw std::runtime_error("C3 cut triangles differ in topology or material.");
          triangles.erase(partner);
        }
      }

      template<class Field>
      void project(Field& field) const
      {
        const auto& space = field.getFiniteElementSpace();
        auto read = [&](Index vertex) {
          const auto dofs = space.getDOFs(0, vertex);
          return Math::SpatialVector<Real>{
            field.getData()(dofs(0)), field.getData()(dofs(1)), field.getData()(dofs(2))};
        };
        auto write = [&](Index vertex, const Math::SpatialVector<Real>& value) {
          const auto dofs = space.getDOFs(0, vertex);
          for (size_t j = 0; j < 3; ++j)
            field.getData()(dofs(j)) = value(j);
        };
        const auto axis = axialDirection();
        for (auto [slave, master] : m_pairs)
        {
          if (slave == master)
          {
            write(slave, axis.dot(read(slave)) * axis);
            continue;
          }
          Math::SpatialVector<Real> average = (read(slave) + m_rotation.transpose() * read(master)) * Real(0.5);
          average(0) = average(2) = (average(0) + average(2)) * Real(0.5);
          write(slave, average);
          write(master, m_rotation * average);
        }
        for (Index i : m_wall)
          write(i, Math::SpatialVector<Real>{0, 0, 0});
      }

    private:
      Math::SpatialMatrix<Real> m_rotation = RotationPairs.front().rotation;
      std::map<Index, Index> m_pairs;
      std::set<Index> m_wall;
  };
}
#endif
