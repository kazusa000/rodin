#ifndef KELVIN_BALL_PERIODIC_CUTS_H
#define KELVIN_BALL_PERIODIC_CUTS_H

#include <map>
#include <numeric>
#include <set>
#include <Eigen/IterativeLinearSolvers>
#include "Common.h"
#include "SewedOutput.h"

namespace KelvinBall
{
  /// Reconstruct semantic Gamma labels from the actual material adjacency.
  /// MMG's volume conversion exports external boundary triangles, not these
  /// internal faces, when face/cell incidence is present. This changes labels
  /// only: it neither moves coordinates nor adds/removes any material face.
  inline void labelMaterialInterface(Mesh& mesh)
  {
    mesh.getConnectivity().compute(2,3);
    const auto& owners = mesh.getConnectivity().getIncidence(2,3);
    for (auto face=mesh.getPolytope(2);face;++face)
    {
      bool solid=false,fluid=false;
      for (Index owner:owners.at(face->getIndex()))
      {
        const auto material=mesh.getPolytope(3,owner)->getAttribute();
        solid=solid || material==Obstacle;
        fluid=fluid || material==Fluid;
      }
      if (solid && fluid) mesh.setAttribute({2,face->getIndex()},Gamma);
      else if (face->getAttribute()==Gamma) mesh.setAttribute({2,face->getIndex()},{});
    }
  }

  /** Exact P1 scalar identification on congruent rotated cuts, Local Eigen only.
   * No point is moved. Noncongruent triangulations are rejected, not welded.
   * The tolerance is for floating-point coordinate matching, not seam repair.
   */
  class PeriodicCuts
  {
    using Key = std::array<long long, 3>;
    using FaceKey = std::array<Key, 3>;
    public:
      explicit PeriodicCuts(const Mesh& mesh) : m_mesh(mesh), m_roots(mesh.getVertexCount())
      {
        std::iota(m_roots.begin(), m_roots.end(), Index(0));
        Real scale = 1;
        for (Index i = 0; i < mesh.getVertexCount(); ++i)
          scale = std::max(scale, mesh.getVertexCoordinates(i).norm());
        m_tolerance = Real(1e-11) * scale;
        for (const auto& pair : RotationPairs)
        {
          std::map<Key, Index> master;
          std::set<FaceKey> masterFaces, slaveFaces;
          for (auto face = mesh.getPolytope(2); face; ++face)
          {
            if (face->getAttribute() != pair.master) continue;
            FaceKey triangle;
            for (size_t k = 0; k < 3; ++k)
            {
              const Index i = face->getVertices()[k];
              const Key point = key(mesh.getVertexCoordinates(i));
              master.emplace(point, i);
              triangle[k] = point;
            }
            std::sort(triangle.begin(), triangle.end());
            masterFaces.insert(triangle);
          }
          for (auto face = mesh.getPolytope(2); face; ++face)
          {
            if (face->getAttribute() != pair.slave) continue;
            FaceKey triangle;
            for (size_t k = 0; k < 3; ++k)
            {
              const Index i = face->getVertices()[k];
              const Math::SpatialPoint x = pair.rotation * mesh.getVertexCoordinates(i);
              const Key point = key(x);
              const auto found = master.find(point);
              if (found == master.end() ||
                  (mesh.getVertexCoordinates(found->second) - x).norm() > m_tolerance)
                throw std::runtime_error("Periodic cut has an unpaired vertex.");
              const Index a = root(i), b = root(found->second);
              m_roots[std::max(a, b)] = std::min(a, b);
              triangle[k] = point;
            }
            std::sort(triangle.begin(), triangle.end());
            slaveFaces.insert(triangle);
          }
          if (masterFaces.empty() || slaveFaces != masterFaces)
            throw std::runtime_error("Periodic cut triangulations are not rotationally congruent: master=" +
              std::to_string(pair.master) + " faces=" + std::to_string(masterFaces.size()) +
              ", slave=" + std::to_string(pair.slave) + " faces=" + std::to_string(slaveFaces.size()));
        }
        for (Index i = 0; i < m_roots.size(); ++i) m_roots[i] = root(i);
      }

      template <class Field>
      void project(Field& field) const
      {
        const auto& space = field.getFiniteElementSpace();
        std::vector<Real> sum(m_roots.size(), 0);
        std::vector<size_t> count(m_roots.size(), 0);
        for (Index i = 0; i < m_roots.size(); ++i)
        {
          sum[m_roots[i]] += field.getData()(space.getDOFs(0, i)(0));
          ++count[m_roots[i]];
        }
        for (Index i = 0; i < m_roots.size(); ++i)
          field.getData()(space.getDOFs(0, i)(0)) = sum[m_roots[i]] / count[m_roots[i]];
      }

      template <class Space, class System>
      Real solve(const Space& space, System& system) const
      {
        using Sparse = Eigen::SparseMatrix<Real>;
        std::map<Index, Index> columns;
        std::vector<Eigen::Triplet<Real>> entries;
        for (Index i = 0; i < m_roots.size(); ++i)
        {
          const auto [it, inserted] = columns.emplace(m_roots[i], columns.size());
          entries.emplace_back(space.getDOFs(0, i)(0), it->second, Real(1));
        }
        Sparse expansion(space.getSize(), columns.size());
        expansion.setFromTriplets(entries.begin(), entries.end());
        Sparse reduced = expansion.transpose() * system.getOperator() * expansion;
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> rhs = expansion.transpose() * system.getVector();
        Eigen::ConjugateGradient<Sparse, Eigen::Lower | Eigen::Upper> solver;
        solver.setTolerance(Real(1e-13));
        solver.setMaxIterations(10000);
        solver.compute(reduced);
        const Eigen::Matrix<Real, Eigen::Dynamic, 1> solution = solver.solve(rhs);
        if (solver.info() != Eigen::Success ||
            (reduced * solution - rhs).norm() / std::max(rhs.norm(), Real(1)) > Real(1e-10))
          throw std::runtime_error("Exact periodic scalar solve failed its reduced residual.");
        system.getSolution() = expansion * solution;
        return (reduced * solution - rhs).norm() / std::max(rhs.norm(), Real(1));
      }

      void checkClosedInterface() const
      {
        std::map<std::array<Key, 2>, size_t> edges;
        std::set<FaceKey> triangles;
        for (const auto& rotation : SewedOutput::getCubeRotations())
          for (auto face = m_mesh.getPolytope(2); face; ++face)
          {
            if (face->getAttribute() != Gamma) continue;
            FaceKey triangle;
            for (size_t k = 0; k < 3; ++k)
              triangle[k] = key(rotation * m_mesh.getVertexCoordinates(face->getVertices()[k]));
            std::sort(triangle.begin(), triangle.end());
            if (triangle[0] == triangle[1] || triangle[1] == triangle[2] ||
                !triangles.insert(triangle).second)
              throw std::runtime_error("Periodic interface has a degenerate or duplicate face.");
            for (size_t k = 0; k < 3; ++k)
            {
              std::array<Key, 2> edge{triangle[k], triangle[(k + 1) % 3]};
              std::sort(edge.begin(), edge.end());
              ++edges[edge];
            }
          }
        if (edges.empty()) throw std::runtime_error("Periodic interface is empty.");
        for (const auto& [edge, count] : edges)
          if (count != 2)
            throw std::runtime_error("Periodic interface is open or nonmanifold.");
      }

    private:
      Key key(const Math::SpatialPoint& x) const
      {
        return {std::llround(x(0) / m_tolerance), std::llround(x(1) / m_tolerance),
                std::llround(x(2) / m_tolerance)};
      }
      Index root(Index i) const
      {
        while (m_roots[i] != i) i = m_roots[i];
        return i;
      }
      const Mesh& m_mesh;
      std::vector<Index> m_roots;
      Real m_tolerance;
  };
}
#endif
