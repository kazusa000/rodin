/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_CHAMBER_SYMMETRY_H
#define KELVIN_BALL_CHAMBER_SYMMETRY_H

#include "Common.h"
#include "SewedOutput.h"
#include <Eigen/Eigenvalues>
#include <map>
#include <tuple>

namespace KelvinBall
{
  /** O24 orbit identifications on chamber coordinates, without a full-volume mesh.
   * Architecture: canonicalize each vertex orbit, project its stabilizer, then
   * transfer scalar/vector traces or form the admissible tensor-state basis.
   * Coordinate quantization is the existing sewn-topology tolerance, not a repair.
   */
  class ChamberSymmetry
  {
    public:
      using Key = std::tuple<long long, long long, long long>;
      static constexpr Real tolerance = 1e-10;
      struct Source { Index vertex; Math::SpatialMatrix<Real> rotation; };
      struct Orbit { Math::SpatialPoint point; std::vector<Source> sources; };

      /// Canonical coordinate key; x is a physical point.
      static Key key(const Math::SpatialPoint& x)
      {
        return {std::llround(x(0)/tolerance),std::llround(x(1)/tolerance),
          std::llround(x(2)/tolerance)};
      }

      /// Group coordinates of the supplied chamber mesh; no mesh copies are made.
      explicit ChamberSymmetry(const Mesh& mesh)
      {
        for (Index v=0; v<mesh.getVertexCount(); ++v)
        {
          const auto& x=mesh.getVertexCoordinates(v);
          auto best=key(x);
          Math::SpatialMatrix<Real> transform(3,3); transform.setIdentity();
          Math::SpatialPoint canonical=x;
          for (const auto& g:SewedOutput::getCubeRotations())
          {
            const Math::SpatialPoint y=g*x;
            if (key(y)<best) { best=key(y); transform=g; canonical=y; }
          }
          auto [it,inserted]=m_orbits.try_emplace(best,Orbit{canonical,{}});
          it->second.sources.push_back({v,transform});
        }
      }

      /// Projects a chamber scalar trace onto coincident orbit representatives.
      template <class Field> void scalar(Field& field) const
      {
        const auto& space=field.getFiniteElementSpace();
        for (const auto& [key,orbit]:m_orbits)
        {
          Real value=0;
          for (const auto& s:orbit.sources) value+=field.getData()(space.getDOFs(0,s.vertex)(0));
          value/=orbit.sources.size();
          for (const auto& s:orbit.sources) field.getData()(space.getDOFs(0,s.vertex)(0))=value;
        }
      }

      /// Projects a chamber vector trace, including axes with nontrivial stabilizers.
      template <class Field> void vector(Field& field) const
      {
        const auto& space=field.getFiniteElementSpace();
        for (const auto& [key,orbit]:m_orbits)
        {
          Math::SpatialVector<Real> value(3); value.setZero();
          for (const auto& s:orbit.sources)
          {
            Math::SpatialVector<Real> local(3); const auto& dofs=space.getDOFs(0,s.vertex);
            for (size_t j=0;j<3;++j) local(j)=field.getData()(dofs(j));
            value+=s.rotation*local;
          }
          value/=orbit.sources.size();
          Math::SpatialVector<Real> projected(3); projected.setZero(); size_t count=0;
          for (const auto& g:SewedOutput::getCubeRotations())
            if (ChamberSymmetry::key(Math::SpatialPoint(g*orbit.point))==key)
            { projected+=g*value; ++count; }
          projected/=count;
          for (const auto& s:orbit.sources)
          {
            const auto local=s.rotation.transpose()*projected;
            const auto& dofs=space.getDOFs(0,s.vertex);
            for (size_t j=0;j<3;++j) field.getData()(dofs(j))=local(j);
          }
        }
      }

      /// Basis E for F(gx)=g F(x) g^T on a matrix-valued chamber P1 space.
      template <class Space> Math::Matrix<Real> tensorBasis(const Space& space) const
      {
        Math::Matrix<Real> result= Math::Matrix<Real>::Zero(space.getSize(),space.getSize());
        size_t column=0;
        for (const auto& [key,orbit]:m_orbits)
        {
          Math::Matrix<Real> projector=Math::Matrix<Real>::Zero(9,9); size_t count=0;
          for (const auto& g:SewedOutput::getCubeRotations())
            if (ChamberSymmetry::key(Math::SpatialPoint(g*orbit.point))==key)
            {
              for (size_t a=0;a<3;++a) for (size_t b=0;b<3;++b)
                for (size_t c=0;c<3;++c) for (size_t d=0;d<3;++d)
                  projector(3*a+b,3*c+d)+=g(a,c)*g(b,d);
              ++count;
            }
          projector/=count;
          Eigen::SelfAdjointEigenSolver<Math::Matrix<Real>> eig(projector);
          if (eig.info()!=Eigen::Success) throw std::runtime_error("O24 stabilizer basis failed.");
          for (size_t k=0;k<9;++k)
          {
            if (eig.eigenvalues()(k)<Real(.5)) continue;
            for (const auto& s:orbit.sources)
            {
              const auto& dofs=space.getDOFs(0,s.vertex);
              for (size_t a=0;a<3;++a) for (size_t b=0;b<3;++b)
                for (size_t c=0;c<3;++c) for (size_t d=0;d<3;++d)
                  result(dofs(3*a+b),column)+=s.rotation(c,a)*s.rotation(d,b)*
                    eig.eigenvectors()(3*c+d,k);
            }
            ++column;
          }
        }
        return result.leftCols(column);
      }
    private:
      std::map<Key,Orbit> m_orbits;
  };
}
#endif
