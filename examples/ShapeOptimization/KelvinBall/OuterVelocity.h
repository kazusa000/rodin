/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_OUTER_VELOCITY_H
#define KELVIN_BALL_OUTER_VELOCITY_H

#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include "Common.h"

namespace KelvinBall
{
  /** Six prescribed outer traces for a fixed P1 mesh, not an optimizer policy.
   * Rows contain x,y,z then T0,T1,T2,R0,R1,R2 (three components each).
   * Coordinate lookup fails on missing/re-meshed nodes; no nearest extrapolation.
   */
  class OuterVelocity
  {
    public:
      explicit OuterVelocity(const std::string& file = "")
        : m_sphere(file == "analytic-sphere"), m_zero(file.empty())
      {
        if (m_sphere || m_zero)
          return;
        std::ifstream input(file);
        size_t count = 0;
        if (!(input >> count) || count == 0)
          throw std::runtime_error("Invalid outer velocity table header.");
        for (size_t i = 0; i < count; ++i)
        {
          Row row;
          for (Real& value : row)
            if (!(input >> value) || !std::isfinite(value))
              throw std::runtime_error("Invalid outer velocity table row.");
          const Key k = key(row[0], row[1], row[2]);
          if (!m_rows.emplace(k, row).second)
            throw std::runtime_error("Duplicate outer velocity coordinate.");
        }
        std::string extra;
        if (input >> extra)
          throw std::runtime_error("Trailing outer velocity table data.");
      }

      Math::SpatialVector<Real> value(const Geometry::Point& p, size_t load) const
      {
        return value(p.getPhysicalCoordinates(), load);
      }

      Math::SpatialVector<Real> value(const Math::SpatialPoint& x, size_t load) const
      {
        if (load >= 6)
          throw std::runtime_error("Invalid rigid load index.");
        Math::SpatialVector<Real> result(3);
        result.setZero();
        if (m_zero)
          return result;
        if (m_sphere)
        {
          const Real r = x.norm();
          if (r < 1 - 1e-10)
            throw std::runtime_error("Sphere exterior trace requested inside sphere.");
          Math::SpatialVector<Real> e(3);
          e.setZero();
          e(load % 3) = 1;
          if (load < 3)
            result = (Real(3)/(4*r) + Real(1)/(4*r*r*r))*e
              + (Real(3)/(4*r) - Real(3)/(4*r*r*r))*x*(x(load)/(r*r));
          else
          {
            result(0) = e(1)*x(2) - e(2)*x(1);
            result(1) = e(2)*x(0) - e(0)*x(2);
            result(2) = e(0)*x(1) - e(1)*x(0);
            result /= r*r*r;
          }
          return result;
        }
        const auto found = m_rows.find(key(x(0), x(1), x(2)));
        if (found == m_rows.end())
          throw std::runtime_error("Outer velocity table has no matching mesh node.");
        const auto& row = found->second;
        for (size_t j = 0; j < 3; ++j)
        {
          if (std::abs(row[j] - x(j)) > 1e-9)
            throw std::runtime_error("Outer velocity coordinate mismatch.");
          result(j) = row[3 + 3*load + j];
        }
        return result;
      }

    private:
      using Key = std::array<long long, 3>;
      using Row = std::array<Real, 21>;
      static Key key(Real x, Real y, Real z)
      {
        return {std::llround(x*1e9), std::llround(y*1e9), std::llround(z*1e9)};
      }
      bool m_sphere;
      bool m_zero;
      std::map<Key, Row> m_rows;
  };
}
#endif
