#ifndef RODIN_VARIATIONAL_REGULARIZEDSTOKESLET_H
#define RODIN_VARIATIONAL_REGULARIZEDSTOKESLET_H

#include <cmath>
#include <stdexcept>
#include "Rodin/Geometry/Point.h"
#include "Rodin/Math/SpatialMatrix.h"

namespace Rodin::Variational
{
  /** Smooth three-dimensional single-layer kernel, including 1/(8 pi mu).
   * The declared smoothness opts Potential into tensor surface quadrature.
   * This is a regularized model, not a singular-kernel quadrature shortcut.
   */
  class RegularizedStokeslet
  {
    public:
      static constexpr bool SmoothKernel = true;

      RegularizedStokeslet(Real epsilon, Real viscosity)
        : m_epsilon(epsilon), m_viscosity(viscosity)
      {
        if (!std::isfinite(epsilon) || epsilon <= 0 ||
            !std::isfinite(viscosity) || viscosity <= 0)
          throw std::invalid_argument("Positive finite Stokeslet epsilon and viscosity required.");
      }

      void operator()(Math::SpatialMatrix<Real>& out,
        const Geometry::Point& x, const Geometry::Point& y) const
      {
        evaluate(out, x.vector() - y.vector());
      }

      void evaluate(Math::SpatialMatrix<Real>& out,
        const Math::SpatialVector<Real>& r) const
      {
        if (r.size() != 3) throw std::invalid_argument("Stokeslet requires spatial dimension three.");
        const Real r2 = r.squaredNorm(), e2 = m_epsilon * m_epsilon;
        const Real pi = std::acos(Real(-1));
        out.resize(3, 3);
        out = ((r2 + 2 * e2) * Math::SpatialMatrix<Real>::Identity(3, 3) +
               r * r.transpose()) / (8 * pi * m_viscosity * std::pow(r2 + e2, Real(1.5)));
      }

      Real getEpsilon() const { return m_epsilon; }
      Real getViscosity() const { return m_viscosity; }

      /** Pull back H:dS to separation and epsilon, with no coordinate Jacobian stored. */
      void differential(const Math::SpatialVector<Real>& r,
        const Math::SpatialMatrix<Real>& H, Math::SpatialVector<Real>& dr, Real& de) const
      {
        const Real e = m_epsilon, r2 = r.squaredNorm(), d = r2 + e * e;
        const Real inv3 = std::pow(d, Real(-1.5)), inv5 = inv3 / d;
        const Real factor = 1 / (8 * std::acos(Real(-1)) * m_viscosity);
        const Real trace = H.trace();
        const Real contraction = (r2 + 2 * e * e) * trace + r.dot(H * r);
        dr = factor * ((2 * trace * r + H * r + H.transpose() * r) * inv3 -
                       3 * contraction * r * inv5);
        de = factor * e * (4 * trace * inv3 - 3 * contraction * inv5);
      }

    private:
      Real m_epsilon, m_viscosity;
  };
}
#endif
