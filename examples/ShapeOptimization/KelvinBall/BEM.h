#ifndef KELVIN_BALL_BEM_H
#define KELVIN_BALL_BEM_H

#include "Common.h"
#include "SewedOutput.h"
#include <Rodin/Variational/RegularizedStokeslet.h>
#include <array>

namespace KelvinBall
{
  /** Native P1 Galerkin resistance and exact fixed-topology coordinate derivative.
   * Local CPU/dense support only. No panel selection or geometry repair.
   */
  class BEM
  {
    public:
      /// Opt-in wall timings and raw-basis one-norm condition diagnostics.
      struct Diagnostics
      {
        std::array<Real,2> assemblySeconds{}, reductionSeconds{}, factorSeconds{}, solveSeconds{};
        Real basisSeconds=0, conditionSeconds=0, conditionOne=0, conditionOneEstimate=0;
        Real operatorNormOne=0, inverseNormOne=0, inverseResidual=0;
      };

      struct Result
      {
        Real k = 0, c = 0, q = 0, rho = 0, epsilon = 0;
        Real residual = 0, reciprocity = 0;
        size_t panels = 0, openEdges = 0, stateDOFs = 0, independentDOFs = 0;
        Math::Matrix<Real> resistance, gradient;
      };

      BEM(size_t order = 3, Real epsilonFactor = 0.6, Real viscosity = Mu)
        : m_order(order), m_epsilonFactor(epsilonFactor), m_viscosity(viscosity) {}

      Result evaluate(const Mesh& surface, bool derivative = true) const;
      Result evaluateQuotient(const Mesh& surface, bool derivative = true,
                              Diagnostics* diagnostics = nullptr) const;
      Result evaluateChamber(const Mesh& chamber, const VelocitySpace& shapeSpace,
                             Math::Vector<Real>& load, bool derivative = true) const;

    private:
      size_t m_order;
      Real m_epsilonFactor, m_viscosity;
  };
}
#endif
