/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Metrics.h"

namespace KelvinBall
{
  Values Metrics::evaluateAxial(const VelocitySpace& Vh, const PressureSpace& Qh,
    const RotatedNitscheIntegrator* coupling, GridFunction<VelocitySpace, Math::Vector<Real>>& translation,
    GridFunction<VelocitySpace, Math::Vector<Real>>& rotation, GridFunction<PressureSpace, Math::Vector<Real>>& pT,
    GridFunction<PressureSpace, Math::Vector<Real>>& pR) const
  {
    if (coupling)
      checkPairs(Vh.getMesh(), *coupling);
    TrialFunction u(Vh);
    TrialFunction p(Qh);
    TestFunction v(Vh);
    TestFunction q(Qh);
    P0g gauge(Vh.getMesh());
    TrialFunction lambda(gauge);
    TestFunction eta(gauge);
    const auto Du = Real(0.5) * (Jacobian(u) + Jacobian(u).T());
    const auto Dv = Real(0.5) * (Jacobian(v) + Jacobian(v).T());
    const auto tau = m_parameters.cellStabilization();
    const auto axis = axialDirection();
    const auto solve = [&](const auto& rigid, auto& state, auto& pressure) {
      Problem stokes(u, p, lambda, v, q, eta);
      const auto trace = coupling
        ? coupling->stokesTrace(Tuple(std::cref(u), std::cref(p)),
            Tuple(std::cref(v), std::cref(q)), Mu, m_parameters.nitschePenalty,
            -m_parameters.pressureStabilization())
        : ProblemBody<Math::SparseMatrix<Real>, Math::Vector<Real>, Real>{};
      stokes = Integral(2 * Mu * Du, Dv) + trace.getBFs() - Integral(p, Div(v)) -
        Integral(Div(u), q) - Integral(tau * Grad(p), Grad(q)) +
        Integral(lambda, q) + Integral(p, eta) + DirichletBC(u, rigid).on(Gamma) +
        DirichletBC(u, VectorFunction{0, 0, 0}).on(Outer);
      stokes.assemble();
      solveDirect(stokes);
      const auto& system = stokes.getLinearSystem();
      const Real residual =
        (system.getOperator() * system.getSolution() - system.getVector()).norm() /
        std::max(Real(1), system.getVector().norm());
      if (!std::isfinite(residual) || residual > LinearResidualTolerance)
        throw std::runtime_error("Axial Stokes residual exceeds tolerance.");
      state = u.getSolution();
      pressure = p.getSolution();
    };
    const Real a = axis(0);
    solve(VectorFunction{a,a,a}, translation, pT);
    solve(VectorFunction{a * (F::z-F::y), a * (F::x-F::z), a * (F::y-F::x)}, rotation, pR);
    BilinearForm energy(u,v);
    energy = Integral(2 * Mu * Du,Dv);
    energy.assemble();
    const Real copies = coupling ? ChamberMultiplicity : 1;
    Values values;
    values.k = copies * energy(translation,translation);
    values.c = copies * energy(translation,rotation);
    values.q = copies * energy(rotation,rotation);
    if (!(values.k > 0 && values.q > 0 && std::isfinite(values.c)))
      throw std::runtime_error("Invalid axial resistance.");
    const Real scale = std::sqrt(values.k * values.q);
    values.rho = std::abs(values.c) / scale;
    values.couplingSymmetry = std::abs(values.c - copies * energy(rotation, translation)) / scale;
    if (coupling)
      values.nitscheJump = std::max(coupling->vectorJump(translation),coupling->vectorJump(rotation));
    return values;
  }

  Values Metrics::evaluateChamber(Mesh& fluid) const
  {
    prepare(fluid);
    VelocitySpace velocitySpace(fluid, 3);
    PressureSpace pressureSpace(fluid);
    RotatedNitscheIntegrator coupling(fluid, MasterCuts, 0.01, 0.25);
    GridFunction uT0(velocitySpace), uT1(velocitySpace), uT2(velocitySpace),
      uR0(velocitySpace), uR1(velocitySpace), uR2(velocitySpace);
    GridFunction pT0(pressureSpace), pT1(pressureSpace), pT2(pressureSpace),
      pR0(pressureSpace), pR1(pressureSpace), pR2(pressureSpace);
    if constexpr (AxialC3)
      return evaluateAxial(velocitySpace, pressureSpace, &coupling, uT0, uR0, pT0, pR0);
    return evaluateChamber(velocitySpace, pressureSpace, coupling, uT0, uT1, uT2, uR0,
      uR1, uR2, pT0, pT1, pT2, pR0, pR1, pR2);
  }

  Values Metrics::evaluateSewed(Mesh& fluid) const
  {
    prepare(fluid);
    VelocitySpace velocitySpace(fluid, 3);
    PressureSpace pressureSpace(fluid);
    GridFunction uT0(velocitySpace), uT1(velocitySpace), uT2(velocitySpace),
      uR0(velocitySpace), uR1(velocitySpace), uR2(velocitySpace);
    GridFunction pT0(pressureSpace), pT1(pressureSpace), pT2(pressureSpace),
      pR0(pressureSpace), pR1(pressureSpace), pR2(pressureSpace);
    if constexpr (AxialC3)
      return evaluateAxial(velocitySpace, pressureSpace, nullptr, uT0, uR0, pT0, pR0);
    solveFullFamily(velocitySpace, pressureSpace, VectorFunction{1, 0, 0},
      VectorFunction{0, 1, 0}, VectorFunction{0, 0, 1}, uT0, uT1, uT2, pT0, pT1, pT2);
    solveFullFamily(velocitySpace, pressureSpace, VectorFunction{0, -F::z, F::y},
      VectorFunction{F::z, 0, -F::x}, VectorFunction{-F::y, F::x, 0}, uR0, uR1, uR2, pR0,
      pR1, pR2);
    return resistance(velocitySpace, uT0, uT1, uT2, uR0, uR1, uR2, 1);
  }
}
