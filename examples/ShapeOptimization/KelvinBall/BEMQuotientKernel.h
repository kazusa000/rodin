/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef KELVIN_BALL_BEM_QUOTIENT_KERNEL_H
#define KELVIN_BALL_BEM_QUOTIENT_KERNEL_H
#include "SewedOutput.h"
#include <Rodin/Variational/RegularizedStokeslet.h>
#include <Rodin/Math/SpatialTensor.h>
namespace KelvinBall
{
  /** Matrix-response single layer on one O24 chamber.
   * K_Q[F](x,y)=sum_g S_epsilon(x-g y) g F(y) g^T.
   * Both physical-vector and rigid-load indices rotate. Only image points,
   * never image meshes or independent image DOFs, are created.
   */
  class BEMQuotientKernel
  {
    public:
      static constexpr bool SmoothKernel=true;
      /// Construct with the fixed regularization and viscosity of the discrete model.
      BEMQuotientKernel(Real epsilon,Real viscosity):m_kernel(epsilon,viscosity) {}
      /// Fourth-order kernel acting on a 3x3 response at physical points x,y.
      void operator()(Math::SpatialTensor<Real,4>& out,const Point& x,const Point& y) const
      {
        out.resize({3,3,3,3}); out.setZero();
        for (const auto& g:SewedOutput::getCubeRotations())
        {
          Math::SpatialMatrix<Real> K(3,3);
          m_kernel.evaluate(K,Math::SpatialVector<Real>(x.vector()-g*y.vector()));
          const Math::SpatialMatrix<Real> Kg=K*g;
          for (size_t a=0;a<3;++a) for (size_t b=0;b<3;++b)
            for (size_t c=0;c<3;++c) for (size_t d=0;d<3;++d)
              out(a,b,c,d)+=Kg(a,c)*g(b,d);
        }
      }
      /// Same positive symmetric quadrature as the full-surface reference model.
      static const QF::QuadratureFormulaBase& getQuadrature(size_t order,Polytope::Type geometry)
      { return RegularizedStokeslet::getQuadrature(order,geometry); }
    private:
      RegularizedStokeslet m_kernel;
  };
}
#endif
