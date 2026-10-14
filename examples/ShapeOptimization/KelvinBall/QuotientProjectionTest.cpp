#include "Sphere.h"
#include "SewedOutput.h"
#include "RotatedNitscheIntegrator.h"
#include <iostream>

using namespace KelvinBall;

int main()
{
  try
  {
    Configuration configuration;
    configuration.points = 5;
    configuration.finalize();
    auto mesh = Sphere(configuration).makeUniformChamber();
    P1 scalarSpace(mesh);
    VelocitySpace vectorSpace(mesh, 3);
    GridFunction scalar(scalarSpace), vector(vectorSpace);
    for (Eigen::Index i = 0; i < scalar.getData().size(); ++i)
      scalar.getData()(i) = std::sin(Real(i + 1) * .73);
    for (Eigen::Index i = 0; i < vector.getData().size(); ++i)
      vector.getData()(i) = std::sin(Real(i + 1) * .51);
    const SewedOutput quotient(mesh);
    VelocitySpace fullSpace(quotient.getMesh(), 3);
    GridFunction fullVector(fullSpace);
    quotient.setVector(fullVector, vector);
    Math::Matrix<Real> fullGradient(quotient.getMesh().getVertexCount(), 3);
    for (Index i = 0; i < quotient.getMesh().getVertexCount(); ++i)
      for (size_t j = 0; j < 3; ++j) fullGradient(i, j) = std::cos(Real(3*i+j+1)*.37);
    const auto load = quotient.pullBackCoordinateGradient(fullGradient, vectorSpace);
    const Real left = load.dot(vector.getData());
    Real right = 0;
    for (Index i = 0; i < quotient.getMesh().getVertexCount(); ++i)
    {
      const auto& dofs = fullSpace.getDOFs(0, i);
      for (size_t j = 0; j < 3; ++j) right += fullGradient(i, j)*fullVector.getData()(dofs(j));
    }
    if (std::abs(left-right) > 1e-11*std::max(Real(1),std::abs(right)))
      throw std::runtime_error("BEM coordinate pullback is not the adjoint of averaged transport.");
    quotient.enforceScalarRepresentatives(scalar);
    quotient.enforceVectorRepresentatives(vector);
    const Math::Vector<Real> firstScalar = scalar.getData(), firstVector = vector.getData();
    quotient.enforceScalarRepresentatives(scalar);
    quotient.enforceVectorRepresentatives(vector);
    if ((firstScalar-scalar.getData()).norm() > 1e-12 ||
        (firstVector-vector.getData()).norm() > 1e-12)
      throw std::runtime_error("Quotient trace projection is not idempotent.");
    RotatedNitscheIntegrator coupling(mesh,
      FlatSet<Attribute>{SigmaPlus, SigmaMinus, SigmaXYPlus, SigmaXYMinus}, .01, .25);
    const Real scalarJump = coupling.scalarJump(scalar), vectorJump = coupling.vectorJump(vector);
    if (scalarJump > 1e-11 || vectorJump > 1e-11)
      throw std::runtime_error("Projected uniform-chamber traces do not match.");
    std::cout << "adjoint_error=" << std::abs(left-right) << " scalar_jump=" << scalarJump
      << " vector_jump=" << vectorJump << '\n';
    return 0;
  }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
