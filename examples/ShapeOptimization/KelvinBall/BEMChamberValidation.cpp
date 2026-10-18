/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "BEM.h"
#include "Sphere.h"
#include <Rodin/MMG/MeshLoader.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace KelvinBall;
int main(int argc,char** argv)
{
  try
  {
    Configuration configuration;
    configuration.points=7;configuration.hminFactor=.8;configuration.hmaxFactor=1.2;
    configuration.seedAmplitude=.6;configuration.finalize();
    MMG::Mesh mesh;
    if (argc==2)
    {
      std::ifstream input(argv[1]);if (!input) throw std::runtime_error("Cannot read chamber mesh.");
      MMG::MeshLoader(mesh).load(input);
    }
    else mesh=Sphere(configuration).discretize(true,0,true).mesh;
    mesh.getConnectivity().compute(2,3);
    VelocitySpace space(mesh,3);Math::Vector<Real> load;
    const auto result=BEM().evaluateChamber(mesh,space,load,true);
    GridFunction direction(space);
    direction=VectorFunction(size_t(3),[](const Point& p) {
      Math::SpatialVector<Real> x=p.vector();
      const Real r=x.norm();if (r==0) return x;
      const Real a=x(0)/r,b=x(1)/r,c=x(2)/r;
      return Math::SpatialVector<Real>(x*(27*a*b*c*(a*a-b*b)*(b*b-c*c)*(c*c-a*a)));
    });
    const Real analytic=load.dot(direction.getData());
    std::cout<<std::setprecision(17)<<"state rho="<<result.rho<<" residual="<<result.residual
      <<" open_edges="<<result.openEdges<<" chamber_panels="<<result.panels
      <<" state_dofs="<<result.stateDOFs<<" independent_dofs="<<result.independentDOFs<<std::endl;
    Real finalError=1;
    for (Real step:{Real(1e-3),Real(5e-4),Real(2.5e-4)})
    {
      Real values[2];size_t side=0;
      for (Real sign:{Real(-1),Real(1)})
      {
        MMG::Mesh perturbed(mesh);
        for (Index v=0;v<mesh.getVertexCount();++v)
        {
          Math::SpatialVector<Real> movement(3);const auto& dofs=space.getDOFs(0,v);
          for (size_t j=0;j<3;++j) movement(j)=direction.getData()(dofs(j));
          perturbed.setVertexCoordinates(v,mesh.getVertexCoordinates(v)+sign*step*movement);
        }
        VelocitySpace movedSpace(perturbed,3);Math::Vector<Real> unused;
        values[side++]=BEM().evaluateChamber(perturbed,movedSpace,unused,false).rho;
      }
      const Real fd=(values[1]-values[0])/(2*step);
      finalError=std::abs(fd-analytic)/std::max(std::abs(analytic),Real(1e-20));
      std::cout<<"finite_difference step="<<step<<" analytic="<<analytic<<" fd="<<fd
        <<" relative_error="<<finalError<<std::endl;
    }
    std::ofstream record("bem-chamber-validation.json");
    const bool passed=result.residual<=1e-9 && result.openEdges==0 && finalError<.002;
    record<<std::setprecision(17)<<"{\"passed\":"<<(passed ? "true":"false")
      <<",\"design_domain\":\"O24-chamber\",\"bem_domain\":\"infinite-exterior\""
      <<",\"rho\":"<<result.rho<<",\"residual\":"<<result.residual
      <<",\"chamber_panels\":"<<result.panels<<",\"state_dofs\":"<<result.stateDOFs
      <<",\"independent_dofs\":"<<result.independentDOFs<<",\"open_edges\":"<<result.openEdges
      <<",\"analytic_derivative\":"<<analytic<<",\"fd_relative_error\":"<<finalError<<"}\n";
    if (!record) throw std::runtime_error("Cannot save chamber derivative validation.");
    return passed ? 0 : 2;
  }
  catch (const std::exception& error) { std::cerr<<error.what()<<std::endl;return 1; }
}
