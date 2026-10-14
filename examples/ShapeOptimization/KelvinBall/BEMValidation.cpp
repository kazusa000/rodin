#include "BEM.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

using namespace KelvinBall;

namespace
{
  using Vertices = std::vector<Math::SpatialPoint>;
  using Faces = std::vector<std::array<Index, 3>>;

  KelvinBall::Mesh surface(const Vertices& vertices, const Faces& faces)
  {
    KelvinBall::Mesh::Builder builder;
    builder.initialize(3).nodes(vertices.size());
    for (const auto& x : vertices) builder.vertex(x);
    for (const auto& face : faces)
    {
      IndexArray indices(3);
      for (size_t i = 0; i < 3; ++i) indices(i) = face[i];
      Index index;
      builder.polytope(Geometry::Polytope::Type::Triangle, std::move(indices), index);
    }
    return builder.finalize();
  }

  std::pair<Vertices, Faces> sphere(size_t refinement)
  {
    Vertices vertices;
    for (size_t axis = 0; axis < 3; ++axis)
      for (Real sign : {Real(1), Real(-1)})
      {
        Math::SpatialPoint x(3); x.setZero(); x(axis) = sign;
        vertices.push_back(x);
      }
    Faces faces;
    for (Index x : {0, 1}) for (Index y : {2, 3}) for (Index z : {4, 5})
    {
      std::array<Index, 3> f{x, y, z};
      if (vertices[x].dot(vertices[y].cross(vertices[z])) < 0) std::swap(f[1], f[2]);
      faces.push_back(f);
    }
    for (size_t level = 0; level < refinement; ++level)
    {
      std::map<std::pair<Index, Index>, Index> midpoints;
      auto midpoint = [&](Index a, Index b) {
        auto [position, inserted] = midpoints.emplace(std::minmax(a, b), vertices.size());
        if (inserted)
        {
          Math::SpatialPoint x = vertices[a] + vertices[b]; x.normalize();
          vertices.push_back(x);
        }
        return position->second;
      };
      Faces refined;
      for (const auto& f : faces)
      {
        const Index a = midpoint(f[0], f[1]), b = midpoint(f[1], f[2]), c = midpoint(f[2], f[0]);
        refined.push_back({f[0], a, c}); refined.push_back({a, f[1], b});
        refined.push_back({c, b, f[2]}); refined.push_back({a, b, c});
      }
      faces = std::move(refined);
    }
    return {vertices, faces};
  }

  void require(bool condition, const char* message)
  {
    if (!condition) throw std::runtime_error(message);
  }

  void writeSurface(const Vertices& vertices, const Faces& faces, const char* path)
  {
    std::ofstream file(path); file << std::setprecision(17);
    file << "MeshVersionFormatted 2\nDimension 3\nVertices\n" << vertices.size() << '\n';
    for (const auto& x : vertices) file << x.transpose() << " 0\n";
    file << "Triangles\n" << faces.size() << '\n';
    for (const auto& f : faces) file << f[0]+1 << ' ' << f[1]+1 << ' ' << f[2]+1 << " 13\n";
    file << "End\n";
    require(bool(file), "Cannot save actual BEM surface.");
  }

  Real volume(const Vertices& vertices, const Faces& faces)
  {
    Real result = 0;
    for (const auto& f : faces) result += vertices[f[0]].dot(vertices[f[1]].cross(vertices[f[2]]))/6;
    return result;
  }
}

int main()
{
  try
  {
    std::cout << std::setprecision(17);
    BEM bem;
    auto [vertices, faces] = sphere(3);
    const auto spherical = bem.evaluate(surface(vertices, faces), false);
    const Real pi = std::acos(Real(-1));
    const Real kError = std::abs(spherical.k/(6*pi)-1), qError = std::abs(spherical.q/(8*pi)-1);
    require(spherical.openEdges == 0 && spherical.panels == 512, "Sphere is not closed.");
    require(kError < .2 && qError < .2 && spherical.rho < 1e-10, "Sphere resistance validation failed.");
    std::cout << "sphere k_error=" << kError << " q_error=" << qError
      << " residual=" << spherical.residual << " reciprocity=" << spherical.reciprocity << '\n';
    for (auto& x : vertices)
    {
      const Real a=x(0), b=x(1), c=x(2);
      x *= 1 + 40*a*b*c*(a*a-b*b)*(b*b-c*c)*(c*c-a*a);
    }
    const auto initial = bem.evaluate(surface(vertices, faces));
    require(initial.openEdges == 0 && initial.rho > 1e-9, "Chiral seed coupling is unresolved.");
    Vertices direction(vertices.size());
    Real analytic = 0, dilation = 0;
    for (Index i=0; i<vertices.size(); ++i)
    {
      direction[i].resize(3);
      for (size_t j=0; j<3; ++j) direction[i](j)=std::sin(Real(3*i+j+1)*.731);
      analytic += initial.gradient.row(i).dot(direction[i].getData());
      dilation += initial.gradient.row(i).dot(vertices[i].getData());
    }
    require(std::abs(dilation) < 1e-9, "Area-dependent epsilon dilation derivative failed.");
    Real smallestError=1;
    for (Real h : {Real(1e-3), Real(5e-4), Real(2.5e-4)})
    {
      auto plus=vertices, minus=vertices;
      for (Index i=0; i<vertices.size(); ++i) { plus[i] += h*direction[i]; minus[i] -= h*direction[i]; }
      const Real fd=(bem.evaluate(surface(plus,faces),false).rho-bem.evaluate(surface(minus,faces),false).rho)/(2*h);
      smallestError=std::abs(fd-analytic)/std::max(std::abs(analytic),Real(1e-10));
      std::cout << "finite_difference h=" << h << " analytic=" << analytic << " fd=" << fd << " relative_error=" << smallestError << '\n';
    }
    require(smallestError < 2e-3, "Native BEM coordinate gradient failed finite differences.");
    // A real, closed-surface, volume-preserving update followed by a fresh six-state solve.
    const Real originalVolume=volume(vertices,faces);
    Real gradientNorm=initial.gradient.norm();
    require(gradientNorm > 0, "Zero BEM update direction.");
    bool accepted=false;
    Real step=.02;
    BEM::Result updated;
    Vertices candidate;
    for (size_t attempt=0; attempt<12; ++attempt, step/=2)
    {
      candidate=vertices;
      for (Index i=0; i<vertices.size(); ++i) candidate[i] += step/gradientNorm*initial.gradient.row(i).transpose();
      const Real scale=std::cbrt(originalVolume/volume(candidate,faces));
      for (auto& x : candidate) x *= scale;
      updated=bem.evaluate(surface(candidate,faces),false);
      if (updated.rho > initial.rho) { accepted=true; break; }
    }
    require(accepted && updated.openEdges==0, "No improving real BEM update was accepted.");
    writeSurface(vertices,faces,"bem-initial.mesh");
    writeSurface(candidate,faces,"bem-updated.mesh");
    std::ofstream record("bem-validation.json"); record << std::setprecision(17);
    record << "{\"passed\":true,\"sphere_k_relative_error\":" << kError
      << ",\"sphere_q_relative_error\":" << qError << ",\"gradient_relative_error\":" << smallestError
      << ",\"dilation_derivative\":" << dilation << ",\"initial_rho\":" << initial.rho
      << ",\"updated_rho\":" << updated.rho << ",\"step\":" << step
      << ",\"relative_volume_change\":" << volume(candidate,faces)/originalVolume-1
      << ",\"residual\":" << updated.residual << ",\"open_edges\":" << updated.openEdges << "}\n";
    require(bool(record), "Cannot save validation evidence.");
    std::cout << "accepted_update rho_before=" << initial.rho << " rho_after=" << updated.rho << " step=" << step << '\n';
    return 0;
  }
  catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
