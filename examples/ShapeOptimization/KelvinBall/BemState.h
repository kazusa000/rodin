/* Distributed under the Boost Software License, Version 1.0. */
#ifndef KELVIN_BALL_BEM_STATE_H
#define KELVIN_BALL_BEM_STATE_H

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "Metrics.h"
#include <Rodin/IO/MEDIT.h>

namespace KelvinBall
{
  // Example-specific interchange: the callback differentiates every original
  // Gamma panel. The vector is an integrated nodal covector, not a density.
  struct BemState
  {
    Values metrics;
    std::vector<Math::SpatialVector<Real>> gradient;
    Real residual = 0;

    static BemState evaluate(Mesh& mesh, size_t iteration)
    {
      const char* callback = std::getenv("KELVIN_BEM_COMMAND");
      if (!callback || !*callback)
        throw std::runtime_error("BEM mode requires KELVIN_BEM_COMMAND.");
      const std::string stem = "bem-state-" + std::to_string(iteration);
      for (const auto* extension : {".mesh", ".json", ".nodes"})
        if (std::filesystem::exists(stem + extension))
          throw std::runtime_error("Refusing to overwrite a BEM state artifact.");
      {
        // Mesh::save currently uses the stream's default six digits. BEM must
        // evaluate the actual native coordinates, not a rounded surrogate.
        std::ofstream output(stem + ".mesh");
        if (!output) throw std::runtime_error("Cannot write original BEM input.");
        output << std::setprecision(std::numeric_limits<Real>::max_digits10);
        IO::MeshPrinter<IO::FileFormat::MEDIT, Context::Local>(mesh).print(output);
        output.flush();
        if (!output) throw std::runtime_error("BEM input write failed.");
      }
      // Names above are generated, shell-safe relative filenames. The command
      // itself is an explicit trusted runner setting, never mesh content.
      const std::string command = std::string(callback) + " " + stem +
        ".mesh --json-output " + stem + ".json --nodal-output " + stem + ".nodes";
      if (std::system(command.c_str()) != 0)
        throw std::runtime_error("BEM state/derivative callback failed.");
      return read(mesh, stem + ".nodes");
    }

    static BemState read(const Mesh& mesh, const std::string& filename)
    {
      std::ifstream input(filename);
      std::string token, value;
      if (!(input >> token) || token != "KELVIN_BEM_ETA_GRADIENT_V2")
        throw std::runtime_error("Invalid BEM nodal protocol.");
      std::map<std::string, std::string> header;
      while (input >> token && token != "nodes")
      {
        if (!(input >> value) || !header.emplace(token, value).second)
          throw std::runtime_error("Malformed/duplicate BEM metadata.");
      }
      auto number = [&](const std::string& key) {
        const auto found = header.find(key);
        if (found == header.end()) throw std::runtime_error("Missing BEM key: " + key);
        size_t consumed = 0;
        const Real result = std::stod(found->second, &consumed);
        if (consumed != found->second.size() || !std::isfinite(result))
          throw std::runtime_error("Nonfinite/invalid BEM key: " + key);
        return result;
      };
      if (token != "nodes" || header.at("symmetry") != (Tetrahedral ? "terra" : "cubic") ||
          number("symmetry_order") != ChamberMultiplicity ||
          number("vertex_count") != mesh.getVertexCount() || number("viscosity") != Mu)
        throw std::runtime_error("BEM symmetry, viscosity or mesh size mismatch.");
      BemState result;
      result.metrics.k = number("k");
      result.metrics.c = number("c");
      result.metrics.q = number("q");
      result.metrics.rho = number("eta");
      result.residual = number("linear_relative_residual");
      if (!(result.metrics.k > 0 && result.metrics.q > 0) || result.residual < 0 ||
          result.residual > 1e-10 ||
          std::abs(result.metrics.rho - std::abs(result.metrics.c) /
            std::sqrt(result.metrics.k * result.metrics.q)) >
            1e-12 * std::max(Real(1), result.metrics.rho))
        throw std::runtime_error("BEM resistance/residual validation failed.");
      // These are FEM state diagnostics, unavailable rather than fabricated.
      result.metrics.nitscheJump = std::numeric_limits<Real>::quiet_NaN();
      result.metrics.couplingSymmetry = std::numeric_limits<Real>::quiet_NaN();
      result.gradient.resize(mesh.getVertexCount());
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
      {
        size_t id;
        Math::SpatialVector<Real> x(3), normal(3), g(3);
        if (!(input >> id) || id != vertex + 1)
          throw std::runtime_error("BEM original vertex ordering mismatch.");
        for (auto* vector : {&x, &g, &normal})
          for (Eigen::Index component = 0; component < 3; ++component)
            if (!(input >> (*vector)(component)) || !std::isfinite((*vector)(component)))
              throw std::runtime_error("Malformed/nonfinite BEM node.");
        if ((x - mesh.getVertexCoordinates(vertex)).norm() > 1e-11 * std::max(Real(1), x.norm()))
          throw std::runtime_error("Stale BEM vertex coordinates.");
        result.gradient[vertex] = g;
      }
      if (!(input >> token) || token != "triangles")
        throw std::runtime_error("BEM Gamma connectivity is missing.");
      size_t count = 0;
      std::vector<bool> onGamma(mesh.getVertexCount(), false);
      for (auto face = mesh.getPolytope(mesh.getDimension() - 1); face; ++face)
        if (face->getAttribute() == Gamma)
        {
          const auto& vertices = face->getVertices();
          if (vertices.size() != 3) throw std::runtime_error("BEM needs triangular Gamma.");
          for (Index vertex : vertices)
          {
            size_t id;
            if (!(input >> id) || id != vertex + 1)
              throw std::runtime_error("BEM original Gamma connectivity mismatch.");
            onGamma[vertex] = true;
          }
          ++count;
        }
      if (count != number("gamma_triangle_count") || !(input >> token) || token != "end" || input >> token)
        throw std::runtime_error("BEM panel count or payload length mismatch.");
      for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
        if (!onGamma[vertex] && result.gradient[vertex].norm() != 0)
          throw std::runtime_error("BEM derivative outside Gamma.");
      return result;
    }

    template <class Space>
    Math::Vector<Real> nodalLoad(const Space& space) const
    {
      Math::Vector<Real> load = Math::Vector<Real>::Zero(space.getSize());
      for (Index vertex = 0; vertex < gradient.size(); ++vertex)
      {
        const auto dofs = space.getDOFs(0, vertex);
        if (dofs.size() != 3) throw std::runtime_error("BEM bridge requires vector P1.");
        for (Eigen::Index component = 0; component < 3; ++component)
          load(dofs(component)) += gradient[vertex](component);
      }
      return load;
    }
  };
}
#endif
