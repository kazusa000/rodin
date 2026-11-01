/*
 * Distributed under the Boost Software License, Version 1.0.
 * See accompanying file LICENSE or https://www.boost.org/LICENSE_1_0.txt.
 */
#ifndef KELVIN_BALL_CHECKPOINT_H
#define KELVIN_BALL_CHECKPOINT_H

#include <algorithm>
#include <array>
#include <vector>
#include <stdexcept>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <Rodin/IO/HDF5.h>
#include <Rodin/MMG.h>

namespace KelvinBall
{
  /** Full native topology/coordinates plus MMG tag persistence for SWIFT restarts. */
  class Checkpoint
  {
    public:
      /** Write a native mesh and the MMG tags omitted by generic mesh persistence. */
      static void saveMesh(const std::string& path, const Rodin::MMG::Mesh& mesh)
      {
        using namespace Rodin;
        mesh.save(path, IO::FileFormat::HDF5);
        const auto file = IO::HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
        if (!file)
          throw std::runtime_error("Cannot append checkpoint MMG tags.");
        const auto group = IO::HDF5::Group(H5Gcreate2(file.get(), "/KelvinMMG",
          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        if (!group)
          throw std::runtime_error("Cannot create checkpoint MMG group.");
        const std::array<const IndexSet*, 6> sets{&mesh.getCorners(), &mesh.getRidges(),
          &mesh.getRequiredVertices(), &mesh.getRequiredEdges(),
          &mesh.getRequiredTriangles(), &mesh.getRequiredTetrahedra()};
        for (size_t i = 0; i < sets.size(); ++i)
        {
          std::vector<IO::HDF5::U64> values(sets[i]->begin(), sets[i]->end());
          std::sort(values.begin(), values.end());
          IO::HDF5::writeVectorDataset(file.get(), "/KelvinMMG/Tags" + std::to_string(i), values);
        }
      }

      /** Load the same native topology and tag indices without remeshing or relabeling. */
      static Rodin::MMG::Mesh loadMesh(const std::string& path)
      {
        using namespace Rodin;
        MMG::Mesh mesh;
        mesh.load(path, IO::FileFormat::HDF5);
        const auto file = IO::HDF5::File(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
        if (!file)
          throw std::runtime_error("Cannot read checkpoint MMG tags.");
        const std::array<IndexSet*, 6> sets{&mesh.getCorners(), &mesh.getRidges(),
          &mesh.getRequiredVertices(), &mesh.getRequiredEdges(),
          &mesh.getRequiredTriangles(), &mesh.getRequiredTetrahedra()};
        const std::array<size_t, 6> dimensions{0, 1, 0, 1, 2, 3};
        for (size_t i = 0; i < sets.size(); ++i)
        {
          const auto values = IO::HDF5::readVectorDataset<IO::HDF5::U64>(
            file.get(), "/KelvinMMG/Tags" + std::to_string(i));
          sets[i]->clear();
          for (const auto index : values)
          {
            if (index >= mesh.getConnectivity().getCount(dimensions[i]))
              throw std::runtime_error("Checkpoint MMG tag index is outside its topology.");
            sets[i]->insert(index);
          }
        }
        return mesh;
      }

      /** Persist both independently evolving SWIFT meshes and binary scalar state. */
      static void save(size_t step, const Rodin::MMG::Mesh& mesh,
        const Rodin::MMG::Mesh& background, const std::vector<Rodin::Real>& state)
      {
        std::ostringstream name;
        name << "checkpoint-" << std::setw(6) << std::setfill('0') << step;
        const std::filesystem::path directory(name.str());
        if (!std::filesystem::create_directory(directory))
          throw std::runtime_error("Refuse to overwrite an existing checkpoint.");
        saveMesh((directory / "chamber.h5").string(), mesh);
        saveMesh((directory / "background.h5").string(), background);
        const auto file = Rodin::IO::HDF5::File(H5Fcreate(
          (directory / "state.h5").c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT));
        if (!file)
          throw std::runtime_error("Cannot create checkpoint scalar state.");
        Rodin::IO::HDF5::writeVectorDataset(file.get(), "/State", state);
        if (H5Fflush(file.get(), H5F_SCOPE_GLOBAL) < 0)
          throw std::runtime_error("Cannot flush checkpoint scalar state.");
        std::ofstream ready(directory / "READY");
        ready << "Complete native SWIFT checkpoint; not proof of completed state evaluation.\n";
        if (!ready)
          throw std::runtime_error("Cannot publish checkpoint completion marker.");
      }
  };
}
#endif
