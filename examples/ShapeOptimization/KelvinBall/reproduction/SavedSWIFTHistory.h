/*
 * Distributed under the Boost Software License, Version 1.0.
 * See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt.
 */
#ifndef KELVIN_BALL_SAVED_SWIFT_HISTORY_H
#define KELVIN_BALL_SAVED_SWIFT_HISTORY_H

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <Rodin/IO/HDF5.h>
#include <Rodin/MMG.h>

namespace KelvinBall
{
  /// Inputs emitted before SWIFT fitting, indexed on the independent background.
  struct SavedSWIFTHistory
  {
    std::filesystem::path directory;
    std::vector<double> backgroundSpacing;
    std::map<std::string, std::string> last;

    static std::vector<std::string> columns(const std::string& line)
    {
      std::vector<std::string> fields;
      std::istringstream stream(line);
      for (std::string field; std::getline(stream, field, ','); )
        fields.push_back(field);
      return fields;
    }

    explicit SavedSWIFTHistory(const std::filesystem::path& input) : directory(input)
    {
      std::ifstream history(directory / "kelvin-ball.csv");
      std::string line;
      if (!std::getline(history, line))
        throw std::runtime_error("Missing parent KelvinBall history.");
      const auto header = columns(line);
      while (std::getline(history, line))
      {
        const auto row = columns(line);
        if (row.size() != header.size())
          throw std::runtime_error("Incomplete input history row.");
        for (size_t j = 0; j < row.size(); ++j)
          last[header[j]] = row[j];
        if (std::stoul(last.at("iteration")) != backgroundSpacing.size())
          throw std::runtime_error("Recovery requires contiguous source states starting at zero.");
        backgroundSpacing.push_back(std::stod(last.at("background_h")));
      }
      if (backgroundSpacing.empty() || !std::isfinite(value("rho")))
        throw std::runtime_error("No scored state to continue.");
    }

    double value(const std::string& name) const { return std::stod(last.at(name)); }

    std::string snapshot(const std::string& field, size_t index) const
    {
      std::ostringstream name;
      name << "KelvinBall.Chamber." << field << '.' << std::setw(6) << std::setfill('0')
           << index << ".h5";
      return (directory / name.str()).string();
    }

    template <class Field>
    void readTarget(size_t step, const Rodin::MMG::Mesh& background, Field& target) const
    {
      using namespace Rodin;
      namespace H5 = IO::HDF5;
      const auto chamber = H5::File(H5Fopen(snapshot("mesh", step).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      const auto field = H5::File(H5Fopen(snapshot("Advected", step).c_str(), H5F_ACC_RDONLY, H5P_DEFAULT));
      if (!chamber || !field)
        throw std::runtime_error("Missing chamber mesh or transported field.");
      const auto topology = H5::readVectorDataset<H5::U64>(chamber.get(), "/Mesh/XDMF/Topology");
      const auto values = H5::readVectorDataset<H5::F64>(field.get(), "/GridFunction/Values/Data");
      if (values.size() != background.getVertexCount() || topology.size() != 5 * background.getCellCount())
        throw std::runtime_error("Recovered background has different dimensions at step " + std::to_string(step));
      for (size_t i = 0; i < background.getCellCount(); ++i)
      {
        const auto& vertices = background.getConnectivity().getPolytope(3, i);
        if (topology[5 * i] != 6 || vertices.size() != 4 ||
            !std::equal(vertices.begin(), vertices.end(), topology.begin() + 5 * i + 1))
          throw std::runtime_error("Recovered tetrahedron numbering differs at step " + std::to_string(step));
      }
      for (size_t i = 0; i < values.size(); ++i)
      {
        if (!std::isfinite(values[i]))
          throw std::runtime_error("Nonfinite saved transport value.");
        target.getData()(target.getFiniteElementSpace().getDOFs(0, i)[0]) = values[i];
      }
    }
  };
}
#endif
