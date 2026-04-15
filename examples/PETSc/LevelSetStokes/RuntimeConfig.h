#pragma once

#include <cerrno>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace LevelSetStokes::Runtime
{
  inline const char* envCString(const char* key, const char* fallback)
  {
    const char* value = std::getenv(key);
    return (value && value[0] != '\0') ? value : fallback;
  }

  inline int envInt(const char* key, int fallback)
  {
    const char* value = std::getenv(key);
    if (!value || value[0] == '\0')
      return fallback;

    char* end = nullptr;
    errno = 0;
    const long parsed = std::strtol(value, &end, 10);
    if (errno != 0 || !end || *end != '\0'
        || parsed < std::numeric_limits<int>::min()
        || parsed > std::numeric_limits<int>::max())
    {
      std::ostringstream oss;
      oss << "Invalid integer value for " << key << ": " << value;
      throw std::runtime_error(oss.str());
    }

    return static_cast<int>(parsed);
  }

  inline size_t envSizeT(const char* key, size_t fallback)
  {
    const char* value = std::getenv(key);
    if (!value || value[0] == '\0')
      return fallback;

    char* end = nullptr;
    errno = 0;
    const unsigned long long parsed = std::strtoull(value, &end, 10);
    if (errno != 0 || !end || *end != '\0'
        || parsed > std::numeric_limits<size_t>::max())
    {
      std::ostringstream oss;
      oss << "Invalid size_t value for " << key << ": " << value;
      throw std::runtime_error(oss.str());
    }

    return static_cast<size_t>(parsed);
  }

  inline double envDouble(const char* key, double fallback)
  {
    const char* value = std::getenv(key);
    if (!value || value[0] == '\0')
      return fallback;

    char* end = nullptr;
    errno = 0;
    const double parsed = std::strtod(value, &end);
    if (errno != 0 || !end || *end != '\0')
    {
      std::ostringstream oss;
      oss << "Invalid floating-point value for " << key << ": " << value;
      throw std::runtime_error(oss.str());
    }

    return parsed;
  }
}
