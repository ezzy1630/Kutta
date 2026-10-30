#pragma once

#include <filesystem>
#include <string>

namespace kutta {

// UTF-8 bytes of a path, for C APIs such as HDF5 that take char* file names.
// path::c_str() is wchar_t on Windows, and path::string() converts through the
// ANSI code page there, which loses non-ASCII names; HDF5 expects UTF-8.
inline std::string utf8_path(const std::filesystem::path& path) {
  const auto bytes = path.u8string();
  return {bytes.begin(), bytes.end()};
}

}  // namespace kutta
