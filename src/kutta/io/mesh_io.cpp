#include "kutta/io/mesh_io.hpp"

#include <algorithm>
#include <cctype>
namespace kutta {
PolyMesh read_mesh(const std::filesystem::path& path) {
  auto name = path.filename().string();
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  if (name.ends_with(".msh.h5") || name.ends_with(".cas.h5"))
    return read_fluent_cff(path);
  if (name.ends_with(".msh") || name.ends_with(".cas"))
    return read_fluent_msh(path);
  if (name.ends_with(".gz"))
    throw Error(
        "Compressed mesh input is unsupported: gunzip the file before reading");
  if (name.ends_with(".p2dfmt") || name.ends_with(".p3dfmt") ||
      name.ends_with(".p3d") || name.ends_with(".xyz"))
    return read_plot3d(path);
  throw Error("Unsupported mesh extension: " + path.string());
}
}  // namespace kutta
