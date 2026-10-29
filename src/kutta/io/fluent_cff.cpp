#include "kutta/io/mesh_io.hpp"
namespace kutta {
PolyMesh read_fluent_cff(const std::filesystem::path&, double) {
  throw Error("Fluent CFF reading requires HDF5 (KUTTA_HAS_HDF5=0). Install HDF5 and configure with -DKUTTA_ENABLE_HDF5=ON.");
}
} // namespace kutta
