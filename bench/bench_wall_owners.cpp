// Offline wall-face / owner-centroid map for exported field diagnostics.
// No solver, Device allocation or field modification is involved.
#include <iomanip>
#include <iostream>

#include "kutta/io/mesh_io.hpp"

int main(int argc, char** argv) {
  using namespace kutta;
  try {
    KUTTA_REQUIRE(argc == 3, "usage: bench_wall_owners MESH WALL_PATCH > owners.csv");
    const auto mesh = read_fluent_msh(argv[1]);
    const auto geometry = compute_geometry(mesh);
    bool found = false;
    std::cout << "x,y,z,owner_x,owner_y,owner_z\n" << std::setprecision(17);
    for (const auto& patch : mesh.patches) {
      if (patch.name != argv[2]) continue;
      KUTTA_REQUIRE(patch.kind == PatchKind::Wall, "selected patch must be a wall");
      found = true;
      for (Index face = patch.start; face < patch.start + patch.size; ++face) {
        const auto p = geometry.face_centroid[face];
        const auto c = geometry.cell_centroid[mesh.owner[face]];
        std::cout << p.x << ',' << p.y << ',' << p.z << ',' << c.x << ',' << c.y << ',' << c.z << '\n';
      }
    }
    KUTTA_REQUIRE(found, "wall patch not found");
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
