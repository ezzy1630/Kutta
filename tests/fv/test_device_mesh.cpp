#include "../backend_helpers.hpp"
#include "kutta/fv/device_mesh.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
using namespace kutta::detail;

TEST_CASE("FV upload preserves gather geometry and exact box wall distances") {
  for (bool tet : {false, true}) {
    auto mesh = tet ? make_box_tet(4, 3, 2, {1e8, -2e8, 3e8}, {1e8+2, -2e8+3, 3e8+4}) :
                      make_box_hex(4, 3, 2, {1e8, -2e8, 3e8}, {1e8+2, -2e8+3, 3e8+4});
    mesh.patches[5].kind = PatchKind::Symmetry;
    const auto g = compute_geometry(mesh);
    KUTTA_FOR_EACH_BACKEND(dev) {
      auto m = upload_mesh(*dev, mesh, g);
      const auto wd = dev->download(m.wall_distance), dy = dev->download(m.wall_normal_distance);
      const auto volume = dev->download(m.volume), dx = dev->download(m.owner_to_neighbour.x);
      for (Index c = 0; c < mesh.n_cells; ++c) {
        const auto p = g.cell_centroid[c];
        const double expected = std::min({p.x-1e8, 1e8+2-p.x, p.y+2e8, -2e8+3-p.y, p.z-3e8});
        CHECK(wd[c] == doctest::Approx(expected).epsilon(2e-6));
        CHECK(volume[c] == doctest::Approx(g.cell_volume[c]).epsilon(2e-6));
      }
      for (const auto& p : mesh.patches) for (Index f=p.start; f<p.start+p.size; ++f) {
        const auto r = g.face_centroid[f] - g.cell_centroid[mesh.owner[f]];
        const double expected = p.kind == PatchKind::Wall ? std::abs(dot(r,g.face_area[f]))/norm(g.face_area[f]) : 0;
        CHECK(dy[f-m.n_internal_faces] == doctest::Approx(expected).epsilon(2e-6));
        CHECK(dx[f] == doctest::Approx(r.x).epsilon(2e-6));
      }
      CHECK(dev->download(m.owner) == mesh.owner);
      CHECK(dev->download(m.neighbour) == mesh.neighbour);
    }
  }
}

TEST_CASE("FV wall distance is to the finite surface rather than its infinite plane") {
  auto mesh = make_box_hex(2,2,2,{0,0,0},{2,2,2});
  split_patch(mesh, 0, "wall_tile", [](Vec3d c) { return c.y < 1 && c.z < 1; });
  for (auto& p : mesh.patches) if (p.name != "wall_tile") p.kind = PatchKind::Symmetry;
  const auto g = compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = upload_mesh(*dev,mesh,g);
    auto wd = dev->download(m.wall_distance);
    for (Index c=0;c<mesh.n_cells;++c) {
      auto p=g.cell_centroid[c];
      CHECK(wd[c] == doctest::Approx(std::sqrt(p.x*p.x + std::pow(std::max(0.,p.y-1),2) +
                                               std::pow(std::max(0.,p.z-1),2))).epsilon(2e-6));
    }
  }
  for (auto& p : mesh.patches) p.kind = PatchKind::Symmetry;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m = upload_mesh(*dev,mesh,g);
    for (auto d : dev->download(m.wall_distance)) CHECK(std::isinf(d));
  }
}
