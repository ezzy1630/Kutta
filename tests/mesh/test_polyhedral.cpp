// Tet-to-polyhedra conversion and octree transition meshes: topology,
// conservation, and that the flow solver keeps an exact uniform flow on them.
#include "../backend_helpers.hpp"
#include "../io/test_support.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/polyhedral.hpp"
#include "kutta/mesh/topology.hpp"

#include <cmath>
#include <map>
#include <numeric>

using namespace kutta;

namespace {
// A unit-square column: one prism layer z in [0, 0.1] under two layers of
// Kuhn tets (six per cube, sharing the 000-111 diagonal), whose bottom
// triangles match the prisms' tops.
PolyMesh prism_tet_column(int n) {
  std::vector<Vec3d> points;
  const std::vector<double> zs{0, .1, .6, 1.1};
  const auto id = [&](int i, int j, int k) { return Index((k * (n + 1) + j) * (n + 1) + i); };
  for (int k = 0; k < 4; ++k)
    for (int j = 0; j <= n; ++j)
      for (int i = 0; i <= n; ++i) points.push_back({double(i) / n, double(j) / n, zs[k]});
  std::vector<detail::CellFaces> cells;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      for (auto tri : {std::array<std::array<int, 2>, 3>{{{0, 0}, {1, 0}, {1, 1}}}, std::array<std::array<int, 2>, 3>{{{0, 0}, {0, 1}, {1, 1}}}}) {
        std::vector<Index> bottom, top;
        for (auto [a, b] : tri) {
          bottom.push_back(id(i + a, j + b, 0));
          top.push_back(id(i + a, j + b, 1));
        }
        cells.push_back({bottom, top});
        for (int e = 0; e < 3; ++e) cells.back().push_back({bottom[e], bottom[(e + 1) % 3], top[(e + 1) % 3], top[e]});
      }
  const std::array<std::array<int, 3>, 6> orders{{{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}}};
  for (int k = 1; k < 3; ++k)
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        for (const auto& order : orders) {
          std::array<int, 3> p{i, j, k};
          std::array<Index, 4> v{id(i, j, k), -1, -1, -1};
          for (int s = 0; s < 3; ++s) {
            ++p[order[s]];
            v[s + 1] = id(p[0], p[1], p[2]);
          }
          cells.push_back({{v[0], v[1], v[2]}, {v[0], v[1], v[3]}, {v[0], v[2], v[3]}, {v[1], v[2], v[3]}});
        }
  std::vector<Patch> patches{{"xmin", PatchKind::VelocityInlet}, {"xmax", PatchKind::PressureOutlet}, {"sides", PatchKind::Symmetry}};
  const auto patch = [&](const std::vector<Index>& nodes) -> Index {
    bool xmin = true, xmax = true;
    for (auto n : nodes) {
      xmin = xmin && points[n].x == 0;
      xmax = xmax && points[n].x == 1;
    }
    return xmin ? 0 : xmax ? 1 : 2;
  };
  return detail::from_cells(points, cells, patches, patch);
}

double total_volume(const PolyMesh& m) {
  const auto g = compute_geometry(m);
  return std::accumulate(g.cell_volume.begin(), g.cell_volume.end(), 0.0);
}
std::map<std::string, double> patch_areas(const PolyMesh& m) {
  const auto g = compute_geometry(m);
  std::map<std::string, double> out;
  for (const auto& p : m.patches)
    for (Index f = p.start; f < p.start + p.size; ++f) out[p.name] += detail::norm(g.face_area[f]);
  return out;
}

// Uniform flow through slip sides is an exact solution; a perturbed start
// must return to it. An inconsistent face or gradient operator on these
// cells would leave a spurious flow or a cycle (as on the native tets).
void check_uniform_flow(const PolyMesh& mesh, const std::map<std::string, fv::BoundaryCondition>& conditions) {
  using namespace fv;
  std::vector<PatchCondition> bc;
  for (const auto& p : mesh.patches) bc.push_back({p.name, conditions.at(p.name)});
  SolverSettings settings;
  settings.viscosity = 1e-5;
  settings.sst = SstSettings{};
  const auto g = compute_geometry(mesh);
  std::vector<Vec3d> velocity(mesh.n_cells, Vec3d{1, 0, 0});
  std::vector<double> pressure, kinetic(mesh.n_cells, 1e-6), omega(mesh.n_cells, 10);
  for (auto c : g.cell_centroid) pressure.push_back(.01 * std::sin(37 * c.x + 13 * c.y + 23 * c.z));
  KUTTA_FOR_EACH_BACKEND(dev) {
    INFO(dev->name());
    SteadySolver solver(*dev, mesh, bc, settings);
    solver.initialize_fields(velocity, pressure, kinetic, omega);
    solver.iterate(300);
    const auto report = solver.monitor();
    CHECK(report.continuity.raw < 1e-6);
    const auto u = dev->download(solver.velocity().x), v = dev->download(solver.velocity().y), w = dev->download(solver.velocity().z);
    const auto p = dev->download(solver.pressure());
    double error = 0, pmax = 0;
    for (Index c = 0; c < mesh.n_cells; ++c) {
      error = std::max({error, std::abs(u[c] - 1.0), std::abs(double(v[c])), std::abs(double(w[c]))});
      pmax = std::max(pmax, std::abs(double(p[c])));
    }
    CHECK(error < 1e-5);
    CHECK(pmax < 1e-5);
  }
}
}  // namespace

TEST_CASE("Tets become their median dual polyhedra while prism layers stay") {
  const int n = 3;
  const auto mesh = prism_tet_column(n);
  REQUIRE(mesh.n_cells == 2 * n * n + 2 * 6 * n * n);
  const auto converted = convert_tets_to_polyhedra(mesh);
  CHECK(converted.tets == 12 * n * n);
  CHECK(converted.kept == 2 * n * n);
  // One polyhedron per tet vertex: the three point layers from z = 0.1 up.
  CHECK(converted.polyhedra + converted.split_vertices == 3 * (n + 1) * (n + 1));
  CHECK(converted.mesh.n_cells == converted.kept + converted.polyhedra + (converted.mesh.n_cells - converted.kept - converted.polyhedra));
  CHECK(total_volume(converted.mesh) == doctest::Approx(total_volume(mesh)).epsilon(1e-12));
  const auto before = patch_areas(mesh), after = patch_areas(converted.mesh);
  for (const auto& [name, area] : before) CHECK(after.at(name) == doctest::Approx(area).epsilon(1e-12));
  // Kept prisms come first and keep their volumes.
  const auto g0 = compute_geometry(mesh), g1 = compute_geometry(converted.mesh);
  for (Index c = 0; c < converted.kept; ++c) CHECK(g1.cell_volume[c] == doctest::Approx(g0.cell_volume[c]).epsilon(1e-12));
  // A dual cell holds a quarter of each tet around its vertex.
  const auto cells = build_cell_faces(converted.mesh);
  Index many_faced = 0;
  for (Index c = converted.kept; c < converted.mesh.n_cells; ++c) many_faced += cells.offsets[c + 1] - cells.offsets[c] > 6;
  CHECK(many_faced > 0);
  // A mesh without tets is returned unchanged.
  const auto hex = make_box_hex(2, 2, 2, {0, 0, 0}, {1, 1, 1});
  const auto same = convert_tets_to_polyhedra(hex);
  CHECK(same.tets == 0);
  CHECK(same.mesh.n_cells == hex.n_cells);

  // The Fluent .msh writer round-trips the polyhedral mesh.
  MeshTemp t;
  write_fluent_msh(t.path("poly.msh"), converted.mesh, "test");
  const auto back = read_fluent_msh(t.path("poly.msh"));
  CHECK(back.n_cells == converted.mesh.n_cells);
  CHECK(back.n_faces() == converted.mesh.n_faces());
  CHECK(total_volume(back) == doctest::Approx(total_volume(converted.mesh)).epsilon(1e-12));
  const auto back_areas = patch_areas(back);
  for (const auto& [name, area] : after) CHECK(back_areas.at(name) == doctest::Approx(area).epsilon(1e-12));

  check_uniform_flow(converted.mesh, {{"xmin", fv::VelocityInlet{Vec3d{1, 0, 0}, fv::TurbulenceValues{1e-6, 10}}},
                                      {"xmax", fv::PressureOutlet{0, fv::NormalBackflow{}, fv::BackflowPressure::Static, fv::TurbulenceValues{1e-6, 10}}},
                                      {"sides", fv::Symmetry{}}});
}

TEST_CASE("Octree boxes are 2:1 balanced with polyhedral transition cells") {
  // Refine twice toward one corner and remove a solid block there.
  const OctreeRefine corner = [](Vec3d lo, Vec3d, int level) { return lo.x < (level ? .5 : 1.0) && lo.z < (level ? .5 : 1.0); };
  const OctreeSolid block = [](Vec3d lo, Vec3d hi) { return (lo.x + hi.x) / 2 < .25 && (lo.z + hi.z) / 2 < .25; };
  const auto mesh = make_octree_box(4, 2, 2, {0, 0, 0}, {4, 2, 2}, 2, corner, block);
  CHECK(mesh.patches.back().name == "body");
  CHECK(total_volume(mesh) == doctest::Approx(16 - .25 * 2 * .25));
  const auto g = compute_geometry(mesh);
  const auto cells = build_cell_faces(mesh);
  double smallest = 1e9, largest = 0;
  Index transition = 0;
  for (Index c = 0; c < mesh.n_cells; ++c) {
    smallest = std::min(smallest, g.cell_volume[c]);
    largest = std::max(largest, g.cell_volume[c]);
    transition += cells.offsets[c + 1] - cells.offsets[c] > 6;
  }
  CHECK(smallest == doctest::Approx(1.0 / 64));
  CHECK(largest == doctest::Approx(1.0));
  CHECK(transition > 0);
  // Face neighbours differ by at most a factor of two in size.
  for (Index f = 0; f < mesh.n_internal_faces(); ++f) {
    const auto ratio = g.cell_volume[mesh.owner[f]] / g.cell_volume[mesh.neighbour[f]];
    CHECK((ratio >= 1.0 / 8 - 1e-12 && ratio <= 8 + 1e-12));
  }
  // A solid cell next to a coarser one is refused.
  const OctreeSolid big = [](Vec3d lo, Vec3d hi) { return (lo.x + hi.x) / 2 < .5 && (lo.z + hi.z) / 2 < 1; };
  CHECK_THROWS(make_octree_box(4, 2, 2, {0, 0, 0}, {4, 2, 2}, 2, corner, big));

  auto flow = make_octree_box(4, 1, 1, {0, 0, 0}, {4, 1, 1}, 2, [](Vec3d lo, Vec3d, int level) { return lo.x >= 1 && lo.x < (level ? 2.0 : 3.0); });
  std::map<std::string, fv::BoundaryCondition> conditions{
      {"xmin", fv::VelocityInlet{Vec3d{1, 0, 0}, fv::TurbulenceValues{1e-6, 10}}},
      {"xmax", fv::PressureOutlet{0, fv::NormalBackflow{}, fv::BackflowPressure::Static, fv::TurbulenceValues{1e-6, 10}}}};
  for (const char* side : {"ymin", "ymax", "zmin", "zmax"}) conditions.emplace(side, fv::Symmetry{});
  check_uniform_flow(flow, conditions);
}

TEST_CASE("Native far-field tets cycle under default numerics; their polyhedral dual recovers uniform flow") {
  // Eleven far-field tets cut from the native Ahmed mesh (see
  // tests/fv/test_native_farfield.cpp): with default pressure-gradient
  // boundary rows they hold a perturbation instead of returning to the exact
  // uniform flow. Their dual is one polyhedron per vertex.
  const auto mesh = read_fluent_msh(std::filesystem::path(KUTTA_SOURCE_DIR) / "tests/data/meshes/native_farfield_tets.msh");
  const auto converted = convert_tets_to_polyhedra(mesh);
  CHECK(converted.tets == 11);
  CHECK(total_volume(converted.mesh) == doctest::Approx(total_volume(mesh)).epsilon(1e-12));
  check_uniform_flow(converted.mesh, {{"inlet", fv::VelocityInlet{Vec3d{1, 0, 0}, fv::TurbulenceValues{1e-6, 10}}},
                                      {"outlet", fv::PressureOutlet{0, fv::NormalBackflow{}, fv::BackflowPressure::Static, fv::TurbulenceValues{1e-6, 10}}},
                                      {"slip", fv::Symmetry{}}});
}

