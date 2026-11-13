// Axis-aligned 2:1 octree meshes with polyhedral transition cells.
#include <array>
#include <cstdint>
#include <unordered_map>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta {
using namespace detail;
namespace {
struct Cell {
  std::array<int, 3> lo;  // in finest lattice units
  int size;
};
}  // namespace

PolyMesh make_octree_box(int nx, int ny, int nz, Vec3d lo, Vec3d hi, int levels, const OctreeRefine& refine,
                         const OctreeSolid& solid) {
  KUTTA_REQUIRE(nx > 0 && ny > 0 && nz > 0 && levels >= 0 && levels <= 8, "octree box: invalid dimensions or levels");
  KUTTA_REQUIRE(finite(lo) && finite(hi) && hi.x > lo.x && hi.y > lo.y && hi.z > lo.z, "octree box requires finite lo < hi");
  const int unit = 1 << levels;
  const std::array<int, 3> n{nx * unit, ny * unit, nz * unit};
  KUTTA_REQUIRE(std::int64_t(n[0]) * n[1] * n[2] < (std::int64_t(1) << 31), "octree box: lattice too large");
  const Vec3d h{(hi.x - lo.x) / n[0], (hi.y - lo.y) / n[1], (hi.z - lo.z) / n[2]};
  const auto point = [&](std::array<int, 3> p) { return Vec3d{lo.x + p[0] * h.x, lo.y + p[1] * h.y, lo.z + p[2] * h.z}; };
  const auto bounds = [&](const Cell& c) {
    return std::pair{point(c.lo), point({c.lo[0] + c.size, c.lo[1] + c.size, c.lo[2] + c.size})};
  };
  std::vector<Cell> cells;
  for (int k = 0; k < nz; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i) cells.push_back({{i * unit, j * unit, k * unit}, unit});
  std::vector<Index> owner_of(std::size_t(n[0]) * n[1] * n[2]);
  const auto at = [&](int i, int j, int k) -> Index& { return owner_of[(std::size_t(k) * n[1] + j) * n[0] + i]; };
  const auto paint = [&] {
    for (Index c = 0; c < index(cells.size()); ++c) {
      const auto& q = cells[c];
      for (int k = 0; k < q.size; ++k)
        for (int j = 0; j < q.size; ++j)
          for (int i = 0; i < q.size; ++i) at(q.lo[0] + i, q.lo[1] + j, q.lo[2] + k) = c;
    }
  };
  const auto split = [&](std::vector<char>& marked) {
    std::vector<Cell> next;
    for (std::size_t c = 0; c < cells.size(); ++c) {
      if (!marked[c]) {
        next.push_back(cells[c]);
        continue;
      }
      const int s = cells[c].size / 2;
      for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
          for (int i = 0; i < 2; ++i) next.push_back({{cells[c].lo[0] + i * s, cells[c].lo[1] + j * s, cells[c].lo[2] + k * s}, s});
    }
    cells = std::move(next);
  };
  for (int level = 0; level < levels; ++level) {
    std::vector<char> marked(cells.size(), 0);
    bool any = false;
    for (std::size_t c = 0; c < cells.size(); ++c) {
      const auto [a, b] = bounds(cells[c]);
      marked[c] = cells[c].size == unit >> level && refine(a, b, level);
      any |= marked[c];
    }
    if (any) split(marked);
  }
  // 2:1 balance across faces.
  for (;;) {
    paint();
    std::vector<char> marked(cells.size(), 0);
    bool any = false;
    for (std::size_t c = 0; c < cells.size() && true; ++c) {
      const auto& q = cells[c];
      for (int axis = 0; axis < 3 && !marked[c]; ++axis)
        for (int side = 0; side < 2 && !marked[c]; ++side) {
          const int plane = side ? q.lo[axis] + q.size : q.lo[axis] - 1;
          if (plane < 0 || plane >= n[axis]) continue;
          const int u = (axis + 1) % 3, v = (axis + 2) % 3;
          for (int a = 0; a < q.size && !marked[c]; ++a)
            for (int b = 0; b < q.size; ++b) {
              std::array<int, 3> p{};
              p[axis] = plane;
              p[u] = q.lo[u] + a;
              p[v] = q.lo[v] + b;
              if (cells[at(p[0], p[1], p[2])].size * 2 < q.size) {
                marked[c] = 1;
                break;
              }
            }
        }
      any |= marked[c];
    }
    if (!any) break;
    split(marked);
  }
  // Solid cells leave the mesh; their faces become the "body" wall.
  std::vector<char> removed(cells.size(), 0);
  bool body = false;
  if (solid)
    for (std::size_t c = 0; c < cells.size(); ++c) {
      const auto [a, b] = bounds(cells[c]);
      removed[c] = solid(a, b);
      body |= removed[c];
    }
  std::vector<Index> id(cells.size(), -1);
  Index n_cells = 0;
  for (std::size_t c = 0; c < cells.size(); ++c)
    if (!removed[c]) id[c] = n_cells++;

  std::unordered_map<std::uint64_t, Index> vertex;
  std::vector<Vec3d> points;
  const auto key = [&](std::array<int, 3> p) {
    return (std::uint64_t(p[0]) << 42) | (std::uint64_t(p[1]) << 21) | std::uint64_t(p[2]);
  };
  for (std::size_t c = 0; c < cells.size(); ++c) {
    if (removed[c]) continue;
    for (int corner = 0; corner < 8; ++corner) {
      std::array<int, 3> p{cells[c].lo[0] + (corner & 1) * cells[c].size, cells[c].lo[1] + (corner >> 1 & 1) * cells[c].size,
                           cells[c].lo[2] + (corner >> 2 & 1) * cells[c].size};
      if (vertex.try_emplace(key(p), index(points.size())).second) points.push_back(point(p));
    }
  }
  std::vector<Patch> patches;
  for (const char* name : {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}) patches.push_back({name, PatchKind::Wall});
  if (body) patches.push_back({"body", PatchKind::Wall});
  std::vector<Face> faces;
  // A rectangle in the plane x[axis]=plane, from (a0,b0) to (a1,b1) in the
  // two other axes, wound so its normal points along +axis or -axis. Every
  // lattice vertex on its boundary is a node, so hanging nodes conform.
  const auto rectangle = [&](int axis, int plane, int a0, int b0, int a1, int b1, bool positive) {
    const int u = (axis + 1) % 3, v = (axis + 2) % 3;
    std::vector<Index> loop;
    const auto visit = [&](int a, int b) {
      std::array<int, 3> p{};
      p[axis] = plane;
      p[u] = a;
      p[v] = b;
      if (const auto it = vertex.find(key(p)); it != vertex.end()) loop.push_back(it->second);
    };
    for (int a = a0; a < a1; ++a) visit(a, b0);
    for (int b = b0; b < b1; ++b) visit(a1, b);
    for (int a = a1; a > a0; --a) visit(a, b1);
    for (int b = b1; b > b0; --b) visit(a0, b);
    if (!positive) std::reverse(loop.begin(), loop.end());
    return loop;
  };
  for (std::size_t c = 0; c < cells.size(); ++c) {
    if (removed[c]) continue;
    const auto& q = cells[c];
    for (int axis = 0; axis < 3; ++axis)
      for (int side = 0; side < 2; ++side) {
        const int u = (axis + 1) % 3, v = (axis + 2) % 3;
        const int plane = side ? q.lo[axis] + q.size : q.lo[axis];
        const int across = side ? plane : plane - 1;
        const auto make = [&](Index neighbour, Index patch) {
          faces.push_back({rectangle(axis, plane, q.lo[u], q.lo[v], q.lo[u] + q.size, q.lo[v] + q.size, side), id[c], neighbour, patch});
        };
        if (across < 0 || across >= n[axis]) {
          make(-1, 2 * axis + side);
          continue;
        }
        std::array<int, 3> p{};
        p[axis] = across;
        p[u] = q.lo[u];
        p[v] = q.lo[v];
        const auto other = at(p[0], p[1], p[2]);
        // The smaller cell (or the lower index of two equal ones) owns a face.
        if (cells[other].size < q.size) {
          for (int a = 0; a < q.size; ++a)
            for (int b = 0; b < q.size; ++b) {
              p[u] = q.lo[u] + a;
              p[v] = q.lo[v] + b;
              KUTTA_REQUIRE(!removed[at(p[0], p[1], p[2])], "octree box: a solid cell must only touch cells of its own size; refine around it");
            }
          continue;
        }
        if (removed[other]) {
          make(-1, 6);
          continue;
        }
        if (cells[other].size == q.size && id[other] < id[c]) continue;
        make(id[other], -1);
      }
  }
  auto mesh = pack(std::move(points), n_cells, std::move(faces), std::move(patches));
  validate(mesh);
  return mesh;
}

}  // namespace kutta
