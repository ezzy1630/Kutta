#include "kutta/fv/device_mesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <thread>
#include <vector>

#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta {
namespace {
using namespace detail;

double segment_distance2(Vec3d p, Vec3d a, Vec3d b) {
  const auto ab = b - a;
  const double length2 = dot(ab, ab);
  const double t = length2 == 0 ? 0 : std::clamp(dot(p - a, ab) / length2, 0., 1.);
  const auto r = p - (a + ab * t);
  return dot(r, r);
}

Vec3d segment_closest(Vec3d p, Vec3d a, Vec3d b) {
  const auto ab = b - a;
  const double length2 = dot(ab, ab);
  const double t = length2 == 0 ? 0 : std::clamp(dot(p - a, ab) / length2, 0., 1.);
  return a + ab * t;
}

struct Triangle {
  Vec3d a, b, c;
  // Same projection test as distance2, returning the point instead.
  Vec3d closest(Vec3d p) const {
    const auto ab = b - a, ac = c - a, ap = p - a;
    const auto n = cross(ab, ac);
    const double nn = dot(n, n);
    if (nn > 0) {
      const auto q = p - n * (dot(ap, n) / nn);
      if (dot(cross(ab, q - a), n) >= 0 &&
          dot(cross(c - b, q - b), n) >= 0 &&
          dot(cross(a - c, q - c), n) >= 0) return q;
    }
    Vec3d best = segment_closest(p, a, b);
    for (auto q : {segment_closest(p, b, c), segment_closest(p, c, a)})
      if (dot(p - q, p - q) < dot(p - best, p - best)) best = q;
    return best;
  }
  double distance2(Vec3d p) const {
    const auto ab = b - a, ac = c - a, ap = p - a;
    const auto n = cross(ab, ac);
    const double nn = dot(n, n);
    // Project onto the triangle plane; outside it the minimum is on an edge.
    if (nn > 0) {
      const double h = dot(ap, n);
      const auto q = p - n * (h / nn);
      if (dot(cross(ab, q - a), n) >= 0 &&
          dot(cross(c - b, q - b), n) >= 0 &&
          dot(cross(a - c, q - c), n) >= 0) return h * h / nn;
    }
    return std::min({segment_distance2(p, a, b), segment_distance2(p, b, c),
                     segment_distance2(p, c, a)});
  }
};

struct Bounds {
  Vec3d lo{INFINITY, INFINITY, INFINITY}, hi{-INFINITY, -INFINITY, -INFINITY};
  void extend(Vec3d p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  double distance2(Vec3d p) const {
    const Vec3d d{std::max({lo.x - p.x, 0., p.x - hi.x}),
                  std::max({lo.y - p.y, 0., p.y - hi.y}),
                  std::max({lo.z - p.z, 0., p.z - hi.z})};
    return dot(d, d);
  }
};

class WallBvh {
  struct Node { Bounds bounds; std::size_t begin, end; int left = -1, right = -1; };
  std::vector<Triangle> triangles_;
  std::vector<Node> nodes_;

  int build(std::size_t begin, std::size_t end) {
    const int id = static_cast<int>(nodes_.size());
    nodes_.push_back({{}, begin, end});
    for (auto j = begin; j < end; ++j)
      for (auto p : {triangles_[j].a, triangles_[j].b, triangles_[j].c})
        nodes_[id].bounds.extend(p);
    if (end - begin <= 8) return id;
    const auto extent = nodes_[id].bounds.hi - nodes_[id].bounds.lo;
    const int axis = extent.x >= extent.y && extent.x >= extent.z ? 0 :
                     (extent.y >= extent.z ? 1 : 2);
    const auto center = [axis](const Triangle& t) {
      const auto c = (t.a + t.b + t.c) / 3.;
      return axis == 0 ? c.x : (axis == 1 ? c.y : c.z);
    };
    const auto mid = begin + (end - begin) / 2;
    std::nth_element(triangles_.begin() + begin, triangles_.begin() + mid,
                     triangles_.begin() + end,
                     [&](const Triangle& a, const Triangle& b) { return center(a) < center(b); });
    const int left = build(begin, mid), right = build(mid, end);
    nodes_[id].left = left; nodes_[id].right = right;
    return id;
  }
  void search(int id, Vec3d p, double& best) const {
    const auto& n = nodes_[id];
    if (n.bounds.distance2(p) > best) return;
    if (n.left < 0) {
      for (auto j = n.begin; j < n.end; ++j) best = std::min(best, triangles_[j].distance2(p));
      return;
    }
    const bool left_first = nodes_[n.left].bounds.distance2(p) <= nodes_[n.right].bounds.distance2(p);
    search(left_first ? n.left : n.right, p, best);
    search(left_first ? n.right : n.left, p, best);
  }
  void search_closest(int id, Vec3d p, double& best, Vec3d& point) const {
    const auto& n = nodes_[id];
    if (n.bounds.distance2(p) > best) return;
    if (n.left < 0) {
      for (auto j = n.begin; j < n.end; ++j) {
        const auto q = triangles_[j].closest(p);
        const double d2 = dot(p - q, p - q);
        if (d2 < best) { best = d2; point = q; }
      }
      return;
    }
    const bool left_first = nodes_[n.left].bounds.distance2(p) <= nodes_[n.right].bounds.distance2(p);
    search_closest(left_first ? n.left : n.right, p, best, point);
    search_closest(left_first ? n.right : n.left, p, best, point);
  }
 public:
  explicit WallBvh(const PolyMesh& m) {
    for (const auto& patch : m.patches) if (patch.kind == PatchKind::Wall) {
      for (Index f = patch.start; f < patch.start + patch.size; ++f) {
        const auto center = face_center(m, f);
        const auto ns = nodes(m, f);
        for (std::size_t j = 0; j < ns.size(); ++j)
          triangles_.push_back({center, m.points[ns[j]], m.points[ns[(j + 1) % ns.size()]]});
      }
    }
    if (!triangles_.empty()) build(0, triangles_.size());
  }
  double distance(Vec3d p) const {
    double best = std::numeric_limits<double>::infinity();
    if (!nodes_.empty()) search(0, p, best);
    return std::sqrt(best);
  }
  // Unit vector from the nearest wall point to p; zero on a wall or without one.
  Vec3d direction(Vec3d p) const {
    double best = std::numeric_limits<double>::infinity();
    Vec3d point = p;
    if (!nodes_.empty()) search_closest(0, p, best, point);
    const auto r = p - point;
    const double length = std::sqrt(dot(r, r));
    return length > 0 ? r / length : Vec3d{};
  }
};

// Uploads value(0..n) one float component at a time, so neither a Vec3d
// array nor three float copies of it are held on the host.
template <class F>
DeviceVec3 upload_vectors(Device& dev, std::size_t n, const F& value) {
  std::vector<float> c(n);
  DeviceVec3 out;
  for (std::size_t i = 0; i < n; ++i) c[i] = float(value(i).x);
  out.x = dev.upload<float>(c);
  for (std::size_t i = 0; i < n; ++i) c[i] = float(value(i).y);
  out.y = dev.upload<float>(c);
  for (std::size_t i = 0; i < n; ++i) c[i] = float(value(i).z);
  out.z = dev.upload<float>(c);
  return out;
}
template <class F>
void for_cells(std::size_t cells, const F& body) {
  const auto workers = std::clamp<std::size_t>(std::thread::hardware_concurrency(), 1, cells / 4096 + 1);
  // std::thread, not jthread: the CI macOS libc++ lacks jthread.
  std::vector<std::thread> threads;
  for (std::size_t t = 0; t < workers; ++t)
    threads.emplace_back([&, begin = cells * t / workers, end = cells * (t + 1) / workers] {
      for (auto c = begin; c < end; ++c) body(c);
    });
  for (auto& t : threads) t.join();
}
}  // namespace

DeviceVec3 wall_direction(Device& dev, const PolyMesh& mesh, const MeshGeometry& g) {
  KUTTA_REQUIRE(g.cell_centroid.size() == std::size_t(mesh.n_cells), "FV geometry size mismatch");
  const WallBvh walls(mesh);
  std::vector<Vec3d> direction(std::size_t(mesh.n_cells));
  for_cells(direction.size(), [&](std::size_t c) { direction[c] = walls.direction(g.cell_centroid[c]); });
  return upload_vectors(dev, direction.size(), [&](std::size_t c) { return direction[c]; });
}

DeviceMesh upload_mesh(Device& dev, const PolyMesh& mesh, const MeshGeometry& g) {
  validate(mesh);
  KUTTA_REQUIRE(g.cell_centroid.size() == std::size_t(mesh.n_cells) &&
                    g.cell_volume.size() == g.cell_centroid.size() &&
                    g.face_centroid.size() == mesh.owner.size() &&
                    g.face_area.size() == mesh.owner.size(), "FV geometry size mismatch");
  DeviceMesh m;
  m.n_cells = mesh.n_cells; m.n_faces = mesh.n_faces(); m.n_internal_faces = mesh.n_internal_faces();
  m.owner = dev.upload<Index>(mesh.owner); m.neighbour = dev.upload<Index>(mesh.neighbour);
  const auto faces = std::size_t(m.n_faces), internal = std::size_t(m.n_internal_faces);
  m.area = upload_vectors(dev, faces, [&](std::size_t f) { return g.face_area[f]; });
  const auto rp = [&](std::size_t f) { return g.face_centroid[f] - g.cell_centroid[mesh.owner[f]]; };
  const auto rn = [&](std::size_t f) { return g.face_centroid[f] - g.cell_centroid[mesh.neighbour[f]]; };
  std::vector<float> weights(internal), volumes, distance, normal_distance(m.n_boundary_faces(), 0);
  for (std::size_t f = 0; f < internal; ++f) {
    const double owner_normal = dot(rp(f), g.face_area[f]);
    const double neighbour_normal = -dot(rn(f), g.face_area[f]);
    weights[f] = float(neighbour_normal / (owner_normal + neighbour_normal));
  }
  m.owner_to_face = upload_vectors(dev, faces, rp);
  m.neighbour_to_face = upload_vectors(dev, internal, rn);
  m.owner_to_neighbour = upload_vectors(dev, faces, [&](std::size_t f) {
    return f < internal ? g.cell_centroid[mesh.neighbour[f]] - g.cell_centroid[mesh.owner[f]] : rp(f);
  });
  m.weight = dev.upload<float>(weights);
  const auto cf = build_cell_faces(mesh);
  m.cell_face_offsets = dev.upload<Index>(cf.offsets); m.cell_faces = dev.upload<Index>(cf.faces);
  std::vector<float> sign(cf.signs.begin(), cf.signs.end());
  m.cell_face_sign = dev.upload<float>(sign);
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(std::isfinite(g.cell_volume[c]) && g.cell_volume[c] > 0,
                  "FV requires finite positive cell volumes");
    volumes.push_back(float(g.cell_volume[c]));
  }
  // Wall distance dominates large-mesh startup. Each cell's search is
  // independent, so splitting cells across threads gives the serial result.
  const WallBvh walls(mesh);
  distance.resize(m.n_cells);
  const auto cells = std::size_t(m.n_cells);
  const auto workers = std::clamp<std::size_t>(std::thread::hardware_concurrency(), 1, cells / 4096 + 1);
  {
    // std::thread, not jthread: the CI macOS libc++ lacks jthread.
    std::vector<std::thread> threads;
    for (std::size_t t = 0; t < workers; ++t)
      threads.emplace_back([&, begin = cells * t / workers, end = cells * (t + 1) / workers] {
        for (auto c = begin; c < end; ++c) distance[c] = float(walls.distance(g.cell_centroid[c]));
      });
    for (auto& t : threads) t.join();
  }
  for (const auto& p : mesh.patches) {
    m.patches.push_back({p.name, p.kind, p.start, p.size});
    if (p.kind == PatchKind::Wall)
      for (Index f = p.start; f < p.start + p.size; ++f)
        normal_distance[f - m.n_internal_faces] = float(std::abs(dot(rp(std::size_t(f)), g.face_area[f])) / norm(g.face_area[f]));
  }
  m.volume = dev.upload<float>(volumes); m.wall_distance = dev.upload<float>(distance);
  m.wall_normal_distance = dev.upload<float>(normal_distance);
  return m;
}
}  // namespace kutta
