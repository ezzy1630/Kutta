#include "kutta/fv/device_mesh.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

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

struct Triangle {
  Vec3d a, b, c;
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
};

DeviceVec3 upload_vectors(Device& dev, const std::vector<Vec3d>& values) {
  std::vector<float> x, y, z;
  x.reserve(values.size()); y.reserve(values.size()); z.reserve(values.size());
  for (auto v : values) { x.push_back(float(v.x)); y.push_back(float(v.y)); z.push_back(float(v.z)); }
  return {dev.upload<float>(x), dev.upload<float>(y), dev.upload<float>(z)};
}
}  // namespace

DeviceMesh upload_mesh(Device& dev, const PolyMesh& mesh, const MeshGeometry& g) {
  validate(mesh);
  KUTTA_REQUIRE(g.cell_centroid.size() == std::size_t(mesh.n_cells) &&
                    g.cell_volume.size() == g.cell_centroid.size() &&
                    g.face_centroid.size() == mesh.owner.size() &&
                    g.face_area.size() == mesh.owner.size(), "FV geometry size mismatch");
  DeviceMesh m;
  m.n_cells = mesh.n_cells; m.n_faces = mesh.n_faces(); m.n_internal_faces = mesh.n_internal_faces();
  m.owner = dev.upload<Index>(mesh.owner); m.neighbour = dev.upload<Index>(mesh.neighbour);
  m.area = upload_vectors(dev, g.face_area);
  std::vector<Vec3d> rp, rn, d;
  std::vector<float> weights, volumes, distance, normal_distance(m.n_boundary_faces(), 0);
  for (Index f = 0; f < m.n_faces; ++f) {
    rp.push_back(g.face_centroid[f] - g.cell_centroid[mesh.owner[f]]);
    if (f < m.n_internal_faces) {
      rn.push_back(g.face_centroid[f] - g.cell_centroid[mesh.neighbour[f]]);
      d.push_back(g.cell_centroid[mesh.neighbour[f]] - g.cell_centroid[mesh.owner[f]]);
      const double owner_normal = dot(rp.back(), g.face_area[f]);
      const double neighbour_normal = -dot(rn.back(), g.face_area[f]);
      weights.push_back(float(neighbour_normal / (owner_normal + neighbour_normal)));
    } else d.push_back(rp.back());
  }
  m.owner_to_face = upload_vectors(dev, rp); m.neighbour_to_face = upload_vectors(dev, rn);
  m.owner_to_neighbour = upload_vectors(dev, d); m.weight = dev.upload<float>(weights);
  const auto cf = build_cell_faces(mesh);
  m.cell_face_offsets = dev.upload<Index>(cf.offsets); m.cell_faces = dev.upload<Index>(cf.faces);
  std::vector<float> sign(cf.signs.begin(), cf.signs.end());
  m.cell_face_sign = dev.upload<float>(sign);
  const WallBvh walls(mesh);
  for (Index c = 0; c < m.n_cells; ++c) {
    KUTTA_REQUIRE(std::isfinite(g.cell_volume[c]) && g.cell_volume[c] > 0,
                  "FV requires finite positive cell volumes");
    volumes.push_back(float(g.cell_volume[c]));
    distance.push_back(float(walls.distance(g.cell_centroid[c])));
  }
  for (const auto& p : mesh.patches) {
    m.patches.push_back({p.name, p.kind, p.start, p.size});
    if (p.kind == PatchKind::Wall)
      for (Index f = p.start; f < p.start + p.size; ++f)
        normal_distance[f - m.n_internal_faces] = float(std::abs(dot(rp[f], g.face_area[f])) / norm(g.face_area[f]));
  }
  m.volume = dev.upload<float>(volumes); m.wall_distance = dev.upload<float>(distance);
  m.wall_normal_distance = dev.upload<float>(normal_distance);
  return m;
}
}  // namespace kutta
