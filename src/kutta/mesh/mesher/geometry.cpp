#include "internal.hpp"

namespace kutta::mesher::internal {
using namespace detail;

Adjacency cell_faces(const Work& w) {
  Adjacency a;
  a.offsets.assign(w.cells + 1, 0);
  for (const auto& f : w.faces) { ++a.offsets[f.owner + 1]; if (f.neighbour >= 0) ++a.offsets[f.neighbour + 1]; }
  for (Index c = 0; c < w.cells; ++c) a.offsets[c + 1] += a.offsets[c];
  a.faces.resize(a.offsets.back()); a.signs.resize(a.offsets.back());
  auto next = a.offsets;
  for (Index f = 0; f < index(w.faces.size()); ++f) {
    const auto& face = w.faces[f];
    a.faces[next[face.owner]] = f; a.signs[next[face.owner]++] = 1;
    if (face.neighbour >= 0) { a.faces[next[face.neighbour]] = f; a.signs[next[face.neighbour]++] = -1; }
  }
  return a;
}

PointFaces point_faces(const Work& w) {
  PointFaces p;
  p.offsets.assign(w.points.size() + 1, 0);
  for (const auto& f : w.faces) for (auto n : f.nodes) ++p.offsets[n + 1];
  for (std::size_t i = 0; i < w.points.size(); ++i) p.offsets[i + 1] += p.offsets[i];
  p.faces.resize(p.offsets.back());
  auto next = p.offsets;
  for (Index f = 0; f < index(w.faces.size()); ++f) for (auto n : w.faces[f].nodes) p.faces[next[n]++] = f;
  return p;
}

FaceGeom polygon(const std::vector<Vec3d>& pts, const std::vector<Index>& ns) {
  const auto origin = pts[ns[0]];
  Vec3d center{};
  for (auto n : ns) center = center + (pts[n] - origin);
  center = origin + center / double(ns.size());
  Vec3d area_vector{};
  for (std::size_t j = 0; j < ns.size(); ++j)
    area_vector = area_vector + cross(pts[ns[j]] - center, pts[ns[(j + 1) % ns.size()]] - center) * 0.5;
  const auto magnitude = norm(area_vector);
  if (!(magnitude > 0) || !std::isfinite(magnitude)) return {center, {}};
  const auto normal = area_vector / magnitude;
  Vec3d moment{};
  double area = 0;
  for (std::size_t j = 0; j < ns.size(); ++j) {
    const auto a = pts[ns[j]] - center, b = pts[ns[(j + 1) % ns.size()]] - center;
    const auto t = dot(cross(a, b) * 0.5, normal);
    area += t;
    moment = moment + (a + b) * (t / 3);
  }
  if (!(area > 0)) return {center, {}};
  return {center + moment / area, area_vector};
}

namespace {
void cell(const Work& w, const Adjacency& a, Geometry& g, Index c) {
  const auto first = a.offsets[c], last = a.offsets[c + 1];
  const auto anchor = g.face[a.faces[first]].centroid;
  Vec3d estimate{};
  for (auto i = first; i < last; ++i) estimate = estimate + (g.face[a.faces[i]].centroid - anchor);
  estimate = anchor + estimate / double(last - first);
  double volume = 0;
  Vec3d moment{};
  for (auto i = first; i < last; ++i) {
    const auto& f = g.face[a.faces[i]];
    const auto r = f.centroid - estimate;
    const auto v = a.signs[i] * dot(f.area, r) / 3;
    volume += v;
    moment = moment + r * (0.75 * v);
  }
  g.volume[c] = volume;
  g.centroid[c] = volume > 0 ? estimate + moment / volume : estimate;
  (void)w;
}
double quality(const Work& w, const Adjacency& a, const Geometry& g, Index c) {
  if (!(g.volume[c] > 0)) return 0;
  // Faces far below the cell's face scale underflow in the float solver.
  const auto tiny = degenerate_face * std::pow(g.volume[c], 2. / 3);
  double q = 1;
  const auto cosine = [](Vec3d x, Vec3d y) { const auto d = norm(x) * norm(y); return d > 0 ? dot(x, y) / d : -1.; };
  for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i) {
    const auto f = a.faces[i];
    const auto& fg = g.face[f];
    const auto area = fg.area * double(a.signs[i]);
    if (!(norm(area) > tiny)) return 0;
    q = std::min(q, cosine(area, fg.centroid - g.centroid[c]));
    const auto& face = w.faces[f];
    if (face.neighbour >= 0) {
      const auto other = a.signs[i] > 0 ? face.neighbour : face.owner;
      q = std::min(q, cosine(area, g.centroid[other] - g.centroid[c]));
    }
  }
  return std::max(q, 0.);
}
}  // namespace

Geometry geometry(const Work& w, const Adjacency& a) {
  Geometry g;
  g.face.resize(w.faces.size());
  g.centroid.resize(w.cells); g.volume.resize(w.cells); g.quality.resize(w.cells);
  parallel_for(index(w.faces.size()), [&](Index f) { g.face[f] = polygon(w.points, w.faces[f].nodes); });
  parallel_for(w.cells, [&](Index c) { cell(w, a, g, c); });
  parallel_for(w.cells, [&](Index c) { g.quality[c] = quality(w, a, g, c); });
  return g;
}

std::vector<Index> with_neighbours(const Work& w, const Adjacency& a, const std::vector<Index>& cells) {
  std::vector<Index> out = cells;
  for (auto c : cells)
    for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i) {
      const auto& f = w.faces[a.faces[i]];
      if (f.neighbour >= 0) out.push_back(f.owner == c ? f.neighbour : f.owner);
    }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

void update(const Work& w, const Adjacency& a, Geometry& g, const std::vector<Index>& faces, const std::vector<Index>& cells) {
  parallel_for(index(faces.size()), [&](Index i) { g.face[faces[i]] = polygon(w.points, w.faces[faces[i]].nodes); });
  parallel_for(index(cells.size()), [&](Index i) { cell(w, a, g, cells[i]); });
  const auto touched = with_neighbours(w, a, cells);
  parallel_for(index(touched.size()), [&](Index i) { g.quality[touched[i]] = quality(w, a, g, touched[i]); });
}

// Ericson, Real-Time Collision Detection 5.1.5.
Vec3d closest_on_triangle(Vec3d p, Vec3d a, Vec3d b, Vec3d c) {
  const auto ab = b - a, ac = c - a, ap = p - a;
  const auto d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0 && d2 <= 0) return a;
  const auto bp = p - b;
  const auto d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) return b;
  const auto vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + ab * (d1 / (d1 - d3));
  const auto cp = p - c;
  const auto d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) return c;
  const auto vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + ac * (d2 / (d2 - d6));
  const auto va = d3 * d6 - d5 * d4;
  if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
  const auto denom = 1 / (va + vb + vc);
  return a + ab * (vb * denom) + ac * (vc * denom);
}
}  // namespace kutta::mesher::internal
