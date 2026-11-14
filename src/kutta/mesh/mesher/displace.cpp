#include "internal.hpp"

namespace kutta::mesher::internal {
using namespace detail;
namespace {
constexpr int band_depth = displace_band_depth, smooth_iterations = 12, zero_after = 10;
double& at(Vec3d& v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }
}  // namespace

DisplaceResult displace(Work& w, const Adjacency& adjacency, const PointFaces& pf, const std::vector<Index>& sources,
                        const std::vector<Vec3d>& target, std::vector<double>& scale, double min_quality, const Progress& progress, std::string_view stage) {
  const auto np = index(w.points.size());
  std::vector<int> depth(np, -1);
  for (auto p : sources) depth[p] = 0;
  std::vector<Index> band;
  {
    std::vector<Index> front = sources;
    for (int d = 1; d <= band_depth; ++d) {
      std::vector<Index> next;
      for (auto p : front) point_neighbours(w, pf, p, [&](Index q) { if (depth[q] < 0) { depth[q] = d; next.push_back(q); } });
      band.insert(band.end(), next.begin(), next.end());
      front = std::move(next);
    }
  }
  std::vector<Index> moving_faces, moving_cells;
  {
    std::vector<char> face_flag(w.faces.size()), cell_flag(w.cells);
    for (Index p = 0; p < np; ++p)
      if (depth[p] >= 0)
        for (auto i = pf.offsets[p]; i < pf.offsets[p + 1]; ++i) face_flag[pf.faces[i]] = 1;
    for (Index f = 0; f < index(w.faces.size()); ++f)
      if (face_flag[f]) {
        moving_faces.push_back(f);
        cell_flag[w.faces[f].owner] = 1;
        if (w.faces[f].neighbour >= 0) cell_flag[w.faces[f].neighbour] = 1;
      }
    for (Index c = 0; c < w.cells; ++c) if (cell_flag[c]) moving_cells.push_back(c);
  }
  const auto origin = w.points;
  std::vector<Vec3d> displacement(np), next(np);
  auto apply = [&] {
    for (auto p : sources) displacement[p] = target[p] * scale[p];
    for (auto p : band) displacement[p] = {};
    for (int k = 0; k < smooth_iterations; ++k) {
      parallel_for(index(band.size()), [&](Index i) {
        const auto p = band[i];
        Vec3d sum{};
        int count = 0;
        point_neighbours(w, pf, p, [&](Index q) { sum = sum + displacement[q]; ++count; });
        auto d = count ? sum / double(count) : Vec3d{};
        for (int axis = 0; axis < 3; ++axis) if (w.planes[p] & (3 << (2 * axis))) at(d, axis) = 0;
        next[p] = d;
      });
      for (auto p : band) displacement[p] = next[p];
    }
    for (auto p : sources) w.points[p] = origin[p] + displacement[p];
    for (auto p : band) w.points[p] = origin[p] + displacement[p];
  };
  apply();
  auto g = geometry(w, adjacency);
  const auto check = with_neighbours(w, adjacency, moving_cells);
  DisplaceResult result;
  std::vector<int> seen(np, -1);
  for (;; ++result.rounds) {
    std::vector<Index> bad;
    for (auto c : check) if (g.quality[c] < min_quality) bad.push_back(c);
    if (result.rounds == 0) result.first_bad = index(bad.size());
    if (bad.empty()) break;
    if (result.rounds > zero_after + 2) { result.stuck = std::move(bad); break; }
    // Sources within band_depth point hops of a bad cell retreat.
    std::vector<Index> front;
    for (auto c : bad)
      for (auto i = adjacency.offsets[c]; i < adjacency.offsets[c + 1]; ++i)
        for (auto p : w.faces[adjacency.faces[i]].nodes)
          if (depth[p] >= 0 && seen[p] != result.rounds) { seen[p] = result.rounds; front.push_back(p); }
    // The retreat starts at the bad cells' own points and widens if needed.
    const int reach = std::min(band_depth, result.rounds / 2);
    for (int d = 0; d <= reach && !front.empty(); ++d) {
      std::vector<Index> grow;
      for (auto p : front) {
        if (depth[p] == 0) scale[p] = result.rounds + 1 >= zero_after ? 0 : scale[p] * .5;
        if (d < reach)
          point_neighbours(w, pf, p, [&](Index q) { if (depth[q] >= 0 && depth[q] < depth[p] && seen[q] != result.rounds) { seen[q] = result.rounds; grow.push_back(q); } });
      }
      front = std::move(grow);
    }
    apply();
    update(w, adjacency, g, moving_faces, moving_cells);
    if (progress) progress(stage, stage == "layers" ? .78 : .7, "Retreating where cell quality would fail");
  }
  return result;
}
}  // namespace kutta::mesher::internal
