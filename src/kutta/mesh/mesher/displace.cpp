#include "internal.hpp"

namespace kutta::mesher::internal {
using namespace detail;
namespace {
constexpr int band_depth = 3, smooth_iterations = 12, zero_after = 10;
// A later run restarts the scales of sources this many wall-point hops from
// a changed target (twice the retreat reach), so retreat near a dropped
// stack is redone while the rest of the mesh keeps its converged scales.
constexpr int reset_hops = 2 * band_depth;
double& at(Vec3d& v, int a) { return a == 0 ? v.x : a == 1 ? v.y : v.z; }
}  // namespace

Displacement::Displacement(Work& w, const Adjacency& adjacency, const PointFaces& pf, std::vector<Index> sources, double min_quality)
    : w_(w), adjacency_(adjacency), pf_(pf), sources_(std::move(sources)), min_quality_(min_quality) {
  const auto np = index(w.points.size());
  scale.assign(np, 1.);
  depth_.assign(np, -1);
  for (auto p : sources_) depth_[p] = 0;
  std::vector<Index> front = sources_;
  for (int d = 1; d <= band_depth; ++d) {
    std::vector<Index> next;
    for (auto p : front) point_neighbours(w, pf, p, [&](Index q) { if (depth_[q] < 0) { depth_[q] = d; next.push_back(q); } });
    band_.insert(band_.end(), next.begin(), next.end());
    front = std::move(next);
  }
  std::vector<char> face_flag(w.faces.size()), cell_flag(w.cells);
  for (Index p = 0; p < np; ++p)
    if (depth_[p] >= 0)
      for (auto i = pf.offsets[p]; i < pf.offsets[p + 1]; ++i) face_flag[pf.faces[i]] = 1;
  for (Index f = 0; f < index(w.faces.size()); ++f)
    if (face_flag[f]) {
      moving_faces_.push_back(f);
      cell_flag[w.faces[f].owner] = 1;
      if (w.faces[f].neighbour >= 0) cell_flag[w.faces[f].neighbour] = 1;
    }
  for (Index c = 0; c < w.cells; ++c) if (cell_flag[c]) moving_cells_.push_back(c);
  origin_ = w.points;
  displacement_.assign(np, {});
  value_.assign(np, {});
  target_.assign(np, {});
  mark_.assign(np, 0);
  reset_mark_.assign(np, 0);
  slot_.assign(np, -1);
  face_mark_.assign(w.faces.size(), 0);
  cell_mark_.assign(w.cells, 0);
}

// Every source at target*scale, the band smoothed from zero (the first run).
void Displacement::apply_all() {
  const auto np = index(w_.points.size());
  std::vector<Vec3d> next(np);
  for (auto p : sources_) displacement_[p] = value_[p] = target_[p] * scale[p];
  for (auto p : band_) displacement_[p] = {};
  for (int k = 0; k < smooth_iterations; ++k) {
    parallel_for(index(band_.size()), [&](Index i) {
      const auto p = band_[i];
      Vec3d sum{};
      int count = 0;
      point_neighbours(w_, pf_, p, [&](Index q) { sum = sum + displacement_[q]; ++count; });
      auto d = count ? sum / double(count) : Vec3d{};
      for (int axis = 0; axis < 3; ++axis) if (w_.planes[p] & (3 << (2 * axis))) at(d, axis) = 0;
      next[p] = d;
    });
    for (auto p : band_) displacement_[p] = next[p];
  }
  for (auto p : sources_) w_.points[p] = origin_[p] + displacement_[p];
  for (auto p : band_) w_.points[p] = origin_[p] + displacement_[p];
}

// The band smoothing is linear in the source values with zero band start and
// zero outside the band, so a change at a few sources adds the smoothing of
// their change alone. After k sweeps that is nonzero only within k band hops
// of the changed sources, so restricting it there is exact. Returns the
// points moved.
std::vector<Index> Displacement::apply_changed() {
  ++stamp_;
  std::vector<Index> changed;
  for (auto p : sources_) {
    const auto v = target_[p] * scale[p];
    if (v.x == value_[p].x && v.y == value_[p].y && v.z == value_[p].z) continue;
    changed.push_back(p);
  }
  if (changed.empty()) return {};
  std::vector<Vec3d> delta(changed.size());
  for (std::size_t i = 0; i < changed.size(); ++i) {
    const auto p = changed[i];
    const auto v = target_[p] * scale[p];
    delta[i] = v - value_[p];
    value_[p] = v;
    mark_[p] = stamp_;
    slot_[p] = index(i);
  }
  // Band points by hop level from the changed sources.
  std::vector<Index> region;
  std::vector<std::size_t> level_end;
  {
    std::vector<Index> front = changed;
    for (int k = 0; k < smooth_iterations; ++k) {
      std::vector<Index> next;
      for (auto p : front)
        point_neighbours(w_, pf_, p, [&](Index q) {
          if (depth_[q] > 0 && mark_[q] != stamp_) { mark_[q] = stamp_; slot_[q] = index(region.size()); region.push_back(q); next.push_back(q); }
        });
      level_end.push_back(region.size());
      front = std::move(next);
    }
  }
  std::vector<Vec3d> current(region.size()), next(region.size());
  const auto value = [&](Index q) -> Vec3d {
    if (mark_[q] != stamp_) return {};
    return depth_[q] == 0 ? delta[std::size_t(slot_[q])] : current[std::size_t(slot_[q])];
  };
  for (int k = 0; k < smooth_iterations; ++k) {
    const auto n = index(level_end[std::size_t(k)]);
    parallel_for(n, [&](Index i) {
      const auto p = region[std::size_t(i)];
      Vec3d sum{};
      int count = 0;
      point_neighbours(w_, pf_, p, [&](Index q) { sum = sum + value(q); ++count; });
      auto d = count ? sum / double(count) : Vec3d{};
      for (int axis = 0; axis < 3; ++axis) if (w_.planes[p] & (3 << (2 * axis))) at(d, axis) = 0;
      next[std::size_t(i)] = d;
    });
    std::copy(next.begin(), next.begin() + std::ptrdiff_t(n), current.begin());
  }
  for (std::size_t i = 0; i < changed.size(); ++i) {
    const auto p = changed[i];
    displacement_[p] = value_[p];
    w_.points[p] = origin_[p] + displacement_[p];
  }
  for (std::size_t i = 0; i < region.size(); ++i) {
    const auto p = region[i];
    displacement_[p] = displacement_[p] + current[i];
    w_.points[p] = origin_[p] + displacement_[p];
  }
  changed.insert(changed.end(), region.begin(), region.end());
  return changed;
}

// Recomputes the geometry of every face and cell on the moved points and
// adds the cells (and their neighbours, whose quality also changes) to `checked`.
void Displacement::refresh(const std::vector<Index>& points, std::vector<Index>& checked) {
  ++stamp_;
  std::vector<Index> faces, cells;
  for (auto p : points)
    for (auto i = pf_.offsets[p]; i < pf_.offsets[p + 1]; ++i) {
      const auto f = pf_.faces[i];
      if (face_mark_[f] == stamp_) continue;
      face_mark_[f] = stamp_;
      faces.push_back(f);
      for (auto c : {w_.faces[f].owner, w_.faces[f].neighbour})
        if (c >= 0 && cell_mark_[c] != stamp_) { cell_mark_[c] = stamp_; cells.push_back(c); }
    }
  update(w_, adjacency_, geometry_, faces, cells);
  const auto touched = with_neighbours(w_, adjacency_, cells);
  checked.insert(checked.end(), touched.begin(), touched.end());
  std::sort(checked.begin(), checked.end());
  checked.erase(std::unique(checked.begin(), checked.end()), checked.end());
}

// Marks the sources within reset_hops wall-point hops of `front` (not yet
// marked this run) and appends them to `reset`.
void Displacement::grow(std::vector<Index> front, std::vector<Index>& reset) {
  for (auto& p : front) if (reset_mark_[p] != run_) { reset_mark_[p] = run_; reset.push_back(p); } else p = -1;
  std::erase(front, -1);
  for (int d = 0; d < reset_hops && !front.empty(); ++d) {
    std::vector<Index> next;
    for (auto p : front)
      point_neighbours(w_, pf_, p, [&](Index q) { if (depth_[q] == 0 && reset_mark_[q] != run_) { reset_mark_[q] = run_; reset.push_back(q); next.push_back(q); } });
    front = std::move(next);
  }
}

// Halves the scales near bad cells (zero after zero_after rounds), as
// displace does. When `allowed` is set, a retreat reaching a source outside
// it stops and returns that source in `outside` instead.
DisplaceResult Displacement::retreat(std::vector<Index>& check, const std::vector<std::uint32_t>* allowed, std::vector<Index>& outside,
                                     const Progress& progress, std::string_view stage) {
  DisplaceResult result;
  std::vector<int> seen(w_.points.size(), -1);
  for (;; ++result.rounds) {
    std::vector<Index> bad;
    for (auto c : check) if (geometry_.quality[c] < min_quality_) bad.push_back(c);
    if (result.rounds == 0) result.first_bad = index(bad.size());
    if (bad.empty()) break;
    if (result.rounds > zero_after + 2) { result.stuck = std::move(bad); break; }
    // Sources within band_depth point hops of a bad cell retreat.
    std::vector<Index> front;
    for (auto c : bad)
      for (auto i = adjacency_.offsets[c]; i < adjacency_.offsets[c + 1]; ++i)
        for (auto p : w_.faces[adjacency_.faces[i]].nodes)
          if (depth_[p] >= 0 && seen[p] != result.rounds) { seen[p] = result.rounds; front.push_back(p); }
    // The retreat starts at the bad cells' own points and widens if needed.
    const int reach = std::min(band_depth, result.rounds / 2);
    std::vector<Index> halve;
    for (int d = 0; d <= reach && !front.empty(); ++d) {
      std::vector<Index> grow;
      for (auto p : front) {
        if (depth_[p] == 0) {
          if (allowed && (*allowed)[p] != run_) outside.push_back(p);
          else halve.push_back(p);
        }
        if (d < reach)
          point_neighbours(w_, pf_, p, [&](Index q) { if (depth_[q] >= 0 && depth_[q] < depth_[p] && seen[q] != result.rounds) { seen[q] = result.rounds; grow.push_back(q); } });
      }
      front = std::move(grow);
    }
    if (!outside.empty()) return result;
    for (auto p : halve) scale[p] = result.rounds + 1 >= zero_after ? 0 : scale[p] * .5;
    if (whole_) {
      apply_all();
      update(w_, adjacency_, geometry_, moving_faces_, moving_cells_);
    } else {
      refresh(apply_changed(), check);
    }
    if (progress) progress(stage, stage == "layers" ? .78 : .7, "Retreating where cell quality would fail");
  }
  return result;
}

DisplaceResult Displacement::run(const std::vector<Vec3d>& target, const Progress& progress, std::string_view stage) {
  std::vector<Index> check, outside;
  if (first_) {
    // Whole mesh, as a single displace call.
    first_ = false;
    for (auto p : sources_) target_[p] = target[p];
    apply_all();
    geometry_ = geometry(w_, adjacency_);
    check = with_neighbours(w_, adjacency_, moving_cells_);
    const auto result = retreat(check, nullptr, outside, progress, stage);
    whole_ = false;
    return result;
  }
  // Sources near a changed target restart at scale 1 and retreat afresh.
  // Where that retreat reaches a source beyond them, the restart widens to
  // include it and repeats, so no source is halved on top of its earlier
  // retreat (which would compound into a collapse cascade).
  ++run_;
  std::vector<Index> changed, reset;
  for (auto p : sources_)
    if (!(target[p].x == target_[p].x && target[p].y == target_[p].y && target[p].z == target_[p].z)) { target_[p] = target[p]; changed.push_back(p); }
  grow(changed, reset);
  for (;;) {
    for (auto p : reset) scale[p] = 1;
    refresh(apply_changed(), check);
    outside.clear();
    const auto result = retreat(check, &reset_mark_, outside, progress, stage);
    if (outside.empty()) return result;
    grow(outside, reset);
  }
}

DisplaceResult displace(Work& w, const Adjacency& adjacency, const PointFaces& pf, const std::vector<Index>& sources,
                        const std::vector<Vec3d>& target, std::vector<double>& scale, double min_quality, const Progress& progress, std::string_view stage) {
  Displacement d(w, adjacency, pf, sources, min_quality);
  d.scale = scale;
  const auto result = d.run(target, progress, stage);
  scale = d.scale;
  return result;
}
}  // namespace kutta::mesher::internal
