#pragma once
#include <algorithm>
#include <cstdint>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::detail {
// Records keyed by positive source id. Files declare ids in contiguous ranges
// (sections, zones), so the table holds ranges, not one heap node per id:
// that cost about 1 KB per cell at read time. Iteration is in ascending id,
// and rank() is an id's position in that order.
template <class T>
class IdRanges {
 public:
  // Adds ids [first, first + values.size()); false if any of them is present.
  bool insert(Index first, std::vector<T> values) {
    if (values.empty()) return true;
    const auto last = std::int64_t(first) + std::int64_t(values.size()) - 1;
    auto at = std::lower_bound(ranges_.begin(), ranges_.end(), first,
                               [](const Range& r, Index id) { return r.first < id; });
    if ((at != ranges_.end() && at->first <= last) ||
        (at != ranges_.begin() && std::prev(at)->last() >= first))
      return false;
    ranges_.insert(at, Range{first, 0, std::move(values)});
    std::size_t base = 0;
    for (auto& r : ranges_) {
      r.base = base;
      base += r.values.size();
    }
    size_ = base;
    return true;
  }
  const T* find(Index id) const {
    const auto* r = range(id);
    return r ? &r->values[std::size_t(id - r->first)] : nullptr;
  }
  T* find(Index id) { return const_cast<T*>(std::as_const(*this).find(id)); }
  // Ascending position of a present id.
  std::size_t rank(Index id) const {
    const auto* r = range(id);
    KUTTA_REQUIRE(r, "Fluent: unknown id " + std::to_string(id));
    return r->base + std::size_t(id - r->first);
  }
  std::size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }
  Index max_id() const { return ranges_.empty() ? 0 : ranges_.back().last(); }
  template <class F>
  void each(F&& f) const {
    for (const auto& r : ranges_)
      for (std::size_t i = 0; i < r.values.size(); ++i) f(Index(r.first + Index(i)), r.values[i]);
  }

 private:
  struct Range {
    Index first;
    std::size_t base;  // rank of `first`
    std::vector<T> values;
    Index last() const { return Index(first + Index(values.size()) - 1); }
  };
  const Range* range(Index id) const {
    auto at = std::upper_bound(ranges_.begin(), ranges_.end(), id,
                               [](Index i, const Range& r) { return i < r.first; });
    if (at == ranges_.begin()) return nullptr;
    --at;
    return id <= at->last() ? &*at : nullptr;
  }
  std::vector<Range> ranges_;  // sorted by first, disjoint
  std::size_t size_ = 0;
};
struct FluentFace {
  Index c0 = 0, c1 = 0;
  int zone = -1;
  Index first_node = 0, node_count = 0;  // into FluentRaw::face_nodes
};
struct FluentZone {
  std::string name, type;
};
struct FluentRaw {
  int dimension = 3;
  IdRanges<Vec3d> nodes;
  IdRanges<int> cells;  // active source cell -> source zone
  IdRanges<FluentFace> faces;
  std::vector<Index> face_nodes;  // source node ids
  std::map<int, FluentZone> zones;
  std::vector<PeriodicFacePair> periodic;  // source face ids until assemble
  std::span<const Index> nodes_of(const FluentFace& f) const {
    return std::span<const Index>(face_nodes).subspan(std::size_t(f.first_node), std::size_t(f.node_count));
  }
};
std::string fluent_type(int id);
PatchKind patch_kind(const std::string& type);
FluentMeshData assemble_fluent(const FluentRaw& raw, double thickness);
}  // namespace kutta::detail
