#pragma once
#include <map>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::detail {
struct FluentFace {
  std::vector<Index> nodes;
  Index c0 = 0, c1 = 0;
  int zone = -1;
};
struct FluentZone {
  std::string name, type;
};
struct FluentRaw {
  int dimension = 3;
  std::map<Index, Vec3d> nodes;
  std::map<Index, int> cells;  // active source cell -> source zone
  std::map<Index, FluentFace> faces;
  std::map<int, FluentZone> zones;
  std::vector<PeriodicFacePair> periodic;  // source face ids until assemble
};
std::string fluent_type(int id);
PatchKind patch_kind(const std::string& type);
FluentMeshData assemble_fluent(const FluentRaw& raw, double thickness);
}  // namespace kutta::detail
