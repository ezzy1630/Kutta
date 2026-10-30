#pragma once

#include "kutta/core/json_value.hpp"
#include "kutta/io/view.hpp"
#include "kutta/io/view_geometry.hpp"
#include "kutta/mesh/topology.hpp"
#include <map>
#include <memory>
#include <variant>

namespace kutta::post {
struct Rake {
  Vec3d start, end;
  int count = 2;
};
struct PlaneGrid {
  Vec3d origin, u, v;
  int nu = 2, nv = 2;
};
struct SurfaceSeeds {
  std::vector<std::string> zones;
  int count = 100;
  bool random = false;
  std::uint64_t random_seed = 0;
  double offset = 0;
};
using Seeds = std::variant<Rake, PlaneGrid, SurfaceSeeds>;
struct Frame {
  double time = 0;
  std::filesystem::path run;
};
struct Trace {
  std::string name = "streamlines";
  Seeds seeds = Rake{};
  std::vector<std::string> fields{"velocity_magnitude", "pressure"};
  bool pathline = false;
  int direction = 1; // 0 = both
  int max_steps = 10000;
  double max_length = 1000, max_time = 1e6, start_time = 0;
  double step_length = 0, relative_tolerance = 1e-5, absolute_tolerance = 0;
  std::vector<Frame> frames;
};
struct Iso {
  std::string name, field;
  double value = 0;
  int component = -1;
  std::vector<std::string> fields{"velocity_magnitude", "pressure"};
};
struct LineProbe {
  std::string name;
  Rake line;
  std::vector<std::string> fields;
};
struct SectionProbe {
  std::string name;
  std::vector<std::string> zones;
  view::Plane plane;
  std::vector<std::string> fields;
};
using Probe = std::variant<LineProbe, SectionProbe>;
struct Integral {
  std::string name = "surface-integrals";
  std::vector<std::string> zones, fields{"pressure"};
  Vec3d origin{}, drag{1, 0, 0}, lift{0, 1, 0}, side{0, 0, 1};
  std::optional<view::Plane> pressure_plane;
};
struct Request {
  std::vector<Trace> traces;
  std::vector<Iso> isosurfaces;
  std::vector<Probe> probes;
  std::vector<Integral> integrals;
  std::vector<std::string> fields; // extra fields on ordinary surfaces/slices
  std::map<std::string, double> reference;
};
Request read_request(const std::filesystem::path &,
                     bool streamlines_only = false);

struct Sample {
  std::array<double, 3> value{};
  std::uint32_t components = 1;
};
struct FieldInfo {
  std::string name, unit;
  std::uint32_t components = 1;
  bool boundary_only = false;
};
struct Walk {
  Index cell = -1, face = -1;
  Vec3d point{};
  bool boundary = false;
};

// All storage is host-side export data. No calls enter the solver's hot path.
// Unlimited inverse-distance LSQ is exact for affine manufactured fields.
class Volume {
public:
  explicit Volume(const view::Snapshot &,
                  const Volume *same_geometry = nullptr);
  ~Volume();
  Volume(const Volume &) = delete;
  Volume &operator=(const Volume &) = delete;
  const view::Snapshot &snapshot() const;
  const MeshGeometry &geometry() const;
  const kutta::CellFaces &cell_faces() const;
  double scale() const;
  double tolerance() const;
  std::vector<std::uint32_t> face_triangles(Index face) const;
  Vec3d decomposition_center(Index cell) const;
  FieldInfo info(std::string_view) const;
  std::vector<FieldInfo> available() const;
  Sample sample(std::string_view, Index cell, Vec3d point) const;
  Sample boundary(std::string_view, Index face, Vec3d point) const;
  std::array<Vec3d, 3> velocity_gradient(Index cell,
                                         std::string_view suffix = {}) const;
  Sample node(std::string_view, Index node) const;
  Index locate(Vec3d point) const;
  bool contains(Index cell, Vec3d point) const;
  Walk walk(Vec3d from, Index cell, Vec3d to) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
std::vector<Vec3d> seeds(const Volume &, const Seeds &);
struct Curves {
  view::Geometry geometry;
  Json::Array lines;
};
Curves trace(const Volume &, const Trace &,
             std::span<const Volume *const> frames = {});
view::Geometry isosurface(const Volume &, const Iso &);
Curves probe(const Volume &, const Probe &);
Json integrate(const Volume &, const Integral &);
void enrich(view::Snapshot &, const Request &);
struct Outputs {
  Json::Object manifest;
  std::vector<view::FieldRange> fields;
};
Outputs write_outputs(const std::filesystem::path &stage,
                      const view::Snapshot &, const Request &);
} // namespace kutta::post
