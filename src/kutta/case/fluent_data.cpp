#include "kutta/case/fluent_data.hpp"
#include "kutta/case/hdf5.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace kutta::fluent {
std::span<const double> FieldSection::at(std::uint64_t id) const {
  if (id < first || id > last)
    return {};
  return std::span<const double>(values).subspan(static_cast<std::size_t>(id - first) * components,
                                                 components);
}
std::span<const double> Field::at(std::uint64_t id) const {
  auto it = std::upper_bound(sections.begin(), sections.end(), id,
                             [](auto n, const FieldSection &s) { return n < s.first; });
  return it == sections.begin() ? std::span<const double>{} : std::prev(it)->at(id);
}
std::size_t Field::size() const {
  std::size_t n = 0;
  for (auto &s : sections)
    n += s.values.size() / s.components;
  return n;
}
std::vector<ZoneFieldView> FluentData::zone_field(FieldLocation location, std::int64_t id,
                                                  std::string_view field) const {
  const auto &zones = location == FieldLocation::cell ? mesh.cell_zones : mesh.face_zones;
  auto z = std::find_if(zones.begin(), zones.end(), [id](auto &z) { return z.zone_id == id; });
  if (z == zones.end())
    throw std::out_of_range("unknown Fluent zone ID");
  const auto &fields = location == FieldLocation::cell ? cell_fields : face_fields;
  auto f = fields.find(std::string(field));
  if (f == fields.end())
    throw std::out_of_range("unknown Fluent field");
  std::vector<ZoneFieldView> views;
  for (auto &s : f->second.sections) {
    auto first = std::max(s.first, z->first), last = std::min(s.last, z->last);
    if (first <= last)
      views.push_back({first, last, s.components,
                       std::span<const double>(s.values).subspan(
                           static_cast<std::size_t>(first - s.first) * s.components,
                           static_cast<std::size_t>(last - first + 1) * s.components)});
  }
  return views;
}
#if KUTTA_HAS_HDF5
namespace {
std::uint64_t required_id(hid_t f, const std::string &path, const char *name) {
  auto n = h5::attribute(f, path, name);
  if (!n || !*n)
    throw std::runtime_error("missing/zero ID attribute: " + path + "/" + name);
  return *n;
}
void contained(std::uint64_t first, std::uint64_t last, const std::vector<ZoneRange> &zones) {
  auto next = first;
  for (auto &z : zones) {
    if (z.last < next)
      continue;
    if (z.first > next)
      break;
    if (z.last >= last)
      return;
    next = z.last + 1;
  }
  throw std::runtime_error("data IDs lie outside matching case topology");
}
std::map<std::string, Field> fields(hid_t f, const std::string &base,
                                    const std::vector<ZoneRange> &zones) {
  std::map<std::string, Field> out;
  if (!h5::exists(f, base))
    return out;
  for (auto &name : h5::children(f, base)) {
    if (name == "fields")
      continue;
    const auto group = base + "/" + name;
    Field field;
    field.name = name;
    for (auto &part : h5::children(f, group)) {
      auto path = group + "/" + part;
      FieldSection s;
      s.first = required_id(f, path, "minId");
      s.last = required_id(f, path, "maxId");
      if (s.last < s.first)
        throw std::runtime_error("reversed field ID interval");
      contained(s.first, s.last, zones);
      std::vector<hsize_t> dims;
      s.values = h5::numbers<double>(f, path, &dims);
      if (dims.empty() || dims[0] != s.last - s.first + 1 || (dims.size() == 2 && dims[1] == 0))
        throw std::runtime_error("field shape disagrees with ID interval: " + path);
      s.components = dims.size() == 2 ? static_cast<std::size_t>(dims[1]) : 1;
      if (name == "SV_WALL_SHEAR" && s.components != 3)
        throw std::runtime_error("wall shear must have three row-major components");
      for (auto scalar : {"SV_P", "SV_U", "SV_V", "SV_W", "SV_K", "SV_O", "SV_D", "SV_MU_T",
                          "SV_DENSITY", "SV_WALL_YPLUS", "SV_FLUX"})
        if (name == scalar && s.components != 1)
          throw std::runtime_error("scalar field has multiple components: " + name);
      if (!field.sections.empty() && s.components != field.sections.front().components)
        throw std::runtime_error("inconsistent field component count");
      field.sections.push_back(std::move(s));
    }
    std::sort(field.sections.begin(), field.sections.end(),
              [](auto &a, auto &b) { return a.first < b.first; });
    for (std::size_t i = 1; i < field.sections.size(); ++i)
      if (field.sections[i - 1].last >= field.sections[i].first)
        throw std::runtime_error("overlapping field chunks: " + name);
    if (field.sections.empty())
      throw std::runtime_error("field group has no chunks: " + name);
    out.emplace(name, std::move(field));
  }
  return out;
}
void complete_primary_fields(const FluentData &out, const FluentCase &c) {
  for (auto field_name : {"SV_P", "SV_U", "SV_V", "SV_W"}) {
    auto field = out.cell_fields.find(field_name);
    if (field == out.cell_fields.end())
      continue;
    for (const auto &range : c.mesh.cell_zones) {
      if (!c.cell_zones.empty()) {
        auto zone = std::find_if(c.cell_zones.begin(), c.cell_zones.end(),
                                 [&](auto &z) { return z.id == range.zone_id; });
        if (zone == c.cell_zones.end())
          throw std::runtime_error("case topology has no matching cell-zone settings");
        if (zone->type != ZoneType::fluid)
          continue;
      }
      auto next = range.first;
      bool covered = false;
      for (const auto &chunk : field->second.sections) {
        if (chunk.last < next)
          continue;
        if (chunk.first > next)
          break;
        if (chunk.last >= range.last) {
          covered = true;
          break;
        }
        next = chunk.last + 1;
      }
      if (!covered)
        throw std::runtime_error("primary field " + std::string(field_name) +
                                 " does not cover expected fluid cell IDs in zone " +
                                 std::to_string(range.zone_id));
    }
  }
}
void read_residuals(hid_t f, FluentData &out) {
  const std::string root = "/results/residuals/phase-1";
  if (!h5::exists(f, root)) {
    out.issues.push_back({"residual-history", "absent"});
    return;
  }
  for (auto &eq : h5::children(f, root)) {
    auto path = root + "/" + eq;
    ResidualHistory r;
    std::vector<hsize_t> shape;
    auto iterations = h5::numbers<double>(f, path + "/iterations", &shape);
    if (shape.size() != 1)
      throw std::runtime_error("residual iteration rank must be one");
    auto values = h5::numbers<double>(f, path + "/data", &shape);
    if (shape.size() != 2 || shape[0] != iterations.size() || shape[1] != 4)
      throw std::runtime_error("residual shape mismatch");
    for (std::size_t i = 0; i < iterations.size(); ++i) {
      auto n = iterations[i];
      if (!std::isfinite(n) || n < 0 ||
          n >= static_cast<double>(std::numeric_limits<std::uint64_t>::max()) || std::floor(n) != n)
        throw std::runtime_error("invalid residual iteration");
      auto it = static_cast<std::uint64_t>(n);
      if (!r.iterations.empty() && it <= r.iterations.back())
        throw std::runtime_error("non-increasing residual iterations");
      r.iterations.push_back(it);
      r.rows.push_back({values[4 * i], values[4 * i + 1], values[4 * i + 2], values[4 * i + 3]});
      out.last_residual_iteration = std::max(out.last_residual_iteration.value_or(0), it);
    }
    out.residuals.emplace(eq, std::move(r));
  }
}
} // namespace
#endif
FluentData read_data(const std::filesystem::path &path, const FluentCase &c) {
#if KUTTA_HAS_HDF5
  if (c.mesh.cell_zones.empty() || c.mesh.face_zones.empty())
    throw std::runtime_error("data reader needs matching case topology metadata");
  auto f = h5::open(path);
  FluentData out;
  out.mesh = c.mesh;
  for (auto &key : h5::children(f, "/settings"))
    out.raw_settings.emplace(key, h5::string(f, "/settings/" + key));
  if (auto it = out.raw_settings.find("Data Variables"); it != out.raw_settings.end()) {
    auto s = parse_sexpr(it->second);
    if (s.elements.size() != 2 || s.elements[0].kind != SExpr::Kind::integer ||
        s.elements[0].integer != 37)
      throw std::runtime_error("malformed Data Variables");
    out.variables.alist = std::move(s.elements[1]);
    if (auto n = out.variables.value("number-of-iterations");
        n && n->kind == SExpr::Kind::integer && n->integer >= 0) {
      out.iteration_count = static_cast<std::uint64_t>(n->integer);
      out.variables.interpreted.insert("number-of-iterations");
    }
  }
  if (!out.iteration_count)
    out.issues.push_back(
        {"iteration-count", "absent or unparsed; last residual iteration exposed separately"});
  const std::string phase = "/results/1/phase-1";
  if (!h5::exists(f, phase))
    throw std::runtime_error("missing results/1/phase-1");
  for (auto &result : h5::children(f, "/results"))
    if (result != "1" && result != "residuals")
      out.issues.push_back({"results/" + result, "unparsed additional result set"});
  for (auto &p : h5::children(f, "/results/1"))
    if (p != "phase-1")
      out.issues.push_back({"results/1/" + p, "unparsed additional phase"});
  out.cell_fields = fields(f, phase + "/cells", c.mesh.cell_zones);
  out.face_fields = fields(f, phase + "/faces", c.mesh.face_zones);
  complete_primary_fields(out, c);
  read_residuals(f, out);
  return out;
#else
  (void)path;
  (void)c;
  throw std::runtime_error("HDF5 support is not built");
#endif
}
} // namespace kutta::fluent
