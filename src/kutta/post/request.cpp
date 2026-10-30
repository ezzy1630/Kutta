#include "kutta/core/json.hpp"
#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>

namespace kutta::post {
namespace {
using Object = Json::Object;
const Object &object(const Json &j, std::string_view path) {
  const auto *p = std::get_if<Object>(&j.value);
  KUTTA_REQUIRE(p, std::string(path) + ": expected object");
  return *p;
}
const Json *get(const Object &o, std::string_view key) {
  const auto it = o.find(std::string(key));
  return it == o.end() ? nullptr : &it->second;
}
const Json &required(const Object &o, std::string_view key) {
  const auto *p = get(o, key);
  KUTTA_REQUIRE(p, "missing post-processing property: " + std::string(key));
  return *p;
}
void keys(const Object &o, std::initializer_list<std::string_view> allowed) {
  for (const auto &[key, value] : o) {
    (void)value;
    KUTTA_REQUIRE(std::find(allowed.begin(), allowed.end(), key) !=
                      allowed.end(),
                  "unknown post-processing property: " + key);
  }
}
std::string string(const Json &j) {
  const auto *p = std::get_if<std::string>(&j.value);
  KUTTA_REQUIRE(p && !p->empty(), "expected nonempty string");
  return *p;
}
std::string string(const Object &o, std::string_view key,
                   std::string fallback) {
  const auto *p = get(o, key);
  return p ? string(*p) : std::move(fallback);
}
double number(const Json &j) {
  const auto *p = std::get_if<double>(&j.value);
  KUTTA_REQUIRE(p && std::isfinite(*p), "expected finite number");
  return *p;
}
double number(const Object &o, std::string_view key, double fallback) {
  const auto *p = get(o, key);
  return p ? number(*p) : fallback;
}
int integer(const Object &o, std::string_view key, int fallback,
            int max = 1000000) {
  const auto n = number(o, key, fallback);
  KUTTA_REQUIRE(n >= 1 && n <= max && std::floor(n) == n,
                "invalid positive integer: " + std::string(key));
  return int(n);
}
// The cache identity: the operation exactly as parsed, plus its effective
// (possibly defaulted, ordinal-dependent) name. Json objects are ordered maps,
// so the dump is canonical.
std::string identity(Object o, const std::string &name) {
  o.insert_or_assign("name", name);
  return dump_json(o);
}
const Json::Array &array(const Json &j) {
  const auto *p = std::get_if<Json::Array>(&j.value);
  KUTTA_REQUIRE(p, "expected array");
  return *p;
}
Vec3d vector(const Json &j) {
  const auto &a = array(j);
  KUTTA_REQUIRE(a.size() == 3, "expected vec3");
  return {number(a[0]), number(a[1]), number(a[2])};
}
Vec3d vector(const Object &o, std::string_view key, Vec3d fallback) {
  const auto *p = get(o, key);
  return p ? vector(*p) : fallback;
}
std::vector<std::string> strings(const Object &o, std::string_view key,
                                 std::vector<std::string> fallback = {}) {
  const auto *p = get(o, key);
  if (!p)
    return fallback;
  std::vector<std::string> out;
  std::set<std::string> seen;
  for (const auto &v : array(*p)) {
    auto s = string(v);
    KUTTA_REQUIRE(seen.insert(s).second, "duplicate name: " + s);
    out.push_back(std::move(s));
  }
  return out;
}
Rake line(const Object &o) {
  return {vector(required(o, "start")), vector(required(o, "end")),
          integer(o, "count", 100)};
}
Seeds seed_spec(const Json &j) {
  const auto &o = object(j, "seeds");
  const auto kind = string(required(o, "type"));
  if (kind == "rake") {
    keys(o, {"type", "start", "end", "count"});
    return line(o);
  }
  if (kind == "plane") {
    keys(o, {"type", "origin", "u", "v", "nu", "nv"});
    PlaneGrid p{vector(required(o, "origin")), vector(required(o, "u")),
                vector(required(o, "v")), integer(o, "nu", 10),
                integer(o, "nv", 10)};
    KUTTA_REQUIRE(detail::norm(detail::cross(p.u, p.v)) > 0 &&
                      std::uint64_t(p.nu) * p.nv <= 1000000,
                  "degenerate or oversized seed plane");
    return p;
  }
  KUTTA_REQUIRE(kind == "surface", "seed type must be rake, plane or surface");
  keys(o, {"type", "zones", "count", "distribution", "random_seed", "offset"});
  SurfaceSeeds s;
  s.zones = strings(o, "zones");
  KUTTA_REQUIRE(!s.zones.empty(), "surface seeds require exact zone names");
  s.count = integer(o, "count", 100);
  const auto distribution = string(o, "distribution", "uniform");
  KUTTA_REQUIRE(distribution == "uniform" || distribution == "random",
                "surface distribution must be uniform or random");
  s.random = distribution == "random";
  const auto random_seed = number(o, "random_seed", 0);
  KUTTA_REQUIRE(random_seed >= 0 && random_seed <= 9007199254740991.0 &&
                    std::floor(random_seed) == random_seed,
                "invalid random_seed");
  s.random_seed = std::uint64_t(random_seed);
  s.offset = number(o, "offset", 0);
  KUTTA_REQUIRE(s.offset >= 0,
                "surface seed offset must be nonnegative inward distance");
  return s;
}
Trace trace_spec(const Json &j, const std::filesystem::path &base,
                 std::size_t i) {
  const auto &o = object(j, "streamlines");
  keys(o, {"name", "mode", "seeds", "fields", "direction", "max_steps",
           "max_length", "max_time", "start_time", "step_length",
           "relative_tolerance", "absolute_tolerance", "frames"});
  Trace t;
  t.name = string(o, "name", "streamlines-" + std::to_string(i));
  t.identity = identity(o, t.name);
  t.seeds = seed_spec(required(o, "seeds"));
  t.fields = strings(o, "fields", t.fields);
  const auto mode = string(o, "mode", "streamline");
  KUTTA_REQUIRE(mode == "streamline" || mode == "pathline",
                "mode must be streamline or pathline");
  t.pathline = mode == "pathline";
  const auto direction = string(o, "direction", "forward");
  KUTTA_REQUIRE(direction == "forward" || direction == "backward" ||
                    direction == "both",
                "invalid trace direction");
  t.direction = direction == "both" ? 0 : direction == "backward" ? -1 : 1;
  t.max_steps = integer(o, "max_steps", 10000);
  t.max_length = number(o, "max_length", t.max_length);
  t.max_time = number(o, "max_time", t.max_time);
  t.start_time = number(o, "start_time", 0);
  t.step_length = number(o, "step_length", 0);
  t.relative_tolerance = number(o, "relative_tolerance", 1e-5);
  t.absolute_tolerance = number(o, "absolute_tolerance", 0);
  KUTTA_REQUIRE(t.max_length > 0 && t.max_time > 0 && t.step_length >= 0 &&
                    t.relative_tolerance > 0 && t.relative_tolerance <= .1 &&
                    t.absolute_tolerance >= 0,
                "invalid trace integration bounds/tolerances");
  if (const auto *frames = get(o, "frames")) {
    KUTTA_REQUIRE(t.pathline, "time frames are only valid for pathlines");
    for (const auto &frame : array(*frames)) {
      const auto &f = object(frame, "frame");
      keys(f, {"time", "run"});
      Frame out{number(required(f, "time")), [&] {
                  const auto s = string(required(f, "run"));
                  return std::filesystem::path(
                      std::u8string(s.begin(), s.end()));
                }()};
      if (out.run.is_relative())
        out.run = base / out.run;
      KUTTA_REQUIRE(t.frames.empty() || out.time > t.frames.back().time,
                    "frame times must be strictly increasing");
      t.frames.push_back(std::move(out));
    }
    KUTTA_REQUIRE(t.frames.size() >= 2 &&
                      t.start_time >= t.frames.front().time &&
                      t.start_time <= t.frames.back().time,
                  "pathlines need two frames bracketing start_time");
  }
  return t;
}
} // namespace
Request read_request(const std::filesystem::path &path, bool streamlines_only) {
  KUTTA_REQUIRE(std::filesystem::file_size(path) <= 1024 * 1024,
                "post request exceeds 1 MiB");
  std::ifstream f(path, std::ios::binary);
  KUTTA_REQUIRE(f, "cannot open post request");
  const std::string text{std::istreambuf_iterator<char>(f), {}};
  const auto json = parse_json(text);
  const auto &root = object(json, "post request");
  KUTTA_REQUIRE(required(root, "v") == Json(1), "post request requires v:1");
  Request out;
  if (streamlines_only && get(root, "seeds")) {
    auto entry = root;
    entry.erase("v");
    out.traces.push_back(trace_spec(entry, path.parent_path(), 0));
    return out;
  }
  keys(root, {"v", "fields", "fields_policy", "reference", "streamlines",
              "isosurfaces", "probes", "integrals"});
  out.fields = strings(root, "fields");
  const auto policy = string(root, "fields_policy", "strict");
  KUTTA_REQUIRE(policy == "strict" || policy == "available",
                "fields_policy must be strict or available");
  out.available_fields = policy == "available";
  if (const auto *ref = get(root, "reference"))
    for (const auto &[key, value] : object(*ref, "reference")) {
      KUTTA_REQUIRE(key == "density" || key == "viscosity" ||
                        key == "reference_density" || key == "velocity" ||
                        key == "pressure" || key == "area" || key == "length",
                    "unknown reference property: " + key);
      const auto n = number(value);
      KUTTA_REQUIRE(key == "pressure" || n > 0,
                    "reference property must be positive: " + key);
      out.reference.emplace(key, n);
    }
  if (const auto *traces = get(root, "streamlines"))
    for (const auto &j : array(*traces))
      out.traces.push_back(
          trace_spec(j, path.parent_path(), out.traces.size()));
  if (const auto *isos = get(root, "isosurfaces"))
    for (const auto &j : array(*isos)) {
      const auto &o = object(j, "isosurface");
      keys(o, {"name", "field", "value", "component", "fields"});
      Iso iso;
      iso.name =
          string(o, "name", "iso-" + std::to_string(out.isosurfaces.size()));
      iso.identity = identity(o, iso.name);
      iso.field = string(required(o, "field"));
      iso.value = number(required(o, "value"));
      iso.fields = strings(o, "fields", iso.fields);
      const auto component = number(o, "component", -1);
      KUTTA_REQUIRE(component >= -1 && component <= 2 &&
                        std::floor(component) == component,
                    "component must be -1 (magnitude), 0, 1 or 2");
      iso.component = int(component);
      out.isosurfaces.push_back(std::move(iso));
    }
  if (const auto *probes = get(root, "probes"))
    for (const auto &j : array(*probes)) {
      const auto &o = object(j, "probe");
      const auto type = string(o, "type", "line"),
                 name = string(o, "name",
                               "probe-" + std::to_string(out.probes.size()));
      const auto fields =
          strings(o, "fields", {"pressure", "velocity_magnitude"});
      if (type == "line") {
        keys(o, {"type", "name", "start", "end", "count", "fields"});
        out.probes.emplace_back(
            LineProbe{name, line(o), fields, identity(o, name)});
      } else {
        KUTTA_REQUIRE(type == "section", "probe type must be line or section");
        keys(o, {"type", "name", "zones", "origin", "normal", "fields"});
        SectionProbe p{
            name,
            strings(o, "zones"),
            {vector(required(o, "origin")), vector(required(o, "normal"))},
            fields, identity(o, name)};
        KUTTA_REQUIRE(!p.zones.empty() && detail::norm(p.plane.normal) > 0,
                      "section requires zones and nonzero plane normal");
        p.plane.normal =
            detail::operator/(p.plane.normal, detail::norm(p.plane.normal));
        out.probes.emplace_back(std::move(p));
      }
    }
  if (const auto *integrals = get(root, "integrals"))
    for (const auto &j : array(*integrals)) {
      const auto &o = object(j, "integral");
      keys(o, {"name", "zones", "fields", "origin", "drag", "lift", "side",
               "pressure_plane"});
      Integral r;
      r.name = string(o, "name",
                      "integrals-" + std::to_string(out.integrals.size()));
      r.identity = identity(o, r.name);
      r.zones = strings(o, "zones");
      r.fields = strings(o, "fields", r.fields);
      r.origin = vector(o, "origin", {});
      r.drag = vector(o, "drag", r.drag);
      r.lift = vector(o, "lift", r.lift);
      r.side = vector(o, "side", r.side);
      if (const auto *plane = get(o, "pressure_plane")) {
        const auto &p = object(*plane, "pressure_plane");
        keys(p, {"origin", "normal"});
        r.pressure_plane = view::Plane{vector(required(p, "origin")),
                                       vector(required(p, "normal"))};
        KUTTA_REQUIRE(detail::norm(r.pressure_plane->normal) > 0,
                      "pressure-plane normal must be nonzero");
        r.pressure_plane->normal = detail::operator/(
            r.pressure_plane->normal, detail::norm(r.pressure_plane->normal));
      }
      for (auto *d : {&r.drag, &r.lift, &r.side}) {
        KUTTA_REQUIRE(detail::norm(*d) > 0, "report direction must be nonzero");
        *d = detail::operator/(*d, detail::norm(*d));
      }
      out.integrals.push_back(std::move(r));
    }
  KUTTA_REQUIRE(out.traces.size() + out.isosurfaces.size() + out.probes.size() +
                        out.integrals.size() <=
                    256,
                "too many post-processing operations");
  std::set<std::string> names;
  const auto add = [&](const auto &r) {
    KUTTA_REQUIRE(names.insert(r.name).second,
                  "duplicate operation name: " + r.name);
  };
  for (const auto &t : out.traces)
    add(t);
  for (const auto &i : out.isosurfaces)
    add(i);
  for (const auto &p : out.probes)
    std::visit(add, p);
  for (const auto &i : out.integrals)
    add(i);
  return out;
}
} // namespace kutta::post
