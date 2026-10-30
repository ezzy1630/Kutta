#include "aero.hpp"
#include "kutta/mesh/detail.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>

namespace kutta::post {
using namespace detail;
namespace {
using Object = Json::Object;
Json vector(Vec3d v) { return Json::Array{v.x, v.y, v.z}; }
Json optional(std::optional<double> v) { return v ? Json(*v) : Json(nullptr); }
bool positive(std::optional<double> v) { return v && std::isfinite(*v) && *v > 0; }

// loads.json is written by solve, so a malformed file is a broken run (exit 1).
const Json &field(const Object &o, const std::string &key) {
  const auto it = o.find(key);
  KUTTA_REQUIRE(it != o.end(), "loads.json: missing " + key);
  return it->second;
}
template <class T> const T &as(const Json &j, const std::string &what) {
  const auto *p = std::get_if<T>(&j.value);
  KUTTA_REQUIRE(p, "loads.json: invalid " + what);
  return *p;
}
double finite(const Json &j, const std::string &what) {
  const auto v = as<double>(j, what);
  KUTTA_REQUIRE(std::isfinite(v), "loads.json: invalid " + what);
  return v;
}
Vec3d vec3(const Json &j, const std::string &what) {
  const auto &a = as<Json::Array>(j, what);
  KUTTA_REQUIRE(a.size() == 3, "loads.json: invalid " + what);
  return {finite(a[0], what), finite(a[1], what), finite(a[2], what)};
}

// Spec parsing: every failure names its property path.
[[noreturn]] void invalid(const std::string &path, const std::string &why) {
  throw InvalidAeroSpec("aero spec " + path + ": " + why);
}
const Object &spec_object(const Json &j, const std::string &path,
                          std::initializer_list<std::string_view> allowed) {
  const auto *o = std::get_if<Object>(&j.value);
  if (!o)
    invalid(path, "expected an object");
  for (const auto &[key, value] : *o) {
    (void)value;
    if (std::find(allowed.begin(), allowed.end(), key) == allowed.end())
      invalid(path == "" ? key : path + "." + key, "unknown property");
  }
  return *o;
}
double spec_number(const Json &j, const std::string &path) {
  const auto *n = std::get_if<double>(&j.value);
  if (!n || !std::isfinite(*n))
    invalid(path, "expected a finite number");
  return *n;
}
Vec3d spec_vec3(const Json &j, const std::string &path) {
  const auto *a = std::get_if<Json::Array>(&j.value);
  if (!a || a->size() != 3)
    invalid(path, "expected [x,y,z]");
  return {spec_number((*a)[0], path), spec_number((*a)[1], path),
          spec_number((*a)[2], path)};
}
std::vector<std::string> spec_names(const Json &j, const std::string &path) {
  const auto *a = std::get_if<Json::Array>(&j.value);
  if (!a || a->empty())
    invalid(path, "expected a nonempty array of zone names");
  std::vector<std::string> out;
  for (const auto &v : *a) {
    const auto *s = std::get_if<std::string>(&v.value);
    if (!s || s->empty())
      invalid(path, "expected zone names");
    if (std::find(out.begin(), out.end(), *s) != out.end())
      invalid(path, "zone listed twice: " + *s);
    out.push_back(*s);
  }
  return out;
}
const Json *get(const Object &o, const std::string &key) {
  const auto it = o.find(key);
  return it == o.end() ? nullptr : &it->second;
}

struct Sum {
  Vec3d pressure, viscous, pressure_moment, viscous_moment;
  void add(const ZoneLoads &z) {
    pressure = pressure + z.pressure;
    viscous = viscous + z.viscous;
    pressure_moment = pressure_moment + z.pressure_moment;
    viscous_moment = viscous_moment + z.viscous_moment;
  }
  Vec3d force() const { return pressure + viscous; }
};
struct Frame {
  Vec3d drag, lift, side, origin, centre;
  // Dynamic pressure and the force (q A) and moment (q A L) denominators.
  std::optional<double> q, force_scale, moment_scale;
  Object projected(Vec3d f) const {
    return {{"drag", dot(f, drag)}, {"lift", dot(f, lift)}, {"side", dot(f, side)}};
  }
  // Moment about the moment centre from one about the loads origin.
  Vec3d moment(Vec3d about_origin, Vec3d force) const {
    return about_origin + cross(origin - centre, force);
  }
};
Json ratio(double value, std::optional<double> denominator) {
  return denominator ? Json(value / *denominator) : Json(nullptr);
}
} // namespace

Json loads_json(const Loads &l) {
  Json::Array zones, report_zones(l.report_zones.begin(), l.report_zones.end());
  for (const auto &z : l.zones)
    zones.push_back(Object{{"name", z.name},
                           {"zone_id", double(z.zone_id)},
                           {"area", z.area},
                           {"pressure", vector(z.pressure)},
                           {"viscous", vector(z.viscous)},
                           {"pressure_moment", vector(z.pressure_moment)},
                           {"viscous_moment", vector(z.viscous_moment)}});
  const auto &r = l.reference;
  return Object{{"v", 1},
                {"origin", vector(l.origin)},
                {"units", Object{{"force", "N"}, {"moment", "N m"}, {"area", "m^2"}}},
                {"reference", Object{{"area", optional(r.area)},
                                     {"density", optional(r.density)},
                                     {"velocity", optional(r.velocity)},
                                     {"length", optional(r.length)},
                                     {"pressure", optional(r.pressure)}}},
                {"report_zones", std::move(report_zones)},
                {"zones", std::move(zones)},
                {"provenance", l.provenance},
                {"approximate", l.approximate},
                {"modified", l.modified},
                {"iterations", l.iterations}};
}
Loads parse_loads(const Json &j) {
  const auto &o = as<Object>(j, "root");
  KUTTA_REQUIRE(field(o, "v") == Json(1), "loads.json: expected version 1");
  Loads l;
  l.origin = vec3(field(o, "origin"), "origin");
  const auto &r = as<Object>(field(o, "reference"), "reference");
  for (auto [key, target] : {std::pair{"area", &l.reference.area}, {"density", &l.reference.density},
                             {"velocity", &l.reference.velocity}, {"length", &l.reference.length},
                             {"pressure", &l.reference.pressure}}) {
    const auto &v = field(r, key);
    if (!std::holds_alternative<std::nullptr_t>(v.value))
      *target = finite(v, std::string("reference.") + key);
  }
  for (const auto &z : as<Json::Array>(field(o, "zones"), "zones")) {
    const auto &zo = as<Object>(z, "zone");
    ZoneLoads zone{as<std::string>(field(zo, "name"), "zone name"),
                   std::int64_t(finite(field(zo, "zone_id"), "zone_id")),
                   finite(field(zo, "area"), "zone area"),
                   vec3(field(zo, "pressure"), "pressure"),
                   vec3(field(zo, "viscous"), "viscous"),
                   vec3(field(zo, "pressure_moment"), "pressure_moment"),
                   vec3(field(zo, "viscous_moment"), "viscous_moment")};
    KUTTA_REQUIRE(std::none_of(l.zones.begin(), l.zones.end(), [&](const auto &p) { return p.name == zone.name; }),
                  "loads.json: duplicate zone " + zone.name);
    l.zones.push_back(std::move(zone));
  }
  for (const auto &z : as<Json::Array>(field(o, "report_zones"), "report_zones")) {
    const auto &name = as<std::string>(z, "report zone");
    KUTTA_REQUIRE(std::any_of(l.zones.begin(), l.zones.end(), [&](const auto &p) { return p.name == name; }),
                  "loads.json: report zone without loads: " + name);
    l.report_zones.push_back(name);
  }
  l.provenance = as<std::string>(field(o, "provenance"), "provenance");
  l.approximate = as<bool>(field(o, "approximate"), "approximate");
  l.modified = as<bool>(field(o, "modified"), "modified");
  l.iterations = int(finite(field(o, "iterations"), "iterations"));
  return l;
}

AeroSpec parse_aero_spec(const Json &j) {
  const auto &o = spec_object(j, "", {"v", "axes", "moment_centre", "car_zones", "groups", "vehicle"});
  AeroSpec s;
  if (const auto *v = get(o, "v"); v && *v != Json(1))
    invalid("v", "expected protocol version 1");
  if (const auto *axes = get(o, "axes")) {
    const auto &a = spec_object(*axes, "axes", {"drag", "lift"});
    for (auto [key, target] : {std::pair{"drag", &s.drag}, {"lift", &s.lift}})
      if (const auto *v = get(a, key)) {
        const auto axis = spec_vec3(*v, std::string("axes.") + key);
        if (!(norm(axis) > 0))
          invalid(std::string("axes.") + key, "must be nonzero");
        *target = axis / norm(axis);
      }
    if (std::abs(dot(s.drag, s.lift)) > 1e-6)
      invalid("axes", "drag and lift must be orthogonal");
  }
  if (const auto *c = get(o, "moment_centre"))
    s.moment_centre = spec_vec3(*c, "moment_centre");
  if (const auto *c = get(o, "car_zones"))
    s.car_zones = spec_names(*c, "car_zones");
  if (const auto *g = get(o, "groups")) {
    const auto *groups = std::get_if<Json::Array>(&g->value);
    if (!groups)
      invalid("groups", "expected an array");
    std::set<std::string> zones;
    for (std::size_t i = 0; i < groups->size(); ++i) {
      const auto path = "groups[" + std::to_string(i) + "]";
      const auto &group = spec_object((*groups)[i], path, {"name", "zones"});
      const auto *name = get(group, "name");
      const auto *text = name ? std::get_if<std::string>(&name->value) : nullptr;
      if (!text || text->empty())
        invalid(path + ".name", "expected a nonempty name");
      if (std::any_of(s.groups.begin(), s.groups.end(), [&](const auto &p) { return p.name == *text; }))
        invalid(path + ".name", "duplicate group name: " + *text);
      const auto *members = get(group, "zones");
      if (!members)
        invalid(path + ".zones", "missing");
      auto names = spec_names(*members, path + ".zones");
      for (const auto &z : names)
        if (!zones.insert(z).second)
          invalid(path + ".zones", "zone is in more than one group: " + z);
      s.groups.push_back({*text, std::move(names)});
    }
  }
  if (const auto *v = get(o, "vehicle")) {
    const auto &vehicle = spec_object(*v, "vehicle", {"front_axle", "wheelbase"});
    const auto *axle = get(vehicle, "front_axle"), *wheelbase = get(vehicle, "wheelbase");
    if (!axle || !wheelbase)
      invalid("vehicle", "needs front_axle and wheelbase");
    s.vehicle = AeroSpec::Vehicle{spec_vec3(*axle, "vehicle.front_axle"), spec_number(*wheelbase, "vehicle.wheelbase")};
    if (!(s.vehicle->wheelbase > 0))
      invalid("vehicle.wheelbase", "must be positive");
  }
  return s;
}
AeroSpec read_aero_spec(const std::filesystem::path &path) {
  std::error_code ec;
  const auto size = std::filesystem::file_size(path, ec);
  if (ec)
    throw InvalidAeroSpec("cannot read aero spec: " + path.string());
  if (size > 1024 * 1024)
    throw InvalidAeroSpec("aero spec exceeds 1 MiB");
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw InvalidAeroSpec("cannot read aero spec: " + path.string());
  const std::string text{std::istreambuf_iterator<char>(in), {}};
  Json parsed;
  try {
    parsed = parse_json(text);
  } catch (const std::exception &e) {
    throw InvalidAeroSpec(std::string("aero spec: ") + e.what());
  }
  return parse_aero_spec(parsed);
}

Json aero(const Loads &loads, const AeroSpec &spec) {
  const auto find = [&](const std::string &name) {
    return std::find_if(loads.zones.begin(), loads.zones.end(), [&](const auto &z) { return z.name == name; });
  };
  // Defaults match the case's own force reports, else the whole wall set.
  std::vector<std::string> car = spec.car_zones ? *spec.car_zones : loads.report_zones;
  if (!spec.car_zones && car.empty())
    for (const auto &z : loads.zones)
      car.push_back(z.name);
  for (const auto &name : car)
    if (find(name) == loads.zones.end())
      invalid("car_zones", "not a wall zone of this run: " + name);
  std::set<std::string> grouped;
  for (const auto &g : spec.groups)
    for (const auto &z : g.zones) {
      if (std::find(car.begin(), car.end(), z) == car.end())
        invalid("groups", "zone is not a car zone: " + z);
      grouped.insert(z);
    }

  const auto &ref = loads.reference;
  Frame frame{spec.drag, spec.lift, cross(spec.lift, spec.drag), loads.origin, spec.moment_centre, {}, {}, {}};
  if (positive(ref.density) && positive(ref.velocity))
    frame.q = .5 * *ref.density * *ref.velocity * *ref.velocity;
  if (frame.q && positive(ref.area))
    frame.force_scale = *frame.q * *ref.area;
  if (frame.force_scale && positive(ref.length))
    frame.moment_scale = *frame.force_scale * *ref.length;

  Sum total;
  for (const auto &z : loads.zones)
    if (std::find(car.begin(), car.end(), z.name) != car.end())
      total.add(z);
  const double total_downforce = -dot(total.force(), frame.lift);
  const auto row = [&](const std::string &name, const char *kind, const std::vector<std::string> &zones) {
    Sum s;
    for (const auto &z : zones)
      s.add(*find(z));
    const auto force = s.force();
    const auto moment = frame.moment(s.pressure_moment + s.viscous_moment, force);
    const double drag = dot(force, frame.drag), lift = dot(force, frame.lift);
    const double pitch = dot(moment, frame.side), roll = dot(moment, frame.drag), yaw = dot(moment, frame.lift);
    Json::Array names(zones.begin(), zones.end());
    return Object{
        {"name", name},
        {"kind", kind},
        {"zones", std::move(names)},
        {"force", Object{{"pressure", frame.projected(s.pressure)},
                         {"viscous", frame.projected(s.viscous)},
                         {"total", frame.projected(force)}}},
        {"downforce", -lift},
        {"moment", Object{{"pitch", pitch}, {"roll", roll}, {"yaw", yaw}}},
        {"coefficients",
         Object{{"cd", ratio(drag, frame.force_scale)},
                {"cl", ratio(lift, frame.force_scale)},
                {"cs", ratio(dot(force, frame.side), frame.force_scale)},
                {"cd_pressure", ratio(dot(s.pressure, frame.drag), frame.force_scale)},
                {"cd_viscous", ratio(dot(s.viscous, frame.drag), frame.force_scale)},
                {"cl_pressure", ratio(dot(s.pressure, frame.lift), frame.force_scale)},
                {"cl_viscous", ratio(dot(s.viscous, frame.lift), frame.force_scale)},
                // cl*A and cd*A need no reference area.
                {"cla", ratio(lift, frame.q)},
                {"cda", ratio(drag, frame.q)},
                {"cm_pitch", ratio(pitch, frame.moment_scale)},
                {"cm_roll", ratio(roll, frame.moment_scale)},
                {"cm_yaw", ratio(yaw, frame.moment_scale)}}},
        {"downforce_share", total_downforce != 0 ? Json(-lift / total_downforce) : Json(nullptr)}};
  };
  Json::Array rows, ungrouped;
  for (const auto &z : loads.zones)
    if (std::find(car.begin(), car.end(), z.name) != car.end()) {
      rows.push_back(row(z.name, "zone", {z.name}));
      if (!grouped.contains(z.name))
        ungrouped.push_back(z.name);
    }
  for (const auto &g : spec.groups)
    rows.push_back(row(g.name, "group", g.zones));
  std::vector<std::string> ordered;
  for (const auto &z : loads.zones)
    if (std::find(car.begin(), car.end(), z.name) != car.end())
      ordered.push_back(z.name);
  rows.push_back(row("Total", "total", ordered));

  const auto force = total.force();
  const double total_drag = dot(force, frame.drag);
  Object balance{{"status", "needs_vehicle"}, {"front_axle", nullptr}, {"wheelbase", nullptr},
                 {"front_load", nullptr}, {"rear_load", nullptr}, {"front_percent", nullptr},
                 {"cop_from_front_axle", nullptr}};
  if (spec.vehicle) {
    // Moment about the front contact point; its side component is carried by
    // the rear tyres at one wheelbase.
    const auto &v = *spec.vehicle;
    const auto about_front = frame.moment(total.pressure_moment + total.viscous_moment, force) +
                             cross(spec.moment_centre - v.front_axle, force);
    const double rear = dot(about_front, frame.side) / v.wheelbase, front = total_downforce - rear;
    const bool defined = total_downforce > 0;
    balance = Object{{"status", defined ? "defined" : "no_downforce"},
                     {"front_axle", vector(v.front_axle)},
                     {"wheelbase", v.wheelbase},
                     {"front_load", front},
                     {"rear_load", rear},
                     {"front_percent", defined ? Json(100 * front / total_downforce) : Json(nullptr)},
                     {"cop_from_front_axle", defined ? Json(v.wheelbase * rear / total_downforce) : Json(nullptr)}};
  }
  return Object{{"v", 1},
                {"provenance", loads.provenance},
                {"approximate", loads.approximate},
                {"modified", loads.modified},
                {"iterations", loads.iterations},
                {"reference", Object{{"area", optional(ref.area)},
                                     {"density", optional(ref.density)},
                                     {"velocity", optional(ref.velocity)},
                                     {"length", optional(ref.length)},
                                     {"dynamic_pressure", optional(frame.q)}}},
                {"axes", Object{{"drag", vector(frame.drag)}, {"lift", vector(frame.lift)}, {"side", vector(frame.side)}}},
                {"moment_centre", vector(spec.moment_centre)},
                {"rows", std::move(rows)},
                {"ungrouped", std::move(ungrouped)},
                {"lift_to_drag", total_drag != 0 ? Json(total_downforce / total_drag) : Json(nullptr)},
                {"balance", std::move(balance)}};
}
} // namespace kutta::post
