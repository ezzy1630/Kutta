#include "kutta/setup/kutta_case.hpp"

#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "kutta/setup/convergence.hpp"
#include "kutta/setup/fluent.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>

namespace kutta::setup {
namespace {
using fluent::SettingIssue;
using Object = Json::Object;
constexpr std::string_view kutta_case_rejected = "Kutta case rejected:";

// Fluent's defaults for a new steady, pressure-based, single-phase case, as
// lowering writes them. Sources (docs/kutta-case.md has the full table):
// [UG]/[TG] the Fluent 2026 R1 User's/Theory Guide section; [UG fig] a value
// shown in a UG dialog figure (docs/research/fluent-numerics-and-files.md);
// [O] the value Fluent 22.2 and 23.1 write into new cases where the guides
// state none (the mixing elbow, elbow_param and nozzle examples agree).
namespace defaults {
constexpr double operating_pressure = 101325;              // UG 9.15
constexpr double air_density = 1.225, air_viscosity = 1.7894e-05;  // default fluid "air" [O]
constexpr bool production_limiter = true;                  // TG 4.21, UG 16.16.6
constexpr double production_limiter_clip = 10;             // TG 4.21
constexpr bool low_re = false, curvature = false, compressibility = false;  // UG 16.16.11, TG 4.19, UG 16.6.4
constexpr bool kato_launder = false, vorticity_production = false;           // TG 4.21; [O]
constexpr double min_k = 1e-14, min_omega = 1e-20, min_epsilon = 1e-20;      // [O]
constexpr double max_viscosity_ratio = 1e5;                // [O]
constexpr double courant = 200;                            // UG 37.3.1.5
constexpr double explicit_relaxation = 0.5;                // UG 37.3.1.5 (Coupled, pseudo time off)
constexpr double time_scale_factor = 1;                    // UG 37.14.2.2
constexpr double pseudo_time_step = 1;                     // [O]; used only when not automatic
constexpr int hybrid_iterations = 10;                      // [O]
constexpr double reference_area = 1, reference_length = 1, reference_density = 1.225, reference_velocity = 1,
                 reference_pressure = 0, reference_viscosity = 1.7894e-05, reference_depth = 1;  // UG 42.13 fields, [O] values
constexpr double intensity = 0.05, viscosity_ratio = 10;   // UG 8.4 (Intensity and Viscosity Ratio)
constexpr double length_scale = 1, hydraulic_diameter = 1, k = 1, omega = 1, epsilon = 1, intermittency = 1;  // [O]
constexpr double residual_criterion = 1e-3;                // UG 37.17
// Pseudo-time explicit relaxation [UG fig, O]; classic under-relaxation [UG fig, O].
constexpr RelaxationInput pseudo{0.5, 0.5, 1, 1, 0.75, 0.75, 0.75, 1, 0.75};
constexpr RelaxationInput classic{0.3, 0.7, 1, 1, 0.8, 0.8, 0.8, 1, 0.8};
} // namespace defaults

// ---- numbers and text ----------------------------------------------------------------------

// Shortest text that reads back as the same double.
std::string number_text(double v) {
  char buf[64];
  const auto r = std::to_chars(buf, buf + sizeof buf, v);
  return std::string(buf, r.ptr);
}
// A Scheme real: Fluent's reader types "1" as an integer, "1.0" as a real.
std::string real_text(double v) {
  auto s = number_text(v);
  if (s.find_first_of(".en") == std::string::npos) s += ".0";
  return s;
}
std::string flag_text(bool v) { return v ? "#t" : "#f"; }
std::string scheme_string(std::string_view s) {
  std::string out = "\"";
  for (char c : s) {
    if (c == '"' || c == '\\') out += '\\';
    out += c;
  }
  return out + '"';
}
bool plain_symbol(std::string_view s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), [](unsigned char c) {
    return std::isalnum(c) || std::string_view("-_.:+*/<>=!?").find(char(c)) != std::string_view::npos;
  });
}
std::string symbol_text(std::string_view s) { return plain_symbol(s) ? std::string(s) : scheme_string(s); }

std::size_t edit_distance(std::string_view a, std::string_view b) {
  std::vector<std::size_t> row(b.size() + 1);
  for (std::size_t j = 0; j <= b.size(); ++j) row[j] = j;
  for (std::size_t i = 1; i <= a.size(); ++i) {
    std::size_t prev = row[0];
    row[0] = i;
    for (std::size_t j = 1; j <= b.size(); ++j) {
      const auto cur = row[j];
      row[j] = std::min({row[j] + 1, row[j - 1] + 1, prev + (a[i - 1] == b[j - 1] ? 0 : 1)});
      prev = cur;
    }
  }
  return row[b.size()];
}
std::string suggestion(std::string_view key, std::span<const std::string_view> known) {
  std::string_view best;
  std::size_t distance = 3;
  for (auto k : known)
    if (const auto d = edit_distance(key, k); d < distance) distance = d, best = k;
  return best.empty() ? std::string{} : " (did you mean \"" + std::string(best) + "\"?)";
}

// ---- strict JSON reading -------------------------------------------------------------------

class Errors {
 public:
  void add(std::string path, std::string reason) { list.push_back({std::move(path), std::move(reason)}); }
  bool empty() const { return list.empty(); }
  [[noreturn]] void raise() { throw RejectedSettings(std::move(list), kutta_case_rejected); }
  void finish() {
    if (!list.empty()) raise();
  }
  std::vector<SettingIssue> list;
};

// One JSON object: each read marks its key used; finish() rejects the rest.
class Reader {
 public:
  Reader(const Json& value, std::string path, Errors& errors) : path_(std::move(path)), errors_(errors) {
    object_ = std::get_if<Object>(&value.value);
    if (!object_) errors_.add(where(), "expected an object");
  }
  bool ok() const { return object_ != nullptr; }
  bool has(std::string_view key) const { return object_ && object_->contains(std::string(key)); }
  const Json* raw(std::string_view key) {
    known_.push_back(key);
    if (!object_) return nullptr;
    const auto it = object_->find(std::string(key));
    return it == object_->end() ? nullptr : &it->second;
  }
  std::string at(std::string_view key) const { return path_.empty() ? std::string(key) : path_ + "/" + std::string(key); }
  std::optional<double> number(std::string_view key, std::optional<double> min = {}, std::optional<double> max = {},
                               bool exclusive_min = false) {
    const auto* v = raw(key);
    if (!v) return {};
    const auto* n = std::get_if<double>(&v->value);
    if (!n) {
      errors_.add(at(key), "expected a number");
      return {};
    }
    if (min && (exclusive_min ? *n <= *min : *n < *min)) {
      errors_.add(at(key), std::string("must be ") + (exclusive_min ? "greater than " : "at least ") + number_text(*min));
      return {};
    }
    if (max && *n > *max) {
      errors_.add(at(key), "must be at most " + number_text(*max));
      return {};
    }
    return *n;
  }
  std::optional<double> positive(std::string_view key) { return number(key, 0.0, {}, true); }
  std::optional<int> integer(std::string_view key, int min, int max) {
    const auto n = number(key, double(min), double(max));
    if (n && *n != std::floor(*n)) {
      errors_.add(at(key), "must be a whole number");
      return {};
    }
    return n ? std::optional<int>(int(*n)) : std::nullopt;
  }
  std::optional<bool> boolean(std::string_view key) {
    const auto* v = raw(key);
    if (!v) return {};
    if (const auto* b = std::get_if<bool>(&v->value)) return *b;
    errors_.add(at(key), "expected true or false");
    return {};
  }
  std::optional<std::string> string(std::string_view key) {
    const auto* v = raw(key);
    if (!v) return {};
    if (const auto* s = std::get_if<std::string>(&v->value)) return *s;
    errors_.add(at(key), "expected a string");
    return {};
  }
  template <class E>
  std::optional<E> choice(std::string_view key, std::initializer_list<std::pair<std::string_view, E>> options) {
    const auto s = string(key);
    if (!s) return {};
    std::vector<std::string_view> names;
    for (const auto& [name, value] : options) {
      if (name == *s) return value;
      names.push_back(name);
    }
    std::string list;
    for (auto n : names) list += (list.empty() ? "" : ", ") + std::string(n);
    errors_.add(at(key), "\"" + *s + "\" is not one of: " + list + suggestion(*s, names));
    return {};
  }
  std::optional<Vec3> vector(std::string_view key, bool nonzero = false) {
    const auto* v = raw(key);
    if (!v) return {};
    const auto* a = std::get_if<Json::Array>(&v->value);
    if (!a || a->size() != 3 || !std::all_of(a->begin(), a->end(), [](const Json& x) { return std::holds_alternative<double>(x.value); })) {
      errors_.add(at(key), "expected an array of three numbers [x, y, z]");
      return {};
    }
    Vec3 out{std::get<double>((*a)[0].value), std::get<double>((*a)[1].value), std::get<double>((*a)[2].value)};
    if (nonzero && out == Vec3{}) {
      errors_.add(at(key), "must be a nonzero vector");
      return {};
    }
    return out;
  }
  // Rejects every property no reader asked for.
  void finish() {
    if (!object_) return;
    for (const auto& [key, value] : *object_) {
      (void)value;
      if (std::find(known_.begin(), known_.end(), key) == known_.end())
        errors_.add(at(key), "unknown property" + suggestion(key, known_));
    }
  }
  std::string where() const { return path_.empty() ? "case" : path_; }
  Errors& errors() { return errors_; }

 private:
  const Object* object_ = nullptr;
  std::string path_;
  Errors& errors_;
  std::vector<std::string_view> known_;
};

constexpr std::initializer_list<std::pair<std::string_view, Upwind>> upwind_names{
    {"first-order-upwind", Upwind::first_order}, {"second-order-upwind", Upwind::second_order}};

TurbulenceInput read_turbulence(Reader& parent, std::string_view key) {
  TurbulenceInput t;
  const auto* v = parent.raw(key);
  if (!v) return t;
  Reader r(*v, parent.at(key), parent.errors());
  using M = TurbulenceInput::Method;
  t.method = r.choice<M>("method", {{"intensity-viscosity-ratio", M::intensity_viscosity_ratio},
                                    {"intensity-length-scale", M::intensity_length_scale},
                                    {"intensity-hydraulic-diameter", M::intensity_hydraulic_diameter},
                                    {"k-omega", M::k_omega},
                                    {"k-epsilon", M::k_epsilon}});
  const auto method = t.method.value_or(M::intensity_viscosity_ratio);
  // Only the values the selected method uses may appear.
  constexpr std::string_view method_names[] = {"intensity-viscosity-ratio", "intensity-length-scale",
                                               "intensity-hydraulic-diameter", "k-omega", "k-epsilon"};
  const auto use = [&](std::string_view name, bool used, auto read) {
    if (r.has(name) && !used) {
      (void)r.raw(name);
      r.errors().add(r.at(name), "not used by turbulence method " + std::string(method_names[int(method)]) +
                                     (t.method ? "" : " (the default)"));
      return decltype(read()){};
    }
    return read();
  };
  const bool intensity = method != M::k_omega && method != M::k_epsilon;
  t.intensity = use("intensity", intensity, [&] { return r.number("intensity", 0.0, 1.0); });
  t.viscosity_ratio = use("viscosity_ratio", method == M::intensity_viscosity_ratio, [&] { return r.positive("viscosity_ratio"); });
  t.length_scale = use("length_scale", method == M::intensity_length_scale, [&] { return r.positive("length_scale"); });
  t.hydraulic_diameter = use("hydraulic_diameter", method == M::intensity_hydraulic_diameter, [&] { return r.positive("hydraulic_diameter"); });
  t.k = use("k", method == M::k_omega || method == M::k_epsilon, [&] { return r.number("k", 0.0); });
  t.omega = use("omega", method == M::k_omega, [&] { return r.positive("omega"); });
  t.epsilon = use("epsilon", method == M::k_epsilon, [&] { return r.positive("epsilon"); });
  r.finish();
  return t;
}

BoundaryInput read_boundary(const Json& value, const std::string& path, Errors& errors) {
  Reader r(value, path, errors);
  enum class Type { wall, velocity_inlet, pressure_outlet, symmetry, interior };
  const auto type = r.choice<Type>("type", {{"wall", Type::wall}, {"velocity-inlet", Type::velocity_inlet},
                                            {"pressure-outlet", Type::pressure_outlet}, {"symmetry", Type::symmetry},
                                            {"interior", Type::interior}});
  if (!type) {
    if (r.ok() && !r.has("type")) errors.add(r.at("type"), "required: wall, velocity-inlet, pressure-outlet, symmetry or interior");
    return SymmetryInput{};
  }
  BoundaryInput out;
  switch (*type) {
    case Type::wall: {
      WallInput w;
      using M = WallInput::Motion;
      w.motion = r.choice<M>("motion", {{"stationary", M::stationary}, {"moving", M::moving}, {"rotating", M::rotating}});
      const auto motion = w.motion.value_or(M::stationary);
      const auto only = [&](std::string_view key, M needed) {
        if (r.has(key) && motion != needed) {
          (void)r.raw(key);
          errors.add(r.at(key), std::string("applies only to motion \"") + (needed == M::moving ? "moving" : "rotating") + "\"");
          return false;
        }
        return true;
      };
      if (only("speed", M::moving)) w.speed = r.number("speed");
      if (only("direction", M::moving)) w.direction = r.vector("direction", true);
      if (only("velocity", M::moving)) w.velocity = r.vector("velocity");
      if (only("angular_speed", M::rotating)) w.angular_speed = r.number("angular_speed");
      if (only("axis_origin", M::rotating)) w.axis_origin = r.vector("axis_origin");
      if (only("axis_direction", M::rotating)) w.axis_direction = r.vector("axis_direction", true);
      if (r.has("relative") && motion == M::stationary) {
        (void)r.raw("relative");
        errors.add(r.at("relative"), "applies only to moving and rotating walls");
      } else
        w.relative = r.boolean("relative");
      if (w.velocity && (w.speed || w.direction)) errors.add(r.at("velocity"), "give either velocity components or speed and direction, not both");
      out = w;
      break;
    }
    case Type::velocity_inlet: {
      VelocityInletInput v;
      v.speed = r.number("speed", 0.0);
      v.direction = r.vector("direction", true);
      v.velocity = r.vector("velocity");
      v.gauge_pressure = r.number("gauge_pressure");
      v.turbulence = read_turbulence(r, "turbulence");
      if (!r.has("turbulence")) v.turbulence.reset();
      v.intermittency = r.number("intermittency", 0.0, 1.0);
      if (v.velocity && (v.speed || v.direction)) errors.add(r.at("velocity"), "give either velocity components or speed (and direction), not both");
      if (v.direction && !v.speed) errors.add(r.at("direction"), "requires speed");
      out = v;
      break;
    }
    case Type::pressure_outlet: {
      PressureOutletInput p;
      using B = PressureOutletInput::Backflow;
      p.gauge_pressure = r.number("gauge_pressure");
      p.backflow_pressure = r.choice<B>("backflow_pressure", {{"total", B::total}, {"static", B::static_pressure}});
      p.prevent_reverse_flow = r.boolean("prevent_reverse_flow");
      p.turbulence = read_turbulence(r, "turbulence");
      if (!r.has("turbulence")) p.turbulence.reset();
      p.intermittency = r.number("intermittency", 0.0, 1.0);
      out = p;
      break;
    }
    case Type::symmetry: out = SymmetryInput{}; break;
    case Type::interior: out = InteriorInput{}; break;
  }
  r.finish();
  return out;
}

RelaxationInput read_relaxation(Reader& parent, std::string_view key) {
  RelaxationInput out;
  const auto* v = parent.raw(key);
  if (!v) return out;
  Reader r(*v, parent.at(key), parent.errors());
  const auto factor = [&](std::string_view name) { return r.number(name, 0.0, 1.0, true); };
  out.pressure = factor("pressure");
  out.momentum = factor("momentum");
  out.density = factor("density");
  out.body_force = factor("body_force");
  out.k = factor("k");
  out.omega = factor("omega");
  out.epsilon = factor("epsilon");
  out.turbulent_viscosity = factor("turbulent_viscosity");
  out.intermittency = factor("intermittency");
  r.finish();
  return out;
}
} // namespace

bool is_kutta_case(const std::filesystem::path& p) {
  auto ext = p.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return char(std::tolower(c)); });
  return ext == ".json";
}

KuttaCase parse_kutta_case(std::string_view text) {
  Json root;
  try {
    root = parse_json(text);
  } catch (const std::exception& e) {
    throw RejectedSettings({{"case", std::string("invalid JSON: ") + e.what()}}, kutta_case_rejected);
  }
  Errors errors;
  KuttaCase c;
  Reader r(root, "", errors);
  if (!r.ok()) errors.raise();
  if (const auto version = r.number("kutta_case"); version != double(kutta_case_version))
    errors.add("kutta_case", r.has("kutta_case") ? "unsupported version; this Kutta reads version 1"
                                                  : "required: 1 (the Kutta case format version)");
  if (const auto mesh = r.string("mesh"); mesh && !mesh->empty())
    c.mesh = std::filesystem::path(std::u8string(mesh->begin(), mesh->end()));
  else
    errors.add("mesh", "required: the mesh file, relative to this file (.msh, .msh.h5, .cas or .cas.h5)");

  if (const auto* m = r.raw("models")) {
    Reader models(*m, "models", errors);
    const auto viscous = models.choice<ViscousModel>("viscous", {{"laminar", ViscousModel::laminar},
                                                                 {"k-omega-sst", ViscousModel::k_omega_sst},
                                                                 {"k-epsilon-standard", ViscousModel::k_epsilon_standard}});
    if (viscous) c.models.viscous = *viscous;
    else if (models.ok() && !models.has("viscous")) errors.add("models/viscous", "required: laminar, k-omega-sst or k-epsilon-standard");
    if (const auto t = models.choice<bool>("transition", {{"none", false}, {"intermittency", true}})) c.models.intermittency_transition = *t;
    if (const auto* sst = models.raw("k_omega_sst")) {
      Reader s(*sst, "models/k_omega_sst", errors);
      auto& o = c.models.k_omega_sst;
      o.production_limiter = s.boolean("production_limiter");
      o.production_limiter_clip_factor = s.positive("production_limiter_clip_factor");
      o.low_re_corrections = s.boolean("low_re_corrections");
      o.curvature_correction = s.boolean("curvature_correction");
      o.compressibility_effects = s.boolean("compressibility_effects");
      o.kato_launder = s.boolean("kato_launder");
      o.vorticity_based_production = s.boolean("vorticity_based_production");
      s.finish();
      if (c.models.viscous != ViscousModel::k_omega_sst) errors.add("models/k_omega_sst", "applies only to viscous k-omega-sst");
    }
    if (const auto* ke = models.raw("k_epsilon_standard")) {
      Reader s(*ke, "models/k_epsilon_standard", errors);
      c.models.k_epsilon.compressibility_effects = s.boolean("compressibility_effects");
      s.finish();
      if (c.models.viscous != ViscousModel::k_epsilon_standard) errors.add("models/k_epsilon_standard", "applies only to viscous k-epsilon-standard");
    }
    models.finish();
    if (c.models.intermittency_transition == true && c.models.viscous != ViscousModel::k_omega_sst)
      errors.add("models/transition", "the intermittency transition model runs with k-omega-sst only");
  } else
    errors.add("models", "required: at least {\"viscous\": ...}");

  if (const auto* m = r.raw("material")) {
    Reader s(*m, "material", errors);
    c.material.name = s.string("name");
    c.material.density = s.positive("density");
    c.material.viscosity = s.positive("viscosity");
    s.finish();
  }
  c.operating_pressure = r.number("operating_pressure", 0.0);

  if (const auto* b = r.raw("boundaries")) {
    if (const auto* zones = std::get_if<Object>(&b->value))
      for (const auto& [name, value] : *zones) c.boundaries.emplace(name, read_boundary(value, "boundaries/" + name, errors));
    else
      errors.add("boundaries", "expected an object keyed by mesh zone name");
  } else
    errors.add("boundaries", "required: a condition for every boundary zone of the mesh");

  if (const auto* m = r.raw("reference_values")) {
    Reader s(*m, "reference_values", errors);
    auto& v = c.reference;
    v.area = s.positive("area");
    v.length = s.positive("length");
    v.density = s.positive("density");
    v.velocity = s.positive("velocity");
    v.pressure = s.number("pressure");
    v.viscosity = s.positive("viscosity");
    v.depth = s.positive("depth");
    s.finish();
  }
  if (const auto* m = r.raw("methods")) {
    Reader s(*m, "methods", errors);
    auto& v = c.methods;
    v.coupling = s.choice<CouplingScheme>("pressure_velocity_coupling", {{"simple", CouplingScheme::simple},
                                                                          {"simplec", CouplingScheme::simplec},
                                                                          {"coupled", CouplingScheme::coupled}});
    v.pseudo_transient = s.boolean("pseudo_transient");
    v.gradient = s.choice<GradientMethod>("gradient", {{"least-squares-cell-based", GradientMethod::least_squares_cell},
                                                       {"green-gauss-cell-based", GradientMethod::green_gauss_cell}});
    v.pressure = s.choice<PressureInterpolation>("pressure", {{"standard", PressureInterpolation::standard},
                                                              {"second-order", PressureInterpolation::second_order}});
    v.momentum = s.choice<Upwind>("momentum", upwind_names);
    v.turbulent_kinetic_energy = s.choice<Upwind>("turbulent_kinetic_energy", upwind_names);
    v.specific_dissipation_rate = s.choice<Upwind>("specific_dissipation_rate", upwind_names);
    v.turbulent_dissipation_rate = s.choice<Upwind>("turbulent_dissipation_rate", upwind_names);
    v.intermittency = s.choice<Upwind>("intermittency", upwind_names);
    v.flux_type = s.choice<FluxType>("flux_type", {{"auto", FluxType::automatic}, {"rhie-chow-momentum", FluxType::momentum},
                                                    {"rhie-chow-distance", FluxType::distance}});
    v.gradient_limiter = s.choice<GradientLimiter>("gradient_limiter", {{"standard", GradientLimiter::standard},
                                                                         {"differentiable", GradientLimiter::differentiable},
                                                                         {"none", GradientLimiter::none}});
    v.warped_face_gradient_correction = s.boolean("warped_face_gradient_correction");
    v.high_order_term_relaxation = s.boolean("high_order_term_relaxation");
    s.finish();
  }
  if (const auto* m = r.raw("controls")) {
    Reader s(*m, "controls", errors);
    auto& v = c.controls;
    v.flow_courant_number = s.positive("flow_courant_number");
    v.pseudo_time_scale_factor = s.positive("pseudo_time_scale_factor");
    if (const auto* step = s.raw("pseudo_time_step"); step && std::get_if<std::string>(&step->value)) {
      if (std::get<std::string>(step->value) != "automatic") errors.add("controls/pseudo_time_step", "expected \"automatic\" or a time step in seconds");
    } else if (step)
      v.pseudo_time_step = s.positive("pseudo_time_step");
    v.pseudo_relaxation = read_relaxation(s, "pseudo_relaxation");
    v.relaxation = read_relaxation(s, "relaxation");
    v.explicit_relaxation = read_relaxation(s, "explicit_relaxation");
    if (v.explicit_relaxation != RelaxationInput{.pressure = v.explicit_relaxation.pressure, .momentum = v.explicit_relaxation.momentum})
      errors.add("controls/explicit_relaxation", "only pressure and momentum have Coupled explicit relaxation");
    if (const auto* l = s.raw("limits")) {
      Reader limits(*l, "controls/limits", errors);
      v.min_turbulent_kinetic_energy = limits.positive("min_turbulent_kinetic_energy");
      v.min_specific_dissipation_rate = limits.positive("min_specific_dissipation_rate");
      v.min_turbulent_dissipation_rate = limits.positive("min_turbulent_dissipation_rate");
      v.max_turbulent_viscosity_ratio = limits.positive("max_turbulent_viscosity_ratio");
      limits.finish();
    }
    s.finish();
  }
  if (const auto* m = r.raw("initialization")) {
    Reader s(*m, "initialization", errors);
    auto& v = c.initialization;
    using M = KuttaCase::Initialization::Method;
    v.method = s.choice<M>("method", {{"hybrid", M::hybrid}, {"standard", M::standard}});
    const bool standard = v.method == M::standard;
    if (s.has("hybrid_iterations") && standard) {
      (void)s.raw("hybrid_iterations");
      errors.add("initialization/hybrid_iterations", "applies only to hybrid initialization");
    } else
      v.hybrid_iterations = s.integer("hybrid_iterations", 1, 100000);
    for (auto key : {"velocity", "gauge_pressure", "k", "omega", "epsilon", "intermittency"})
      if (s.has(key) && !standard) {
        (void)s.raw(key);
        errors.add(s.at(key), "applies only to standard initialization (\"method\": \"standard\")");
      }
    if (standard) {
      v.velocity = s.vector("velocity");
      v.gauge_pressure = s.number("gauge_pressure");
      v.k = s.number("k", 0.0);
      v.omega = s.positive("omega");
      v.epsilon = s.positive("epsilon");
      v.intermittency = s.number("intermittency", 0.0, 1.0);
    }
    s.finish();
  }
  if (const auto* m = r.raw("reports")) {
    if (const auto* reports = std::get_if<Object>(&m->value)) {
      for (const auto& [name, value] : *reports) {
        Reader s(value, "reports/" + name, errors);
        ForceReportInput f;
        using K = ForceReportInput::Kind;
        const auto kind = s.choice<K>("type", {{"force", K::force}, {"drag", K::drag}, {"lift", K::lift}});
        if (kind) f.kind = *kind;
        else if (s.ok() && !s.has("type")) errors.add(s.at("type"), "required: force, drag or lift");
        if (const auto* z = s.raw("zones"); z && std::get_if<Json::Array>(&z->value)) {
          for (const auto& x : std::get<Json::Array>(z->value))
            if (const auto* str = std::get_if<std::string>(&x.value)) f.zones.push_back(*str);
            else errors.add(s.at("zones"), "expected zone names");
        } else
          errors.add(s.at("zones"), "required: the wall zones whose force is reported");
        if (s.ok() && f.zones.empty() && s.has("zones")) errors.add(s.at("zones"), "must name at least one wall zone");
        f.direction = s.vector("direction", true);
        f.coefficient = s.boolean("coefficient");
        s.finish();
        if (s.ok()) c.reports.emplace(name, std::move(f));
      }
    } else
      errors.add("reports", "expected an object keyed by report name");
  }
  if (const auto* m = r.raw("monitors")) {
    Reader s(*m, "monitors", errors);
    if (const auto* res = s.raw("residuals")) {
      if (const auto* rows = std::get_if<Object>(&res->value))
        for (const auto& [name, value] : *rows) {
          Reader row(value, "monitors/residuals/" + name, errors);
          ResidualInput in;
          in.check = row.boolean("check");
          in.absolute_criterion = row.number("absolute_criterion", 0.0);
          row.finish();
          c.residuals.emplace(name, in);
        }
      else
        errors.add("monitors/residuals", "expected an object keyed by residual name");
    }
    s.finish();
  }
  r.finish();
  errors.finish();
  return c;
}

KuttaCase read_kutta_case_file(const std::filesystem::path& path) {
  std::error_code ec;
  const auto size = std::filesystem::file_size(path, ec);
  if (ec) throw Error("cannot read Kutta case " + path.string() + ": " + ec.message());
  if (size > 1024 * 1024) throw RejectedSettings({{"case", "Kutta case file exceeds 1 MiB"}}, kutta_case_rejected);
  std::ifstream in(path, std::ios::binary);
  if (!in) throw Error("cannot open Kutta case " + path.string());
  return parse_kutta_case(std::string{std::istreambuf_iterator<char>(in), {}});
}

// ---- writing -------------------------------------------------------------------------------
namespace {
// An ordered JSON value for the file writer (Json::Object sorts its keys).
struct Out {
  enum class Kind { value, object, array } kind = Kind::value;
  std::string text;  // scalar JSON text
  std::vector<std::pair<std::string, Out>> members;
  std::vector<Out> items;
  static Out num(double v) { return {Kind::value, number_text(v), {}, {}}; }
  static Out str(std::string_view s) {
    std::ostringstream o;
    o << '"';
    for (unsigned char ch : s) {
      if (ch == '"' || ch == '\\') o << '\\' << ch;
      else if (ch < 32) {
        char buf[8];
        std::snprintf(buf, sizeof buf, "\\u%04x", ch);
        o << buf;
      } else o << ch;
    }
    o << '"';
    return {Kind::value, o.str(), {}, {}};
  }
  static Out boolean(bool b) { return {Kind::value, b ? "true" : "false", {}, {}}; }
  static Out vec(const Vec3& v) {
    Out a{Kind::array, {}, {}, {}};
    for (double x : v) a.items.push_back(num(x));
    return a;
  }
  static Out object() { return {Kind::object, {}, {}, {}}; }
  Out& add(std::string key, Out v) {
    members.emplace_back(std::move(key), std::move(v));
    return members.back().second;
  }
  template <class T>
  void opt(std::string key, const std::optional<T>& v) {
    if (!v) return;
    if constexpr (std::is_same_v<T, double>) add(std::move(key), num(*v));
    else if constexpr (std::is_same_v<T, int>) add(std::move(key), num(*v));
    else if constexpr (std::is_same_v<T, bool>) add(std::move(key), boolean(*v));
    else if constexpr (std::is_same_v<T, Vec3>) add(std::move(key), vec(*v));
    else add(std::move(key), str(*v));
  }
  bool empty() const { return kind == Kind::object ? members.empty() : kind == Kind::array ? items.empty() : false; }
  // Arrays of scalars and objects with only scalar members stay on one line.
  bool flat() const {
    if (kind == Kind::value) return true;
    if (kind == Kind::array) return std::all_of(items.begin(), items.end(), [](const Out& o) { return o.kind == Kind::value; });
    return members.size() <= 4 && std::all_of(members.begin(), members.end(), [](const auto& m) {
             return m.second.kind == Kind::value || (m.second.kind == Kind::array && m.second.flat());
           });
  }
  void write(std::ostream& s, int indent) const {
    if (kind == Kind::value) {
      s << text;
      return;
    }
    const bool one_line = flat();
    const char open = kind == Kind::object ? '{' : '[', close = kind == Kind::object ? '}' : ']';
    s << open;
    const auto n = kind == Kind::object ? members.size() : items.size();
    for (std::size_t i = 0; i < n; ++i) {
      if (i) s << ',';
      if (one_line) s << (i ? " " : "");
      else s << '\n' << std::string(std::size_t(indent + 2), ' ');
      if (kind == Kind::object) {
        s << str(members[i].first).text << ": ";
        members[i].second.write(s, indent + 2);
      } else
        items[i].write(s, indent + 2);
    }
    if (!one_line && n) s << '\n' << std::string(std::size_t(indent), ' ');
    s << close;
  }
};
std::string_view upwind_text(Upwind u) { return u == Upwind::first_order ? "first-order-upwind" : "second-order-upwind"; }
Out turbulence_out(const TurbulenceInput& t) {
  using M = TurbulenceInput::Method;
  auto o = Out::object();
  if (t.method) {
    constexpr std::string_view names[] = {"intensity-viscosity-ratio", "intensity-length-scale", "intensity-hydraulic-diameter", "k-omega", "k-epsilon"};
    o.add("method", Out::str(names[int(*t.method)]));
  }
  (void)M::k_omega;
  o.opt("intensity", t.intensity);
  o.opt("viscosity_ratio", t.viscosity_ratio);
  o.opt("length_scale", t.length_scale);
  o.opt("hydraulic_diameter", t.hydraulic_diameter);
  o.opt("k", t.k);
  o.opt("omega", t.omega);
  o.opt("epsilon", t.epsilon);
  return o;
}
Out relaxation_out(const RelaxationInput& r) {
  auto o = Out::object();
  o.opt("pressure", r.pressure);
  o.opt("momentum", r.momentum);
  o.opt("density", r.density);
  o.opt("body_force", r.body_force);
  o.opt("k", r.k);
  o.opt("omega", r.omega);
  o.opt("epsilon", r.epsilon);
  o.opt("turbulent_viscosity", r.turbulent_viscosity);
  o.opt("intermittency", r.intermittency);
  return o;
}
} // namespace

std::string kutta_case_text(const KuttaCase& c) {
  auto root = Out::object();
  root.add("kutta_case", Out::num(kutta_case_version));
  const auto mesh = c.mesh.generic_u8string();
  root.add("mesh", Out::str(std::string(mesh.begin(), mesh.end())));
  {
    auto& m = root.add("models", Out::object());
    constexpr std::string_view viscous[] = {"laminar", "k-omega-sst", "k-epsilon-standard"};
    m.add("viscous", Out::str(viscous[int(c.models.viscous)]));
    if (c.models.intermittency_transition) m.add("transition", Out::str(*c.models.intermittency_transition ? "intermittency" : "none"));
    auto sst = Out::object();
    const auto& s = c.models.k_omega_sst;
    sst.opt("production_limiter", s.production_limiter);
    sst.opt("production_limiter_clip_factor", s.production_limiter_clip_factor);
    sst.opt("low_re_corrections", s.low_re_corrections);
    sst.opt("curvature_correction", s.curvature_correction);
    sst.opt("compressibility_effects", s.compressibility_effects);
    sst.opt("kato_launder", s.kato_launder);
    sst.opt("vorticity_based_production", s.vorticity_based_production);
    if (!sst.empty()) m.add("k_omega_sst", std::move(sst));
    auto ke = Out::object();
    ke.opt("compressibility_effects", c.models.k_epsilon.compressibility_effects);
    if (!ke.empty()) m.add("k_epsilon_standard", std::move(ke));
  }
  {
    auto m = Out::object();
    m.opt("name", c.material.name);
    m.opt("density", c.material.density);
    m.opt("viscosity", c.material.viscosity);
    if (!m.empty()) root.add("material", std::move(m));
  }
  root.opt("operating_pressure", c.operating_pressure);
  {
    auto& b = root.add("boundaries", Out::object());
    for (const auto& [name, condition] : c.boundaries) {
      auto z = Out::object();
      std::visit([&](const auto& v) {
        using T = std::decay_t<decltype(v)>;
        if constexpr (std::is_same_v<T, WallInput>) {
          z.add("type", Out::str("wall"));
          if (v.motion) {
            constexpr std::string_view motions[] = {"stationary", "moving", "rotating"};
            z.add("motion", Out::str(motions[int(*v.motion)]));
          }
          z.opt("relative", v.relative);
          z.opt("speed", v.speed);
          z.opt("direction", v.direction);
          z.opt("velocity", v.velocity);
          z.opt("angular_speed", v.angular_speed);
          z.opt("axis_origin", v.axis_origin);
          z.opt("axis_direction", v.axis_direction);
        } else if constexpr (std::is_same_v<T, VelocityInletInput>) {
          z.add("type", Out::str("velocity-inlet"));
          z.opt("speed", v.speed);
          z.opt("direction", v.direction);
          z.opt("velocity", v.velocity);
          z.opt("gauge_pressure", v.gauge_pressure);
          if (v.turbulence) z.add("turbulence", turbulence_out(*v.turbulence));
          z.opt("intermittency", v.intermittency);
        } else if constexpr (std::is_same_v<T, PressureOutletInput>) {
          z.add("type", Out::str("pressure-outlet"));
          z.opt("gauge_pressure", v.gauge_pressure);
          if (v.backflow_pressure) z.add("backflow_pressure", Out::str(*v.backflow_pressure == PressureOutletInput::Backflow::total ? "total" : "static"));
          z.opt("prevent_reverse_flow", v.prevent_reverse_flow);
          if (v.turbulence) z.add("turbulence", turbulence_out(*v.turbulence));
          z.opt("intermittency", v.intermittency);
        } else if constexpr (std::is_same_v<T, SymmetryInput>) z.add("type", Out::str("symmetry"));
        else z.add("type", Out::str("interior"));
      }, condition);
      b.add(name, std::move(z));
    }
  }
  {
    auto r = Out::object();
    const auto& v = c.reference;
    r.opt("area", v.area);
    r.opt("length", v.length);
    r.opt("density", v.density);
    r.opt("velocity", v.velocity);
    r.opt("pressure", v.pressure);
    r.opt("viscosity", v.viscosity);
    r.opt("depth", v.depth);
    if (!r.empty()) root.add("reference_values", std::move(r));
  }
  {
    auto m = Out::object();
    const auto& v = c.methods;
    if (v.coupling) {
      constexpr std::string_view names[] = {"simple", "simplec", "coupled"};
      m.add("pressure_velocity_coupling", Out::str(names[int(*v.coupling)]));
    }
    m.opt("pseudo_transient", v.pseudo_transient);
    if (v.gradient) m.add("gradient", Out::str(*v.gradient == GradientMethod::least_squares_cell ? "least-squares-cell-based" : "green-gauss-cell-based"));
    if (v.pressure) m.add("pressure", Out::str(*v.pressure == PressureInterpolation::standard ? "standard" : "second-order"));
    for (auto [key, value] : {std::pair{"momentum", v.momentum}, {"turbulent_kinetic_energy", v.turbulent_kinetic_energy},
                              {"specific_dissipation_rate", v.specific_dissipation_rate},
                              {"turbulent_dissipation_rate", v.turbulent_dissipation_rate}, {"intermittency", v.intermittency}})
      if (value) m.add(key, Out::str(upwind_text(*value)));
    if (v.flux_type) {
      constexpr std::string_view names[] = {"auto", "rhie-chow-momentum", "rhie-chow-distance"};
      m.add("flux_type", Out::str(names[int(*v.flux_type)]));
    }
    if (v.gradient_limiter) {
      constexpr std::string_view names[] = {"standard", "differentiable", "none"};
      m.add("gradient_limiter", Out::str(names[int(*v.gradient_limiter)]));
    }
    m.opt("warped_face_gradient_correction", v.warped_face_gradient_correction);
    m.opt("high_order_term_relaxation", v.high_order_term_relaxation);
    if (!m.empty()) root.add("methods", std::move(m));
  }
  {
    auto m = Out::object();
    const auto& v = c.controls;
    m.opt("flow_courant_number", v.flow_courant_number);
    m.opt("pseudo_time_scale_factor", v.pseudo_time_scale_factor);
    m.opt("pseudo_time_step", v.pseudo_time_step);
    for (auto [key, value] : {std::pair{"pseudo_relaxation", &v.pseudo_relaxation}, {"relaxation", &v.relaxation},
                              {"explicit_relaxation", &v.explicit_relaxation}})
      if (auto r = relaxation_out(*value); !r.empty()) m.add(key, std::move(r));
    auto limits = Out::object();
    limits.opt("min_turbulent_kinetic_energy", v.min_turbulent_kinetic_energy);
    limits.opt("min_specific_dissipation_rate", v.min_specific_dissipation_rate);
    limits.opt("min_turbulent_dissipation_rate", v.min_turbulent_dissipation_rate);
    limits.opt("max_turbulent_viscosity_ratio", v.max_turbulent_viscosity_ratio);
    if (!limits.empty()) m.add("limits", std::move(limits));
    if (!m.empty()) root.add("controls", std::move(m));
  }
  {
    auto m = Out::object();
    const auto& v = c.initialization;
    if (v.method) m.add("method", Out::str(*v.method == KuttaCase::Initialization::Method::hybrid ? "hybrid" : "standard"));
    m.opt("hybrid_iterations", v.hybrid_iterations);
    m.opt("velocity", v.velocity);
    m.opt("gauge_pressure", v.gauge_pressure);
    m.opt("k", v.k);
    m.opt("omega", v.omega);
    m.opt("epsilon", v.epsilon);
    m.opt("intermittency", v.intermittency);
    if (!m.empty()) root.add("initialization", std::move(m));
  }
  if (!c.reports.empty()) {
    auto& m = root.add("reports", Out::object());
    for (const auto& [name, r] : c.reports) {
      auto o = Out::object();
      constexpr std::string_view kinds[] = {"force", "drag", "lift"};
      o.add("type", Out::str(kinds[int(r.kind)]));
      auto zones = Out{Out::Kind::array, {}, {}, {}};
      for (const auto& z : r.zones) zones.items.push_back(Out::str(z));
      o.add("zones", std::move(zones));
      o.opt("direction", r.direction);
      o.opt("coefficient", r.coefficient);
      m.add(name, std::move(o));
    }
  }
  if (!c.residuals.empty()) {
    auto& m = root.add("monitors", Out::object());
    auto& rows = m.add("residuals", Out::object());
    for (const auto& [name, r] : c.residuals) {
      auto o = Out::object();
      o.opt("check", r.check);
      o.opt("absolute_criterion", r.absolute_criterion);
      rows.add(name, std::move(o));
    }
  }
  std::ostringstream s;
  root.write(s, 0);
  s << '\n';
  return s.str();
}

// ---- mesh facts ----------------------------------------------------------------------------
MeshInfo mesh_info(const PolyMesh& mesh) {
  MeshInfo out;
  out.cells = std::size_t(mesh.n_cells);
  out.faces = std::size_t(mesh.n_faces());
  out.nodes = std::size_t(mesh.n_points());
  // read_mesh extrudes a 2D mesh one cell thick between two zone-less
  // symmetry caps named front and back (kutta/io/mesh_io.hpp).
  int caps = 0;
  for (const auto& p : mesh.patches)
    if (p.source_zone_id < 0 && p.kind == PatchKind::Symmetry && p.size == mesh.n_cells && (p.name == "front" || p.name == "back")) ++caps;
  out.dimension = caps == 2 ? 2 : 3;
  if (out.dimension == 3) {
    const auto cf = build_cell_faces(mesh);
    for (Index cell = 0; cell < mesh.n_cells && !out.polyhedral; ++cell) {
      int tri = 0, quad = 0, other = 0;
      for (auto i = cf.offsets[cell]; i < cf.offsets[cell + 1]; ++i) {
        const auto n = detail::nodes(mesh, cf.faces[i]).size();
        (n == 3 ? tri : n == 4 ? quad : other)++;
      }
      const bool tet = tri == 4 && quad == 0, pyramid = tri == 4 && quad == 1, wedge = tri == 2 && quad == 3, hex = tri == 0 && quad == 6;
      out.polyhedral = other > 0 || !(tet || pyramid || wedge || hex);
    }
  }
  // Readers name a cell zone's type only when the file does; otherwise they keep
  // Fluent's element-zone code (section 12: 1 fluid, 17 solid).
  const auto cell_type = [](const std::string& t) {
    return t == "fluent-cell-zone-type-1" ? std::string("fluid") : t == "fluent-cell-zone-type-17" ? std::string("solid") : t;
  };
  for (const auto& z : mesh.cell_zones)
    out.zones.push_back({z.source_zone_id, z.name, cell_type(z.source_type), MeshZoneInfo::Kind::cell, z.cells.size()});
  for (const auto& z : mesh.internal_face_zones)
    out.zones.push_back({z.source_zone_id, z.name, z.source_type, MeshZoneInfo::Kind::interior, z.faces.size()});
  for (const auto& p : mesh.patches)
    if (p.source_zone_id >= 0) out.zones.push_back({p.source_zone_id, p.name, p.source_type, MeshZoneInfo::Kind::boundary, std::size_t(p.size)});
  std::sort(out.zones.begin(), out.zones.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
  return out;
}

// ---- lowering ------------------------------------------------------------------------------
namespace {
fluent::ViscousModel fluent_model(ViscousModel m) {
  return m == ViscousModel::laminar ? fluent::ViscousModel::laminar
         : m == ViscousModel::k_omega_sst ? fluent::ViscousModel::kw_sst : fluent::ViscousModel::ke_standard;
}
// Alist text in Fluent's own syntax, parsed by the same reader as a .cas.
class Alist {
 public:
  void atom(std::string_view key, std::string value) { s_ << "\n (" << key << " . " << value << ')'; }
  void list(std::string_view key, std::string value) { s_ << "\n (" << key << ' ' << value << ')'; }
  void real(std::string_view key, double v) { atom(key, real_text(v)); }
  void integer(std::string_view key, long long v) { atom(key, std::to_string(v)); }
  void flag(std::string_view key, bool v) { atom(key, flag_text(v)); }
  void vector(std::array<std::string_view, 3> keys, const Vec3& v) {
    for (int i = 0; i < 3; ++i) real(keys[std::size_t(i)], v[std::size_t(i)]);
  }
  std::string text() const { return '(' + s_.str() + ')'; }

 private:
  std::ostringstream s_;
};
struct Solved {
  bool laminar, sst, ke, transition;
};
void turbulence_entries(Alist& a, const std::optional<TurbulenceInput>& in, const Solved& m) {
  using M = TurbulenceInput::Method;
  const TurbulenceInput t = in.value_or(TurbulenceInput{});
  const auto method = t.method.value_or(M::intensity_viscosity_ratio);
  // ke-spec: 0 k and omega/epsilon, 1 intensity and length scale,
  // 2 intensity and viscosity ratio, 3 intensity and hydraulic diameter.
  const int spec = method == M::k_omega || method == M::k_epsilon ? 0 : method == M::intensity_length_scale ? 1
                   : method == M::intensity_viscosity_ratio ? 2 : 3;
  (void)m;
  a.integer("ke-spec", spec);
  a.real("turb-intensity", t.intensity.value_or(defaults::intensity));
  a.real("turb-viscosity-ratio", t.viscosity_ratio.value_or(defaults::viscosity_ratio));
  a.real("turb-length-scale", t.length_scale.value_or(defaults::length_scale));
  a.real("turb-hydraulic-diam", t.hydraulic_diameter.value_or(defaults::hydraulic_diameter));
  a.real("k", t.k.value_or(defaults::k));
  a.real("o", t.omega.value_or(defaults::omega));
  a.real("e", t.epsilon.value_or(defaults::epsilon));
}
double relax(const RelaxationInput& in, const RelaxationInput& fallback, std::optional<double> RelaxationInput::*member) {
  return (in.*member).value_or(*(fallback.*member));
}
}  // namespace

fluent::FluentCase lower_kutta_case(const KuttaCase& c, const MeshInfo& mesh) {
  Errors errors;
  const bool laminar = c.models.viscous == ViscousModel::laminar, sst = c.models.viscous == ViscousModel::k_omega_sst,
             ke = c.models.viscous == ViscousModel::k_epsilon_standard;
  const bool transition = c.models.intermittency_transition.value_or(false);
  const Solved solved{laminar, sst, ke, transition};
  const bool planar = mesh.dimension == 2;
  const auto coupling = c.methods.coupling.value_or(CouplingScheme::coupled);
  // UG 37.14: Global Time Step is the default with the Coupled scheme.
  const bool pseudo = c.methods.pseudo_transient.value_or(coupling == CouplingScheme::coupled);

  // Zones: every boundary zone of the mesh needs a condition; interior zones
  // are interior; every cell zone is fluid with the case's material.
  std::map<std::string, const MeshZoneInfo*> by_name;
  for (const auto& z : mesh.zones) by_name.emplace(z.name, &z);
  std::string known;
  for (const auto& z : mesh.zones)
    if (z.kind == MeshZoneInfo::Kind::boundary) known += (known.empty() ? "" : ", ") + z.name;
  for (const auto& [name, b] : c.boundaries) {
    const auto it = by_name.find(name);
    if (it == by_name.end()) {
      std::vector<std::string_view> names;
      for (const auto& z : mesh.zones) names.push_back(z.name);
      errors.add("boundaries/" + name, "the mesh has no zone named \"" + name + "\"" + suggestion(name, names) + "; boundary zones: " + known);
      continue;
    }
    const auto kind = it->second->kind;
    if (kind == MeshZoneInfo::Kind::cell) errors.add("boundaries/" + name, "is a cell zone; cell zones are fluid and take no boundary condition");
    else if (kind == MeshZoneInfo::Kind::interior && !std::holds_alternative<InteriorInput>(b))
      errors.add("boundaries/" + name, "is an interior face zone of the mesh; only \"interior\" applies");
    else if (kind == MeshZoneInfo::Kind::boundary && std::holds_alternative<InteriorInput>(b))
      errors.add("boundaries/" + name, "is a boundary of the mesh; interior applies only to interior face zones");
  }
  for (const auto& z : mesh.zones) {
    if (z.kind == MeshZoneInfo::Kind::boundary && !c.boundaries.contains(z.name))
      errors.add("boundaries/" + z.name, "required: mesh boundary zone \"" + z.name + "\" (" + z.type + ", " + std::to_string(z.count) +
                                             " faces) has no condition");
    if (z.kind == MeshZoneInfo::Kind::cell && z.type != "fluid")
      errors.add("mesh/" + z.name, "cell zone type " + z.type + " is unsupported; a Kutta case solves fluid zones only");
  }
  const auto z_zero = [&](const std::optional<Vec3>& v, const std::string& path) {
    if (planar && v && (*v)[2] != 0) errors.add(path, "a 2D mesh is solved in its x-y plane; the z component must be 0");
  };
  const auto check_turbulence = [&](const std::optional<TurbulenceInput>& t, const std::string& path) {
    if (!t) return;
    if (laminar) errors.add(path, "laminar flow has no turbulence inputs");
    if (t->method == TurbulenceInput::Method::k_omega && !sst) errors.add(path + "/method", "k-omega applies only to viscous k-omega-sst");
    if (t->method == TurbulenceInput::Method::k_epsilon && !ke) errors.add(path + "/method", "k-epsilon applies only to viscous k-epsilon-standard");
  };

  // ---- Rampant variables (section 37) ----
  Alist rp, config;
  config.flag("rp-seg?", true);
  config.flag("rp-visc?", true);
  config.flag("rp-unsteady?", false);
  config.flag("rp-des?", false);
  config.flag("rp-sa?", false);
  config.flag("rp-3d?", !planar);
  config.flag("rp-double?", true);
  config.flag("rf-energy?", false);
  config.flag("rp-lam?", laminar);
  config.flag("rp-ke?", ke);
  config.flag("rp-kw?", sst);
  config.flag("rp-turb?", !laminar);
  config.flag("rp-trb-scl?", !laminar);
  config.flag("rp-trans-one-eqn?", transition);
  config.flag("rp-trans-sst?", false);
  config.atom("rp-transition-model", transition ? "gamma-transport-eqn" : "none");
  rp.list("case-config", config.text());
  rp.real("operating-pressure", c.operating_pressure.value_or(defaults::operating_pressure));
  rp.flag("gravity?", false);
  {
    const auto material = c.material.name.value_or("air");
    if (!plain_symbol(material)) errors.add("material/name", "use letters, digits and -_. only");
    rp.list("materials", "((" + material + " fluid (density . " + real_text(c.material.density.value_or(defaults::air_density)) +
                             ") (viscosity . " + real_text(c.material.viscosity.value_or(defaults::air_viscosity)) + ")))");
  }
  // Model selection and options.
  if (sst) {
    rp.flag("kw-sst-on?", true);
    rp.flag("kw-bsl-on?", false);
    rp.flag("kw-std-on?", false);
    for (auto key : {"kw-geko-on?", "kw-easm-on?", "kw-wj-bsl-earsm-on?"}) rp.flag(key, false);
  }
  if (ke) {
    rp.flag("rng-ke-on?", false);
    rp.flag("realizable-epsilon?", false);
  }
  if (!laminar) {
    const auto& o = c.models.k_omega_sst;
    rp.flag(ke ? "ke-production-limiter?" : "kw-not-trans-sst-production-limiter?", ke ? false : o.production_limiter.value_or(defaults::production_limiter));
    rp.real(ke ? "ke-pk-factor" : "sst-pk-factor", o.production_limiter_clip_factor.value_or(defaults::production_limiter_clip));
    rp.flag(ke ? "ke-low-re?" : "kw-low-re-mod?", o.low_re_corrections.value_or(defaults::low_re));
    rp.flag("curvature-correction?", o.curvature_correction.value_or(defaults::curvature));
    rp.real("curvature-correction-ccurv", 1);
    rp.flag("kw-sst-strain?", false);  // [O]; the adapter requires it off
    rp.flag(ke ? "ke-kato-launder?" : "kw-not-trans-sst-kato-launder?", o.kato_launder.value_or(defaults::kato_launder));
    rp.flag("kw-vorticity-pk?", o.vorticity_based_production.value_or(defaults::vorticity_production));
    rp.flag("kw-shear-mod?", true);  // [O]
    rp.flag("turb-compress-mod?", ke ? c.models.k_epsilon.compressibility_effects.value_or(defaults::compressibility)
                                     : o.compressibility_effects.value_or(defaults::compressibility));
    rp.real("omega/min-limit", c.controls.min_specific_dissipation_rate.value_or(defaults::min_omega));
    rp.atom("turbulence/wall-omega-treatment", "correlation");  // TG 4.18.3 "current default"
    rp.real("k/min-limit", c.controls.min_turbulent_kinetic_energy.value_or(defaults::min_k));
    rp.real("turb-viscosity/max-limit", c.controls.max_turbulent_viscosity_ratio.value_or(defaults::max_viscosity_ratio));
  }
  if (sst) {
    rp.real("sst-limit-cross-diff-to-n-times-w-destruct", 1e6);  // [O]
    rp.flag("sst-linearize-neg-cross-diff?", true);              // [O]
  }
  if (ke) {
    // TG 4.3.1.2 constants; standard wall functions (UG 16.16); the two
    // undocumented solver switches Fluent writes on in every k-epsilon case
    // [O], so Kutta runs them only under --approximate ke-options, as for a
    // Fluent case.
    for (auto [key, value] : {std::pair{"kecmu", 0.09}, {"kec1", 1.44}, {"kec2", 1.92}, {"kesigk", 1.0}, {"kesige", 1.3}}) rp.real(key, value);
    rp.flag("ke-std-wall?", true);
    for (auto key : {"ke-enh-wall?", "ke-scalable-wall?", "ke-non-eqm", "ke-nwt-ml?", "ke-udf-wall?"}) rp.flag(key, false);
    rp.flag("ke-source-mod", true);
    rp.flag("ke-time-marching?", true);
    rp.real("epsilon/min-limit", c.controls.min_turbulent_dissipation_rate.value_or(defaults::min_epsilon));
  }
  if (laminar && (c.controls.min_turbulent_kinetic_energy || c.controls.min_specific_dissipation_rate ||
                  c.controls.min_turbulent_dissipation_rate || c.controls.max_turbulent_viscosity_ratio))
    errors.add("controls/limits", "turbulence limits do not apply to laminar flow");
  if (c.controls.min_specific_dissipation_rate && !sst) errors.add("controls/limits/min_specific_dissipation_rate", "applies only to k-omega-sst");
  if (c.controls.min_turbulent_dissipation_rate && !ke) errors.add("controls/limits/min_turbulent_dissipation_rate", "applies only to k-epsilon-standard");
  if (transition) {
    // TG 4.8 model constants as Fluent 22.2 writes them [O].
    for (auto [key, value] : {std::pair{"trans-one-eqn-ctu1", 100.0}, {"trans-one-eqn-ctu2", 1000.0}, {"trans-one-eqn-ctu3", 1.0},
                              {"trans-one-eqn-cpg1", 1.0}, {"trans-one-eqn-cpg2", -0.5}, {"trans-one-eqn-cpg3", 0.0},
                              {"trans-one-eqn-lim-cpg1", 1.5}, {"trans-one-eqn-lim-cpg2", 3.0}, {"trans-one-eqn-csep", 1.0},
                              {"trans-one-eqn-tke-term-coef", 1.0}, {"trans-one-eqn-crossflow-crsf", 1.0}, {"intermit-sig", 1.0},
                              {"trans-intermit/min-limit", 1e-7}, {"trans-intermit/max-limit", 1.0}})
      rp.real(key, value);
    rp.flag("trans-one-eqn-crossflow?", false);
    rp.flag("trans-intermit-smoothing?", false);
    rp.integer("trans-one-eqn-dudn-method", 1);
  }
  // Methods (UG 37.3, research table A3).
  const auto upwind = [](std::optional<Upwind> u, Upwind fallback) { return u.value_or(fallback) == Upwind::first_order ? 0 : 1; };
  rp.integer("flow/scheme", coupling == CouplingScheme::simple ? 20 : coupling == CouplingScheme::simplec ? 21 : 24);
  rp.integer("pressure/scheme", c.methods.pressure.value_or(PressureInterpolation::second_order) == PressureInterpolation::standard ? 10 : 12);
  rp.integer("mom/scheme", upwind(c.methods.momentum, Upwind::second_order));
  if (!laminar) rp.integer("k/scheme", upwind(c.methods.turbulent_kinetic_energy, Upwind::first_order));
  if (sst) rp.integer("omega/scheme", upwind(c.methods.specific_dissipation_rate, Upwind::first_order));
  if (ke) rp.integer("epsilon/scheme", upwind(c.methods.turbulent_dissipation_rate, Upwind::first_order));
  if (transition) rp.integer("intermit/scheme", upwind(c.methods.intermittency, Upwind::first_order));
  if (laminar && c.methods.turbulent_kinetic_energy) errors.add("methods/turbulent_kinetic_energy", "laminar flow solves no turbulence equation");
  if (!sst && c.methods.specific_dissipation_rate) errors.add("methods/specific_dissipation_rate", "applies only to k-omega-sst");
  if (!ke && c.methods.turbulent_dissipation_rate) errors.add("methods/turbulent_dissipation_rate", "applies only to k-epsilon-standard");
  if (!transition && c.methods.intermittency) errors.add("methods/intermittency", "applies only with the intermittency transition model");
  const auto gradient = c.methods.gradient.value_or(GradientMethod::least_squares_cell);
  rp.flag("recon/cell-lsq?", gradient == GradientMethod::least_squares_cell);
  rp.flag("recon/node-lsq?", false);
  const auto limiter = c.methods.gradient_limiter.value_or(GradientLimiter::standard);
  rp.flag("recon/limit?", limiter != GradientLimiter::none);
  rp.integer("recon/lim/type", limiter == GradientLimiter::differentiable ? 2 : 0);
  rp.flag("recon/lim/dx-to-face?", true);
  rp.flag("recon/lim/filter/on?", false);
  // UG: warped-face gradient correction is on by default for 3D meshes with polyhedral cells.
  rp.flag("recon/cell-lsf?", c.methods.warped_face_gradient_correction.value_or(!planar && mesh.polyhedral));
  rp.flag("recon/relax/relax?", c.methods.high_order_term_relaxation.value_or(false));
  rp.real("recon/relax/steady-urf", 0.25);
  const auto flux = c.methods.flux_type.value_or(FluxType::automatic);
  rp.flag("pbs/flux-auto-select?", flux == FluxType::automatic);
  rp.integer("pbs/flux-index", flux == FluxType::distance ? 0 : 1);
  rp.integer("wall-distance/method", 2);
  rp.flag("wall-distance/quad-correction?", true);
  // Controls.
  rp.flag("pseudo-transient-formulation?", pseudo);
  rp.integer("pseudo-time-method/coupled-pbns/dt-method", pseudo && coupling == CouplingScheme::coupled ? 1 : 0);
  rp.integer("pseudo-time-method/segregated-pbns/dt-method", pseudo && coupling != CouplingScheme::coupled ? 1 : 0);
  rp.integer("pseudo-transient-mode", 0);
  rp.integer("pseudo/autotime-lscale-option", 1);  // Conservative, UG 37.14.2.3
  rp.flag("pseudo-auto-time-step?", !c.controls.pseudo_time_step);
  rp.real("pseudo-time-step", c.controls.pseudo_time_step.value_or(defaults::pseudo_time_step));
  rp.real("pseudo-auto-time-step", 1);  // Fluent's last automatic step, display only [O]
  rp.real("pseudo-auto-time-step-scale-factor", c.controls.pseudo_time_scale_factor.value_or(defaults::time_scale_factor));
  rp.real("pressure-coupled/courant-number", c.controls.flow_courant_number.value_or(defaults::courant));
  if (!pseudo && (c.controls.pseudo_time_step || c.controls.pseudo_time_scale_factor || c.controls.pseudo_relaxation != RelaxationInput{}))
    errors.add("controls", "pseudo-time controls apply only with the pseudo-transient formulation");
  const auto& pr = c.controls.pseudo_relaxation;
  const auto& ur = c.controls.relaxation;
  const auto& er = c.controls.explicit_relaxation;
  rp.real("pressure-coupled/pressure/explicit-relax", er.pressure.value_or(defaults::explicit_relaxation));
  rp.real("pressure-coupled/mom/explicit-relax", er.momentum.value_or(defaults::explicit_relaxation));
  if ((er.pressure || er.momentum) && (pseudo || coupling != CouplingScheme::coupled))
    errors.add("controls/explicit_relaxation", "applies only to Coupled without the pseudo-transient formulation");
  if (ur != RelaxationInput{} && pseudo) errors.add("controls/relaxation", "classic under-relaxation applies only without the pseudo-transient formulation");
  rp.real("pressure-coupled/pressure/pseudo-explicit-relax", relax(pr, defaults::pseudo, &RelaxationInput::pressure));
  rp.real("pressure-coupled/mom/pseudo-explicit-relax", relax(pr, defaults::pseudo, &RelaxationInput::momentum));
  std::vector<std::pair<std::string_view, std::optional<double> RelaxationInput::*>> equations{
      {"pressure", &RelaxationInput::pressure}, {"mom", &RelaxationInput::momentum}, {"density", &RelaxationInput::density},
      {"body-force", &RelaxationInput::body_force}};
  if (sst) equations.insert(equations.end(), {{"k", &RelaxationInput::k}, {"omega", &RelaxationInput::omega}});
  if (ke) equations.insert(equations.end(), {{"k", &RelaxationInput::k}, {"epsilon", &RelaxationInput::epsilon}});
  if (transition) equations.push_back({"intermit", &RelaxationInput::intermittency});
  if (!laminar) equations.push_back({"turb-viscosity", &RelaxationInput::turbulent_viscosity});
  for (const auto& [eq, member] : equations) {
    rp.real(std::string(eq) + "/relax", relax(ur, defaults::classic, member));
    rp.real(std::string(eq) + "/pseudo-relax", relax(pr, defaults::pseudo, member));
  }
  for (const auto* set : {&pr, &ur}) {
    const auto path = set == &pr ? "controls/pseudo_relaxation" : "controls/relaxation";
    if (laminar && (set->k || set->omega || set->epsilon || set->turbulent_viscosity)) errors.add(path, "laminar flow relaxes no turbulence equation");
    if (!sst && set->omega) errors.add(std::string(path) + "/omega", "applies only to k-omega-sst");
    if (!ke && set->epsilon) errors.add(std::string(path) + "/epsilon", "applies only to k-epsilon-standard");
    if (!transition && set->intermittency) errors.add(std::string(path) + "/intermittency", "applies only with the transition model");
  }
  // Reference values (UG 42.13).
  const auto& ref = c.reference;
  rp.real("reference-area", ref.area.value_or(defaults::reference_area));
  rp.real("reference-length", ref.length.value_or(defaults::reference_length));
  rp.real("reference-density", ref.density.value_or(defaults::reference_density));
  rp.real("reference-velocity", ref.velocity.value_or(defaults::reference_velocity));
  rp.real("reference-pressure", ref.pressure.value_or(defaults::reference_pressure));
  rp.real("reference-viscosity", ref.viscosity.value_or(defaults::reference_viscosity));
  rp.real("reference-depth", ref.depth.value_or(defaults::reference_depth));
  if (planar && ref.depth && *ref.depth != 1)
    errors.add("reference_values/depth", "a 2D mesh is read as a one-cell extrusion 1 m deep; the reference depth must be 1");
  for (auto axis : {"x", "y", "z"}) rp.real(std::string("pressure-reference/") + axis, 0);
  // Initialization.
  using Init = KuttaCase::Initialization::Method;
  const auto& init = c.initialization;
  rp.flag("hyb-init?", init.method.value_or(Init::hybrid) == Init::hybrid);
  rp.flag("fmg-init/enable?", false);
  rp.real("hyb-init/niter", init.hybrid_iterations.value_or(defaults::hybrid_iterations));
  rp.list("hyb-init/explicit-relax", "(1.0 1.0)");
  rp.list("hyb-init/relax", "(1.0 1.0)");
  rp.list("hyb-init/residual-tolerance", "(1e-06 1e-06)");
  rp.list("hyb-init/solve?", "(#t #t)");
  rp.list("hyb-init/patch", "(0.0 0.0)");
  rp.flag("hyb-init/turb-domain-averaged?", true);
  rp.flag("hyb-init/patch-ads?", false);
  rp.flag("hyb-init/constant-vel-vector?", false);
  rp.real("hyb-init/turb-intens", 0.05);
  rp.real("hyb-init/turb-vr", 10);
  rp.flag("hyb-init/pressure-from-bc?", false);
  const auto v0 = init.velocity.value_or(Vec3{});
  z_zero(init.velocity, "initialization/velocity");
  rp.real("x-velocity/default", v0[0]);
  rp.real("y-velocity/default", v0[1]);
  rp.real("z-velocity/default", v0[2]);
  rp.real("pressure/default", init.gauge_pressure.value_or(0));
  rp.real("k/default", init.k.value_or(defaults::k));
  rp.real("omega/default", init.omega.value_or(defaults::omega));
  rp.real("epsilon/default", init.epsilon.value_or(defaults::epsilon));
  rp.real("intermit/default", init.intermittency.value_or(defaults::intermittency));
  if (laminar && (init.k || init.omega || init.epsilon)) errors.add("initialization", "laminar flow has no turbulence initial values");
  if (init.omega && !sst) errors.add("initialization/omega", "applies only to k-omega-sst");
  if (init.epsilon && !ke) errors.add("initialization/epsilon", "applies only to k-epsilon-standard");
  if (init.intermittency && !transition) errors.add("initialization/intermittency", "applies only with the transition model");
  // Residual monitors: UG 37.17 default criteria, every solved equation checked.
  {
    const auto names = solved_residuals(fluent_model(c.models.viscous), transition, planar);
    std::string rows;
    for (const auto& [name, r] : c.residuals)
      if (std::find(names.begin(), names.end(), name) == names.end()) {
        std::string list;
        for (auto n : names) list += (list.empty() ? "" : ", ") + std::string(n);
        errors.add("monitors/residuals/" + name, "not a residual of this case; residuals: " + list);
      }
    for (auto name : names) {
      const auto it = c.residuals.find(std::string(name));
      const auto row = it == c.residuals.end() ? ResidualInput{} : it->second;
      rows += "(" + std::string(name) + " #t 0. " + flag_text(row.check.value_or(true)) + ' ' +
              real_text(row.absolute_criterion.value_or(defaults::residual_criterion)) + " 0. 0.05) ";
    }
    rp.list("residuals/settings", "(" + rows + ")");
  }
  // Force reports (TG 25.2.1).
  {
    std::string reports;
    for (const auto& [name, r] : c.reports) {
      using K = ForceReportInput::Kind;
      for (const auto& z : r.zones) {
        const auto it = c.boundaries.find(z);
        if (it == c.boundaries.end() || !std::holds_alternative<WallInput>(it->second))
          errors.add("reports/" + name + "/zones", "\"" + z + "\" is not a wall zone of this case");
      }
      const Vec3 direction = r.direction.value_or(r.kind == K::lift ? Vec3{0, 1, 0} : Vec3{1, 0, 0});
      z_zero(r.direction, "reports/" + name + "/direction");
      std::string zones;
      for (const auto& z : r.zones) zones += ' ' + scheme_string(z);
      reports += "((name . " + scheme_string(name) + ") (report-definition " +
                 std::string(r.kind == K::drag ? "force-drag" : r.kind == K::lift ? "force-lift" : "force-force") +
                 " (zone-names" + zones + ") (force-vector " + real_text(direction[0]) + ' ' + real_text(direction[1]) + ' ' +
                 real_text(direction[2]) + ") (scaled? . " + flag_text(r.coefficient.value_or(r.kind != K::force)) + "))) ";
    }
    rp.list("monitor/report-definitions", "(" + reports + ")");
  }

  // ---- Thread variables (section 39) ----
  std::vector<fluent::SExpr> threads;
  const auto thread = [&](const MeshZoneInfo& z, std::string_view type, const Alist& a) {
    threads.push_back(fluent::parse_sexpr("(39 (" + std::to_string(z.id) + ' ' + std::string(type) + ' ' + symbol_text(z.name) + " 1)" + a.text() + ')'));
  };
  for (const auto& z : mesh.zones) {
    if (z.kind == MeshZoneInfo::Kind::cell) {
      Alist a;
      a.atom("material", c.material.name.value_or("air"));
      a.flag("mrf-motion?", false);
      a.real("mrf-omega", 0);
      a.integer("mrf-relative-to-thread", -1);
      a.integer("motion-spec", 0);
      a.vector({"mrf-grid-x-vel", "mrf-grid-y-vel", "mrf-grid-z-vel"}, {0, 0, 0});
      a.vector({"mrf-x-origin", "mrf-y-origin", "mrf-z-origin"}, {0, 0, 0});
      a.vector({"mrf-ai", "mrf-aj", "mrf-ak"}, {0, 0, 1});
      a.flag("porous?", false);
      a.flag("laminar?", false);
      a.real("porosity", 1);
      a.flag("rel-vel-resistance?", true);
      thread(z, "fluid", a);
      continue;
    }
    const auto it = c.boundaries.find(z.name);
    if (z.kind == MeshZoneInfo::Kind::interior || it == c.boundaries.end()) {
      if (z.kind == MeshZoneInfo::Kind::interior) thread(z, "interior", Alist{});
      continue;
    }
    const auto path = "boundaries/" + z.name;
    std::visit([&](const auto& b) {
      using T = std::decay_t<decltype(b)>;
      Alist a;
      if constexpr (std::is_same_v<T, WallInput>) {
        using M = WallInput::Motion;
        const auto motion = b.motion.value_or(M::stationary);
        a.integer("motion-bc", motion == M::stationary ? 0 : 1);
        a.flag("relative?", b.relative.value_or(true));
        if (motion == M::moving) {
          a.flag("rotating?", false);
          a.flag("components?", b.velocity.has_value());
          if (b.velocity) a.vector({"u", "v", "w"}, *b.velocity);
          a.real("vmag", b.speed.value_or(0));
          a.vector({"ni", "nj", "nk"}, b.direction.value_or(Vec3{1, 0, 0}));
          z_zero(b.velocity, path + "/velocity");
          z_zero(b.direction, path + "/direction");
        }
        if (motion == M::rotating) {
          a.flag("rotating?", true);
          a.real("omega", b.angular_speed.value_or(0));
          a.vector({"x-origin", "y-origin", "z-origin"}, b.axis_origin.value_or(Vec3{}));
          a.vector({"ai", "aj", "ak"}, b.axis_direction.value_or(Vec3{0, 0, 1}));
          if (planar && b.axis_direction && ((*b.axis_direction)[0] != 0 || (*b.axis_direction)[1] != 0))
            errors.add(path + "/axis_direction", "a 2D mesh rotates about the z axis only");
        }
        a.flag("specified-shear?", false);
        a.integer("shear-bc", 0);
        a.integer("rough-bc", 0);
        a.vector({"shear-x", "shear-y", "shear-z"}, {0, 0, 0});
        a.real("roughness-height", 0);
        a.real("roughness-const", 0.5);
        thread(z, "wall", a);
      } else if constexpr (std::is_same_v<T, VelocityInletInput>) {
        const int spec = b.velocity ? 1 : b.direction ? 0 : 2;  // components, magnitude and direction, normal
        a.integer("velocity-spec", spec);
        a.real("vmag", b.speed.value_or(0));
        if (b.direction) a.vector({"ni", "nj", "nk"}, *b.direction);
        if (b.velocity) a.vector({"u", "v", "w"}, *b.velocity);
        a.integer("frame-of-reference", 0);
        a.real("p", b.gauge_pressure.value_or(0));
        turbulence_entries(a, b.turbulence, solved);
        check_turbulence(b.turbulence, path + "/turbulence");
        if (transition) a.real("intermit", b.intermittency.value_or(defaults::intermittency));
        else if (b.intermittency) errors.add(path + "/intermittency", "applies only with the intermittency transition model");
        z_zero(b.velocity, path + "/velocity");
        z_zero(b.direction, path + "/direction");
        thread(z, "velocity-inlet", a);
      } else if constexpr (std::is_same_v<T, PressureOutletInput>) {
        a.real("p", b.gauge_pressure.value_or(0));
        a.integer("p-backflow-spec", b.backflow_pressure.value_or(PressureOutletInput::Backflow::total) == PressureOutletInput::Backflow::total ? 0 : 1);
        a.flag("prevent-reverse-flow?", b.prevent_reverse_flow.value_or(false));
        a.integer("direction-spec", 1);  // backflow normal to the boundary
        a.real("p-profile-multiplier", 1);
        turbulence_entries(a, b.turbulence, solved);
        check_turbulence(b.turbulence, path + "/turbulence");
        if (transition) a.real("intermit", b.intermittency.value_or(defaults::intermittency));
        else if (b.intermittency) errors.add(path + "/intermittency", "applies only with the intermittency transition model");
        thread(z, "pressure-outlet", a);
      } else if constexpr (std::is_same_v<T, SymmetryInput>) thread(z, "symmetry", a);
    }, it->second);
  }
  errors.finish();

  std::map<std::string, std::vector<fluent::SExpr>> sections;
  sections["Rampant Variables"].push_back(fluent::parse_sexpr("(37 " + rp.text() + ")"));
  sections["Thread Variables"] = std::move(threads);
  auto out = fluent::extract_case(std::move(sections), std::nullopt, "Kutta case v" + std::to_string(kutta_case_version));
  // No Fluent writer produced these settings; its version is not a gap.
  std::erase_if(out.issues, [](const auto& i) { return i.path == "solver-version"; });
  out.mesh.dimension = mesh.dimension;
  out.dimension = mesh.dimension;
  out.mesh.cells = mesh.cells;
  out.mesh.faces = mesh.faces;
  out.mesh.nodes = mesh.nodes;
  // Ranges are positions in Kutta's mesh order, not Fluent cell or face IDs:
  // a Kutta case has no Fluent data file to seed from.
  std::uint64_t cell = 1, face = 1;
  for (const auto& z : mesh.zones) {
    if (z.count == 0) continue;
    auto& target = z.kind == MeshZoneInfo::Kind::cell ? out.mesh.cell_zones : out.mesh.face_zones;
    auto& next = z.kind == MeshZoneInfo::Kind::cell ? cell : face;
    target.push_back({z.id, z.name, next, next + z.count - 1, {}});
    next += z.count;
  }
  return out;
}

void apply_boundary_kinds(const KuttaCase& c, PolyMesh& mesh) {
  for (auto& p : mesh.patches) {
    const auto it = c.boundaries.find(p.name);
    if (p.source_zone_id < 0 || it == c.boundaries.end()) continue;
    std::visit([&](const auto& b) {
      using T = std::decay_t<decltype(b)>;
      if constexpr (std::is_same_v<T, WallInput>) p.kind = PatchKind::Wall, p.source_type = "wall";
      else if constexpr (std::is_same_v<T, VelocityInletInput>) p.kind = PatchKind::VelocityInlet, p.source_type = "velocity-inlet";
      else if constexpr (std::is_same_v<T, PressureOutletInput>) p.kind = PatchKind::PressureOutlet, p.source_type = "pressure-outlet";
      else if constexpr (std::is_same_v<T, SymmetryInput>) p.kind = PatchKind::Symmetry, p.source_type = "symmetry";
    }, it->second);
  }
}

LoadedKuttaCase load_kutta_case(const std::filesystem::path& setup_file) {
  LoadedKuttaCase out;
  out.setup_path = std::filesystem::absolute(setup_file).lexically_normal();
  out.source = read_kutta_case_file(setup_file);
  out.mesh_path = out.source.mesh.is_absolute() ? out.source.mesh : (out.setup_path.parent_path() / out.source.mesh).lexically_normal();
  if (!std::filesystem::is_regular_file(out.mesh_path))
    throw RejectedSettings({{"mesh", "mesh file not found: " + out.mesh_path.string()}}, kutta_case_rejected);
  auto name = out.mesh_path.filename().string();
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char ch) { return char(std::tolower(ch)); });
  if (!name.ends_with(".msh") && !name.ends_with(".msh.h5") && !name.ends_with(".cas") && !name.ends_with(".cas.h5"))
    throw RejectedSettings({{"mesh", "a Kutta case reads Fluent meshes (.msh, .msh.h5, or the mesh of a .cas/.cas.h5): " + out.mesh_path.string()}},
                           kutta_case_rejected);
  out.mesh = read_mesh(out.mesh_path);
  out.settings = lower_kutta_case(out.source, mesh_info(out.mesh));
  out.settings.source_path = out.setup_path;
  apply_boundary_kinds(out.source, out.mesh);
  return out;
}

std::string kutta_case_identity(std::string_view setup_sha256, std::string_view mesh_sha256) {
  return sha256(std::string(setup_sha256) + "\n" + std::string(mesh_sha256));
}


// ---- Fluent case -> Kutta case -------------------------------------------------------------
ConvertedCase kutta_case_from_fluent(const fluent::FluentCase& c, const MeshInfo& mesh, std::filesystem::path mesh_path) {
  using namespace fluent;
  ConvertedCase out;
  auto& k = out.kutta;
  k.mesh = std::move(mesh_path);
  // A change to the solved problem blocks the conversion; anything else is
  // reported as informational.
  const auto omit = [&](std::string path, std::string reason, bool changes_solution = true) {
    out.omitted.push_back({std::move(path), std::move(reason), changes_solution ? IssueClass::unsupported : IssueClass::informational});
  };
  const auto value = [&](const std::optional<ScalarValue>& v, const std::string& path) -> std::optional<double> {
    if (!v) return {};
    if (v->constant) return v->constant;
    if (v->resolved) {
      omit(path, "expression " + (v->expression_name.empty() ? v->definition : v->expression_name) + " written as its value " + number_text(*v->resolved), false);
      return v->resolved;
    }
    omit(path, "profile or unresolved expression " + v->binding + " cannot be written");
    return {};
  };
  const auto set = [](std::optional<double>& target, std::optional<double> v, double fallback) {
    if (v && *v != fallback) target = v;
  };
  const auto set_flag = [](std::optional<bool>& target, std::optional<bool> v, bool fallback) {
    if (v && *v != fallback) target = v;
  };
  const auto raw_number = [&](const RawSettings& raw, std::string_view key) -> std::optional<double> {
    const auto* v = raw.value(key);
    return v ? v->number() : std::nullopt;
  };
  const auto raw_flag = [&](const RawSettings& raw, std::string_view key) -> std::optional<bool> {
    const auto* v = raw.value(key);
    return v && v->kind == SExpr::Kind::boolean ? std::optional<bool>(v->boolean) : std::nullopt;
  };
  const auto expect = [&](std::optional<double> v, double expected, std::string path) {
    if (v && *v != expected) omit(std::move(path), "is " + number_text(*v) + "; a Kutta case uses " + number_text(expected));
  };
  const auto expect_flag = [&](std::optional<bool> v, bool expected, std::string path) {
    if (v && *v != expected) omit(std::move(path), std::string("is ") + flag_text(*v) + "; a Kutta case uses " + flag_text(expected));
  };
  const auto vec = [&](const ScalarVector& v, const std::string& path) -> std::optional<Vec3> {
    Vec3 out3{};
    for (int i = 0; i < 3; ++i) {
      const auto x = value(v[std::size_t(i)], path + "/" + "xyz"[i]);
      if (!x) return {};
      out3[std::size_t(i)] = *x;
    }
    return out3;
  };

  const auto& t = c.viscous;
  if (t.model != fluent::ViscousModel::laminar && t.model != fluent::ViscousModel::kw_sst && t.model != fluent::ViscousModel::ke_standard)
    throw RejectedSettings({{"viscous-model", "a Kutta case expresses laminar, k-omega SST and standard k-epsilon; this case uses " +
                                                  std::string(t.model ? name(*t.model) : "an unparsed model")}},
                           "Not expressible in a Kutta case (nothing written):");
  const bool laminar = t.model == fluent::ViscousModel::laminar, sst = t.model == fluent::ViscousModel::kw_sst, ke = t.model == fluent::ViscousModel::ke_standard;
  k.models.viscous = laminar ? setup::ViscousModel::laminar : sst ? setup::ViscousModel::k_omega_sst : setup::ViscousModel::k_epsilon_standard;
  const bool transition = t.transition == TransitionModel::intermittency;
  if (t.transition && *t.transition != TransitionModel::none && !transition)
    omit("case-config/rp-transition-model", "only the intermittency transition model is expressible");
  if (transition) k.models.intermittency_transition = true;
  if (transition) {
    const std::map<std::string, double> constants{{"trans-one-eqn-ctu1", 100}, {"trans-one-eqn-ctu2", 1000}, {"trans-one-eqn-ctu3", 1},
        {"trans-one-eqn-cpg1", 1}, {"trans-one-eqn-cpg2", -0.5}, {"trans-one-eqn-cpg3", 0}, {"trans-one-eqn-lim-cpg1", 1.5},
        {"trans-one-eqn-lim-cpg2", 3}, {"trans-one-eqn-csep", 1}, {"trans-one-eqn-tke-term-coef", 1}, {"trans-one-eqn-crossflow-crsf", 1},
        {"intermit-sig", 1}, {"trans-intermit/min-limit", 1e-7}, {"trans-intermit/max-limit", 1}, {"intermit/default", 1}};
    for (const auto& [key, v] : t.transition_coefficients)
      if (const auto it = constants.find(key); it != constants.end() && key != "intermit/default") expect(v, it->second, "rp/" + key);
    expect_flag(t.transition_crossflow, false, "rp/trans-one-eqn-crossflow?");
    expect_flag(t.transition_smoothing, false, "rp/trans-intermit-smoothing?");
    if (t.transition_normal_derivative && *t.transition_normal_derivative != 1) omit("rp/trans-one-eqn-dudn-method", "only method 1 is expressible");
  }
  if (sst) {
    auto& o = k.models.k_omega_sst;
    set_flag(o.production_limiter, t.production_limiter, defaults::production_limiter);
    set(o.production_limiter_clip_factor, t.production_factor, defaults::production_limiter_clip);
    set_flag(o.low_re_corrections, t.low_re, defaults::low_re);
    set_flag(o.curvature_correction, t.curvature_correction, defaults::curvature);
    set_flag(o.compressibility_effects, t.compressibility, defaults::compressibility);
    set_flag(o.kato_launder, t.kato_launder, defaults::kato_launder);
    set_flag(o.vorticity_based_production, t.vorticity_production, defaults::vorticity_production);
    if (t.curvature_correction == true) expect(t.curvature_factor, 1, "rp/curvature-correction-ccurv");
    expect(raw_number(c.variables, "sst-limit-cross-diff-to-n-times-w-destruct"), 1e6, "rp/sst-limit-cross-diff-to-n-times-w-destruct");
    expect_flag(raw_flag(c.variables, "sst-linearize-neg-cross-diff?"), true, "rp/sst-linearize-neg-cross-diff?");
    if (c.numerics.wall_distance_method) expect(double(*c.numerics.wall_distance_method), 2, "rp/wall-distance/method");
    expect_flag(c.numerics.wall_distance_quad_correction, true, "rp/wall-distance/quad-correction?");
    if (t.wall_omega == WallOmega::tabulated) omit("rp/turbulence/wall-omega-treatment", "only the correlation is expressible");
    const std::map<std::string, double> coefficients{{"sst-a1", .31}, {"sst-beta-i1", .075}, {"sst-beta-i2", .0828}, {"sst-sig-k1", 1.176},
        {"sst-sig-k2", 1}, {"sst-sig-w1", 2}, {"sst-sig-w2", 1.168}, {"sst-kappa", .41}, {"kw-beta-star-inf", .09}, {"kw-alpha-star-inf", 1},
        {"kw-beta-i", .072}, {"kw-alpha-0", .11111}, {"kw-r-w", 2.95}, {"kw-r-k", 6}, {"kw-r-beta", 8}, {"kw-mt0", .25}, {"kw-zeta-star", 1.5},
        {"kw-alpha-inf", .52}, {"kw-sig-k", 2}, {"kw-sig-w", 2}};
    for (const auto& [key, v] : t.coefficients)
      if (const auto it = coefficients.find(key); it == coefficients.end() || std::abs(v - it->second) > 1e-10) omit("rp/" + key, "model coefficients are Fluent's defaults in a Kutta case");
  }
  if (ke) {
    set_flag(k.models.k_epsilon.compressibility_effects, t.compressibility, defaults::compressibility);
    expect_flag(t.production_limiter, false, "rp/ke-production-limiter?");
    expect_flag(t.low_re, false, "rp/ke-low-re?");
    expect_flag(t.kato_launder, false, "rp/ke-kato-launder?");
    expect_flag(t.curvature_correction, false, "rp/curvature-correction?");
    for (auto [key, fallback] : {std::pair{"ke-source-mod", true}, {"ke-time-marching?", true}, {"ke-std-wall?", true}, {"ke-enh-wall?", false},
                                 {"ke-scalable-wall?", false}, {"ke-non-eqm", false}, {"ke-nwt-ml?", false}, {"ke-udf-wall?", false}})
      expect_flag(raw_flag(c.variables, key), fallback, std::string("rp/") + key);
    for (auto [key, fallback] : {std::pair{"kecmu", 0.09}, {"kec1", 1.44}, {"kec2", 1.92}, {"kesigk", 1.0}, {"kesige", 1.3}})
      if (const auto v = raw_number(c.variables, key); v && std::abs(*v - fallback) > 2e-7 * fallback) omit(std::string("rp/") + key, "model constants are Fluent's defaults in a Kutta case");
  }
  if (!laminar) {
    expect_flag(t.sst_strain, false, "rp/kw-sst-strain?");
    expect_flag(t.shear_modification, true, "rp/kw-shear-mod?");
    set(k.controls.min_turbulent_kinetic_energy, raw_number(c.variables, "k/min-limit"), defaults::min_k);
    set(k.controls.max_turbulent_viscosity_ratio, raw_number(c.variables, "turb-viscosity/max-limit"), defaults::max_viscosity_ratio);
    if (sst) set(k.controls.min_specific_dissipation_rate, t.omega_floor, defaults::min_omega);
    if (ke) set(k.controls.min_turbulent_dissipation_rate, raw_number(c.variables, "epsilon/min-limit"), defaults::min_epsilon);
  }
  if (c.gravity == true) omit("rp/gravity?", "gravity is not expressible");
  if (c.energy == true) omit("case-config/rf-energy?", "the energy equation is not expressible; with constant density and viscosity it does not change the flow", false);
  if (c.dimension == 2 && c.config.value("rp-axi?") && raw_flag(c.config, "rp-axi?") == true) omit("case-config/rp-axi?", "axisymmetric cases are not expressible");
  set(k.operating_pressure, c.operating_pressure, defaults::operating_pressure);

  // Material and cell zones.
  std::optional<std::string> material;
  for (const auto& z : c.cell_zones) {
    const auto prefix = "zone/" + z.name + "/";
    if (z.type != ZoneType::fluid) omit(prefix + "type", "only fluid cell zones are expressible");
    if (z.frame_motion == true) omit(prefix + "mrf-motion?", "moving reference frames are not expressible in a version 1 Kutta case");
    if (z.porous == true) omit(prefix + "porous?", "porous zones are not expressible in a version 1 Kutta case");
    if (z.laminar == true) omit(prefix + "laminar?", "laminar zones are not expressible");
    if (z.type == ZoneType::fluid && z.material) {
      if (material && *material != *z.material) omit(prefix + "material", "a Kutta case has one fluid material");
      material = z.material;
    }
  }
  if (material) {
    const auto m = std::find_if(c.materials.begin(), c.materials.end(), [&](const auto& x) { return x.name == *material; });
    if (m != c.materials.end()) {
      if (m->name != "air") k.material.name = m->name;
      set(k.material.density, value(m->density, "material/" + m->name + "/density"), defaults::air_density);
      set(k.material.viscosity, value(m->viscosity, "material/" + m->name + "/viscosity"), defaults::air_viscosity);
    }
  }

  // Boundary conditions.
  using M = TurbulenceInput::Method;
  const auto turbulence = [&](const BoundaryZone& z, const std::string& prefix) -> std::optional<TurbulenceInput> {
    if (laminar || !z.turbulence_spec) return {};
    TurbulenceInput in;
    const auto spec = *z.turbulence_spec;
    const auto method = spec == TurbulenceSpec::k_dissipation ? (sst ? M::k_omega : M::k_epsilon)
                        : spec == TurbulenceSpec::intensity_length ? M::intensity_length_scale
                        : spec == TurbulenceSpec::intensity_viscosity_ratio ? M::intensity_viscosity_ratio : M::intensity_hydraulic_diameter;
    if (method != M::intensity_viscosity_ratio) in.method = method;
    const auto get = [&](const std::string& key) -> std::optional<double> {
      const auto it = z.turbulence.find(key);
      return it == z.turbulence.end() ? std::nullopt : value(it->second, prefix + key);
    };
    if (method == M::k_omega || method == M::k_epsilon) {
      set(in.k, get("k"), defaults::k);
      if (sst) set(in.omega, get("omega"), defaults::omega);
      else set(in.epsilon, get("epsilon"), defaults::epsilon);
    } else {
      set(in.intensity, get("turb-intensity"), defaults::intensity);
      if (method == M::intensity_viscosity_ratio) set(in.viscosity_ratio, get("turb-viscosity-ratio"), defaults::viscosity_ratio);
      if (method == M::intensity_length_scale) set(in.length_scale, get("turb-length-scale"), defaults::length_scale);
      if (method == M::intensity_hydraulic_diameter) set(in.hydraulic_diameter, get("turb-hydraulic-diam"), defaults::hydraulic_diameter);
    }
    if (in == TurbulenceInput{}) return {};
    return in;
  };
  for (const auto& z : c.boundaries) {
    const auto prefix = "zone/" + z.name + "/";
    switch (z.type) {
      case ZoneType::interior: break;
      case ZoneType::symmetry: k.boundaries.emplace(z.name, SymmetryInput{}); break;
      case ZoneType::wall: {
        WallInput w;
        using Motion = WallInput::Motion;
        if (z.wall_motion == WallMotion::translating) {
          w.motion = Motion::moving;
          if (z.velocity_spec == VelocitySpec::components) w.velocity = vec(z.velocity, prefix + "velocity");
          else {
            w.speed = value(z.velocity_magnitude, prefix + "vmag");
            if (const auto d = vec(z.direction, prefix + "direction"); d && *d != Vec3{1, 0, 0}) w.direction = d;
          }
        } else if (z.wall_motion == WallMotion::rotating) {
          w.motion = Motion::rotating;
          w.angular_speed = value(z.angular_velocity, prefix + "omega");
          if (const auto o = vec(z.axis_origin, prefix + "axis-origin"); o && *o != Vec3{}) w.axis_origin = o;
          if (const auto a = vec(z.axis_direction, prefix + "axis-direction"); a && *a != Vec3{0, 0, 1}) w.axis_direction = a;
        } else if (z.wall_motion != WallMotion::stationary)
          omit(prefix + "motion-bc", "unparsed wall motion");
        if (w.motion) set_flag(w.relative, z.relative_motion, true);
        if (z.specified_shear == true || (z.shear_bc && *z.shear_bc != 0) || (z.rough_bc && *z.rough_bc != 0))
          omit(prefix + "shear-bc", "specified shear and roughness laws are not expressible");
        if (const auto h = value(z.roughness_height, prefix + "roughness-height"); h && *h != 0) omit(prefix + "roughness-height", "wall roughness is not expressible");
        k.boundaries.emplace(z.name, w);
        break;
      }
      case ZoneType::velocity_inlet: {
        VelocityInletInput v;
        if (z.velocity_spec == VelocitySpec::components) v.velocity = vec(z.velocity, prefix + "velocity");
        else {
          v.speed = value(z.velocity_magnitude, prefix + "vmag");
          if (v.speed == 0.0) v.speed.reset();
          if (z.velocity_spec == VelocitySpec::magnitude_direction) {
            v.direction = vec(z.direction, prefix + "direction");
            if (!v.speed) v.speed = 0.0;  // a direction needs its speed
          }
        }
        if (z.frame_of_reference && *z.frame_of_reference != 0) omit(prefix + "frame-of-reference", "only absolute inlet velocities are expressible");
        set(v.gauge_pressure, value(z.gauge_pressure, prefix + "p"), 0);
        v.turbulence = turbulence(z, prefix);
        if (transition) set(v.intermittency, value(z.intermittency, prefix + "intermit"), defaults::intermittency);
        k.boundaries.emplace(z.name, v);
        break;
      }
      case ZoneType::pressure_outlet: {
        PressureOutletInput p;
        using B = PressureOutletInput::Backflow;
        set(p.gauge_pressure, value(z.gauge_pressure, prefix + "p"), 0);
        if (z.backflow_pressure_spec == 1) p.backflow_pressure = B::static_pressure;
        else if (z.backflow_pressure_spec && *z.backflow_pressure_spec != 0) omit(prefix + "p-backflow-spec", "only total (0) and static (1) are expressible");
        if (!z.raw.entry("prevent-reverse-flow?"))
          omit(prefix + "prevent-reverse-flow?", "absent (Fluent " + c.solver_version.value_or("unknown") +
                   " writer); the Kutta case states Fluent's default, off, which is how Kutta's reverse-flow approximation runs it", false);
        set_flag(p.prevent_reverse_flow, z.prevent_reverse_flow, false);
        expect(raw_number(z.raw, "direction-spec"), 1, prefix + "direction-spec");
        if (const auto m = raw_number(z.raw, "p-profile-multiplier")) expect(m, 1, prefix + "p-profile-multiplier");
        p.turbulence = turbulence(z, prefix);
        if (transition) set(p.intermittency, value(z.intermittency, prefix + "intermit"), defaults::intermittency);
        k.boundaries.emplace(z.name, p);
        break;
      }
      default: omit(prefix + "type", "boundary type " + z.fluent_type + " is not expressible in a version 1 Kutta case");
    }
  }

  // Methods.
  const auto& n = c.numerics;
  if (const auto scheme = n.run_coupling ? n.run_coupling : n.coupling) {
    if (*scheme == Coupling::piso) omit("rp/flow/scheme", "PISO is not expressible");
    else if (*scheme != Coupling::coupled) k.methods.coupling = *scheme == Coupling::simple ? CouplingScheme::simple : CouplingScheme::simplec;
  }
  const auto coupling = k.methods.coupling.value_or(CouplingScheme::coupled);
  const bool pseudo = uses_pseudo_time(c);
  set_flag(k.methods.pseudo_transient, pseudo, coupling == CouplingScheme::coupled);
  if (n.gradient == Gradient::green_gauss_cell) k.methods.gradient = GradientMethod::green_gauss_cell;
  else if (n.gradient == Gradient::green_gauss_node) omit("rp/recon/node-lsq?", "node-based gradients are not expressible");
  if (n.pressure == Scheme::standard_pressure) k.methods.pressure = PressureInterpolation::standard;
  else if (n.pressure && *n.pressure != Scheme::second_pressure) omit("rp/pressure/scheme", std::string(name(*n.pressure)) + " is not expressible");
  const auto scheme = [&](std::optional<Scheme> s, Upwind fallback, std::optional<Upwind>& target, const std::string& path) {
    if (!s) return;
    if (*s != Scheme::first_upwind && *s != Scheme::second_upwind) {
      omit(path, std::string(name(*s)) + " is not expressible");
      return;
    }
    const auto u = *s == Scheme::first_upwind ? Upwind::first_order : Upwind::second_order;
    if (u != fallback) target = u;
  };
  scheme(n.momentum, Upwind::second_order, k.methods.momentum, "rp/mom/scheme");
  if (!laminar) scheme(n.k, Upwind::first_order, k.methods.turbulent_kinetic_energy, "rp/k/scheme");
  if (sst) scheme(n.omega, Upwind::first_order, k.methods.specific_dissipation_rate, "rp/omega/scheme");
  if (ke) scheme(n.epsilon, Upwind::first_order, k.methods.turbulent_dissipation_rate, "rp/epsilon/scheme");
  if (transition) scheme(n.intermittency, Upwind::first_order, k.methods.intermittency, "rp/intermit/scheme");
  if (!n.flux_auto_select && !n.flux_index)
    omit("rp/pbs/flux-index", "absent (case predates the flux-type selector); written as Fluent's default, Auto Select (momentum-based Rhie-Chow), "
                               "the momentum-weighted flux Kutta runs for such cases", false);
  else if (n.flux_index == 0) k.methods.flux_type = FluxType::distance;
  else if (n.flux_auto_select == false) k.methods.flux_type = FluxType::momentum;
  if (n.limiting == false) k.methods.gradient_limiter = GradientLimiter::none;
  else if (n.limiter == fluent::Limiter::differentiable) k.methods.gradient_limiter = GradientLimiter::differentiable;
  expect_flag(n.limit_to_face, true, "rp/recon/lim/dx-to-face?");
  expect_flag(n.limiter_filter, false, "rp/recon/lim/filter/on?");
  set_flag(k.methods.warped_face_gradient_correction, n.warped_face_correction, mesh.dimension == 3 && mesh.polyhedral);
  set_flag(k.methods.high_order_term_relaxation, n.high_order_relaxation, false);

  // Controls.
  if (pseudo) {
    set(k.controls.pseudo_time_scale_factor, n.pseudo.scale_factor, defaults::time_scale_factor);
    if (n.pseudo.automatic == false) k.controls.pseudo_time_step = n.pseudo.time_step;
    if (coupling == CouplingScheme::coupled) set(k.controls.flow_courant_number, n.courant, defaults::courant);
    if (n.pseudo.mode && *n.pseudo.mode != 0) omit("rp/pseudo-transient-mode", "only mode 0 is expressible");
    if (n.pseudo.length_scale_method && *n.pseudo.length_scale_method != 1) omit("rp/pseudo/autotime-lscale-option", "only the conservative length scale is expressible");
  }
  const auto relaxation = [&](RelaxationInput& target, const RelaxationInput& fallback, auto lookup) {
    for (auto [eq, member] : std::initializer_list<std::pair<std::string_view, std::optional<double> RelaxationInput::*>>{
             {"pressure", &RelaxationInput::pressure}, {"mom", &RelaxationInput::momentum}, {"density", &RelaxationInput::density},
             {"body-force", &RelaxationInput::body_force}, {"k", &RelaxationInput::k}, {"omega", &RelaxationInput::omega},
             {"epsilon", &RelaxationInput::epsilon}, {"turb-viscosity", &RelaxationInput::turbulent_viscosity},
             {"intermit", &RelaxationInput::intermittency}}) {
      const bool solved = eq == "pressure" || eq == "mom" || eq == "density" || eq == "body-force" ||
                          (!laminar && (eq == "k" || eq == "turb-viscosity")) || (sst && eq == "omega") || (ke && eq == "epsilon") ||
                          (transition && eq == "intermit");
      if (solved) set(target.*member, lookup(std::string(eq)), *(fallback.*member));
    }
  };
  if (pseudo)
    relaxation(k.controls.pseudo_relaxation, defaults::pseudo, [&](const std::string& eq) -> std::optional<double> {
      const auto key = coupling == CouplingScheme::coupled && (eq == "pressure" || eq == "mom") ? "pressure-coupled/" + eq + "/pseudo-explicit-relax" : eq + "/pseudo-relax";
      const auto it = n.pseudo_relaxation.find(key);
      return it == n.pseudo_relaxation.end() ? std::nullopt : std::optional<double>(it->second);
    });
  else {
    relaxation(k.controls.relaxation, defaults::classic, [&](const std::string& eq) -> std::optional<double> {
      if (coupling == CouplingScheme::coupled && (eq == "pressure" || eq == "mom")) return {};
      const auto it = n.relaxation.find(eq);
      return it == n.relaxation.end() ? std::nullopt : std::optional<double>(it->second);
    });
    if (coupling == CouplingScheme::coupled) {
      set(k.controls.explicit_relaxation.pressure, raw_number(c.variables, "pressure-coupled/pressure/explicit-relax"), defaults::explicit_relaxation);
      set(k.controls.explicit_relaxation.momentum, raw_number(c.variables, "pressure-coupled/mom/explicit-relax"), defaults::explicit_relaxation);
    }
  }

  // Reference values.
  set(k.reference.area, c.reference.area, defaults::reference_area);
  set(k.reference.length, c.reference.length, defaults::reference_length);
  set(k.reference.density, c.reference.density, defaults::reference_density);
  set(k.reference.velocity, c.reference.velocity, defaults::reference_velocity);
  set(k.reference.pressure, c.reference.pressure, defaults::reference_pressure);
  set(k.reference.viscosity, c.reference.viscosity, defaults::reference_viscosity);
  set(k.reference.depth, raw_number(c.variables, "reference-depth"), defaults::reference_depth);
  for (auto axis : {"x", "y", "z"}) expect(raw_number(c.variables, std::string("pressure-reference/") + axis), 0, std::string("rp/pressure-reference/") + axis);

  // Initialization.
  using Init = KuttaCase::Initialization::Method;
  const auto hybrid = raw_flag(c.variables, "hyb-init?");
  if (hybrid == false) {
    k.initialization.method = Init::standard;
    if (const auto x = raw_number(c.variables, "x-velocity/default"), y = raw_number(c.variables, "y-velocity/default"),
        zz = raw_number(c.variables, "z-velocity/default");
        x && y && zz && Vec3{*x, *y, *zz} != Vec3{})
      k.initialization.velocity = Vec3{*x, *y, *zz};
    set(k.initialization.gauge_pressure, raw_number(c.variables, "pressure/default"), 0);
    if (!laminar) set(k.initialization.k, raw_number(c.variables, "k/default"), defaults::k);
    if (sst) set(k.initialization.omega, raw_number(c.variables, "omega/default"), defaults::omega);
    if (ke) set(k.initialization.epsilon, raw_number(c.variables, "epsilon/default"), defaults::epsilon);
    if (transition) set(k.initialization.intermittency, raw_number(c.variables, "intermit/default"), defaults::intermittency);
  } else if (hybrid == true) {
    if (const auto n_iter = raw_number(c.variables, "hyb-init/niter"); n_iter && *n_iter != defaults::hybrid_iterations)
      k.initialization.hybrid_iterations = int(*n_iter);
    const auto pair = [&](std::string_view key, std::string_view expected) {
      const auto* v = c.variables.value(key);
      if (v && to_string(*v) != expected) {
        // Compare numerically: writers differ in "1" vs "1.0".
        const auto e = parse_sexpr(expected);
        bool same = v->elements.size() == e.elements.size();
        for (std::size_t i = 0; same && i < e.elements.size(); ++i)
          same = v->elements[i].number() ? v->elements[i].number() == e.elements[i].number() : v->elements[i] == e.elements[i];
        if (!same) omit("rp/" + std::string(key), "hybrid initialization controls other than the iteration count are Fluent's defaults in a Kutta case");
      }
    };
    pair("hyb-init/explicit-relax", "(1.0 1.0)");
    pair("hyb-init/relax", "(1.0 1.0)");
    pair("hyb-init/residual-tolerance", "(1e-06 1e-06)");
    pair("hyb-init/solve?", "(#t #t)");
    pair("hyb-init/patch", "(0.0 0.0)");
    expect_flag(raw_flag(c.variables, "hyb-init/turb-domain-averaged?"), true, "rp/hyb-init/turb-domain-averaged?");
    expect_flag(raw_flag(c.variables, "hyb-init/patch-ads?"), false, "rp/hyb-init/patch-ads?");
    expect_flag(raw_flag(c.variables, "hyb-init/constant-vel-vector?"), false, "rp/hyb-init/constant-vel-vector?");
    expect_flag(raw_flag(c.variables, "hyb-init/pressure-from-bc?"), false, "rp/hyb-init/pressure-from-bc?");
    expect(raw_number(c.variables, "hyb-init/turb-intens"), 0.05, "rp/hyb-init/turb-intens");
    expect(raw_number(c.variables, "hyb-init/turb-vr"), 10, "rp/hyb-init/turb-vr");
  }
  expect_flag(raw_flag(c.variables, "fmg-init/enable?"), false, "rp/fmg-init/enable?");

  // Reports and residual monitors.
  for (const auto& r : c.reports) {
    using RK = ReportDefinition::Kind;
    if (r.kind == RK::moment || r.kind == RK::other) {
      omit("report/" + r.name, "report type " + r.fluent_type + " is not evaluated by Kutta and is not expressible", false);
      continue;
    }
    if (r.kind == RK::downforce || !r.zone_names || !r.direction) {
      omit("report/" + r.name, "only force, drag and lift reports with named zones and a direction are expressible");
      continue;
    }
    ForceReportInput f;
    f.kind = r.kind == RK::drag ? ForceReportInput::Kind::drag : r.kind == RK::lift ? ForceReportInput::Kind::lift : ForceReportInput::Kind::force;
    f.zones = *r.zone_names;
    const Vec3 d{(*r.direction)[0], (*r.direction)[1], (*r.direction)[2]};
    if (d != (f.kind == ForceReportInput::Kind::lift ? Vec3{0, 1, 0} : Vec3{1, 0, 0})) f.direction = d;
    if (r.scaled && *r.scaled != (f.kind != ForceReportInput::Kind::force)) f.coefficient = r.scaled;
    k.reports.emplace(r.name, std::move(f));
  }
  if (const auto* rows = c.variables.value("residuals/settings"); rows && rows->kind == SExpr::Kind::list) {
    const auto names = solved_residuals(t.model, transition, c.dimension == 2);
    for (const auto& row : rows->elements) {
      if (row.elements.size() != 7 || !row.elements[0].atom_text()) continue;
      const std::string eq(*row.elements[0].atom_text());
      if (std::find(names.begin(), names.end(), eq) == names.end()) continue;
      ResidualInput in;
      if (row.elements[3].kind == SExpr::Kind::boolean && !row.elements[3].boolean) in.check = false;
      if (const auto crit = row.elements[4].number(); crit && *crit != defaults::residual_criterion) in.absolute_criterion = crit;
      if (in != ResidualInput{}) k.residuals.emplace(eq, in);
    }
  }
  return out;
}

// ---- starter cases -------------------------------------------------------------------------
std::string_view role_name(ZoneRole r) {
  constexpr std::string_view names[] = {"velocity-inlet", "pressure-outlet", "moving-ground", "wall", "symmetry", "interior"};
  return names[int(r)];
}
std::optional<ZoneRole> parse_role(std::string_view s) {
  for (int i = 0; i <= int(ZoneRole::interior); ++i)
    if (role_name(ZoneRole(i)) == s) return ZoneRole(i);
  return {};
}

namespace {
std::vector<std::string> tokens(std::string_view name) {
  std::vector<std::string> out;
  std::string cur;
  for (unsigned char ch : name) {
    if (std::isalnum(ch)) cur += char(std::tolower(ch));
    else if (!cur.empty()) out.push_back(std::exchange(cur, {}));
  }
  if (!cur.empty()) out.push_back(cur);
  return out;
}
// The app's aero groups (app/src/core/Protocol.cpp suggestAeroGroups), first match wins.
std::string aero_group(std::string_view name) {
  static const std::vector<std::pair<std::string, std::regex>> rules{
      {"Front wing", std::regex(R"((^|[^a-z])(fw|front[-_ ]?wing|nose[-_ ]?wing))", std::regex::icase)},
      {"Rear wing", std::regex(R"((^|[^a-z])(rw|rear[-_ ]?wing|drs))", std::regex::icase)},
      {"Wing", std::regex(R"(wing|flap|element|main[-_ ]?plane|foil|slat|gurney|endplate)", std::regex::icase)},
      {"Undertray / diffuser", std::regex(R"(undertray|under[-_ ]?tray|diffuser|floor|underbody|splitter)", std::regex::icase)},
      {"Wheels", std::regex(R"(wheel|tyre|tire|rim)", std::regex::icase)}};
  const std::string s(name);
  for (const auto& [group, pattern] : rules)
    if (std::regex_search(s, pattern)) return group;
  return "Body";
}
} // namespace

RoleSuggestion suggest_role(const MeshZoneInfo& z) {
  if (z.kind == MeshZoneInfo::Kind::interior) return {ZoneRole::interior, "interior face zone of the mesh", {}};
  const auto words = tokens(z.name);
  const auto has = [&](std::initializer_list<std::string_view> list) -> std::optional<std::string> {
    for (const auto& w : words)
      for (auto l : list)
        if (w == l) return w;
    return {};
  };
  if (const auto w = has({"inlet", "inflow", "upstream", "intake"})) return {ZoneRole::velocity_inlet, "name contains \"" + *w + "\"", {}};
  if (const auto w = has({"outlet", "outflow", "downstream", "exit"})) return {ZoneRole::pressure_outlet, "name contains \"" + *w + "\"", {}};
  if (const auto w = has({"ground", "road", "belt"})) return {ZoneRole::moving_ground, "name contains \"" + *w + "\"; the road moves with the freestream", {}};
  if (const auto w = has({"symmetry", "sym", "symm", "mirror", "midplane", "centreplane", "centerplane", "top", "ceiling", "bottom",
                          "side", "sides", "lateral", "farfield", "far"}))
    return {ZoneRole::symmetry, "name contains \"" + *w + "\"; domain boundaries are slip planes", {}};
  if (z.type == "velocity-inlet" || z.type == "pressure-inlet" || z.type == "mass-flow-inlet")
    return {ZoneRole::velocity_inlet, "mesh zone type " + z.type, {}};
  if (z.type == "pressure-outlet" || z.type == "outflow") return {ZoneRole::pressure_outlet, "mesh zone type " + z.type, {}};
  if (z.type == "symmetry" || z.type == "pressure-far-field") return {ZoneRole::symmetry, "mesh zone type " + z.type, {}};
  return {ZoneRole::wall, "a wall of the car", aero_group(z.name)};
}

KuttaCase fsae_starter(const MeshInfo& mesh, std::filesystem::path mesh_path, const StarterOptions& o) {
  Errors errors;
  if (!(o.speed > 0) || !std::isfinite(o.speed)) errors.add("speed", "must be a positive speed in m/s");
  for (const auto& [name, role] : o.roles)
    if (std::none_of(mesh.zones.begin(), mesh.zones.end(), [&](const auto& z) { return z.name == name && z.kind != MeshZoneInfo::Kind::cell; }))
      errors.add("role/" + name, "the mesh has no face zone named \"" + name + "\"");
  KuttaCase c;
  c.mesh = std::move(mesh_path);
  c.models.viscous = ViscousModel::k_omega_sst;
  c.material = {"air", defaults::air_density, defaults::air_viscosity};
  const bool planar = mesh.dimension == 2;
  const Vec3 flow{1, 0, 0}, down = planar ? Vec3{0, -1, 0} : Vec3{0, 0, -1};
  // Low freestream turbulence, as on the road or in a wind tunnel.
  TurbulenceInput freestream;
  freestream.intensity = 0.01;
  freestream.viscosity_ratio = 10;
  std::vector<std::string> car;
  for (const auto& z : mesh.zones) {
    if (z.kind == MeshZoneInfo::Kind::cell) continue;
    const auto it = o.roles.find(z.name);
    const auto role = it != o.roles.end() ? it->second : suggest_role(z).role;
    if ((role == ZoneRole::interior) != (z.kind == MeshZoneInfo::Kind::interior)) {
      errors.add("role/" + z.name, z.kind == MeshZoneInfo::Kind::interior ? "is an interior face zone; only interior applies"
                                                                          : "is a boundary zone; interior applies only to interior face zones");
      continue;
    }
    switch (role) {
      case ZoneRole::velocity_inlet: {
        VelocityInletInput v;
        v.speed = o.speed;
        v.direction = flow;
        v.turbulence = freestream;
        c.boundaries.emplace(z.name, v);
        break;
      }
      case ZoneRole::pressure_outlet: {
        PressureOutletInput p;
        p.gauge_pressure = 0;
        p.turbulence = freestream;
        c.boundaries.emplace(z.name, p);
        break;
      }
      case ZoneRole::moving_ground: {
        WallInput w;
        w.motion = WallInput::Motion::moving;
        w.speed = o.speed;
        w.direction = flow;
        c.boundaries.emplace(z.name, w);
        break;
      }
      case ZoneRole::wall:
        c.boundaries.emplace(z.name, WallInput{});
        car.push_back(z.name);
        break;
      case ZoneRole::symmetry: c.boundaries.emplace(z.name, SymmetryInput{}); break;
      case ZoneRole::interior: break;  // the default for interior zones
    }
  }
  errors.finish();
  c.reference.velocity = o.speed;
  c.reference.density = defaults::air_density;
  c.reference.area = defaults::reference_area;
  c.reference.length = defaults::reference_length;
  // Second-order turbulence transport, the usual external-aero choice.
  c.methods.turbulent_kinetic_energy = Upwind::second_order;
  c.methods.specific_dissipation_rate = Upwind::second_order;
  if (!car.empty()) {
    c.reports.emplace("drag", ForceReportInput{ForceReportInput::Kind::drag, car, flow, true});
    c.reports.emplace("downforce", ForceReportInput{ForceReportInput::Kind::lift, car, down, true});
  }
  return c;
}

} // namespace kutta::setup
