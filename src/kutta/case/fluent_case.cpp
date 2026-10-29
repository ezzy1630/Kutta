#include "kutta/case/fluent_case.hpp"
#include "kutta/case/hdf5.hpp"
#include <algorithm>
#include <charconv>
#include <fstream>
#include <limits>
#include <sstream>

namespace kutta::fluent {
const SExpr *RawSettings::entry(std::string_view key) const {
  if (auto p = alist_entry(alist, key))
    return p;
  for (const auto &e : alist.elements)
    if (!e.elements.empty() && e.elements[0].atom_text() == key)
      return &e;
  return nullptr;
}
const SExpr *RawSettings::value(std::string_view key) const {
  auto e = entry(key);
  return e ? entry_value(*e) : nullptr;
}
std::vector<std::string> RawSettings::unparsed_names() const {
  std::set<std::string> keys;
  for (auto &e : alist.elements)
    if (!e.elements.empty())
      if (auto s = e.elements[0].atom_text())
        if (!interpreted.contains(std::string(*s)))
          keys.emplace(*s);
  return {keys.begin(), keys.end()};
}
namespace {
using K = SExpr::Kind;
std::optional<int> integer(const SExpr *x) {
  if (x && x->kind == K::integer && x->integer >= std::numeric_limits<int>::min() &&
      x->integer <= std::numeric_limits<int>::max())
    return static_cast<int>(x->integer);
  return {};
}
std::optional<bool> boolean(const SExpr *x) {
  if (x && x->kind == K::boolean)
    return x->boolean;
  return {};
}
std::optional<double> number(const SExpr *x) { return x ? x->number() : std::nullopt; }
std::optional<std::string> text(const SExpr *x) {
  if (x)
    if (auto s = x->atom_text())
      return std::string(*s);
  return {};
}
class Extractor {
  FluentCase &c;
  RawSettings &r;
  std::string prefix;
  void accepted(std::string_view key, bool ok, bool required) {
    if (ok)
      r.interpreted.emplace(key);
    else if (r.entry(key) || required)
      c.issues.push_back({prefix + std::string(key), r.entry(key) ? "unparsed: unsupported value"
                                                                  : "absent: no default assumed"});
  }

public:
  Extractor(FluentCase &out, RawSettings &raw, std::string path)
      : c(out), r(raw), prefix(std::move(path)) {}
  std::optional<double> num(std::string_view key, bool required = true) {
    auto selector = text(r.value(std::string(key) + "-sel"));
    if (selector && *selector != "constant") {
      accepted(key, false, true);
      c.issues.back().reason = "unparsed for solve: selector=" + *selector;
      return {};
    }
    auto v = number(r.value(key));
    accepted(key, v.has_value(), required);
    return v;
  }
  std::optional<int> code(std::string_view key, bool required = true) {
    auto v = integer(r.value(key));
    accepted(key, v.has_value(), required);
    return v;
  }
  std::optional<bool> flag(std::string_view key, bool required = true) {
    auto v = boolean(r.value(key));
    accepted(key, v.has_value(), required);
    return v;
  }
  std::optional<std::string> str(std::string_view key, bool required = true) {
    auto v = text(r.value(key));
    accepted(key, v.has_value(), required);
    return v;
  }
  template <class E>
  std::optional<E> enumeration(std::string_view key, std::initializer_list<E> allowed,
                               bool required = true) {
    auto v = integer(r.value(key));
    if (v)
      for (E e : allowed)
        if (static_cast<int>(e) == *v) {
          accepted(key, true, required);
          return e;
        }
    accepted(key, false, required);
    return {};
  }
  std::optional<ScalarValue> scalar(std::string_view key, bool required = true) {
    auto entry = r.entry(key);
    if (!entry) {
      accepted(key, false, required);
      return {};
    }
    const SExpr *v = entry_value(*entry);
    ScalarValue result;
    auto selector = text(r.value(std::string(key) + "-sel"));
    if (selector && *selector != "constant") {
      result.kind =
          *selector == "expression" ? ScalarValue::Kind::expression : ScalarValue::Kind::unparsed;
      result.binding = text(r.value(std::string(key) + "-expr")).value_or(*selector);
      accepted(key, false, true);
      return result;
    }
    if (auto n = number(v)) {
      result.kind = ScalarValue::Kind::constant;
      result.constant = n;
      accepted(key, true, required);
      return result;
    }
    // BC/property alternatives are ordered; the first selects the active method.
    const SExpr *method = nullptr;
    if (entry->elements.size() > 1 && entry->elements[1].kind == K::list)
      method = &entry->elements[1];
    if (v && v->kind == K::list && !v->elements.empty() && v->elements[0].kind == K::list)
      method = &v->elements[0];
    if (!method || method->elements.empty()) {
      accepted(key, false, required);
      return result;
    }
    auto tag = method->elements[0].atom_text();
    auto mv = entry_value(*method);
    if (tag == "constant" && number(mv)) {
      result.kind = ScalarValue::Kind::constant;
      result.constant = number(mv);
      accepted(key, true, required);
      return result;
    }
    result.binding = to_string(*method);
    if (tag == "parameter") {
      result.kind = ScalarValue::Kind::parameter;
      result.binding = text(mv).value_or(to_string(*method));
      for (std::size_t i = 2; i < entry->elements.size(); ++i) {
        auto &branch = entry->elements[i];
        if (!branch.elements.empty() && branch.elements[0].atom_text() == "constant")
          result.constant = number(entry_value(branch));
      }
    } else if (tag == "expression") {
      result.kind = ScalarValue::Kind::expression;
      result.binding = text(mv).value_or(to_string(*method));
    } else if (tag == "profile") {
      result.kind = ScalarValue::Kind::profile;
      result.binding = to_string(*method);
    }
    accepted(key, false, required);
    c.issues.back().reason = "unparsed for solve: " + to_string(*method);
    return result;
  }
  ScalarVector vector(std::array<std::string, 3> keys, bool required = true) {
    return {scalar(keys[0], required), scalar(keys[1], required), scalar(keys[2], required)};
  }
  void vector_override(std::string_view key, ScalarVector &target) {
    auto e = r.entry(key);
    if (!e)
      return;
    auto v = entry_value(*e);
    if (v && v->kind == K::boolean && !v->boolean) {
      r.interpreted.emplace(key);
      return;
    }
    const std::vector<SExpr> *parts = &e->elements;
    std::size_t begin = 1;
    if (v && v->kind == K::list && v->elements.size() == 3) {
      parts = &v->elements;
      begin = 0;
    }
    target = {};
    if (parts->size() != begin + 3) {
      accepted(key, false, true);
      return;
    }
    RawSettings components;
    for (int i = 0; i < 3; ++i) {
      SExpr component;
      component.elements.push_back(SExpr::symbol(std::to_string(i)));
      const auto &part = (*parts)[begin + i];
      if (part.number())
        component.elements.push_back(part);
      else if (part.kind == K::list) {
        if (!part.elements.empty() && part.elements[0].atom_text())
          component.elements.push_back(part);
        else
          component.elements.insert(component.elements.end(), part.elements.begin(),
                                    part.elements.end());
      } else
        component.elements.push_back(part);
      components.alist.elements.push_back(std::move(component));
    }
    Extractor sub(c, components, prefix + std::string(key) + "/");
    target = sub.vector({"0", "1", "2"});
    bool constant = true;
    for (auto &component : target)
      constant = constant && component && component->kind == ScalarValue::Kind::constant;
    accepted(key, constant, true);
  }
};
ZoneType zone_type(std::string_view s) {
  if (s == "fluid")
    return ZoneType::fluid;
  if (s == "solid")
    return ZoneType::solid;
  if (s == "wall")
    return ZoneType::wall;
  if (s == "symmetry")
    return ZoneType::symmetry;
  if (s == "velocity-inlet")
    return ZoneType::velocity_inlet;
  if (s == "pressure-outlet")
    return ZoneType::pressure_outlet;
  if (s == "interior")
    return ZoneType::interior;
  return ZoneType::other;
}
void model(FluentCase &c) {
  auto &r = c.variables;
  Extractor x(c, r, "rp/");
  Extractor cfg(c, c.config, "case-config/");
  auto three = cfg.flag("rp-3d?");
  if (three)
    c.dimension = *three ? 3 : 2;
  c.double_precision = cfg.flag("rp-double?");
  c.unsteady = cfg.flag("rp-unsteady?");
  c.pressure_based = cfg.flag("rp-seg?");
  auto lam = cfg.flag("rp-lam?"), ke = cfg.flag("rp-ke?"), kw = cfg.flag("rp-kw?");
  auto &t = c.viscous;
  if (lam && *lam && ke && !*ke && kw && !*kw)
    t.model = ViscousModel::laminar;
  else if (lam && !*lam && ke && *ke && kw && !*kw) {
    auto rng = x.flag("rng-ke-on?"), realizable = x.flag("realizable-epsilon?");
    if (rng && realizable && !(*rng && *realizable))
      t.model = *rng ? ViscousModel::ke_rng
                     : (*realizable ? ViscousModel::ke_realizable : ViscousModel::ke_standard);
  } else if (lam && !*lam && kw && *kw && ke && !*ke) {
    auto sst = x.flag("kw-sst-on?"), bsl = x.flag("kw-bsl-on?"), standard = x.flag("kw-std-on?");
    auto geko = x.flag("kw-geko-on?"), easm = x.flag("kw-easm-on?"),
         wj = x.flag("kw-wj-bsl-earsm-on?");
    if (sst && bsl && standard && int(*sst) + int(*bsl) + int(*standard) == 1 &&
        !geko.value_or(false) && !easm.value_or(false) && !wj.value_or(false))
      t.model =
          *sst ? ViscousModel::kw_sst : (*bsl ? ViscousModel::kw_bsl : ViscousModel::kw_standard);
  }
  if (!t.model)
    c.issues.push_back(
        {"viscous-model", "unparsed: absent, conflicting or unsupported case-config/model flags"});
  const bool epsilon = t.model == ViscousModel::ke_standard || t.model == ViscousModel::ke_rng ||
                       t.model == ViscousModel::ke_realizable;
  t.production_limiter =
      x.flag(epsilon ? "ke-production-limiter?" : "kw-not-trans-sst-production-limiter?");
  t.production_factor = x.num(epsilon ? "ke-pk-factor" : "sst-pk-factor");
  t.low_re = x.flag(epsilon ? "ke-low-re?" : "kw-low-re-mod?");
  t.curvature_correction = x.flag("curvature-correction?");
  t.curvature_factor = x.num("curvature-correction-ccurv");
  t.sst_strain = x.flag("kw-sst-strain?");
  t.kato_launder = x.flag(epsilon ? "ke-kato-launder?" : "kw-not-trans-sst-kato-launder?");
  t.vorticity_production = x.flag("kw-vorticity-pk?");
  t.shear_modification = x.flag("kw-shear-mod?");
  t.compressibility = x.flag("turb-compress-mod?");
  t.omega_floor = x.num("omega/min-limit");
  auto wall = x.str("turbulence/wall-omega-treatment");
  if (wall == "correlation")
    t.wall_omega = WallOmega::correlation;
  else if (wall == "tabulated")
    t.wall_omega = WallOmega::tabulated;
  else if (wall)
    c.issues.push_back(
        {"rp/turbulence/wall-omega-treatment", "unparsed: unknown wall-omega method"});
  for (auto key :
       {"kw-sig-w",     "kw-sig-k",          "kw-beta-i",  "kw-beta-star-inf", "kw-alpha-0",
        "kw-alpha-inf", "kw-alpha-star-inf", "kw-r-w",     "kw-r-k",           "kw-r-beta",
        "kw-mt0",       "kw-zeta-star",      "sst-a1",     "sst-beta-i1",      "sst-beta-i2",
        "sst-sig-k1",   "sst-sig-k2",        "sst-sig-w1", "sst-sig-w2",       "sst-kappa"})
    if (auto v = x.num(key, false))
      t.coefficients.emplace(key, *v);
}
void numerics(FluentCase &c) {
  Extractor x(c, c.variables, "rp/");
  auto &n = c.numerics;
  auto convection = {Scheme::first_upwind,   Scheme::second_upwind, Scheme::power_law,
                     Scheme::central,        Scheme::quick,         Scheme::muscl,
                     Scheme::bounded_central};
  n.pressure = x.enumeration("pressure/scheme", {Scheme::standard_pressure, Scheme::linear_pressure,
                                                 Scheme::second_pressure,
                                                 Scheme::body_force_weighted, Scheme::presto});
  n.momentum = x.enumeration("mom/scheme", convection);
  n.k = x.enumeration("k/scheme", convection);
  n.omega = x.enumeration("omega/scheme", convection);
  n.coupling = x.enumeration(
      "flow/scheme", {Coupling::simple, Coupling::simplec, Coupling::piso, Coupling::coupled});
  auto cell = x.flag("recon/cell-lsq?"), node = x.flag("recon/node-lsq?");
  if (cell && node)
    n.gradient = *cell ? Gradient::least_squares_cell
                       : (*node ? Gradient::green_gauss_node : Gradient::green_gauss_cell);
  n.limiter = x.enumeration("recon/lim/type", {Limiter::standard});
  n.limiting = x.flag("recon/limit?");
  n.limit_to_face = x.flag("recon/lim/dx-to-face?");
  n.limiter_filter = x.flag("recon/lim/filter/on?");
  n.warped_face_correction = x.flag("recon/cell-lsf?");
  n.high_order_relaxation = x.flag("recon/relax/relax?");
  n.high_order_steady_relax = x.num("recon/relax/steady-urf");
  n.flux_auto_select = x.flag("pbs/flux-auto-select?");
  n.flux_index = x.code("pbs/flux-index");
  n.wall_distance_method = x.code("wall-distance/method");
  n.wall_distance_quad_correction = x.flag("wall-distance/quad-correction?");
  n.pseudo.formulation = x.flag("pseudo-transient-formulation?");
  n.pseudo.automatic = x.flag("pseudo-auto-time-step?");
  n.pseudo.coupled_method = x.code("pseudo-time-method/coupled-pbns/dt-method");
  n.pseudo.segregated_method = x.code("pseudo-time-method/segregated-pbns/dt-method");
  n.pseudo.mode = x.code("pseudo-transient-mode");
  n.pseudo.scale_factor = x.num("pseudo-auto-time-step-scale-factor");
  n.pseudo.time_step = x.num("pseudo-time-step");
  n.pseudo.automatic_time_step = x.num("pseudo-auto-time-step");
  n.pseudo.length_scale_method = x.code("pseudo/autotime-lscale-option");
  n.courant = x.num("pressure-coupled/courant-number");
  for (auto eq : {"pressure", "mom", "k", "omega", "density", "body-force", "turb-viscosity"}) {
    auto key = std::string(eq) + "/relax";
    if (auto v = x.num(key))
      n.relaxation.emplace(eq, *v);
    key = std::string(eq) + "/pseudo-relax";
    if (auto v = x.num(key, false))
      n.pseudo_relaxation.emplace(key, *v);
  }
  for (auto key : {"pressure-coupled/pressure/pseudo-explicit-relax",
                   "pressure-coupled/mom/pseudo-explicit-relax"})
    if (auto v = x.num(key))
      n.pseudo_relaxation.emplace(key, *v);
}
void materials(FluentCase &c) {
  auto v = c.variables.value("materials");
  if (!v || v->kind != K::list) {
    c.issues.push_back({"rp/materials", "absent or unparsed"});
    return;
  }
  c.variables.interpreted.insert("materials");
  for (auto &e : v->elements) {
    if (e.elements.size() < 2 || !text(&e.elements[0]) || !text(&e.elements[1])) {
      c.issues.push_back({"rp/materials", "unparsed material record"});
      continue;
    }
    Material m;
    m.name = *text(&e.elements[0]);
    m.type = *text(&e.elements[1]);
    m.raw.alist.elements.assign(e.elements.begin() + 2, e.elements.end());
    Extractor x(c, m.raw, "material/" + m.name + "/");
    m.density = x.scalar("density");
    m.viscosity = x.scalar("viscosity", m.type == "fluid");
    if (m.density)
      m.nonconstant_density = m.density->kind != ScalarValue::Kind::constant;
    if (m.viscosity)
      m.nonconstant_viscosity = m.viscosity->kind != ScalarValue::Kind::constant;
    c.materials.push_back(std::move(m));
  }
}
void reference(FluentCase &c) {
  Extractor x(c, c.variables, "rp/");
  c.operating_pressure = x.num("operating-pressure");
  auto v = c.variables.value("reference-values");
  RawSettings nested;
  if (v && v->kind == K::list)
    nested.alist = *v;
  Extractor sub(c, nested, "rp/reference-values/");
  auto get = [&](const char *key) {
    auto top = "reference-" + std::string(key);
    return c.variables.entry(top) ? x.num(top) : sub.num(key);
  };
  c.reference = {get("area"),     get("length"),   get("density"),
                 get("velocity"), get("pressure"), get("viscosity")};
  if (v)
    c.variables.interpreted.insert("reference-values");
}
void thread(FluentCase &c, const SExpr &e) {
  if (e.elements.size() != 3 || e.elements[1].elements.size() < 3 || e.elements[2].kind != K::list)
    throw std::runtime_error("malformed section 39");
  auto &h = e.elements[1].elements;
  auto id = integer(&h[0]);
  auto type = text(&h[1]), label = text(&h[2]);
  if (!id || *id <= 0 || !type || !label)
    throw std::runtime_error("invalid thread identity");
  auto domain = h.size() > 3 ? integer(&h[3]) : std::optional<int>{};
  if (h.size() > 3 && !domain)
    throw std::runtime_error("invalid thread domain");
  if (!domain)
    c.issues.push_back({"zone/" + std::to_string(*id) + "/domain", "absent: domain not assumed"});
  auto ztype = zone_type(*type);
  if (ztype == ZoneType::fluid || ztype == ZoneType::solid) {
    CellZone z;
    z.id = *id;
    z.name = *label;
    z.domain = domain;
    z.fluent_type = *type;
    z.type = ztype;
    z.raw.alist = e.elements[2];
    Extractor x(c, z.raw, "zone/" + std::to_string(z.id) + "/");
    z.material = x.str("material");
    z.frame_motion = x.flag("mrf-motion?");
    z.angular_velocity = x.scalar("mrf-omega");
    z.relative_to_zone = x.code("mrf-relative-to-thread");
    z.motion_spec = x.code("motion-spec");
    z.relative_resistance = x.flag("rel-vel-resistance?", false);
    z.frame_velocity = x.vector({"mrf-grid-x-vel", "mrf-grid-y-vel", "mrf-grid-z-vel"});
    z.axis_origin = x.vector({"mrf-x-origin", "mrf-y-origin", "mrf-z-origin"});
    z.axis_direction = x.vector({"mrf-ai", "mrf-aj", "mrf-ak"});
    x.vector_override("reference-frame-velocity-components", z.frame_velocity);
    x.vector_override("reference-frame-axis-origin-components", z.axis_origin);
    x.vector_override("reference-frame-axis-direction-components", z.axis_direction);
    z.porous = x.flag("porous?");
    z.laminar = x.flag("laminar?");
    z.porosity = x.scalar("porosity");
    z.viscous_resistance =
        x.vector({"porous-r/1", "porous-r/2", "porous-r/3"}, z.porous.value_or(false));
    z.inertial_resistance =
        x.vector({"porous-c/1", "porous-c/2", "porous-c/3"}, z.porous.value_or(false));
    for (int i = 0; i < 3; ++i) {
      auto p = "direction-" + std::to_string(i + 1) + "/";
      z.porous_directions[i] =
          x.vector({p + "x", p + "y", p + "z"}, z.porous.value_or(false) && i < 2);
    }
    c.cell_zones.push_back(std::move(z));
    return;
  }
  BoundaryZone z;
  z.id = *id;
  z.name = *label;
  z.domain = domain;
  z.fluent_type = *type;
  z.type = ztype;
  z.raw.alist = e.elements[2];
  Extractor x(c, z.raw, "zone/" + std::to_string(z.id) + "/");
  bool inlet = ztype == ZoneType::velocity_inlet, outlet = ztype == ZoneType::pressure_outlet,
       wall = ztype == ZoneType::wall;
  if (inlet || wall) {
    z.velocity_spec = x.enumeration(
        "velocity-spec",
        {VelocitySpec::magnitude_direction, VelocitySpec::components, VelocitySpec::normal}, inlet);
    z.velocity_magnitude = x.scalar("vmag");
    z.velocity = x.vector({"u", "v", "w"});
    z.direction = x.vector({"ni", "nj", "nk"});
    z.frame_of_reference = x.code("frame-of-reference", inlet);
    x.vector_override("velocity-components", z.velocity);
    x.vector_override("reference-direction-components", z.direction);
    x.vector_override("flow-direction-component", z.direction);
  }
  if (inlet || outlet) {
    z.gauge_pressure = x.scalar("p");
    z.turbulence_spec = x.enumeration(
        "ke-spec", {TurbulenceSpec::k_dissipation, TurbulenceSpec::intensity_length,
                    TurbulenceSpec::intensity_viscosity_ratio, TurbulenceSpec::intensity_diameter});
    for (auto key : {"turb-intensity", "turb-viscosity-ratio", "turb-length-scale",
                     "turb-hydraulic-diam", "k", "o", "e"}) {
      const bool direct = z.turbulence_spec == TurbulenceSpec::k_dissipation;
      const bool omega_model = c.viscous.model == ViscousModel::kw_sst ||
                               c.viscous.model == ViscousModel::kw_bsl ||
                               c.viscous.model == ViscousModel::kw_standard;
      std::string_view k = key;
      bool required =
          (direct && (k == "k" || k == (omega_model ? "o" : "e"))) ||
          (!direct && z.turbulence_spec &&
           (k == "turb-intensity" ||
            (k == "turb-length-scale" && z.turbulence_spec == TurbulenceSpec::intensity_length) ||
            (k == "turb-viscosity-ratio" &&
             z.turbulence_spec == TurbulenceSpec::intensity_viscosity_ratio) ||
            (k == "turb-hydraulic-diam" &&
             z.turbulence_spec == TurbulenceSpec::intensity_diameter)));
      if (auto v = x.scalar(key, required))
        z.turbulence.emplace(k == "o" ? "omega" : (k == "e" ? "epsilon" : key), *v);
    }
    if (outlet) {
      z.backflow_pressure_spec = x.code("press-spec-gen");
      z.prevent_reverse_flow = x.flag("prevent-reverse-flow", false);
    }
  }
  if (wall) {
    auto motion = x.code("motion-bc");
    auto rotating = x.flag("rotating?");
    if (motion == 0)
      z.wall_motion = WallMotion::stationary;
    else if (motion == 1 && rotating)
      z.wall_motion = *rotating ? WallMotion::rotating : WallMotion::translating;
    else
      c.issues.push_back({"zone/" + std::to_string(z.id) + "/motion-bc", "unparsed wall motion"});
    z.relative_motion = x.flag("relative?");
    z.angular_velocity = x.scalar("omega");
    z.axis_origin = x.vector({"x-origin", "y-origin", "z-origin"});
    z.axis_direction = x.vector({"ai", "aj", "ak"});
    z.specified_shear = x.flag("specified-shear?");
    z.shear_bc = x.code("shear-bc");
    z.rough_bc = x.code("rough-bc");
    z.shear = x.vector({"shear-x", "shear-y", "shear-z"});
    z.roughness_height = x.scalar("roughness-height");
    z.roughness_constant = x.scalar("roughness-const");
    x.vector_override("wall-translation-vector", z.direction);
    x.vector_override("axis-origin-component", z.axis_origin);
    x.vector_override("axis-direction-component", z.axis_direction);
    x.vector_override("position-of-rotation-axis", z.axis_origin);
    x.vector_override("direction-component-of-rotation-axis", z.axis_direction);
    x.vector_override("shear-stress-components", z.shear);
  }
  if (ztype == ZoneType::other)
    c.issues.push_back({"zone/" + std::to_string(z.id), "unparsed boundary type: " + *type});
  c.boundaries.push_back(std::move(z));
}
void reports(FluentCase &c) {
  auto v = c.variables.value("monitor/report-definitions");
  if (!v || v->kind != K::list) {
    c.issues.push_back({"rp/monitor/report-definitions", "absent or unparsed"});
    return;
  }
  c.variables.interpreted.insert("monitor/report-definitions");
  for (auto &item : v->elements) {
    ReportDefinition report;
    report.raw.alist = item;
    Extractor x(c, report.raw, "report/");
    report.name = x.str("name").value_or("<unparsed>");
    auto def = report.raw.entry("report-definition");
    if (!def || def->elements.size() < 2) {
      c.issues.push_back({"report/" + report.name, "unparsed report definition"});
      c.reports.push_back(std::move(report));
      continue;
    }
    report.fluent_type = text(&def->elements[1]).value_or("<unparsed>");
    auto &t = report.fluent_type;
    if (t == "force" || t == "force-force")
      report.kind = ReportDefinition::Kind::force;
    else if (t == "drag" || t == "force-drag")
      report.kind = ReportDefinition::Kind::drag;
    else if (t == "lift" || t == "force-lift")
      report.kind = ReportDefinition::Kind::lift;
    else if (t == "moment" || t == "force-moment")
      report.kind = ReportDefinition::Kind::moment;
    RawSettings props;
    props.alist.elements.assign(def->elements.begin() + 2, def->elements.end());
    auto list = [&](std::string_view key) -> const SExpr * {
      auto e = props.entry(key);
      if (!e)
        return nullptr;
      auto value = entry_value(*e);
      return value && value->kind == K::list ? value : e;
    };
    for (auto key : {"zones", "zone-list", "zone-ids", "thread-ids"})
      if (auto zs = list(key)) {
        std::vector<std::int64_t> ids;
        bool valid = true;
        for (auto &z : zs->elements) {
          if (z.atom_text() == key)
            continue;
          if (z.kind == K::integer && z.integer > 0)
            ids.push_back(z.integer);
          else
            valid = false;
        }
        if (valid)
          report.zone_ids = std::move(ids);
      }
    for (auto key : {"zone-names", "thread-names"})
      if (auto zs = list(key)) {
        std::vector<std::string> names;
        bool valid = true;
        for (auto &z : zs->elements) {
          if (z.atom_text() == key)
            continue;
          if (auto s = text(&z))
            names.push_back(*s);
          else
            valid = false;
        }
        if (valid)
          report.zone_names = std::move(names);
      }
    for (auto key : {"force-vector", "direction-vector", "direction"})
      if (auto ds = list(key)) {
        std::vector<double> values;
        bool valid = true;
        for (auto &z : ds->elements) {
          if (z.atom_text() == key)
            continue;
          if (auto n = z.number())
            values.push_back(*n);
          else
            valid = false;
        }
        if (valid && values.size() == 3)
          report.direction = std::array<double, 3>{values[0], values[1], values[2]};
      }
    if (report.kind != ReportDefinition::Kind::other) {
      if (!report.zone_ids && !report.zone_names)
        c.issues.push_back({"report/" + report.name + "/zones", "unparsed force zones"});
      if (!report.direction && report.kind != ReportDefinition::Kind::moment)
        c.issues.push_back({"report/" + report.name + "/direction", "unparsed force direction"});
    } else
      c.issues.push_back({"report/" + report.name, "unparsed report type: " + report.fluent_type});
    c.reports.push_back(std::move(report));
  }
}
} // namespace
FluentCase extract_case(std::map<std::string, std::vector<SExpr>> sections,
                        std::optional<std::string> version, std::optional<std::string> origin) {
  FluentCase c;
  c.raw_sections = std::move(sections);
  c.solver_version = std::move(version);
  c.origin = std::move(origin);
  auto it = c.raw_sections.find("Rampant Variables");
  if (it == c.raw_sections.end() || it->second.size() != 1 || it->second[0].elements.size() != 2 ||
      integer(&it->second[0].elements[0]) != 37)
    throw std::runtime_error("case requires one valid Rampant Variables section 37");
  c.variables.alist = it->second[0].elements[1];
  if (c.variables.alist.kind != K::list || c.variables.alist.tail)
    throw std::runtime_error("rpvars must be a proper alist");
  if (auto config = c.variables.value("case-config"); config && config->kind == K::list) {
    c.config.alist = *config;
    c.variables.interpreted.insert("case-config");
  }
  model(c);
  numerics(c);
  materials(c);
  reference(c);
  reports(c);
  if (auto domains = c.raw_sections.find("Domain Variables"); domains != c.raw_sections.end()) {
    for (auto &form : domains->second) {
      if (form.elements.size() != 2 || form.elements[1].kind != K::list)
        throw std::runtime_error("malformed domain variables");
      for (auto &domain : form.elements[1].elements) {
        if (domain.elements.empty() || domain.elements[0].kind != K::integer ||
            domain.elements[0].integer <= 0)
          throw std::runtime_error("invalid domain ID");
        RawSettings raw;
        raw.alist.elements.assign(domain.elements.begin() + 1, domain.elements.end());
        auto id = domain.elements[0].integer;
        if (!c.domains.emplace(id, std::move(raw)).second)
          throw std::runtime_error("duplicate domain ID");
        c.issues.push_back(
            {"domain/" + std::to_string(id),
             "unparsed domain-specific variables retained by name; not applied globally"});
      }
    }
  }
  if (auto ts = c.raw_sections.find("Thread Variables"); ts != c.raw_sections.end())
    for (auto &e : ts->second)
      thread(c, e);
  else
    c.issues.push_back({"Thread Variables", "absent: boundary/cell-zone settings unavailable"});
  std::set<std::int64_t> ids;
  for (auto &z : c.cell_zones)
    if (!ids.insert(z.id).second)
      throw std::runtime_error("duplicate thread ID");
  for (auto &z : c.boundaries)
    if (!ids.insert(z.id).second)
      throw std::runtime_error("duplicate thread ID");
  // File writers reorder alists; typed records use identity order, raw sections retain source
  // order.
  std::sort(c.boundaries.begin(), c.boundaries.end(), [](auto &a, auto &b) { return a.id < b.id; });
  std::sort(c.cell_zones.begin(), c.cell_zones.end(), [](auto &a, auto &b) { return a.id < b.id; });
  std::sort(c.materials.begin(), c.materials.end(),
            [](auto &a, auto &b) { return a.name < b.name; });
  std::sort(c.reports.begin(), c.reports.end(), [](auto &a, auto &b) { return a.name < b.name; });
  if (std::adjacent_find(c.reports.begin(), c.reports.end(),
                         [](auto &a, auto &b) { return a.name == b.name; }) != c.reports.end())
    throw std::runtime_error("duplicate report name");
  if (std::adjacent_find(c.materials.begin(), c.materials.end(),
                         [](auto &a, auto &b) { return a.name == b.name; }) != c.materials.end())
    throw std::runtime_error("duplicate material name");
  if (!c.solver_version)
    c.issues.push_back({"solver-version", "absent: no default assumed"});
  return c;
}
std::string_view name(ViscousModel v) {
  switch (v) {
  case ViscousModel::laminar:
    return "laminar";
  case ViscousModel::ke_standard:
    return "k-epsilon standard";
  case ViscousModel::ke_rng:
    return "k-epsilon RNG";
  case ViscousModel::ke_realizable:
    return "k-epsilon realizable";
  case ViscousModel::kw_standard:
    return "k-omega standard";
  case ViscousModel::kw_bsl:
    return "k-omega BSL";
  case ViscousModel::kw_sst:
    return "k-omega SST";
  }
  return "unparsed";
}
std::string_view name(Scheme s) {
  switch (s) {
  case Scheme::first_upwind:
    return "FOU";
  case Scheme::second_upwind:
    return "SOU";
  case Scheme::power_law:
    return "power-law";
  case Scheme::central:
    return "central";
  case Scheme::quick:
    return "QUICK";
  case Scheme::muscl:
    return "MUSCL";
  case Scheme::bounded_central:
    return "bounded-central";
  case Scheme::standard_pressure:
    return "standard";
  case Scheme::linear_pressure:
    return "linear";
  case Scheme::second_pressure:
    return "second-order";
  case Scheme::body_force_weighted:
    return "body-force-weighted";
  case Scheme::presto:
    return "PRESTO!";
  }
  return "unparsed";
}
std::string_view name(Coupling v) {
  switch (v) {
  case Coupling::simple:
    return "SIMPLE";
  case Coupling::simplec:
    return "SIMPLEC";
  case Coupling::piso:
    return "PISO";
  case Coupling::coupled:
    return "Coupled";
  }
  return "unparsed";
}
} // namespace kutta::fluent
