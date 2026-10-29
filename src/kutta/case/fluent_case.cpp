#include "kutta/case/fluent_case.hpp"
#include "kutta/case/hdf5.hpp"
#include "kutta/core/parse.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
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
std::string_view trim(std::string_view s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == s.npos)
    return {};
  return s.substr(first, s.find_last_not_of(" \t\r\n") - first + 1);
}
std::optional<double> constant_expression(std::string_view definition) {
  definition = trim(definition);
  auto open = definition.find('[');
  if (open == definition.npos || definition.back() != ']')
    return {};
  auto numeric = trim(definition.substr(0, open));
  if (numeric.starts_with('+'))
    numeric.remove_prefix(1);
  auto number = parse_double(numeric);
  if (!number || !std::isfinite(*number))
    return {};
  const auto unit = trim(definition.substr(open + 1, definition.size() - open - 2));
  double scale = 1;
  if (unit == "in")
    scale = 0.0254;
  else if (unit == "mm")
    scale = 0.001;
  else if (unit != "m/s" && unit != "m" && unit != "Pa" && unit != "kg/m^3" && unit != "rad/s")
    return {};
  auto value = *number * scale;
  return std::isfinite(value) ? std::optional<double>{value} : std::nullopt;
}
void expressions(FluentCase &c) {
  auto list = c.variables.value("named-expressions");
  if (!list)
    return;
  if (list->kind != K::list || list->tail) {
    c.issues.push_back({"rp/named-expressions", "unsupported: malformed expression list"});
    return;
  }
  for (const auto &item : list->elements) {
    NamedExpression e;
    e.name = text(alist_value(item, "name")).value_or("");
    e.parameter_id = text(alist_value(item, "parameterid")).value_or("");
    e.definition = text(alist_value(item, "definition")).value_or("");
    c.named_expressions.push_back(std::move(e));
  }
  c.variables.interpreted.insert("named-expressions");
}
bool resolve_expression(const FluentCase &c, ScalarValue &v) {
  if (v.kind != ScalarValue::Kind::parameter && v.kind != ScalarValue::Kind::expression)
    return false;
  const NamedExpression *match = nullptr;
  bool ambiguous = false;
  for (const auto &e : c.named_expressions)
    if ((v.kind == ScalarValue::Kind::parameter && e.parameter_id == v.binding) || e.name == v.binding) {
      ambiguous = ambiguous || match != nullptr;
      match = &e;
    }
  if (ambiguous)
    return false;
  if (match) {
    v.expression_name = match->name;
    v.definition = match->definition;
  } else if (v.kind == ScalarValue::Kind::expression)
    v.definition = v.binding;
  else
    return false;
  v.resolved = constant_expression(v.definition);
  return v.resolved.has_value();
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
      const bool resolved = resolve_expression(c, result);
      accepted(key, resolved, true);
      if (!resolved)
        c.issues.back().reason = "unresolved expression: " + result.binding;
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
    } else if (tag == "expression") {
      result.kind = ScalarValue::Kind::expression;
      result.binding = text(mv).value_or(to_string(*method));
    } else if (tag == "profile") {
      result.kind = ScalarValue::Kind::profile;
      result.binding = to_string(*method);
    }
    const bool resolved = resolve_expression(c, result);
    accepted(key, resolved, required);
    if (!resolved)
      c.issues.back().reason = "unresolved expression/profile for solve: " + to_string(*method);
    return result;
  }
  ScalarVector vector(std::array<std::string, 3> keys, bool required = true) {
    return {scalar(keys[0], required), scalar(keys[1], required), scalar(keys[2], required)};
  }
  void vector_override(std::string_view key, ScalarVector &target, bool compare = false) {
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
    const auto scalars = target;
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
      constant = constant && component && (component->constant || component->resolved);
    accepted(key, constant, true);
    if (compare) {
      bool conflict = false;
      for (int i = 0; i < 3; ++i) {
        if (!scalars[i] || !target[i])
          continue;
        const auto &a = *scalars[i], &b = *target[i];
        const auto av = a.constant ? a.constant : a.resolved;
        const auto bv = b.constant ? b.constant : b.resolved;
        if (av && bv)
          conflict = conflict || std::abs(*av - *bv) >
              1e-8 * std::max({1.0, std::abs(*av), std::abs(*bv)});
        else
          conflict = conflict || a != b;
      }
      if (conflict)
        c.issues.push_back({prefix + std::string(key), "unsupported: scalar and vector forms disagree"});
    }
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
bool epsilon_model(const FluentCase &c) {
  return c.viscous.model == ViscousModel::ke_standard || c.viscous.model == ViscousModel::ke_rng ||
         c.viscous.model == ViscousModel::ke_realizable;
}
void model_flag_support(FluentCase &c) {
  for (const auto &key : c.variables.unparsed_names()) {
    if (!key.starts_with("kw-") && !key.starts_with("ke-") && !key.starts_with("sst-"))
      continue;
    auto flag = boolean(c.variables.value(key));
    if (!flag)
      continue;
    c.viscous.unsupported_options.emplace(key, *flag);
    if (!*flag)
      continue;
    const bool inactive = c.viscous.model == ViscousModel::laminar ||
        (key.starts_with("ke-") ? !epsilon_model(c) : epsilon_model(c));
    c.issues.push_back({"rp/" + key, inactive
        ? "inactive model flag retained; does not apply to the selected viscous model"
        : "unsupported: active model flag is not in the model allowlist",
        inactive ? IssueClass::informational : IssueClass::unsupported});
  }
}
// The table describes the supported physical problem, not Fluent's defaults.
// New flags are rejected even when a writer introduces them in a future version.
void config_support(FluentCase &c) {
  std::map<std::string, SExpr> expected;
  auto add = [&](std::string_view keys, bool value) {
    std::istringstream words{std::string(keys)};
    for (std::string key; words >> key;) {
      SExpr v;
      v.kind = K::boolean;
      v.boolean = value;
      expected.emplace(key, v);
    }
  };
  add("rp-seg? rp-3d? rp-visc?", true);
  add("rp-acoustics? rp-optics? rp-atm? rp-axi? rp-des? rp-dpm-cache? rp-dpm-cache-serial? "
      "rp-unsteady? rp-dual-time? rp-amg? rp-hvac? rp-inviscid? rp-kklw? rp-les? "
      "rp-structure? rp-structure-damping? rp-structure-linear? rp-structure-nonlinear? "
      "rp-structure-time? rp-lsf? rp-net? rp-react? rp-sa? rp-sas? rp-sge? rp-spe? "
      "rp-spe-part? rp-potential? rp-lithium? rp-dual-potential? rp-electrolysis? "
      "rp-electro-chem? rp-spe-site? rp-spe-surf? rp-trans-one-eqn? rp-trans-sst? "
      "rp-absorbing-media? rp-v2f? sg-cylindrical? sg-disco? sg-montecarlo? sg-bee-gees? "
      "sg-crev? sg-dpm? sg-dqmom-iem? sg-dtrm? sg-dynmesh? sg-overset? sg-ecfm? sg-geqn? "
      "sg-hg? sg-inert? sg-ignite? sg-network? sg-pcb? sg-node-udm? sg-noniterative? "
      "sg-nox? sg-melt? sg-micro-mix? sg-mphase? sg-p1? rp-sootspecies? sg-par-premix? "
      "sg-pb? sg-pdf? sg-pdf-transport? sg-premixed? sg-pull? sg-rosseland? sg-rsm? "
      "sg-s2s? sg-soot? sg-sox? sg-spark? sg-swirl? sg-twotemp? sg-udm? sg-uds? "
      "sg-addon? sg-addon-loaded? sg-vfr? sg-solar? sg-wetsteam? sg-moistair? sg-wallfilm? "
      "rp-spe-surf-rc? rp-watercontent? sg-vbm? rp-ablation?", false);
  expected.emplace("rp-transition-model", SExpr::symbol("none"));
  // These describe precision, execution, initialization or registered storage;
  // neither boolean value changes the selected physical equations.
  const std::set<std::string> metadata = {"rp-double?", "rp-graphics?", "rp-host?", "rp-thread?",
      "hyb-init?", "dpm-cache?", "rp-wallfilm?"};
  const bool turbulent = c.viscous.model && c.viscous.model != ViscousModel::laminar;
  for (auto key : {"rp-turb?", "rp-trb-scl?"}) {
    SExpr v;
    v.kind = K::boolean;
    v.boolean = turbulent;
    expected.emplace(key, v);
  }
  for (auto key : {"rp-lam?", "rp-ke?", "rp-kw?"}) {
    SExpr v;
    v.kind = K::boolean;
    v.boolean = std::string_view(key) == "rp-lam?" ? c.viscous.model == ViscousModel::laminar
        : (std::string_view(key) == "rp-ke?" ? epsilon_model(c)
          : c.viscous.model == ViscousModel::kw_sst || c.viscous.model == ViscousModel::kw_bsl ||
            c.viscous.model == ViscousModel::kw_standard);
    expected.emplace(key, v);
  }
  for (const auto &e : c.config.alist.elements) {
    if (e.elements.empty() || !e.elements[0].atom_text()) {
      c.issues.push_back({"case-config", "unsupported malformed flag"});
      continue;
    }
    const std::string key(*e.elements[0].atom_text());
    auto value = entry_value(e);
    if (key == "rf-energy?")
      continue; // Classified after the active materials are known.
    if (metadata.contains(key) && boolean(value)) {
      c.config.interpreted.insert(key);
      continue;
    }
    auto it = expected.find(key);
    if (value && it != expected.end() && *value == it->second)
      c.config.interpreted.insert(key);
    else
      c.issues.push_back({"case-config/" + key, it == expected.end()
          ? "unsupported: flag is not in the supported configuration allowlist"
          : "unsupported: expected " + to_string(it->second)});
  }
}
void property_support(FluentCase &c) {
  bool constant = !c.cell_zones.empty();
  for (auto &z : c.cell_zones) {
    if (z.type != ZoneType::fluid) {
      constant = false;
      continue;
    }
    auto m = std::find_if(c.materials.begin(), c.materials.end(),
                           [&](auto &m) { return z.material == m.name; });
    const bool ok = m != c.materials.end() && m->density && m->viscosity &&
        (m->density->constant || m->density->resolved) &&
        (m->viscosity->constant || m->viscosity->resolved);
    constant = constant && ok;
    if (!ok)
      c.issues.push_back({"zone/" + std::to_string(z.id) + "/material",
                         "unsupported: requires constant density and viscosity"});
  }
  if (c.energy == true)
    c.issues.push_back({"case-config/rf-energy?", constant
        ? "energy ignored (decoupled): constant density and viscosity"
        : "unsupported: energy coupling cannot be excluded without constant density and viscosity",
        constant ? IssueClass::informational : IssueClass::unsupported});
  if (c.gravity == true)
    c.issues.push_back({"rp/gravity?", constant
        ? "gravity enabled; constant density has no buoyancy coupling"
        : "unsupported: gravity with nonconstant or unverified density",
        constant ? IssueClass::informational : IssueClass::unsupported});
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
  c.energy = cfg.flag("rf-energy?");
  c.gravity = x.flag("gravity?");
  c.gravity_vector = {x.num("gravity/x", c.gravity == true),
                      x.num("gravity/y", c.gravity == true),
                      x.num("gravity/z", c.gravity == true)};
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
  config_support(c);
  const bool epsilon = epsilon_model(c);
  if (t.model && t.model != ViscousModel::laminar && t.model != ViscousModel::kw_sst)
    c.issues.push_back({"viscous-model", "unsupported model: " + std::string(name(*t.model))});
  for (auto key : {"ke-easm-on?", "kw-geko-on?", "kw-easm-on?", "kw-wj-bsl-earsm-on?"})
    if (auto v = x.flag(key, false); v == true)
      c.issues.push_back({"rp/" + std::string(key), "unsupported model option"});
  if (epsilon) {
    for (auto key : {"ke-std-wall?", "ke-non-eqm", "ke-enh-wall?", "ke-scalable-wall?",
                     "ke-nwt-ml?", "ke-udf-wall?", "ke-rough-wall-treatment-r14?"})
      if (auto v = x.flag(key, false))
        t.unsupported_options.emplace(key, *v);
  }
  if (t.model == ViscousModel::laminar) {
    model_flag_support(c);
    return;
  }
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
  model_flag_support(c);
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
  const bool turbulent = c.viscous.model != ViscousModel::laminar;
  if (turbulent) {
    n.k = x.enumeration("k/scheme", convection);
    if (epsilon_model(c))
      n.epsilon = x.enumeration("epsilon/scheme", convection);
    else
      n.omega = x.enumeration("omega/scheme", convection);
  }
  n.coupling = x.enumeration(
      "flow/scheme", {Coupling::simple, Coupling::simplec, Coupling::piso, Coupling::coupled});
  auto cell = x.flag("recon/cell-lsq?"), node = x.flag("recon/node-lsq?");
  if (cell && node)
    n.gradient = *cell ? Gradient::least_squares_cell
                       : (*node ? Gradient::green_gauss_node : Gradient::green_gauss_cell);
  if (n.gradient == Gradient::green_gauss_node)
    c.issues.push_back({"rp/recon/node-lsq?",
        "unverified mapping: recon/node-lsq? to Green-Gauss node-based"});
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
  std::vector<std::string> equations = {"pressure", "mom", "density", "body-force"};
  if (turbulent) {
    equations.push_back("k");
    equations.push_back(epsilon_model(c) ? "epsilon" : "omega");
    equations.push_back("turb-viscosity");
  }
  const bool pseudo = n.pseudo.formulation == true;
  const bool coupled = n.coupling == Coupling::coupled;
  for (const auto &eq : equations) {
    auto key = std::string(eq) + "/relax";
    if (auto v = x.num(key, !pseudo))
      n.relaxation.emplace(eq, *v);
    key = std::string(eq) + "/pseudo-relax";
    if (auto v = x.num(key, pseudo && (!coupled || (eq != "pressure" && eq != "mom"))))
      n.pseudo_relaxation.emplace(key, *v);
  }
  for (auto key : {"pressure-coupled/pressure/pseudo-explicit-relax",
                   "pressure-coupled/mom/pseudo-explicit-relax"})
    if (auto v = x.num(key, pseudo && coupled))
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
      m.nonconstant_density = !m.density->constant && !m.density->resolved;
    if (m.viscosity)
      m.nonconstant_viscosity = !m.viscosity->constant && !m.viscosity->resolved;
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
// Only recognized switches can be neutral or interpreted. Unknown switches are
// held even when false; a future writer must not silently expand solver support.
void zone_support(FluentCase &c, RawSettings &raw, ZoneType type, std::int64_t id) {
  std::map<std::string, SExpr> neutral;
  auto flags = [&](std::string_view keys, bool value) {
    std::istringstream words{std::string(keys)};
    for (std::string key; words >> key;) {
      SExpr v;
      v.kind = K::boolean;
      v.boolean = value;
      neutral.emplace(key, v);
    }
  };
  flags("geom-disable? geom-dir-spec", false);
  if (type == ZoneType::fluid) {
    flags("sources? fixed? cylindrical-fixed-var? mgrid-enable-transient? mgrid-motion? "
          "solid-motion? deactivated? les-embedded? contact-property? glass? conical? cursys? "
          "alt-inertial-form? corey? stone-1? stone-2? rp-cbox-p1 rp-edit-p1 rp-cbox-p2 rp-edit-p2 "
          "rel-perm-tabular-p1? rel-perm-tabular-p2? non-equib-thermal? fanzone? limit-flow-fan "
          "tan-source-term rad-source-term axial-source-term read-fan-curve numerical-beach? "
          "beach-multi-dir? porous-structure? anisotropic-spe-diff? axis-origin-component "
          "axis-direction-component moving-mesh-velocity-components moving-mesh-axis-origin-components "
          "solid-motion-velocity-components solid-motion-axis-origin-components "
          "solid-motion-axis-direction-components none", false);
    // Radiation/chemistry/multiphase switches are inactive in the globally
    // allowlisted single-phase, nonreacting problem. Their model flags are checked.
    flags("radiating? laminar-mut-zero? none? equib-thermal? react? electrolyte? "
          "mp-boiling-zone? beach-damp-relative? active-wetsteam-zone?", true);
  }
  if (type == ZoneType::velocity_inlet || type == ZoneType::pressure_outlet) {
    flags("open-channel-wave-bc? ocw-vel-segregated? moving-object-direction-components "
          "secondary-phase-direction-components primary-phase-direction-components "
          "axis-direction-component axis-origin-component non-equil-boundary? vm-mass-conservation? "
          "volumetric-synthetic-turbulence-generator? species-in-mole-fractions? equ-required? "
          "radial-direction-component solar-direction? solar-irradiation? "
          "fensapice/drop/bccustom? fensapice/drop/dv? fensapice/dpm/surface-injection? "
          "fensapice/drop/cv? mixing-plane-thread? open-channel? radial? avg-press-spec? "
          "targeted-mf-boundary? flow-direction-component velocity-component velocity-components", false);
    flags("vm-streamwise-fluct? parallel-collimated-beam? fensapice/drop/vrh? "
          "solar-fluxes? radiating-s2s-surface?", true);
  }
  if (type == ZoneType::wall) {
    flags("planar-conduction? moving? wall-translation-vector velocity-components int-rad? rough-nasa? rough-shin-et-al? rough-data? inactive "
          "dpm/wall-heat-exchange? dpm/film-condensation? dpm/film-bl-model? dpm/film-stripping? "
          "dpm/initialize-lwf? film-velocity mc-bsource-p? mc-poldfun-p? pold-pair-list-rad "
          "component-of-radiation-direction critical-zone? solar-direction? solar-irradiation? "
          "fsi-interface? react? partially-catalytic? echem-reaction? faradaic-heat? fwd? aero? "
          "cmplx? norm? ablation/flux? film-wall? flux-momentum-components film-relative-vel? "
          "film-bc-imp-press? film-source? momentum-source-components film-phase-change? "
          "film-momentum-coupling? film-splash-wall? film-boundary-separation? film-contact-angle? "
          "film-vof-coupling-high? film-vof-coupling-low? thermal-stabilization? fensapice/ice/hflux? "
          "fensapice/drop/vwet?", false);
    flags("rough-option? dpm/bc-erosion-generic? dpm/bc-erosion-finnie? dpm/bc-erosion-mclaury? "
          "dpm/bc-erosion-oka? dpm/bc-erosion-dnv? dpm/bc-erosion-dnv-ductile? dpm/bc-erosion-shear? "
          "dpm/bc-erosion-shielding? dpm/allow-lwf-to-vof? dpm/allow-vof-to-lwf? "
          "radiating-s2s-surface? parallel-collimated-beam? solar-fluxes?", true);
  }
  if (type == ZoneType::interior)
    flags("is-not-a-rans-les-interface", true);
  auto zero = [&](std::string key) {
    SExpr v;
    v.kind = K::integer;
    neutral.emplace(std::move(key), v);
  };
  if (type == ZoneType::wall)
    zero("method?");
  if (type == ZoneType::velocity_inlet) {
    zero("coordinate-system");
    zero("omega-swirl");
  }
  if (type == ZoneType::pressure_outlet)
    for (auto key : {"gen-nrbc-spec", "p-backflow-spec", "p-backflow-spec-gen", "press-spec-gen"})
      zero(key);
  for (const auto &key : raw.unparsed_names()) {
    const auto *v = raw.value(key);
    auto it = neutral.find(key);
    const bool is_switch = boolean(v).has_value() || key.ends_with('?') || it != neutral.end();
    if (!is_switch)
      continue;
    const bool equal = v && it != neutral.end() &&
        (*v == it->second || (v->number() && it->second.number() && v->number() == it->second.number()));
    if (equal)
      raw.interpreted.insert(key);
    else
      c.issues.push_back({"zone/" + std::to_string(id) + "/" + key,
          it == neutral.end() ? "unsupported: switch is not in this zone type's allowlist"
                              : "unsupported: expected neutral value " + to_string(it->second)});
  }
  // These numeric selectors must be checked even when their typed value was read.
  for (auto key : {"coordinate-system", "omega-swirl", "gen-nrbc-spec", "p-backflow-spec",
                   "p-backflow-spec-gen", "press-spec-gen"}) {
    if (!raw.interpreted.contains(key) || !neutral.contains(key))
      continue;
    auto v = raw.value(key);
    if (!v || v->number() != 0.0)
      c.issues.push_back({"zone/" + std::to_string(id) + "/" + key,
                         "unsupported: expected neutral value 0"});
  }
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
    zone_support(c, z.raw, ztype, z.id);
    if (ztype == ZoneType::solid)
      c.issues.push_back({"zone/" + std::to_string(z.id), "unsupported: solid cell zone"});
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
  if (wall) {
    auto motion = x.code("motion-bc");
    auto rotating = x.flag("rotating?", motion == 1);
    if (motion == 0)
      z.wall_motion = WallMotion::stationary;
    else if (motion == 1 && rotating)
      z.wall_motion = *rotating ? WallMotion::rotating : WallMotion::translating;
    else
      c.issues.push_back({"zone/" + std::to_string(z.id) + "/motion-bc", "unparsed wall motion"});
    z.wall_components = x.flag("components?", z.wall_motion == WallMotion::translating);
    if (z.wall_motion == WallMotion::translating && z.wall_components)
      z.velocity_spec = *z.wall_components ? VelocitySpec::components : VelocitySpec::magnitude_direction;
  }
  if (inlet)
    z.velocity_spec = x.enumeration("velocity-spec",
        {VelocitySpec::magnitude_direction, VelocitySpec::components, VelocitySpec::normal});
  if (inlet || z.wall_motion == WallMotion::translating) {
    if (z.velocity_spec == VelocitySpec::components) {
      const std::string key = inlet ? "velocity-component" : "velocity-components";
      const bool vector_active = z.raw.entry(key) && boolean(z.raw.value(key)) != false;
      z.velocity = x.vector({"u", "v", "w"}, !vector_active);
      x.vector_override(key, z.velocity, true);
    } else if (z.velocity_spec) {
      z.velocity_magnitude = x.scalar("vmag");
      if (z.velocity_spec == VelocitySpec::magnitude_direction) {
        z.direction = x.vector({"ni", "nj", "nk"});
        x.vector_override("reference-direction-components", z.direction, true);
        x.vector_override("flow-direction-component", z.direction, true);
        if (wall)
          x.vector_override("wall-translation-vector", z.direction, true);
      }
    }
    z.frame_of_reference = x.code("frame-of-reference", inlet);
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
      z.backflow_pressure_spec = x.code("p-backflow-spec");
      z.prevent_reverse_flow = x.flag("prevent-reverse-flow?", false);
    }
  }
  if (wall) {
    z.relative_motion = x.flag("relative?");
    if (z.wall_motion == WallMotion::rotating) {
      z.angular_velocity = x.scalar("omega");
      z.axis_origin = x.vector({"x-origin", "y-origin", "z-origin"});
      z.axis_direction = x.vector({"ai", "aj", "ak"});
    }
    z.specified_shear = x.flag("specified-shear?");
    z.shear_bc = x.code("shear-bc");
    z.rough_bc = x.code("rough-bc");
    z.shear = x.vector({"shear-x", "shear-y", "shear-z"});
    z.roughness_height = x.scalar("roughness-height");
    z.roughness_constant = x.scalar("roughness-const");
    x.vector_override("axis-origin-component", z.axis_origin);
    x.vector_override("axis-direction-component", z.axis_direction);
    x.vector_override("position-of-rotation-axis", z.axis_origin);
    x.vector_override("direction-component-of-rotation-axis", z.axis_direction);
    x.vector_override("shear-stress-components", z.shear);
  }
  if (ztype == ZoneType::other)
    c.issues.push_back({"zone/" + std::to_string(z.id), "unparsed boundary type: " + *type});
  zone_support(c, z.raw, ztype, z.id);
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
    for (auto key : {"surface-ids", "zones", "zone-list", "zone-ids", "thread-ids"})
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
        if (valid) {
          std::sort(ids.begin(), ids.end());
          ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
          if (std::string_view(key) == "surface-ids")
            report.surface_ids = std::move(ids);
          else
            report.reported_zone_ids = std::move(ids);
        } else
          c.issues.push_back({"report/" + report.name + "/zones", "unparsed report IDs"});
      }
    for (auto key : {"surface-names", "zone-names", "thread-names"})
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
        if (valid) {
          std::sort(names.begin(), names.end());
          names.erase(std::unique(names.begin(), names.end()), names.end());
          report.zone_names = std::move(names);
        }
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
    Extractor properties(c, props, "report/" + report.name + "/");
    report.scaled = properties.flag("scaled?", report.kind != ReportDefinition::Kind::other);
    if (props.interpreted.contains("scaled?"))
      report.raw.interpreted.insert("scaled?");
    if (report.zone_names && !report.zone_names->empty()) {
      std::vector<std::int64_t> ids;
      bool valid = true;
      for (const auto &name : *report.zone_names) {
        std::vector<std::int64_t> matches;
        for (const auto &z : c.boundaries)
          if (z.name == name)
            matches.push_back(z.id);
        if (matches.size() != 1)
          valid = false;
        else
          ids.push_back(matches[0]);
      }
      if (valid)
        report.zone_ids = std::move(ids);
      else
        c.issues.push_back({"report/" + report.name + "/zones",
                           "unverified: surface name has no unique matching thread"});
    } else if (report.surface_ids || report.reported_zone_ids)
      c.issues.push_back({"report/" + report.name + "/zones",
                         "unverified: ID-only report; surface IDs are not thread IDs"});
    if (report.kind != ReportDefinition::Kind::other) {
      if (!report.zone_ids && !report.zone_names && !report.surface_ids && !report.reported_zone_ids)
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
  expressions(c);
  model(c);
  numerics(c);
  materials(c);
  reference(c);
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
  reports(c);
  property_support(c);
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
