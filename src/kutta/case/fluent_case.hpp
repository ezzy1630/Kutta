#pragma once
#include "kutta/case/sexpr.hpp"
#include <array>
#include <filesystem>
#include <map>
#include <set>

namespace kutta::fluent {
// Missing optionals are never Fluent defaults. issues records absent/unsupported
// required values; all source entries, including duplicate metadata, remain raw.
enum class IssueClass { unsupported, informational };
struct SettingIssue {
  std::string path, reason;
  IssueClass classification = IssueClass::unsupported;
  bool operator==(const SettingIssue &) const = default;
};
struct RawSettings {
  SExpr alist;
  std::set<std::string> interpreted;
  const SExpr *entry(std::string_view name) const;
  const SExpr *value(std::string_view name) const;
  std::vector<std::string> unparsed_names() const;
};
enum class ViscousModel {
  laminar,
  ke_standard,
  ke_rng,
  ke_realizable,
  kw_standard,
  kw_bsl,
  kw_sst
};
enum class Scheme {
  first_upwind = 0,
  second_upwind = 1,
  power_law = 2,
  central = 3,
  quick = 4,
  muscl = 6,
  bounded_central = 7,
  standard_pressure = 10,
  linear_pressure = 11,
  second_pressure = 12,
  body_force_weighted = 13,
  presto = 14
};
enum class Coupling { simple = 20, simplec = 21, piso = 22, coupled = 24 };
enum class Gradient { least_squares_cell, green_gauss_node, green_gauss_cell };
// Only code 0 is verified; other codes remain raw and issue-marked.
enum class Limiter { standard = 0 };
enum class WallOmega { correlation, tabulated };
enum class ZoneType {
  fluid,
  solid,
  wall,
  symmetry,
  velocity_inlet,
  pressure_outlet,
  interior,
  other
};
enum class VelocitySpec { magnitude_direction = 0, components = 1, normal = 2 };
enum class TurbulenceSpec {
  k_dissipation = 0,
  intensity_length = 1,
  intensity_viscosity_ratio = 2,
  intensity_diameter = 3
};
enum class WallMotion { stationary, translating, rotating };
struct ScalarValue {
  enum class Kind { constant, parameter, expression, profile, unparsed };
  Kind kind = Kind::unparsed;
  // Literal constants and resolved expression values have distinct provenance.
  std::optional<double> constant;
  std::string binding, expression_name, definition;
  std::optional<double> resolved;
  bool operator==(const ScalarValue &) const = default;
};
using ScalarVector = std::array<std::optional<ScalarValue>, 3>;
struct Turbulence {
  std::optional<ViscousModel> model;
  std::optional<bool> production_limiter, low_re, curvature_correction, sst_strain, kato_launder,
      vorticity_production, compressibility, shear_modification;
  std::optional<double> production_factor, curvature_factor, omega_floor;
  std::optional<WallOmega> wall_omega;
  std::map<std::string, double> coefficients;
  std::map<std::string, bool> unsupported_options;
  bool operator==(const Turbulence &) const = default;
};
struct Numerics {
  std::optional<Scheme> pressure, momentum, k, omega, epsilon;
  std::optional<Gradient> gradient;
  std::optional<Limiter> limiter;
  std::optional<bool> limiting, limit_to_face, limiter_filter, warped_face_correction,
      high_order_relaxation, flux_auto_select;
  std::optional<int> flux_index, wall_distance_method;
  std::optional<bool> wall_distance_quad_correction;
  std::optional<double> high_order_steady_relax;
  std::optional<Coupling> coupling;
  struct PseudoTime {
    std::optional<bool> formulation, automatic;
    // Raw integer selectors: 0 is off; other meanings are version dependent.
    std::optional<int> coupled_method, segregated_method, mode, length_scale_method;
    std::optional<double> scale_factor, time_step, automatic_time_step;
    bool operator==(const PseudoTime &) const = default;
  } pseudo;
  std::map<std::string, double> relaxation, pseudo_relaxation;
  std::optional<double> courant;
  bool operator==(const Numerics &) const = default;
};
struct Material {
  std::string name, type;
  std::optional<ScalarValue> density, viscosity;
  std::optional<bool> nonconstant_density, nonconstant_viscosity;
  RawSettings raw;
};
struct ReferenceValues {
  std::optional<double> area, length, density, velocity, pressure, viscosity;
  bool operator==(const ReferenceValues &) const = default;
};
struct BoundaryZone {
  std::int64_t id = 0;
  std::optional<std::int64_t> domain;
  std::string name, fluent_type;
  ZoneType type = ZoneType::other;
  std::optional<VelocitySpec> velocity_spec;
  std::optional<ScalarValue> velocity_magnitude, gauge_pressure;
  ScalarVector velocity, direction;
  std::optional<TurbulenceSpec> turbulence_spec;
  std::map<std::string, ScalarValue> turbulence;
  std::optional<WallMotion> wall_motion;
  std::optional<bool> relative_motion, specified_shear, prevent_reverse_flow, wall_components;
  std::optional<int> frame_of_reference, backflow_pressure_spec, shear_bc, rough_bc;
  std::optional<ScalarValue> angular_velocity, roughness_height, roughness_constant;
  ScalarVector axis_origin, axis_direction, shear;
  RawSettings raw;
};
struct CellZone {
  std::int64_t id = 0;
  std::optional<std::int64_t> domain;
  std::string name, fluent_type;
  ZoneType type = ZoneType::other;
  std::optional<std::string> material;
  std::optional<bool> frame_motion, porous, laminar, relative_resistance;
  std::optional<int> relative_to_zone, motion_spec;
  std::optional<ScalarValue> angular_velocity, porosity;
  ScalarVector frame_velocity, axis_origin, axis_direction, viscous_resistance, inertial_resistance;
  // Resistance direction bases; do not assume mesh axes.
  std::array<ScalarVector, 3> porous_directions;
  RawSettings raw;
};
struct ReportDefinition {
  enum class Kind { force, drag, lift, moment, other };
  std::string name, fluent_type;
  Kind kind = Kind::other;
  // zone_ids are resolved thread IDs; source surface IDs use a different ID space.
  std::optional<std::vector<std::int64_t>> zone_ids, surface_ids, reported_zone_ids;
  std::optional<bool> scaled;
  std::optional<std::vector<std::string>> zone_names;
  std::optional<std::array<double, 3>> direction;
  RawSettings raw;
};
// Topology metadata only: no geometry/connectivity is read here. IDs are Fluent
// global, 1-based, inclusive. These ranges bridge to lane B's eventual ID map.
struct ZoneRange {
  std::int64_t zone_id = 0;
  std::string name;
  std::uint64_t first = 0, last = 0;
  std::optional<int> type_code;
  std::uint64_t size() const { return last - first + 1; }
  bool operator==(const ZoneRange &) const = default;
};
struct MeshMetadata {
  std::optional<int> dimension;
  std::optional<std::uint64_t> cells, faces, nodes;
  std::vector<ZoneRange> cell_zones, face_zones;
};
struct NamedExpression {
  std::string name, parameter_id, definition;
};
struct FluentCase {
  std::optional<std::filesystem::path> source_path;
  std::optional<std::string> solver_version, origin;
  std::optional<int> dimension;
  std::optional<bool> double_precision, unsteady, pressure_based, energy, gravity;
  std::array<std::optional<double>, 3> gravity_vector;
  Turbulence viscous;
  Numerics numerics;
  std::optional<double> operating_pressure;
  ReferenceValues reference;
  std::vector<Material> materials;
  std::vector<NamedExpression> named_expressions;
  std::vector<BoundaryZone> boundaries;
  std::vector<CellZone> cell_zones;
  std::vector<ReportDefinition> reports;
  MeshMetadata mesh;
  RawSettings variables, config;
  std::map<std::int64_t, RawSettings> domains;
  std::map<std::string, std::vector<SExpr>> raw_sections;
  std::vector<SettingIssue> issues;
};
// Parse already-isolated sections (37 rpvars; 38 Cortex; 64 domains; 39 threads).
FluentCase extract_case(std::map<std::string, std::vector<SExpr>> sections,
                        std::optional<std::string> version = {},
                        std::optional<std::string> origin = {});
FluentCase read_case(const std::filesystem::path &);
std::string_view name(ViscousModel);
std::string_view name(Scheme);
std::string_view name(Coupling);
} // namespace kutta::fluent
