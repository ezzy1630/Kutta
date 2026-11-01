#pragma once
// A Kutta case: a mesh file plus a small versioned JSON setup (*.kutta.json),
// for users without a Fluent licence. It is lowered to the same in-memory
// FluentCase a .cas file produces, so validation, the settings inventory, the
// app and every solve path are shared. Members mirror the file: an absent
// optional is Fluent's default for a new case (docs/kutta-case.md lists and
// cites each), so a file states only what differs from Fluent's defaults.

#include "kutta/case/fluent_case.hpp"
#include "kutta/core/json_value.hpp"
#include "kutta/mesh/mesh.hpp"

#include <array>
#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace kutta::setup {
inline constexpr int kutta_case_version = 1;
using Vec3 = std::array<double, 3>;

enum class ViscousModel { laminar, k_omega_sst, k_epsilon_standard };
enum class Upwind { first_order, second_order };
enum class PressureInterpolation { standard, second_order };
enum class GradientMethod { least_squares_cell, green_gauss_cell };
enum class CouplingScheme { simple, simplec, coupled };
enum class FluxType { automatic, momentum, distance };
enum class GradientLimiter { standard, differentiable, none };

struct TurbulenceInput {
  enum class Method { intensity_viscosity_ratio, intensity_length_scale, intensity_hydraulic_diameter, k_omega, k_epsilon };
  std::optional<Method> method;
  std::optional<double> intensity, viscosity_ratio, length_scale, hydraulic_diameter, k, omega, epsilon;
  bool operator==(const TurbulenceInput&) const = default;
};
struct WallInput {
  enum class Motion { stationary, moving, rotating };
  std::optional<Motion> motion;
  std::optional<bool> relative;  // relative to the adjacent cell zone
  // moving: speed along direction, or velocity components
  std::optional<double> speed;
  std::optional<Vec3> direction, velocity;
  // rotating: angular speed (rad/s, right-handed about the axis)
  std::optional<double> angular_speed;
  std::optional<Vec3> axis_origin, axis_direction;
  bool operator==(const WallInput&) const = default;
};
// speed alone: magnitude normal to the boundary; speed and direction:
// magnitude and direction; velocity: components.
struct VelocityInletInput {
  std::optional<double> speed, gauge_pressure, intermittency;
  std::optional<Vec3> direction, velocity;
  std::optional<TurbulenceInput> turbulence;
  bool operator==(const VelocityInletInput&) const = default;
};
struct PressureOutletInput {
  enum class Backflow { total, static_pressure };
  std::optional<double> gauge_pressure, intermittency;
  std::optional<Backflow> backflow_pressure;
  std::optional<bool> prevent_reverse_flow;
  std::optional<TurbulenceInput> turbulence;  // backflow turbulence
  bool operator==(const PressureOutletInput&) const = default;
};
struct SymmetryInput { bool operator==(const SymmetryInput&) const = default; };
struct InteriorInput { bool operator==(const InteriorInput&) const = default; };
using BoundaryInput = std::variant<WallInput, VelocityInletInput, PressureOutletInput, SymmetryInput, InteriorInput>;

struct ForceReportInput {
  enum class Kind { force, drag, lift };
  Kind kind = Kind::force;
  std::vector<std::string> zones;
  std::optional<Vec3> direction;
  std::optional<bool> coefficient;
  bool operator==(const ForceReportInput&) const = default;
};
struct ResidualInput {
  std::optional<bool> check;
  std::optional<double> absolute_criterion;
  bool operator==(const ResidualInput&) const = default;
};
// Relaxation factors by equation; which set applies follows Fluent's
// formulation (pseudo time, classic under-relaxation, Coupled explicit).
struct RelaxationInput {
  std::optional<double> pressure, momentum, density, body_force, k, omega, epsilon, turbulent_viscosity, intermittency;
  bool operator==(const RelaxationInput&) const = default;
};

struct KuttaCase {
  std::filesystem::path mesh;  // as written: relative paths are relative to the setup file
  struct Models {
    ViscousModel viscous = ViscousModel::k_omega_sst;  // required in the file
    std::optional<bool> intermittency_transition;       // "transition": "none" | "intermittency"
    struct KOmegaSst {
      std::optional<bool> production_limiter, low_re_corrections, curvature_correction, compressibility_effects,
          kato_launder, vorticity_based_production;
      std::optional<double> production_limiter_clip_factor;
      bool operator==(const KOmegaSst&) const = default;
    } k_omega_sst;
    struct KEpsilon {
      std::optional<bool> compressibility_effects;
      bool operator==(const KEpsilon&) const = default;
    } k_epsilon;
    bool operator==(const Models&) const = default;
  } models;
  struct Material {
    std::optional<std::string> name;
    std::optional<double> density, viscosity;
    bool operator==(const Material&) const = default;
  } material;
  std::optional<double> operating_pressure;
  std::map<std::string, BoundaryInput> boundaries;  // by mesh zone name
  struct Reference {
    std::optional<double> area, length, density, velocity, pressure, viscosity, depth;
    bool operator==(const Reference&) const = default;
  } reference;
  struct Methods {
    std::optional<CouplingScheme> coupling;
    std::optional<bool> pseudo_transient, warped_face_gradient_correction, high_order_term_relaxation;
    std::optional<GradientMethod> gradient;
    std::optional<PressureInterpolation> pressure;
    std::optional<Upwind> momentum, turbulent_kinetic_energy, specific_dissipation_rate, turbulent_dissipation_rate, intermittency;
    std::optional<FluxType> flux_type;
    std::optional<GradientLimiter> gradient_limiter;
    bool operator==(const Methods&) const = default;
  } methods;
  struct Controls {
    std::optional<double> flow_courant_number, pseudo_time_scale_factor, pseudo_time_step;  // step absent: automatic
    RelaxationInput pseudo_relaxation, relaxation, explicit_relaxation;
    std::optional<double> min_turbulent_kinetic_energy, min_specific_dissipation_rate, min_turbulent_dissipation_rate,
        max_turbulent_viscosity_ratio;
    bool operator==(const Controls&) const = default;
  } controls;
  struct Initialization {
    enum class Method { hybrid, standard };
    std::optional<Method> method;
    std::optional<int> hybrid_iterations;
    std::optional<Vec3> velocity;
    std::optional<double> gauge_pressure, k, omega, epsilon, intermittency;
    bool operator==(const Initialization&) const = default;
  } initialization;
  std::map<std::string, ForceReportInput> reports;
  std::map<std::string, ResidualInput> residuals;  // by residual name ("continuity", "x-velocity", ...)
  bool operator==(const KuttaCase&) const = default;
};

// Kutta cases are recognised by extension: *.json (conventionally *.kutta.json).
bool is_kutta_case(const std::filesystem::path&);
// Strict parse: unknown properties, wrong types and invalid combinations are
// all collected and thrown as RejectedSettings with JSON paths.
KuttaCase parse_kutta_case(std::string_view json);
KuttaCase read_kutta_case_file(const std::filesystem::path&);
// The file text: stable key order, shortest round-trip numbers, two-space indent.
std::string kutta_case_text(const KuttaCase&);

// The mesh facts lowering depends on, all read from the PolyMesh.
struct MeshZoneInfo {
  std::int64_t id = 0;
  std::string name, type;  // mesh-file zone type ("wall", "interior", "fluid", ...)
  enum class Kind { boundary, interior, cell } kind = Kind::boundary;
  std::size_t count = 0;  // faces, or cells for cell zones
};
struct MeshInfo {
  int dimension = 3;          // 2: read as a one-cell extrusion with front/back caps
  bool polyhedral = false;    // any cell other than tet, pyramid, wedge or hex
  std::vector<MeshZoneInfo> zones;
  std::size_t cells = 0, faces = 0, nodes = 0;
};
MeshInfo mesh_info(const PolyMesh&);

// The FluentCase a Fluent writer would save for these settings on this mesh,
// with Fluent's defaults for everything the file leaves out. Throws
// RejectedSettings for a setup that does not match the mesh.
fluent::FluentCase lower_kutta_case(const KuttaCase&, const MeshInfo&);
// Patch kinds and source types follow the case's boundary types, not the
// mesh file's hints (meshers often export every zone as a wall).
void apply_boundary_kinds(const KuttaCase&, PolyMesh&);

// A case file loaded with its mesh, ready for the shared setup path.
struct LoadedKuttaCase {
  KuttaCase source;
  std::filesystem::path setup_path, mesh_path;
  fluent::FluentCase settings;
  PolyMesh mesh;
};
LoadedKuttaCase load_kutta_case(const std::filesystem::path& setup_file);
// Case identity for caches and run records: the SHA-256 of the setup file's
// and the mesh file's SHA-256 digests, so a changed mesh is a different case.
std::string kutta_case_identity(std::string_view setup_sha256, std::string_view mesh_sha256);

// The Kutta case expressing a FluentCase's settings: only values that differ
// from Fluent's defaults are written. `omitted` lists settings the format
// cannot express (path: reason); none of them is silently assumed.
struct ConvertedCase {
  KuttaCase kutta;
  std::vector<fluent::SettingIssue> omitted;
};
ConvertedCase kutta_case_from_fluent(const fluent::FluentCase&, const MeshInfo&, std::filesystem::path mesh);

// Starter cases (kutta new-case).
enum class ZoneRole { velocity_inlet, pressure_outlet, moving_ground, wall, symmetry, interior };
std::string_view role_name(ZoneRole);
std::optional<ZoneRole> parse_role(std::string_view);
struct RoleSuggestion {
  ZoneRole role;
  std::string reason;
  std::optional<std::string> group;  // aero group for walls: Front wing, Rear wing, Wing, Undertray / diffuser, Wheels, Body
};
RoleSuggestion suggest_role(const MeshZoneInfo&);
struct StarterOptions {
  double speed = 15;  // m/s, freestream and ground speed
  std::map<std::string, ZoneRole> roles;  // overrides of the suggestions, by zone name
};
// The FSAE external-aero template: flow along +x, up +z (+y in 2D), SST, air,
// velocity inlet, pressure outlet, moving ground, symmetry, stationary car
// walls, drag and downforce on every car wall.
KuttaCase fsae_starter(const MeshInfo&, std::filesystem::path mesh, const StarterOptions&);
} // namespace kutta::setup
