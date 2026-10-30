#pragma once
#include "kutta/core/json_value.hpp"
#include "kutta/core/types.hpp"
#include <filesystem>
#include <optional>
#include <vector>

namespace kutta::post {
// Integrated wall loads of a finished run (DIR/loads.json). Forces are on the
// body in N; moments are about Loads::origin in N m.
struct ZoneLoads {
  std::string name;
  std::int64_t zone_id = 0;
  double area = 0;
  Vec3d pressure{}, viscous{}, pressure_moment{}, viscous_moment{};
};
struct LoadReference {
  std::optional<double> area, density, velocity, length, pressure;
};
struct Loads {
  Vec3d origin;
  LoadReference reference;
  // Union of the case's force-report zones, in report order.
  std::vector<std::string> report_zones;
  std::vector<ZoneLoads> zones;
  std::string provenance;
  bool approximate = false, modified = false;
  int iterations = 0;
};
Json loads_json(const Loads &);
// Unknown properties are ignored so later additive fields stay readable.
Loads parse_loads(const Json &);

struct AeroSpec {
  struct Group {
    std::string name;
    std::vector<std::string> zones;
  };
  struct Vehicle {
    Vec3d front_axle;
    double wheelbase = 0;
  };
  // Unit vectors after parsing.
  Vec3d drag{1, 0, 0}, lift{0, 0, 1}, moment_centre;
  std::optional<std::vector<std::string>> car_zones;
  std::vector<Group> groups;
  std::optional<Vehicle> vehicle;
};
// A spec the user must fix (exit 2), as opposed to an unreadable run (exit 1).
class InvalidAeroSpec : public Error {
 public:
  using Error::Error;
};
AeroSpec parse_aero_spec(const Json &);
AeroSpec read_aero_spec(const std::filesystem::path &);
// The aero.json object. Zone names are checked against the loads here, so an
// unknown or non-car zone also throws InvalidAeroSpec.
Json aero(const Loads &, const AeroSpec &);
} // namespace kutta::post
