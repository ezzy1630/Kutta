#pragma once
#include "kutta/core/json_value.hpp"
#include "kutta/setup/fluent.hpp"
#include <functional>

namespace kutta::setup {
struct RunControls { int iterations=1000,monitor_interval=20; };
struct Setting {
  std::string path,label,group,unit,kind;
  Json value;
  Json::Array options;
  // Fluent's words for options (parallel) and for the current value.
  Json::Array option_labels;
  std::string value_label;
  std::optional<double> min,max;
  bool editable=false;
  // Whether the case as saved uses this setting. With depends_on, active is
  // derived: the controlling enum's current value is one of active_when.
  bool active=true;
  std::optional<std::string> depends_on;
  Json::Array active_when;
  std::string reason="not yet configurable in Kutta";
  // Internal bindings never appear on the wire. Inventory and overrides use
  // this same registry, so an editable field cannot have a missing setter.
  std::function<void(const Json&)> apply;
  std::function<std::string(const Json&)> validate;
};
std::vector<Setting> settings_inventory(fluent::FluentCase&,RunControls&,bool seeded=false);
Json settings_json(std::span<const Setting>);
// Validate every entry first, then apply atomically to the in-memory case.
// The source case file is never changed. Throws RejectedSettings with paths.
Json::Object apply_setup(const std::filesystem::path&,std::span<const Setting>);
} // namespace kutta::setup
