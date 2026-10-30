#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace kutta::setup {
// Permission is per named approximation, never a blanket ignore-unsupported mode.
enum class Approximation { wfgc, compressibility, ke_options };
struct ApproximationInfo {
  Approximation id;
  std::string_view name,setting,fluent_value,kutta_action,effect_documentation,measured_effect;
};
const ApproximationInfo& approximation_info(Approximation id);
std::vector<Approximation> parse_approximations(std::string_view names);
// Shared by the CLI output headers and V6's labelled diagnostic.
std::string approximation_notice(std::span<const Approximation> used);
} // namespace kutta::setup
