#pragma once

#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace kutta::setup {
// Permission is per named approximation, never a blanket ignore-unsupported mode.
enum class Approximation { wfgc, compressibility, ke_options, reverse_flow };
struct ApproximationInfo {
  Approximation id;
  std::string_view name,setting,fluent_value,kutta_action,effect_documentation,measured_effect;
  // The app's issue label and plain-language explanation (docs/app-protocol.md).
  std::string_view label,explanation;
  // A few words naming what the run does differently, for result badges.
  std::string_view short_label;
};
const ApproximationInfo& approximation_info(Approximation id);
// The complete opt-in allowlist, in declaration order.
std::span<const ApproximationInfo> approximations();
std::vector<Approximation> parse_approximations(std::string_view names);
// Shared by the CLI output headers and V6's labelled diagnostic.
std::string approximation_notice(std::span<const Approximation> used);
} // namespace kutta::setup
