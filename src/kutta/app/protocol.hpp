#pragma once
#include "kutta/setup/fluent.hpp"
#include <ostream>

namespace kutta::app {
struct Issue {
  std::string setting, label, classification, approximation, explanation, evidence;
  std::string measured_effect;  // approximable issues: the measured effect, also ending the explanation
  std::string short_label;      // approximable issues: a few words for result badges
};
std::vector<Issue> support_issues(const fluent::FluentCase&, const PolyMesh&, bool seeded);
void write_issues(std::ostream&, std::span<const Issue>);
int inspect_json(int argc, char** argv);
int devices(int argc, char** argv);
void error_event(std::string_view message, std::span<const Issue> issues = {});
} // namespace kutta::app
