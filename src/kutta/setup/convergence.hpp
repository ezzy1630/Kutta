#pragma once
#include "kutta/case/fluent_case.hpp"
#include "kutta/core/json_value.hpp"
#include <array>
#include <span>
#include <optional>
#include <string_view>

namespace kutta::setup {
// The scaled residuals the solver reports, in iteration-event order.
inline constexpr std::array<std::string_view,6> residual_equations{"continuity","x-velocity","y-velocity","z-velocity","k","omega"};
// The residuals reported for a closure: residual_equations for SST, k and
// epsilon for k-epsilon, continuity and momentum only for laminar flow.
std::vector<std::string_view> solved_residuals(std::optional<fluent::ViscousModel>);
struct ResidualCriterion {
  std::string equation;
  bool check=false;
  std::optional<double> absolute{};
};
struct ReportCriterion {
  std::string report;
  bool active=false;
  double stop=1e-3;
  int previous_values=5,ignore_before=100;
};
// Fluent's convergence conditions. Off by default, so a run without a setup
// file still stops at its iteration limit.
struct ConvergenceControls {
  enum class Combine { all, any };
  bool check=false;
  Combine combine=Combine::all;
  std::vector<ResidualCriterion> residuals;
  std::vector<ReportCriterion> reports;
  bool has_conditions() const;
  // The setup-event object: checked residuals and active reports only.
  Json::Object json() const;
};
// A checked residual needs a criterion. Overrides are validated one at a time,
// so this pairing is checked after they are applied, before any solve.
std::vector<fluent::SettingIssue> criteria_issues(const ConvergenceControls&);
struct ConvergenceSample {
  bool met=false;
  // Iteration-event object: {"met":...,"conditions":[...]}.
  Json::Object json;
  // Met conditions as "residuals" or "report:<name>".
  std::vector<std::string> met_by;
};
// Evaluates one monitor sample at a time; keeps only the report history the
// lag residual needs.
class ConvergenceMonitor {
 public:
  explicit ConvergenceMonitor(ConvergenceControls);
  const ConvergenceControls& controls() const { return controls_; }
  // residuals by equation; reports by exact name, projected total force in N.
  ConvergenceSample sample(int iteration,const std::map<std::string,double>& residuals,
                           const std::map<std::string,double>& reports);
 private:
  ConvergenceControls controls_;
  std::vector<std::vector<double>> history_;
};
// Fluent's lag residual over the last previous_values differences, relative to
// the newest value. Null with too few values or a zero newest value.
std::optional<double> lag_residual(std::span<const double> history,int previous_values);
} // namespace kutta::setup
