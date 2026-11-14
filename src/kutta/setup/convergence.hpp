#pragma once
#include "kutta/case/fluent_case.hpp"
#include "kutta/core/json_value.hpp"
#include <array>
#include <deque>
#include <span>
#include <optional>
#include <string_view>

namespace kutta::setup {
// The scaled residuals the solver reports, in iteration-event order.
inline constexpr std::array<std::string_view,6> residual_equations{"continuity","x-velocity","y-velocity","z-velocity","k","omega"};
// The residuals reported for a closure: residual_equations for SST, k and
// epsilon for k-epsilon, continuity and momentum only for laminar flow.
// transition: Fluent's intermittency model is active (adds "intermit").
// planar: a 2D case, whose Fluent residual list has no z-velocity; Kutta's
// one-cell extrusion keeps w at round-off, so it is not reported.
std::vector<std::string_view> solved_residuals(std::optional<fluent::ViscousModel>,bool transition=false,bool planar=false);
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
  // A force, drag, lift or downforce report: only these can be conditions,
  // and they are the reports the forces-settled condition watches.
  bool force=true;
};
// Kutta's forces-settled condition (not a Fluent setting): over the last
// `window` iterations every force report stays within `tolerance` of its mean
// magnitude (max - min <= tolerance * scale, scale = max(|mean|, 1% of the
// largest report's |mean|)). It must hold at two samples one iteration apart,
// so a cycle whose period divides the monitor interval cannot look steady.
// Residuals that stall above their criteria are waived once the forces have
// held for twice the window (docs/app-protocol.md).
struct SettleCriterion {
  bool active=false;
  double tolerance=1e-3;
  int window=300;
};
// Kutta's statistically-steady stop (not a Fluent setting, not convergence).
// Some steady RANS solutions are unstable: the iteration settles into a
// bounded limit cycle (the hex Ahmed body's wake on finer meshes under
// SIMPLEC) and no amount of iterating converges it. Once the forces have failed to settle for twice the
// forces-settled window, every force report is sampled at every iteration.
// At a monitor sample n, the samples in (n-2*window, n-window] and
// (n-window, n] are compared: the run stops, unconverged, when for every
// report the two Hann-weighted window means differ by at most `tolerance`
// times the scale used by forces-settled and the second range is within
// `range` times the first (plus tolerance times the scale), i.e. the cycle is
// neither drifting, growing nor decaying. A decaying cycle is still
// converging, so it runs on to convergence or the iteration limit. The verdict reports each force's
// Hann-weighted mean, min and max over both windows (docs/app-protocol.md,
// Convergence criteria). Fluent's analogue is a report definition's
// "Average Over" running mean.
struct AverageCriterion {
  bool active=false;
  double tolerance=1e-3,range=.1;
  int window=500;
};
// Fluent's convergence conditions plus Kutta's forces-settled condition. The
// settings inventory turns checking and the forces-settled condition on when
// the case has a force report; a default-constructed ConvergenceControls
// checks nothing.
struct ConvergenceControls {
  enum class Combine { all, any };
  bool check=false;
  Combine combine=Combine::all;
  std::vector<ResidualCriterion> residuals;
  std::vector<ReportCriterion> reports;
  SettleCriterion settle;
  AverageCriterion average;
  // Forces-settled is a condition only with at least one force report.
  bool settle_condition() const;
  // The statistically-steady stop needs checking on and a force report.
  bool average_stop() const;
  // The first iteration whose forces are sampled for averaging.
  int average_start() const {return 2*settle.window+1;}
  bool has_conditions() const;
  // The setup-event object: checked residuals and active reports only.
  Json::Object json() const;
};
// A checked residual needs a criterion. Overrides are validated one at a time,
// so this pairing is checked after they are applied, before any solve.
std::vector<fluent::SettingIssue> criteria_issues(const ConvergenceControls&);
struct ConvergenceSample {
  bool met=false;
  // The conditions would be met if the next iteration's sample confirms the
  // settled forces: sample it, off the monitor stride.
  bool confirm=false;
  // Residuals stalled above their criteria were waived (forces settled for
  // twice the window).
  bool residuals_waived=false;
  // Iteration-event object: {"met":...,"conditions":[...]}.
  Json::Object json;
  // Met conditions as "residuals", "report:<name>" or "forces".
  std::vector<std::string> met_by;
  // Not converged, but the forces cycle about a steady mean (AverageCriterion).
  bool statistically_steady=false;
  // The averaging state at this sample, when averaging has started: window,
  // iteration span and per-report mean/min/max/drift/range_change in the report's
  // units. Present on every sample once averaging has collected two windows.
  std::optional<Json::Object> averaged;
};
// Evaluates one monitor sample at a time; keeps only the report history the
// lag residual needs.
class ConvergenceMonitor {
 public:
  explicit ConvergenceMonitor(ConvergenceControls);
  const ConvergenceControls& controls() const { return controls_; }
  // residuals by equation; reports by exact name, projected total force in N.
  // The forces-settled condition reads every force report present in `reports`.
  ConvergenceSample sample(int iteration,const std::map<std::string,double>& residuals,
                           const std::map<std::string,double>& reports);
  // Whether this iteration's forces are wanted for averaging. The caller then
  // passes them to record() (or to sample(), which records them itself).
  bool wants_forces(int iteration) const;
  void record(int iteration,const std::map<std::string,double>& reports);
 private:
  struct Forces {int iteration;std::vector<double> values;};
  // Each force report's relative range over the last `window` iterations;
  // null until the samples cover the whole window.
  std::vector<std::optional<double>> ranges(int iteration,int window) const;
  ConvergenceControls controls_;
  std::vector<std::vector<double>> history_;
  std::vector<std::string> force_reports_;
  std::vector<Forces> forces_;
  int last_iteration_=-1;
  bool last_settled_=false,last_settled_twice_=false;
  // Per-iteration force samples for averaging, newest last, at most two
  // windows: (iteration, values in force_reports_ order).
  std::deque<std::pair<int,std::vector<double>>> averaging_;
  std::optional<Json::Object> averaged(int iteration,bool& steady) const;
};
// Fluent's lag residual over the last previous_values differences, relative to
// the newest value. Null with too few values or a zero newest value.
std::optional<double> lag_residual(std::span<const double> history,int previous_values);
} // namespace kutta::setup
