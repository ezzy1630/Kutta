#pragma once
#include <cstdint>

namespace kutta::app {
// Iteration cadence is a minimum gap. Reserve the previous snapshot cost plus
// a 25% margin before another write. The first write calibrates that estimate.
class LiveCadence {
public:
  explicit LiveCadence(int every) : every_(every), next_(every) {}
  bool due(int iteration, double solve_seconds) const {
    return iteration >= next_ &&
           (count_ == 0 ||
            total_ + 1.25 * last_ <= budget_fraction * solve_seconds);
  }
  void record(int iteration, double seconds) {
    ++count_;
    total_ += seconds;
    last_ = seconds;
    next_ = std::int64_t(iteration) + every_;
  }
  int every() const { return every_; }
  std::uint64_t count() const { return count_; }
  double total_seconds() const { return total_; }
  static constexpr double budget_fraction = .05;
  static constexpr int default_every = 20;

private:
  int every_;
  std::int64_t next_;
  std::uint64_t count_ = 0;
  double total_ = 0, last_ = 0;
};
} // namespace kutta::app
