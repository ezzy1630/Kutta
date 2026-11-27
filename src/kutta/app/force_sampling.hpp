#pragma once

#include "kutta/setup/fluent.hpp"
#include <algorithm>
#include <map>

namespace kutta::app {
class ForceSampling {
public:
  explicit ForceSampling(const std::vector<setup::ForceDefinition> &reports) {
    for (const auto &report : reports) {
      auto patches = report.patches;
      std::sort(patches.begin(), patches.end());
      const auto group =
          groups_.try_emplace(std::move(patches), groups_.size()).first->second;
      by_report_.push_back(group);
    }
    forces_.resize(groups_.size());
  }

  void collect(fv::SteadySolver &solver, double pressure) {
    for (const auto &[patches, group] : groups_)
      forces_[group] = solver.forces(patches, pressure);
  }

  const fv::ForceReport &operator[](std::size_t report) const {
    return forces_[by_report_[report]];
  }

private:
  std::map<std::vector<std::string>, std::size_t> groups_;
  std::vector<std::size_t> by_report_;
  std::vector<fv::ForceReport> forces_;
};
} // namespace kutta::app
