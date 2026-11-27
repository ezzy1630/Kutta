#pragma once
#include "export.hpp"
#include "stop.hpp"

namespace kutta::app {
struct PhysicalRun {
  int steps = 0, max_correctors = 30, mean_from = 0, mean_every = 1;
  double dt = 0, tolerance = 1e-4;
};
void run_transient(Device &, fv::SteadySolver &, const fluent::FluentCase &,
                   const setup::FluentSetup &,
                   const std::vector<setup::ForceDefinition> &,
                   const PhysicalRun &, const std::filesystem::path &,
                   Json::Object metadata, StopRequest &);
} // namespace kutta::app
