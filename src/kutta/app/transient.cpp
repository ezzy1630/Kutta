#include "transient.hpp"
#include "force_sampling.hpp"
#include "kutta/mesh/detail.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace kutta::app {
namespace {
std::ofstream file(const std::filesystem::path &path) {
  std::ofstream stream;
  stream.exceptions(std::ios::badbit | std::ios::failbit);
  stream.open(path);
  stream.imbue(std::locale::classic());
  stream << std::setprecision(17);
  return stream;
}
struct ForceSamples {
  int count = 0;
  double mean = 0, squared_deviations = 0, minimum = 0, maximum = 0;
  void add(double value) {
    if (!count)
      minimum = maximum = value;
    minimum = std::min(minimum, value);
    maximum = std::max(maximum, value);
    const double delta = value - mean;
    mean += delta / ++count;
    squared_deviations += delta * (value - mean);
  }
};

using FieldSums = std::vector<std::vector<double>>;
void accumulate(FieldSums &sums, std::vector<view::Field> &fields, int count,
                bool last) {
  if (sums.empty())
    sums.resize(fields.size());
  for (std::size_t i = 0; i < fields.size(); ++i) {
    auto &values = fields[i].values;
    auto &sum = sums[i];
    if (count == 1)
      sum.assign(values.begin(), values.end());
    else
      for (std::size_t j = 0; j < values.size(); ++j)
        sum[j] += values[j];
    if (last) {
      for (std::size_t j = 0; j < values.size(); ++j)
        values[j] = float(sum[j] / count);
      std::vector<double>().swap(sum);
    }
  }
}
} // namespace
void run_transient(Device &dev, fv::SteadySolver &solver,
                   const fluent::FluentCase &c, const setup::FluentSetup &setup,
                   const std::vector<setup::ForceDefinition> &reports,
                   const PhysicalRun &options, const std::filesystem::path &out,
                   Json::Object metadata, StopRequest &stop) {
  KUTTA_REQUIRE(!std::filesystem::exists(out) || std::filesystem::is_empty(out),
                "transient output directory must be empty");
  std::filesystem::create_directories(out);
  auto history = file(out / "physical-time.csv");
  history << "step,time_s,correctors,max_equation_residual,max_local_predictor_"
             "continuity,max_local_corrected_continuity";
  for (const auto &d : reports)
    if (d.direction) {
      auto label = d.name + " [N]";
      std::string escaped;
      for (char ch : label)
        escaped += ch == '"' ? "\"\"" : std::string(1, ch);
      history << ",\"" << escaped << '"';
    }
  history << '\n';
  const auto &m = solver.device_mesh();
  ForceSampling sampled_forces(reports);
  std::vector<ForceSamples> force_means(reports.size());
  int sampled = 0;
  view::Fields mean;
  FieldSums cell_sums, boundary_sums;
  for (int step = 1; step <= options.steps; ++step) {
    const bool include_in_mean =
        options.mean_from && step >= options.mean_from &&
        (step - options.mean_from) % options.mean_every == 0;
    solver.begin_physical_step(options.dt);
    bool accepted = false;
    double maximum = 0, predictor = 0, corrected = 0;
    for (int i = 1; i <= options.max_correctors; ++i) {
      KUTTA_REQUIRE(!stop.requested(),
                    "physical integration interrupted before completion");
      solver.iterate(1);
      const auto r =
          solver.report().residual_iteration == solver.report().iteration
              ? solver.report()
              : solver.monitor();
      maximum = 0;
      for (auto value :
           {r.u.scaled, r.v.scaled, r.w.scaled, r.k.scaled, r.omega.scaled}) {
        KUTTA_REQUIRE(std::isfinite(value),
                      "nonfinite physical equation residual");
        maximum = std::max(maximum, value);
      }
      predictor = r.physical_predictor_continuity;
      corrected = r.physical_corrected_continuity;
      KUTTA_REQUIRE(std::isfinite(predictor) && std::isfinite(corrected),
                    "nonfinite physical continuity defect");
      if (std::max({maximum, predictor, corrected}) > options.tolerance)
        continue;
      sampled_forces.collect(solver, c.reference.pressure.value_or(0));
      std::ostringstream row;
      row.imbue(std::locale::classic());
      row << std::setprecision(17) << step << ',' << step * options.dt << ','
          << i << ',' << maximum << ',' << predictor << ',' << corrected;
      for (std::size_t j = 0; j < reports.size(); ++j)
        if (reports[j].direction) {
          const auto value =
              detail::dot(sampled_forces[j].force(), *reports[j].direction);
          KUTTA_REQUIRE(std::isfinite(value), "nonfinite transient load");
          if (include_in_mean)
            force_means[j].add(value);
          row << ',' << value;
        }
      history << row.str();
      history << '\n';
      history.flush();
      accepted = true;
      break;
    }
    if (!accepted) {
      const auto worst = solver.worst_time_step_defect();
      const Json::Object detail{
          {"step", step},
          {"time_s", step * options.dt},
          {"max_equation_residual", maximum},
          {"max_correctors", options.max_correctors},
          {"tolerance", options.tolerance},
          {"max_predictor_mass_fraction", predictor},
          {"max_corrected_mass_fraction", corrected},
          {"cell", worst.cell},
          {"cell_volume", worst.volume},
          {"cell_center",
           Json::Array{worst.position.x, worst.position.y, worst.position.z}},
          {"predictor_mass_fraction", worst.predictor_mass_fraction},
          {"corrected_mass_fraction_in_same_cell",
           worst.corrected_mass_fraction}};
      file(out / "failed-step.json") << dump_json(detail) << '\n';
      std::ostringstream message;
      message << "physical step " << step << " did not converge at tolerance "
              << options.tolerance << ": maximum equation residual " << maximum
              << ", predictor continuity " << predictor
              << ", corrected continuity " << corrected;
      throw Error(message.str());
    }
    if (include_in_mean) {
      auto fields = collect_view_fields(dev, solver, c, nullptr, setup,
                                        "Physical-time sample mean", true);
      ++sampled;
      const bool last = options.steps - step < options.mean_every;
      accumulate(cell_sums, fields.cells, sampled, last);
      accumulate(boundary_sums, fields.boundary, sampled, last);
      if (last)
        mean = std::move(fields);
    }
  }
  metadata.emplace("v", 1);
  metadata.emplace("cells", m.n_cells);
  metadata.emplace("backend", backend_name(dev.backend()));
  metadata.emplace("device", dev.name());
  metadata.emplace("dt_s", options.dt);
  metadata.emplace("steps", options.steps);
  metadata.emplace("time_s", options.dt * options.steps);
  metadata.emplace("iterations", solver.report().iteration);
  metadata.emplace("corrector_tolerance", options.tolerance);
  metadata.emplace("max_correctors", options.max_correctors);
  metadata.emplace("scheme", "backward Euler startup, fixed-step BDF2");
  metadata.emplace("verdict", "completed_physical_steps");
  save_view_input(out / "view-input.bin", dev, solver, c, nullptr, setup,
                  "Instantaneous physical-time field");
  if (sampled) {
    std::filesystem::create_directory(out / "mean");
    view::write_snapshot(out / "mean" / "view-input.bin", solver.mesh(), mean);
    auto averaged = metadata;
    averaged.emplace("mean_from_step", options.mean_from);
    averaged.emplace("samples", sampled);
    averaged.emplace("first_sample_time_s", options.mean_from * options.dt);
    averaged.emplace("sample_stride_steps", options.mean_every);
    averaged.emplace("last_sample_time_s",
                     (options.mean_from + (sampled - 1) * options.mean_every) *
                         options.dt);
    Json::Array loads;
    for (std::size_t i = 0; i < reports.size(); ++i)
      if (reports[i].direction) {
        const auto &f = force_means[i];
        loads.emplace_back(Json::Object{
            {"report", reports[i].name},
            {"mean_N", f.mean},
            {"rms_fluctuation_N",
             std::sqrt(std::max(0., f.squared_deviations) / f.count)},
            {"min_N", f.minimum},
            {"max_N", f.maximum}});
      }
    averaged.emplace("loads", std::move(loads));
    averaged.emplace(
        "field_statistic",
        "arithmetic mean of equally spaced accepted physical steps");
    file(out / "mean" / "summary.json") << dump_json(averaged) << '\n';
  }
  file(out / "summary.json") << dump_json(metadata) << '\n';
  std::cout << "Completed " << options.steps << " physical steps at "
            << options.dt << " s; outputs in " << out << '\n';
}
} // namespace kutta::app
