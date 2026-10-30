#pragma once

#include "kutta/case/fluent_case.hpp"
#include "kutta/case/fluent_data.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/setup/approximations.hpp"

namespace kutta::setup {
struct FluentSetup {
  fv::SolverSettings settings;
  std::vector<fv::PatchCondition> boundaries;
  std::vector<std::string> convergence_notes;
  std::vector<Approximation> approximations;
};
// Strict by default: collect every unrepresentable setting path and reason.
// Named opt-ins permit only the setup allowlist; all other failures remain.
// Coupled -> SIMPLEC is a reported convergence-path change.
// Only explicitly allowed, applicable approximations are recorded in the result.
FluentSetup fluent_setup(const fluent::FluentCase& input,const PolyMesh& mesh,
                         std::span<const Approximation> allowed = {});
fv::TurbulenceValues turbulence_values(const fv::TurbulenceCondition& spec,
    double speed,double density,double viscosity);
struct CellFields {
  std::vector<Vec3d> velocity;
  std::vector<double> pressure,k,omega;
};
CellFields fluent_fields(const fluent::FluentCase& input,const fluent::FluentData& data,
                        const PolyMesh& original_mesh);
using Initialization=std::variant<CellFields,fv::HybridInitialization>;
// Parse initialization independently of the selected Device. FMG is unsupported.
Initialization initial_fields(const fluent::FluentCase& input,const PolyMesh& original_mesh);
void initialize(fv::SteadySolver& solver,const Initialization& initialization);
} // namespace kutta::setup
