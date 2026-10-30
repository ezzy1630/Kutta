#pragma once

#include "kutta/case/fluent_case.hpp"
#include "kutta/case/fluent_data.hpp"
#include "kutta/fv/steady_solver.hpp"

namespace kutta::setup {
struct FluentSetup {
  fv::SolverSettings settings;
  std::vector<fv::PatchCondition> boundaries;
  std::vector<std::string> convergence_notes;
};
// Strict flow adapter. Throws with the source setting path when the discrete
// equations cannot be represented. Coupled -> SIMPLEC is a reported path change.
FluentSetup fluent_setup(const fluent::FluentCase& input,const PolyMesh& mesh);
fv::TurbulenceValues turbulence_values(const fv::TurbulenceCondition& spec,
    double speed,double density,double viscosity);
struct CellFields {
  std::vector<Vec3d> velocity;
  std::vector<double> pressure,k,omega;
};
CellFields fluent_fields(const fluent::FluentCase& input,const fluent::FluentData& data,
                        const PolyMesh& original_mesh);
// Standard initialization only. Hybrid/FMG are explicitly rejected.
CellFields initial_fields(const fluent::FluentCase& input,const PolyMesh& original_mesh);
} // namespace kutta::setup
