#pragma once

#include "kutta/setup/fluent.hpp"

namespace kutta::app {
setup::CellFields read_restart(const std::filesystem::path&, std::string_view mesh_hash,
                               Index cells, const fv::SolverSettings&);
void write_restart(const std::filesystem::path&, std::string_view mesh_hash,
                   Device&, const fv::SteadySolver&, const fv::SolverSettings&);
}
