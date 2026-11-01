#pragma once
namespace kutta::app {
// kutta new-case MESH: a Kutta case from a mesh, roles suggested from zone names.
int new_case(int argc, char** argv);
// kutta convert CASE: a Fluent or Kutta case, with optional setup overrides,
// written as a Kutta case.
int convert(int argc, char** argv);
}  // namespace kutta::app
