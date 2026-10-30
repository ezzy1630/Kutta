#include <cstdlib>
#include "../backend_helpers.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/setup/fluent.hpp"

TEST_CASE("FSAE V12 native Fluent polynomial fan import acceptance" * doctest::test_suite("slow")) {
  // A synthetic spelling is not evidence of Fluent's native encoding. Keep
  // this acceptance gate red until an actual saved case supplies that proof.
  const auto* fixture=std::getenv("KUTTA_FSAE_FAN_CASE");
  REQUIRE_MESSAGE((fixture && std::filesystem::is_regular_file(fixture)),
    "Native fan input remains unverified: KUTTA_FSAE_FAN_CASE must identify a saved Fluent case with a polynomial fan thread");
  const auto c=kutta::fluent::read_case(fixture);bool fan=false;
  for(const auto& z:c.boundaries) fan=fan || z.type==kutta::fluent::ZoneType::fan;
  REQUIRE_MESSAGE(fan,"Acceptance fixture must contain an actual fan thread, rather than an MRF rotor");
  const auto mesh=kutta::read_mesh(fixture);
  CHECK_NOTHROW(kutta::setup::fluent_setup(c,mesh));
}
