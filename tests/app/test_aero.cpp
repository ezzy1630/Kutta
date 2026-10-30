#include "../fluent_examples.hpp"
#include "../io/test_support.hpp"
#include <cstdlib>
#include <doctest/doctest.h>
namespace {
std::string quote_aero(const std::filesystem::path& path) {
#ifdef _WIN32
  return '"'+std::filesystem::path(path).make_preferred().string()+'"';
#else
  std::string out="'";for(char c:path.string())out+=c=='\''?"'\\''":std::string(1,c);return out+"'";
#endif
}
}
TEST_CASE("Mixing elbow solve stops on its convergence condition and kutta aero reads its loads") {
#if KUTTA_HAS_HDF5
  const auto fixture=fluent_examples()/"mixing_elbow";
  if(!std::filesystem::exists(fixture/"mixing_elbow.cas.h5")){MESSAGE("SKIP: optional Fluent mixing elbow absent");return;}
  MeshTemp t;
  auto command=quote_aero(KUTTA_PYTHON_PATH)+" "+quote_aero(std::filesystem::path(KUTTA_SOURCE_DIR)/"tests/app/check_aero.py")+" "+quote_aero(KUTTA_CLI_PATH)+" "+quote_aero(fixture)+" "+quote_aero(t.directory);
#ifdef _WIN32
  command='"'+command+'"';
#endif
  CHECK(std::system(command.c_str())==0);
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
