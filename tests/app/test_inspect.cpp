#include <chrono>
#include <cstdlib>
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace {
std::string shell_quote(const std::filesystem::path &path) {
#ifdef _WIN32
  return '"' + std::filesystem::path(path).make_preferred().string() + '"';  // paths cannot contain '"'
#else
  std::string result = "'";
  for (char ch : path.string())
    result += ch == '\'' ? "'\\''" : std::string(1, ch);
  return result + "'";
#endif
}
struct InspectFixture {
  std::filesystem::path path = std::filesystem::temp_directory_path() /
      ("kutta-inspect-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".cas");
  std::filesystem::path output = path.string() + ".out";
  ~InspectFixture() {
    std::error_code ec;
    std::filesystem::remove(path, ec);
    std::filesystem::remove(output, ec);
  }
  std::string inspect(const std::filesystem::path &input) {
    auto command = shell_quote(KUTTA_CLI_PATH) + " inspect " + shell_quote(input) + " > " + shell_quote(output);
#ifdef _WIN32
    command = '"' + command + '"';  // cmd /c strips one outer pair of quotes
#endif
    REQUIRE(std::system(command.c_str()) == 0);
    std::ifstream f(output);
    std::ostringstream text;
    text << f.rdbuf();
    return text.str();
  }
  std::string write_and_inspect(bool coupled) {
    std::ofstream f(path);
    f << R"((0 "fluent24.1.0") (2 3) (37 (
      (case-config ((rp-3d? . #t) (rp-seg? . #t) (rp-double? . #t) (rp-unsteady? . #f)
                    (rp-lam? . #f) (rp-ke? . #t) (rp-kw? . #f) (rf-energy? . #t)))
      (rng-ke-on? #f) (realizable-epsilon? #f)
      (materials ((water fluid (density (constant . 998)) (viscosity (constant . 0.001)))))
      (epsilon/scheme 0) (k/scheme 1) (mom/scheme 1) (pressure/scheme 12)
      (recon/cell-lsq? #f) (recon/node-lsq? #t)
      (pressure/relax 0.3) (mom/relax 0.7) (k/relax 0.8) (epsilon/relax 0.8)
      (pressure-coupled/pressure/pseudo-explicit-relax 0.5)
      (pressure-coupled/mom/pseudo-explicit-relax 0.5)
      (k/pseudo-relax 0.75) (epsilon/pseudo-relax 0.75)
      (flow/scheme )" << (coupled ? 24 : 20) << ") (pseudo-transient-formulation? "
      << (coupled ? "#t" : "#f") << R"() (pseudo-time-method/coupled-pbns/dt-method 1)))
      (39 (2 fluid water 1) ((material . water))))";
    f.close();
    return inspect(path);
  }
};
}
TEST_CASE("Inspect labels active relaxation equations and issue classes") {
  InspectFixture f;
  const auto simple = f.write_and_inspect(false);
  CHECK(simple.find("epsilon=FOU") != std::string::npos);
  const auto schemes_start = simple.find("Schemes:");
  REQUIRE(schemes_start != std::string::npos);
  const auto schemes = simple.substr(schemes_start, simple.find('\n', schemes_start) - schemes_start);
  CHECK(schemes.find("omega=") == std::string::npos);
  CHECK(simple.find("Active under-relaxation (segregated URFs):") != std::string::npos);
  CHECK(simple.find("pressure=0.3") != std::string::npos);
  CHECK(simple.find("GG-node (unverified mapping)") != std::string::npos);
  CHECK(simple.find("Issues: unsupported/required=") != std::string::npos);
  CHECK(simple.find("informational=1") != std::string::npos);
  CHECK(simple.find("[informational] case-config/rf-energy?: energy ignored (decoupled)") != std::string::npos);
  const auto coupled = f.write_and_inspect(true);
  CHECK(coupled.find("Active under-relaxation (coupled pseudo-transient):") != std::string::npos);
  CHECK(coupled.find("pressure=0.5") != std::string::npos);
  CHECK(coupled.find("mom=0.5") != std::string::npos);
  CHECK(coupled.find("epsilon=0.75") != std::string::npos);
}
TEST_CASE("Inspect real static mixer displays resolved expressions instead of saved snapshots") {
  auto root = std::getenv("KUTTA_FLUENT_EXAMPLES");
  auto path = std::filesystem::path(root ? root : KUTTA_SOURCE_DIR "/validation/fluent_examples") /
      "static_mixer/Static_Mixer_Parameters_legacy_text.cas";
  if (!std::filesystem::exists(path)) {
    MESSAGE("optional real-file inspect fixture absent: ", path.string());
    return;
  }
  InspectFixture f;
  auto output = f.inspect(path);
  CHECK(output.find("resolved from expression inlet1_vel") != std::string::npos);
  CHECK(output.find("saved=2") == std::string::npos);
}
