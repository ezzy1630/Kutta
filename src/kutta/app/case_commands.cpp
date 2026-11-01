#include "case_commands.hpp"
#include "protocol.hpp"
#include "kutta/core/json.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/core/path.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/setup/kutta_case.hpp"
#include "kutta/setup/settings.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>

namespace kutta::app {
namespace {
using setup::RejectedSettings;
struct Args {
  std::vector<std::string> positional;
  std::map<std::string, std::vector<std::string>> options;
  std::set<std::string> flags;
};
// --name VALUE options and bare --flags; options may repeat only where allowed.
Args parse_args(int argc, char** argv, const std::set<std::string>& valued, const std::set<std::string>& switches,
                const std::set<std::string>& repeatable) {
  Args a;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (switches.contains(arg)) {
      a.flags.insert(arg);
      continue;
    }
    if (valued.contains(arg)) {
      KUTTA_REQUIRE(i + 1 < argc, arg + " needs a value");
      KUTTA_REQUIRE(repeatable.contains(arg) || !a.options.contains(arg), "duplicate option: " + arg);
      a.options[arg].push_back(argv[++i]);
      continue;
    }
    KUTTA_REQUIRE(!arg.starts_with("--"), "unknown option: " + arg);
    a.positional.push_back(arg);
  }
  return a;
}
std::string option(const Args& a, const std::string& key) {
  const auto it = a.options.find(key);
  return it == a.options.end() ? std::string{} : it->second.front();
}
// The mesh path as the case file should store it: relative to the case's
// directory, so a case folder can move, unless the only common ancestor is
// the filesystem root.
std::filesystem::path relative_mesh(const std::filesystem::path& mesh, const std::filesystem::path& out) {
  const auto base = std::filesystem::absolute(out).parent_path().lexically_normal();
  const auto abs = std::filesystem::absolute(mesh).lexically_normal();
  const auto rel = abs.lexically_relative(base);
  const auto up = std::count(rel.begin(), rel.end(), std::filesystem::path(".."));
  const auto below_root = base.relative_path();
  const auto depth = std::distance(below_root.begin(), below_root.end());
  return rel.empty() || up >= depth ? abs : rel;
}
void write_text(const std::filesystem::path& path, const std::string& text) {
  // Write beside the target and rename, so a failed write never truncates a case.
  const auto tmp = path.string() + ".tmp";
  {
    std::ofstream f(tmp, std::ios::binary);
    KUTTA_REQUIRE(f, "cannot write " + tmp);
    f << text;
    f.close();
    KUTTA_REQUIRE(f, "cannot write " + tmp);
  }
  std::filesystem::rename(tmp, path);
}
Json::Array issues_json(std::span<const fluent::SettingIssue> issues) {
  Json::Array out;
  for (const auto& i : issues)
    out.push_back(Json::Object{{"path", i.path}, {"reason", i.reason},
                               {"class", i.classification == fluent::IssueClass::unsupported ? "unsupported" : "informational"}});
  return out;
}
bool fluent_mesh(const std::filesystem::path& p) {
  auto n = p.filename().string();
  std::transform(n.begin(), n.end(), n.begin(), [](unsigned char c) { return char(std::tolower(c)); });
  return n.ends_with(".msh") || n.ends_with(".msh.h5") || n.ends_with(".cas") || n.ends_with(".cas.h5");
}
int fail(bool json, const std::string& command, const std::exception& e) {
  if (json) {
    std::vector<Issue> issues;
    if (const auto* r = dynamic_cast<const RejectedSettings*>(&e))
      for (const auto& i : r->issues()) issues.push_back({i.path, i.path, "unsupported", {}, i.reason, {}});
    error_event(e.what(), issues);
  }
  std::cerr << "kutta " << command << ": " << e.what() << '\n';
  return 1;
}
}  // namespace

int new_case(int argc, char** argv) {
  bool json = false;
  for (int i = 2; i < argc; ++i) json = json || std::string_view(argv[i]) == "--json";
  constexpr auto usage = "Usage: kutta new-case MESH [--template fsae-external] [--speed M/S] [--role ZONE=ROLE]... [--out CASE.kutta.json [--force]] [--json]";
  Args a;
  try {
    a = parse_args(argc, argv, {"--template", "--speed", "--role", "--out"}, {"--json", "--force"}, {"--role"});
    KUTTA_REQUIRE(a.positional.size() == 1, usage);
  } catch (const std::exception& e) {
    if (json) error_event(e.what());
    std::cerr << e.what() << '\n' << usage << '\n';
    return 2;
  }
  try {
    const std::filesystem::path mesh_path = a.positional[0];
    if (const auto t = option(a, "--template"); !t.empty() && t != "fsae-external")
      throw RejectedSettings({{"template", "unknown template \"" + t + "\"; available: fsae-external"}}, "Invalid options:");
    if (!fluent_mesh(mesh_path))
      throw RejectedSettings({{"mesh", "a Kutta case reads Fluent meshes (.msh, .msh.h5, or the mesh of a .cas/.cas.h5)"}}, "Invalid options:");
    setup::StarterOptions options;
    if (const auto s = option(a, "--speed"); !s.empty()) {
      const auto v = parse_double(s);
      if (!v || !(*v > 0) || !std::isfinite(*v)) throw RejectedSettings({{"speed", "--speed must be a positive number of m/s"}}, "Invalid options:");
      options.speed = *v;
    }
    if (a.options.contains("--role"))
      for (const auto& r : a.options.at("--role")) {
        const auto eq = r.rfind('=');
        const auto role = eq == std::string::npos ? std::nullopt : setup::parse_role(r.substr(eq + 1));
        if (!role)
          throw RejectedSettings({{"role", "--role " + r + ": expected ZONE=ROLE with ROLE one of velocity-inlet, pressure-outlet, "
                                                           "moving-ground, wall, symmetry, interior"}},
                                 "Invalid options:");
        options.roles[r.substr(0, eq)] = *role;
      }
    const std::filesystem::path out = option(a, "--out");
    if (!out.empty()) {
      KUTTA_REQUIRE(setup::is_kutta_case(out), "--out must name a .json file (conventionally CASE.kutta.json)");
      KUTTA_REQUIRE(a.flags.contains("--force") || !std::filesystem::exists(out), out.string() + " exists; pass --force to replace it");
    }
    const auto mesh = read_mesh(mesh_path);
    const auto info = setup::mesh_info(mesh);
    const auto starter = setup::fsae_starter(info, out.empty() ? mesh_path : relative_mesh(mesh_path, out), options);
    // The starter must lower on its own mesh before it is offered.
    (void)setup::lower_kutta_case(starter, info);
    const auto text = setup::kutta_case_text(starter);
    std::vector<std::string> notes;
    notes.push_back("Reference area is Fluent's default, 1 m^2; set your car's frontal area (reference_values/area) for meaningful coefficients.");
    if (info.dimension == 3 && info.polyhedral)
      notes.push_back("The mesh has polyhedral cells, so Fluent's default turns warped-face gradient correction on; Kutta runs it only "
                      "approximately, with WFGC off, after you approve --approximate wfgc.");
    bool inlet = false, outlet = false, wheels = false;
    Json::Array zones;
    for (const auto& z : info.zones) {
      if (z.kind == setup::MeshZoneInfo::Kind::cell) continue;
      const auto s = setup::suggest_role(z);
      const auto it = options.roles.find(z.name);
      const auto role = it != options.roles.end() ? it->second : s.role;
      inlet = inlet || role == setup::ZoneRole::velocity_inlet;
      outlet = outlet || role == setup::ZoneRole::pressure_outlet;
      wheels = wheels || (role == setup::ZoneRole::wall && s.group == "Wheels");
      Json::Object row{{"id", double(z.id)}, {"name", z.name}, {"mesh_type", z.type},
                       {"kind", z.kind == setup::MeshZoneInfo::Kind::interior ? "interior" : "boundary"}, {"faces", double(z.count)},
                       {"role", std::string(setup::role_name(role))}, {"suggested_role", std::string(setup::role_name(s.role))},
                       {"reason", s.reason}, {"group", s.group ? Json(*s.group) : Json(nullptr)}};
      zones.emplace_back(std::move(row));
    }
    if (!inlet) notes.push_back("No zone is a velocity inlet; assign the inflow boundary its role.");
    if (!outlet) notes.push_back("No zone is a pressure outlet; assign the outflow boundary its role.");
    if (wheels) notes.push_back("Wheels are stationary walls in this starter; make them rotating walls (axis, origin, angular speed) if the car's wheels turn.");
    if (!out.empty()) write_text(out, text);
    if (json) {
      Json::Array roles;
      for (int r = 0; r <= int(setup::ZoneRole::interior); ++r) roles.push_back(std::string(setup::role_name(setup::ZoneRole(r))));
      Json::Array note_list(notes.begin(), notes.end());
      std::cout << dump_json(Json::Object{
                       {"v", 1}, {"template", "fsae-external"}, {"speed", options.speed},
                       {"mesh", Json::Object{{"path", utf8_path(std::filesystem::absolute(mesh_path))}, {"dimension", double(info.dimension)},
                                             {"cells", double(info.cells)}, {"faces", double(info.faces)}, {"nodes", double(info.nodes)},
                                             {"polyhedral", info.polyhedral}}},
                       {"zones", std::move(zones)}, {"roles", std::move(roles)}, {"notes", std::move(note_list)},
                       {"case", parse_json(text)}, {"out", out.empty() ? Json(nullptr) : Json(utf8_path(std::filesystem::absolute(out)))}})
                << std::endl;
    } else if (out.empty())
      std::cout << text;
    else {
      std::cout << "Wrote " << out.string() << '\n';
      for (const auto& n : notes) std::cout << "Note: " << n << '\n';
    }
    return 0;
  } catch (const std::exception& e) {
    return fail(json, "new-case", e);
  }
}

int convert(int argc, char** argv) {
  bool json = false;
  for (int i = 2; i < argc; ++i) json = json || std::string_view(argv[i]) == "--json";
  constexpr auto usage = "Usage: kutta convert CASE [--setup overrides.json] --out CASE.kutta.json [--force] [--json]";
  Args a;
  try {
    a = parse_args(argc, argv, {"--setup", "--out"}, {"--json", "--force"}, {});
    KUTTA_REQUIRE(a.positional.size() == 1 && a.options.contains("--out"), usage);
  } catch (const std::exception& e) {
    if (json) error_event(e.what());
    std::cerr << e.what() << '\n' << usage << '\n';
    return 2;
  }
  try {
    const std::filesystem::path input = a.positional[0], out = option(a, "--out");
    KUTTA_REQUIRE(setup::is_kutta_case(out), "--out must name a .json file (conventionally CASE.kutta.json)");
    const bool same = std::filesystem::exists(out) && std::filesystem::equivalent(input, out);
    KUTTA_REQUIRE(a.flags.contains("--force") || !std::filesystem::exists(out), out.string() + " exists; pass --force to replace it");
    std::optional<setup::LoadedKuttaCase> loaded;
    std::filesystem::path mesh_file = input;
    if (setup::is_kutta_case(input)) {
      loaded = setup::load_kutta_case(input);
      mesh_file = loaded->mesh_path;
    }
    auto c = loaded ? std::move(loaded->settings) : fluent::read_case(input);
    const auto mesh = loaded ? std::move(loaded->mesh) : read_mesh(input);
    Json::Object applied;
    if (const auto setup_file = option(a, "--setup"); !setup_file.empty()) {
      setup::RunControls controls;
      auto inventory = setup::settings_inventory(c, controls);
      applied = setup::apply_setup(setup_file, inventory);
    }
    auto converted = setup::kutta_case_from_fluent(c, setup::mesh_info(mesh), relative_mesh(mesh_file, out));
    std::vector<fluent::SettingIssue> blocking;
    for (const auto& i : converted.omitted)
      if (i.classification == fluent::IssueClass::unsupported) blocking.push_back(i);
    if (!blocking.empty()) {
      // Nothing is written: the result would solve a different problem.
      throw RejectedSettings(std::move(blocking), "Not expressible in a Kutta case (nothing written):");
    }
    // The written case must lower on its mesh to the settings it came from.
    (void)setup::lower_kutta_case(converted.kutta, setup::mesh_info(mesh));
    write_text(out, setup::kutta_case_text(converted.kutta));
    Json::Array run_controls;
    for (const auto& [path, v] : applied)
      if (setup::run_control(path)) run_controls.push_back(path);
    if (json)
      std::cout << dump_json(Json::Object{{"v", 1}, {"out", utf8_path(std::filesystem::absolute(out))}, {"replaced_input", same},
                                          {"omitted", issues_json(converted.omitted)}, {"run_controls_not_saved", std::move(run_controls)}})
                << std::endl;
    else {
      for (const auto& i : converted.omitted) std::cout << "Note: " << i.path << ": " << i.reason << '\n';
      for (const auto& p : run_controls) std::cout << "Note: " << std::get<std::string>(p.value) << " is a run control, not a case setting; not saved\n";
      std::cout << "Wrote " << out.string() << '\n';
    }
    return 0;
  } catch (const std::exception& e) {
    return fail(json, "convert", e);
  }
}
}  // namespace kutta::app
