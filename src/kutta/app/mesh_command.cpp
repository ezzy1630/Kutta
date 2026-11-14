#include "mesh_command.hpp"
#include "kutta/mesh/mesher/external_aero.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/core/memory.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/core/path.hpp"
#include "kutta/core/sha256.hpp"
#include "kutta/setup/kutta_case.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>

namespace kutta::app {
namespace {
void write_text(const std::filesystem::path& p,const std::string& text) {
  std::ofstream out(p,std::ios::binary); out << text << '\n'; out.close(); KUTTA_REQUIRE(out,"mesh: cannot write "+p.string());
}
struct Scratch {
  std::filesystem::path path;
  ~Scratch() { std::error_code ec; if(!path.empty()) std::filesystem::remove_all(path,ec); }
};
}
int mesh_command(int argc,char** argv) {
  using namespace mesher;
  bool json=false,progress=false;
  for(int i=2;i<argc;++i) { json=json || std::string_view(argv[i])=="--json"; progress=progress || std::string_view(argv[i])=="--progress"; }
  try {
    std::vector<std::filesystem::path> inputs;
    std::map<std::string,std::string> options; std::set<std::string> flags; std::vector<std::string> wheels;
    for(int i=2;i<argc;++i) {
      const std::string arg=argv[i];
      if(arg=="--case" || arg=="--force" || arg=="--json") { KUTTA_REQUIRE(flags.insert(arg).second,"duplicate option: "+arg); }
      else if(arg=="--setup" || arg=="--out" || arg=="--progress" || arg=="--speed" || arg=="--wheel") {
        KUTTA_REQUIRE(i+1<argc && !std::string_view(argv[i+1]).starts_with("--"),arg+" needs a value");
        const std::string value=argv[++i]; if(arg=="--wheel") wheels.push_back(value);
        else KUTTA_REQUIRE(options.emplace(arg,value).second,"duplicate option: "+arg);
      } else { KUTTA_REQUIRE(!arg.starts_with("--"),"unknown option: "+arg); inputs.emplace_back(arg); }
    }
    KUTTA_REQUIRE(!inputs.empty() && options.contains("--out"),"Usage: kutta mesh GEOM.stl [MORE.stl ...] [--setup mesh.json] --out DIR [--case] [--speed M/S] [--wheel ZONE=R,X,Y,Z,AX,AY,AZ] [--progress jsonl] [--json] [--force]");
    KUTTA_REQUIRE(!options.contains("--progress") || options.at("--progress")=="jsonl","--progress expects jsonl");
    KUTTA_REQUIRE(wheels.empty() || flags.contains("--case"),"--wheel requires --case");
    Settings settings;
    if(options.contains("--setup")) { KUTTA_REQUIRE(std::filesystem::file_size(options.at("--setup"))<=1024*1024,"mesh: setup exceeds 1 MiB"); std::ifstream f(options.at("--setup"),std::ios::binary); KUTTA_REQUIRE(f,"cannot read mesh setup"); const std::string text((std::istreambuf_iterator<char>(f)),{}); settings=parse_settings(text); }
    double speed=20;
    if(options.contains("--speed")) { const auto value=parse_double(options.at("--speed")); KUTTA_REQUIRE(value && *value>0 && std::isfinite(*value),"--speed must be a positive finite speed in m/s"); speed=*value; }
    const auto out=std::filesystem::absolute(options.at("--out")).lexically_normal();
    const bool force=flags.contains("--force");
    for(const auto* name:{"mesh.msh","mesh-report.json","case.kutta.json"}) {
      const auto target=out/name;
      if(std::string_view(name)=="case.kutta.json" && !flags.contains("--case")) continue;
      KUTTA_REQUIRE(force || !std::filesystem::exists(target),target.filename().string()+" exists; pass --force to replace it");
      for(const auto& input:inputs) KUTTA_REQUIRE(!std::filesystem::exists(target) || !std::filesystem::equivalent(input,target),"mesh output must not overwrite an input");
    }
    const auto start=std::chrono::steady_clock::now();
    auto emit=[&](std::string_view stage,double fraction,std::string_view message) {
      if(progress) std::cout << dump_json(Json::Object{{"v",1},{"event","mesh-progress"},{"stage",std::string(stage)},{"fraction",fraction},{"message",std::string(message)}}) << std::endl;
      else if(!json) std::cerr << "mesh " << stage << ": " << message << '\n';
    };
    emit("read",0,"Reading STL geometry"); auto surface=read_stl(inputs,settings);
    emit("repair",.1,"Welded vertices, removed duplicate/degenerate facets, oriented connected shells");
    auto result=generate(surface,settings,speed,emit);
    std::string case_text;
    if(flags.contains("--case")) {
      setup::StarterOptions starter; starter.speed=speed;
      for(const auto& p:result.mesh.patches) starter.roles[p.name]=p.name=="inlet"?setup::ZoneRole::velocity_inlet:p.name=="outlet"?setup::ZoneRole::pressure_outlet:p.name=="ground"?setup::ZoneRole::moving_ground:p.kind==PatchKind::Symmetry?setup::ZoneRole::symmetry:setup::ZoneRole::wall;
      for(const auto& text:wheels) {
        const auto eq=text.rfind('='); KUTTA_REQUIRE(eq!=std::string::npos && eq>0,"--wheel expects ZONE=RADIUS,X,Y,Z,AX,AY,AZ");
        std::istringstream in(text.substr(eq+1)); std::array<double,7> values;
        for(auto& v:values) { std::string part; KUTTA_REQUIRE(bool(std::getline(in,part,',')),"--wheel needs seven values"); const auto d=parse_double(part); KUTTA_REQUIRE(d && std::isfinite(*d),"invalid --wheel number"); v=*d; }
        KUTTA_REQUIRE(in.eof(),"--wheel needs exactly seven values");
        KUTTA_REQUIRE(starter.wheels.emplace(text.substr(0,eq),setup::StarterOptions::Wheel{values[0],{values[1],values[2],values[3]},{values[4],values[5],values[6]}}).second,"duplicate --wheel zone");
      }
      const auto info=setup::mesh_info(result.mesh);
      auto c=setup::fsae_starter(info,"mesh.msh",starter);
      // Snapped faces can be warped, and Kutta does not implement Fluent's
      // warped-face gradient correction: state the scheme it runs rather
      // than inheriting Fluent's default for polyhedral imports.
      c.methods.warped_face_gradient_correction=false;
      (void)setup::lower_kutta_case(c,info); case_text=setup::kutta_case_text(c);
    }
    emit("write",.95,"Writing Fluent ASCII mesh and report");
    std::filesystem::create_directories(out);
    Scratch scratch{out/(".kutta-mesh-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))};
    KUTTA_REQUIRE(std::filesystem::create_directory(scratch.path),"mesh: cannot create staging directory");
    write_fluent_msh(scratch.path/"mesh.msh",result.mesh,"Kutta external aero mesher v1");
    if(!case_text.empty()) write_text(scratch.path/"case.kutta.json",case_text);
    result.report["wall_seconds"]=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    result.report["peak_memory_bytes"]=memory::usage().peak?Json(double(memory::usage().peak)):Json(nullptr);
    Json::Array input_hashes;
    for(const auto& input:inputs) input_hashes.push_back(Json::Object{{"file",utf8_path(input.filename())},{"sha256",sha256_file(input)}});
    result.report["inputs"]=std::move(input_hashes);
    result.report["kutta_commit"]=KUTTA_COMMIT;
    result.report["out"]=utf8_path(out); result.report["case_written"]=!case_text.empty();
    write_text(scratch.path/"mesh-report.json",dump_json(result.report));
    for(const auto* name:{"mesh.msh","mesh-report.json","case.kutta.json"}) if(std::filesystem::exists(scratch.path/name)) {
      if(force) { std::error_code ec; std::filesystem::remove(out/name,ec); }
      std::filesystem::rename(scratch.path/name,out/name);
    }
    result.report["event"]="mesh-done";
    if(json || progress) std::cout << dump_json(result.report) << std::endl;
    else std::cout << "Wrote " << (out/"mesh.msh").string() << " (" << result.mesh.n_cells << " cells); report: " << (out/"mesh-report.json").string() << '\n';
    return 0;
  } catch(const std::exception& e) {
    const std::string advice="Check the named input/setup path; increase max_cells or mesh sizes for budget errors.";
    if(json || progress) std::cout << dump_json(Json::Object{{"v",1},{"event","error"},{"message",e.what()},{"advice",advice}}) << std::endl;
    std::cerr << "kutta mesh: " << e.what() << '\n'; return 1;
  }
}
}
