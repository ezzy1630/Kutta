#include <doctest/doctest.h>
#include "../io/test_support.hpp"
#include "kutta/mesh/mesher/external_aero.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/core/json_value.hpp"
#include <bit>
#include <cstdlib>
#include <numeric>
#include <set>
#include <sstream>

namespace {
using namespace kutta;
using namespace kutta::mesher;
using namespace kutta::detail;
std::string cube(std::string name="body",double scale=1,bool flip=false,Vec3d shift={}) {
  const auto m=make_box_hex(1,1,1,Vec3d{.13,.17,.1}+shift,Vec3d{.73,.77,.7}+shift);
  std::ostringstream out; out << "solid " << name << '\n';
  for(Index f=0;f<m.n_faces();++f) {
    const auto ns=nodes(m,f);
    for(int k=1;k<3;++k) {
      out << "facet normal 0 0 0\nouter loop\n";
      auto tri=std::array{ns[0],ns[k],ns[k+1]}; if(flip && f==0 && k==1) std::swap(tri[0],tri[1]);
      for(auto n:tri) { const auto p=m.points[n]*scale; out << "vertex " << p.x << ' ' << p.y << ' ' << p.z << '\n'; }
      out << "endloop\nendfacet\n";
    }
  }
  return out.str()+"endsolid "+name+"\n";
}
Surface surface(MeshTemp& t,std::string text=cube()) {
  mesh_text(t.path("test.stl"),text); const std::array paths{t.path("test.stl")}; return read_stl(paths,{});
}
Settings tiny(std::string layers=R"({"count":0})") { return parse_settings(R"({"v":1,"base_size":0.4,"surface":{"min_size":0.1,"max_size":0.1},"tunnel":{"upstream":0.3,"downstream":0.3,"side":0.3,"top":0.3},"layers":)"+layers+R"(,"max_cells":200000})"); }
double report(const Json::Object& r,std::initializer_list<std::string_view> path) {
  const Json* j=nullptr; const Json::Object* o=&r;
  for(auto key:path) { j=&o->at(std::string(key)); if(const auto* next=std::get_if<Json::Object>(&j->value)) o=next; }
  return std::get<double>(j->value);
}
double volume(const PolyMesh& m) { const auto g=compute_geometry(m); return std::accumulate(g.cell_volume.begin(),g.cell_volume.end(),0.); }
std::string quote(const std::filesystem::path& p) {
#ifdef _WIN32
  return '"'+p.string()+'"';
#else
  std::string out="'"; for(char c:p.string()) out+=c=='\''?"'\\''":std::string(1,c); return out+"'";
#endif
}
int cli(const std::string& args,const std::filesystem::path& log) {
  auto cmd=quote(KUTTA_CLI_PATH)+" "+args+" > "+quote(log)+" 2>&1";
#ifdef _WIN32
  cmd='"'+cmd+'"';
#endif
  return std::system(cmd.c_str());
}
}
TEST_CASE("mesh setup is strict and path-qualified") {
  CHECK(parse_settings(R"({"v":1})").base_size==.25);
  CHECK_THROWS_WITH_AS(parse_settings(R"({"v":1,"surface":{"mix_size":1}})"),"mesh/surface/mix_size: unknown property",Error);
  CHECK_THROWS(parse_settings(R"({"v":1,"surface":{"min_size":0.5,"max_size":0.1}})"));
  CHECK_THROWS(parse_settings(R"({"v":1,"max_cells":4.2})"));
  CHECK_THROWS(parse_settings(R"({"v":1,"v":1})"));
  CHECK_THROWS(parse_settings(R"({"v":1,"tunnel":{"symmetry":"x"}})"));
  CHECK_THROWS(parse_settings(R"({"v":1,"refinement_boxes":[{"name":"wake","min":[0,0,0],"max":[1,0,1],"size":0.1}]})"));
  CHECK_THROWS_WITH_AS(parse_settings(R"({"v":1,"surface":{"min_size":0.01,"proximity_min_size":0.02}})"),"mesh/surface/proximity_min_size: must not exceed min_size",Error);
  const auto defaults=parse_settings(R"({"v":1})");
  CHECK(defaults.layer_count==4); CHECK(defaults.proximity_floor()==doctest::Approx(defaults.min_size/4));
  CHECK(defaults.last_layer()==doctest::Approx(.002*1.3*1.3*1.3));
}
TEST_CASE("STL ASCII multi-solid repair and inside classification ignore supplied normals") {
  MeshTemp t;
  auto text=cube("front-wing",1,true);
  // Repeat a facet after reversing its winding; the duplicate is removed.
  const auto start=text.find("facet"),end=text.find("endfacet")+8;
  text.insert(text.find("endsolid"),text.substr(start,end-start)+"\n");
  const auto sf=surface(t,text+cube("rear-wing"));
  CHECK(sf.zones==std::vector<std::string>{"front-wing","rear-wing"});
  CHECK(sf.triangles.size()==24);
  CHECK(sf.repairs.at("removed_duplicate_triangles")==Json{1});
  CHECK(sf.repairs.at("open_edges")==Json{0});
  CHECK(std::get<double>(sf.repairs.at("reoriented_triangles").value)>0);
  SurfaceQuery q(sf);
  CHECK(q.sample({.4,.4,.4}).inside);
  CHECK_FALSE(q.sample({-.2,.4,.4}).inside);
  CHECK(q.sample({-.2,.4,.4}).distance==doctest::Approx(.33));
  CHECK_FALSE(q.sample({.4,.4,1.2}).inside);
  auto units=Settings{}; units.units="mm"; mesh_text(t.path("mm.stl"),cube("body",1000)); const std::array paths{t.path("mm.stl")};
  const auto mm=read_stl(paths,units); CHECK(mm.bounds.max.x==doctest::Approx(.73));
  CHECK_THROWS(surface(t,"solid broken\nfacet normal 0 0 1\n"));
}
TEST_CASE("STL small leaks are diagnosed without a flood fill") {
  MeshTemp t; auto text=cube();
  const auto first=text.find("facet"),last=text.find("endfacet")+8;
  text.erase(first,last-first);
  const auto sf=surface(t,text);
  CHECK(std::get<double>(sf.repairs.at("open_edges").value)==3);
  SurfaceQuery query(sf);
  CHECK(query.sample({.4,.4,.4}).inside);
  CHECK_FALSE(query.sample({2,2,2}).inside);
}
TEST_CASE("STL binary header may start solid and truncated binary fails") {
  MeshTemp t; std::string bytes(84,'\0'); bytes.replace(0,5,"solid"); bytes[80]=1;
  auto add_float=[&](float f) { const auto bits=std::bit_cast<std::uint32_t>(f); for(int i=0;i<4;++i) bytes+=char((bits>>(8*i))&255); };
  for(float f:{0.f,0.f,1.f,0.f,0.f,0.f,1.f,0.f,0.f,0.f,1.f,0.f}) add_float(f);
  bytes.append(2,'\0'); mesh_text(t.path("binary.stl"),bytes); const std::array paths{t.path("binary.stl")};
  CHECK(read_stl(paths,{}).triangles.size()==1);
  bytes.pop_back(); mesh_text(t.path("binary.stl"),bytes); CHECK_THROWS(read_stl(paths,{}));
}
TEST_CASE("castellated snapped mesh conserves the body, captures corners and round trips") {
  MeshTemp t; const auto settings=tiny();
  auto tunnel_volume=[](const Result& r) {
    const auto bounds=std::get<Json::Array>(r.report.at("bounds").value);
    const auto lo=std::get<Json::Array>(bounds[0].value),hi=std::get<Json::Array>(bounds[1].value);
    double v=1; for(int a=0;a<3;++a) v*=std::get<double>(hi[a].value)-std::get<double>(lo[a].value); return v;
  };
  auto corner_error=[](const PolyMesh& m,Vec3d shift) {
    double worst=0;
    for(int c=0;c<8;++c) {
      const auto corner=Vec3d{c&1?.73:.13,c&2?.77:.17,c&4?.7:.1}+shift; double nearest=1;
      for(auto p:m.points) nearest=std::min(nearest,norm(p-corner));
      worst=std::max(worst,nearest);
    }
    return worst;
  };
  // On and off the octree lattice every three-plane corner and face is
  // exact; the off-lattice cube needs merges kept apart at one corner.
  const Vec3d shift{.011,.007,.013};
  const auto offset=generate(surface(t,cube("body",1,false,shift)),settings);
  CHECK_NOTHROW(validate(offset.mesh));
  CHECK(volume(offset.mesh)==doctest::Approx(tunnel_volume(offset)-.216).epsilon(1e-9));
  CHECK(corner_error(offset.mesh,shift)<1e-9);
  CHECK(report(offset.report,{"snap","point_distance_m","max"})<1e-9);
  CHECK(report(offset.report,{"snap","stuck_cells_removed"})==0);
  const auto sf=surface(t); const auto result=generate(sf,settings);
  CHECK(result.mesh.n_cells>0); CHECK_NOTHROW(validate(result.mesh));
  CHECK(volume(result.mesh)==doctest::Approx(tunnel_volume(result)-.216).epsilon(1e-9));
  CHECK(corner_error(result.mesh,{})<1e-9);
  CHECK(report(result.report,{"snap","point_distance_m","max"})<1e-9);
  const auto& r=std::get<Json::Object>(result.report.at("quality").value);
  CHECK(report(r,{"orthogonal_min"})>=.15); CHECK(report(r,{"cells_below_0_15"})==0);
  CHECK(report(result.report,{"cell_types","hex_fraction_of_core"})==1);
  write_msh(t.path("test.msh"),result.mesh);
  const auto read=read_mesh(t.path("test.msh")); CHECK_NOTHROW(validate(read));
  CHECK(read.n_cells==result.mesh.n_cells); CHECK(read.n_faces()==result.mesh.n_faces()); CHECK(read.n_points()==result.mesh.n_points());
  REQUIRE(read.patches.size()==result.mesh.patches.size());
  for(std::size_t i=0;i<read.patches.size();++i) { CHECK(read.patches[i].name==result.mesh.patches[i].name); CHECK(read.patches[i].size==result.mesh.patches[i].size); }
  CHECK(volume(read)==doctest::Approx(volume(result.mesh)).epsilon(1e-12));
  auto refined=settings; refined.boxes.push_back({"wake",{{.7,.1,.1},{1.,.8,.7}},.05}); CHECK_NOTHROW(validate(generate(sf,refined).mesh));
  auto limited=settings; limited.max_cells=10; CHECK_THROWS(generate(sf,limited));
  auto symmetry=settings; symmetry.symmetry_y=true; CHECK_NOTHROW(validate(generate(sf,symmetry).mesh));
}
TEST_CASE("prism layers extrude every visible wall face and keep the body volume") {
  MeshTemp t; const auto sf=surface(t);
  const auto plain=generate(sf,tiny());
  const auto result=generate(sf,tiny(R"({"count":3,"first_height":0.004,"growth":1.2})"));
  CHECK_NOTHROW(validate(result.mesh));
  const auto& layers=std::get<Json::Object>(result.report.at("layers").value);
  CHECK(report(layers,{"generated"})==3);
  CHECK(report(layers,{"cells"})==3*report(layers,{"stacks"}));
  CHECK(report(result.report,{"cell_types","prism"})==report(layers,{"cells"}));
  CHECK(report(layers,{"coverage","full"})>.99);
  CHECK(report(layers,{"first_height_achieved","mean"})==doctest::Approx(.004).epsilon(.01));
  CHECK(report(result.report,{"quality","orthogonal_min"})>=.15);
  // Layers fill the gap they push open: the fluid volume is unchanged.
  CHECK(volume(result.mesh)==doctest::Approx(volume(plain.mesh)).epsilon(1e-9));
  auto zones=tiny(R"({"count":2,"first_height":0.004,"growth":1.2,"zones":["missing"]})");
  CHECK_THROWS_WITH_AS(generate(sf,zones),doctest::Contains("mesh/layers/zones: unknown STL zone missing"),Error);
}
TEST_CASE("gap sizing resolves a slot between shells and reports cells across it") {
  MeshTemp t;
  // Two cubes 0.06 m apart in x; the slot needs proximity refinement.
  auto second=cube("rear");
  std::string shifted; std::istringstream in(second); std::string line;
  while(std::getline(in,line)) {
    if(line.rfind("vertex",0)==0) { std::istringstream v(line.substr(7)); double x,y,z; v>>x>>y>>z; std::ostringstream o; o<<"vertex "<<x+.66<<' '<<y<<' '<<z; line=o.str(); }
    shifted+=line+"\n";
  }
  const auto sf=surface(t,cube("front")+shifted);
  auto settings=tiny();
  settings.proximity_cells=2; settings.proximity_min_size=.0125;
  const auto result=generate(sf,settings);
  CHECK_NOTHROW(validate(result.mesh));
  const auto gaps=std::get<Json::Array>(std::get<Json::Object>(result.report.at("castellation").value).at("gaps").value);
  bool found=false;
  for(const auto& g:gaps) {
    const auto& o=std::get<Json::Object>(g.value);
    const auto between=std::get<Json::Array>(o.at("between").value);
    if(between[0]==Json{"front"} && between[1]==Json{"rear"}) { found=true; CHECK(std::get<double>(o.at("min_cells_across").value)>=2); CHECK(std::get<double>(o.at("width_m").value)==doctest::Approx(.06).epsilon(1e-6)); }
  }
  CHECK(found);
}
TEST_CASE("kutta mesh writes JSONL and a case that solves through the real CLI") {
  MeshTemp t; (void)surface(t);
  mesh_text(t.path("mesh.json"),R"({"v":1,"base_size":0.4,"surface":{"min_size":0.1,"max_size":0.1},"tunnel":{"upstream":0.3,"downstream":0.3,"side":0.3,"top":0.3},"layers":{"count":2,"first_height":0.005,"growth":1.2},"max_cells":200000})");
  const auto args="mesh "+quote(t.path("test.stl"))+" --setup "+quote(t.path("mesh.json"))+" --out "+quote(t.path("out"))+" --case --progress jsonl";
  REQUIRE_MESSAGE(cli(args,t.path("mesh.log"))==0,mesh_text(t.path("mesh.log")));
  CHECK(std::filesystem::exists(t.path("out/mesh-report.json")));
  std::istringstream lines(mesh_text(t.path("mesh.log"))); std::string line; int events=0,done=0; double fraction=-1; std::set<std::string> stages;
  while(std::getline(lines,line)) {
    const auto event=std::get<Json::Object>(parse_json(line).value); ++events;
    if(event.at("event")==Json{"mesh-done"}) ++done;
    if(event.at("event")==Json{"mesh-progress"}) { stages.insert(std::get<std::string>(event.at("stage").value)); const auto value=std::get<double>(event.at("fraction").value); CHECK(value>=fraction); CHECK(value<=1); fraction=value; }
  }
  CHECK(events>=8); CHECK(done==1); CHECK(stages.size()==7);
  // A bad setup must produce a terminal error and publish no output.
  mesh_text(t.path("bad.json"),R"({"v":1,"surface":{"min_size":-1}})");
  CHECK(cli("mesh "+quote(t.path("test.stl"))+" --setup "+quote(t.path("bad.json"))+" --out "+quote(t.path("bad"))+" --json",t.path("bad.log"))!=0);
  CHECK_FALSE(std::filesystem::exists(t.path("bad/mesh.msh")));
  CHECK(cli(args,t.path("repeat.log"))!=0);
  REQUIRE_MESSAGE(cli("solve "+quote(t.path("out/case.kutta.json"))+" --iterations 3 --backend cpu --out "+quote(t.path("run")),t.path("solve.log"))==0,mesh_text(t.path("solve.log")));
  const auto summary=std::get<Json::Object>(parse_json(mesh_text(t.path("run/summary.json"))).value);
  CHECK(summary.at("case_format")==Json{"kutta"});
  CHECK(std::filesystem::exists(t.path("run/forces.csv")));
}
