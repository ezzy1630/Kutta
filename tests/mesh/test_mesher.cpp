#include <doctest/doctest.h>
#include "../io/test_support.hpp"
#include "kutta/mesh/mesher/external_aero.hpp"
#include "../../src/kutta/mesh/mesher/internal.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/core/json_value.hpp"
#include <bit>
#include <cstdlib>
#include <numbers>
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
std::string sphere(double r,Vec3d c) {
  // UV sphere, outward winding.
  std::ostringstream out; out << "solid ball\n";
  constexpr int nu=48, nv=24;
  const auto at=[&](int i,int j) { const double u=2*std::numbers::pi*i/nu, v=std::numbers::pi*j/nv; return c+Vec3d{std::sin(v)*std::cos(u),std::sin(v)*std::sin(u),std::cos(v)}*r; };
  const auto tri=[&](Vec3d a,Vec3d b,Vec3d d) { out << "facet normal 0 0 0\nouter loop\n"; for(auto p:{a,b,d}) out << "vertex " << p.x << ' ' << p.y << ' ' << p.z << '\n'; out << "endloop\nendfacet\n"; };
  for(int i=0;i<nu;++i) for(int j=0;j<nv;++j) {
    const auto a=at(i,j), b=at(i+1,j), d=at(i+1,j+1), e=at(i,j+1);
    if(j>0) tri(a,e,b);
    if(j<nv-1) tri(b,e,d);
  }
  return out.str()+"endsolid ball\n";
}
std::string shifted_x(std::string text,double dx) {
  std::string out; std::istringstream in(text); std::string line;
  while(std::getline(in,line)) {
    if(line.rfind("vertex",0)==0) { std::istringstream v(line.substr(7)); double x,y,z; v>>x>>y>>z; std::ostringstream o; o<<"vertex "<<x+dx<<' '<<y<<' '<<z; line=o.str(); }
    out+=line+"\n";
  }
  return out;
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
  CHECK(defaults.layer_count==4); CHECK(defaults.proximity_floor()==doctest::Approx(.002));
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
  // Corner merges leave a few polyhedra; the rest of the core is hex.
  CHECK(report(result.report,{"cell_types","hex_fraction_of_core"})>.9);
  write_fluent_msh(t.path("test.msh"),result.mesh);
  const auto read=read_mesh(t.path("test.msh")); CHECK_NOTHROW(validate(read));
  CHECK(read.n_cells==result.mesh.n_cells); CHECK(read.n_faces()==result.mesh.n_faces()); CHECK(read.n_points()==result.mesh.n_points());
  REQUIRE(read.patches.size()==result.mesh.patches.size());
  for(std::size_t i=0;i<read.patches.size();++i) { CHECK(read.patches[i].name==result.mesh.patches[i].name); CHECK(read.patches[i].size==result.mesh.patches[i].size); }
  CHECK(volume(read)==doctest::Approx(volume(result.mesh)).epsilon(1e-12));
  auto refined=settings; refined.boxes.push_back({"wake",{{.7,.1,.1},{1.,.8,.7}},.05}); CHECK_NOTHROW(validate(generate(sf,refined).mesh));
  auto limited=settings; limited.max_cells=10; CHECK_THROWS(generate(sf,limited));
  auto symmetry=settings; symmetry.symmetry_y=true; CHECK_NOTHROW(validate(generate(sf,symmetry).mesh));
}
TEST_CASE("snap on a curved body keeps quality and surface error under the chord sag") {
  MeshTemp t;
  const auto result=generate(surface(t,sphere(.3,{.43,.47,.4})),parse_settings(R"({"v":1,"base_size":0.4,"surface":{"min_size":0.04,"max_size":0.1},"tunnel":{"upstream":0.3,"downstream":0.3,"side":0.3,"top":0.3},"layers":{"count":0},"max_cells":400000})"));
  CHECK_NOTHROW(validate(result.mesh));
  CHECK(report(result.report,{"quality","orthogonal_min"})>=.15);
  CHECK(report(result.report,{"snap","point_distance_m","max"})<.004);  // facets 0.04 m wide; chord sag 2 mm
  CHECK(report(result.report,{"cell_types","hex"})+report(result.report,{"cell_types","poly"})==result.mesh.n_cells);
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
TEST_CASE("omitted refinement boxes default to body-relative boxes; an explicit list is used as given") {
  CHECK(parse_settings(R"({"v":1})").body_boxes);
  CHECK_FALSE(parse_settings(R"({"v":1,"refinement_boxes":[]})").body_boxes);
  MeshTemp t; const auto sf=surface(t);
  const auto boxes=[](const Result& r) { return std::get<Json::Array>(std::get<Json::Object>(r.report.at("castellation").value).at("refinement_boxes").value); };
  auto settings=tiny();
  const auto defaulted=generate(sf,settings);
  const auto list=boxes(defaulted);
  REQUIRE(list.size()==2);
  const auto& body=std::get<Json::Object>(list[0].value);
  CHECK(body.at("name")==Json{"body"});
  CHECK(std::get<double>(body.at("size").value)==doctest::Approx(2*settings.max_size));
  // The body box encloses the body with 0.1 L margins (0.25 L downstream).
  const auto lo=std::get<Json::Array>(body.at("min").value), hi=std::get<Json::Array>(body.at("max").value);
  CHECK(std::get<double>(lo[0].value)==doctest::Approx(.13-.06)); CHECK(std::get<double>(hi[0].value)==doctest::Approx(.73+.15));
  CHECK(std::get<double>(lo[2].value)==0);  // ground
  CHECK(std::get<Json::Object>(list[1].value).at("name")==Json{"wake"});
  settings.body_boxes=false;
  CHECK(boxes(generate(sf,settings)).empty());
  // Memory estimate: cells times the measured bytes per cell.
  const auto& memory=std::get<Json::Object>(defaulted.report.at("estimated_solver_memory").value);
  CHECK(std::get<double>(memory.at("bytes").value)==doctest::Approx(defaulted.mesh.n_cells*2100.));
  CHECK(std::get<double>(memory.at("coupled_bytes").value)==doctest::Approx(defaulted.mesh.n_cells*2600.));
}
TEST_CASE("incremental displacement matches a fresh displacement of the same target") {
  // A 6x6x6 hex block whose z=0 face is a wall pushed up; the second run
  // changes the target on part of the wall.
  const auto m=make_box_hex(6,6,6,{0,0,0},{1,1,1});
  internal::Work w;
  w.points=m.points; w.cells=m.n_cells; w.kind.assign(std::size_t(m.n_cells),internal::Hex); w.size.assign(std::size_t(m.n_cells),1./6);
  w.planes.assign(m.points.size(),0);
  std::vector<Index> sources;
  for(Index f=0;f<m.n_faces();++f) {
    const auto ns=nodes(m,f);
    internal::Face face{{ns.begin(),ns.end()},m.owner[f],f<m.n_internal_faces()?m.neighbour[f]:-1,-1};
    if(face.neighbour<0) {
      bool bottom=true; for(auto n:ns) bottom=bottom && m.points[n].z<1e-12;
      face.patch=bottom?internal::first_wall:0;
      if(bottom) for(auto n:ns) sources.push_back(n);
    }
    w.faces.push_back(std::move(face));
  }
  std::sort(sources.begin(),sources.end()); sources.erase(std::unique(sources.begin(),sources.end()),sources.end());
  const auto adjacency=internal::cell_faces(w); const auto pf=internal::point_faces(w);
  std::vector<Vec3d> first(w.points.size()),second(w.points.size());
  for(auto p:sources) { first[p]={0,0,.02}; second[p]={0,0,w.points[p].x<.5?.03:.02}; }
  auto incremental=w; internal::Displacement step(incremental,adjacency,pf,sources,.15);
  step.run(first,{},"test"); step.run(second,{},"test");
  auto fresh=w; internal::Displacement once(fresh,adjacency,pf,sources,.15);
  once.run(second,{},"test");
  double worst=0;
  for(std::size_t p=0;p<w.points.size();++p) worst=std::max(worst,norm(incremental.points[p]-fresh.points[p]));
  CHECK(worst<1e-12);
  CHECK(norm(fresh.points[sources[0]]-w.points[sources[0]])>.019);
}
TEST_CASE("mesh report contract: shape counts, skewness range, histogram faces, contacts, progress, zero gap floor") {
  CHECK(parse_settings(R"({"v":1,"surface":{"proximity_min_size":0}})").proximity_floor()==doctest::Approx(.008/4));
  CHECK(parse_settings(R"({"v":1})").proximity_floor()==doctest::Approx(.002));
  CHECK(parse_settings(R"({"v":1,"surface":{"min_size":0.001,"max_size":0.01}})").proximity_floor()==doctest::Approx(.001/4));
  MeshTemp t;
  // The cube stands on the ground: a contact, not a closed slot.
  std::vector<double> fractions;
  const auto result=generate(surface(t,cube("body",1,false,{0,0,-.1})),tiny(R"({"count":2,"first_height":0.01,"growth":1.2})"),20,
                             [&](std::string_view,double f,std::string_view) { fractions.push_back(f); });
  CHECK_NOTHROW(validate(result.mesh));
  REQUIRE(fractions.size()>2);
  for(std::size_t i=1;i<fractions.size();++i) CHECK(fractions[i]>=fractions[i-1]);
  // Shape classes as kutta inspect counts them.
  Index tet=0,hex=0; const auto cf=build_cell_faces(result.mesh);
  for(Index c=0;c<result.mesh.n_cells;++c) {
    bool tri=true,quad=true; for(auto i=cf.offsets[c];i<cf.offsets[c+1];++i) { const auto n=nodes(result.mesh,cf.faces[i]).size(); tri=tri&&n==3; quad=quad&&n==4; }
    tet+=tri&&cf.offsets[c+1]-cf.offsets[c]==4; hex+=quad&&cf.offsets[c+1]-cf.offsets[c]==6;
  }
  CHECK(report(result.report,{"cell_types_by_shape","tet"})==tet);
  CHECK(report(result.report,{"cell_types_by_shape","hex"})==hex);
  CHECK(report(result.report,{"cell_types_by_shape","poly"})==result.mesh.n_cells-tet-hex);
  const auto& quality=std::get<Json::Object>(result.report.at("quality").value);
  CHECK(report(quality,{"skewness_max"})>=0); CHECK(report(quality,{"skewness_max"})<=1);
  CHECK(quality.contains("face_skewness_max"));
  for(const auto& bin:std::get<Json::Array>(std::get<Json::Object>(quality.at("adjacent_volume_ratio_core").value).at("histogram").value)) {
    const auto& o=std::get<Json::Object>(bin.value); CHECK(o.at("faces")==o.at("cells"));
  }
  bool contact=false;
  for(const auto& g:std::get<Json::Array>(std::get<Json::Object>(result.report.at("castellation").value).at("gaps").value)) {
    const auto& o=std::get<Json::Object>(g.value);
    if(std::get<Json::Array>(o.at("between").value)[1]==Json{"ground"}) { contact=std::get<bool>(o.at("contact").value); CHECK(std::get<double>(o.at("width_m").value)==0); }
    else CHECK_FALSE(std::get<bool>(o.at("contact").value));
  }
  CHECK(contact);
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
