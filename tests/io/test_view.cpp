#include "kutta/io/view.hpp"
#include "kutta/io/binary.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/mesher/external_aero.hpp"
#include "kutta/core/json.hpp"
#include <sstream>
#include "test_support.hpp"
#include <doctest/doctest.h>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>

namespace {
using namespace kutta;
std::string quote(const std::filesystem::path& path) {
#ifdef _WIN32
  return '"'+std::filesystem::path(path).make_preferred().string()+'"';
#else
  std::string out="'";for(char c:path.string())out+=c=='\''?"'\\''":std::string(1,c);return out+"'";
#endif
}
}
TEST_CASE("View export preserves known surface and slice geometry, area weights and Fluent deltas") {
  MeshTemp t;view::Snapshot s;s.mesh=make_box_hex(2,1,1,{0,0,0},{2,1,1});s.reference=true;s.provenance="fixture";
  s.mesh.patches[3].kind=PatchKind::Symmetry;
  for(std::size_t i=0;i<s.mesh.patches.size();++i){s.mesh.patches[i].source_zone_id=int(100+i);s.mesh.patches[i].source_type=i==3?"symmetry":"wall";}
  s.cells={{"pressure","Pa",1,{10,20}},{"pressure@fluent","Pa",1,{8,18}},
    {"velocity","m/s",3,{1,2,3,4,5,6}},{"velocity@fluent","m/s",3,{0,1,2,3,4,5}},
    {"k","m^2/s^2",1,{2,3}},{"omega","1/s",1,{5,6}},{"viscosity_ratio","1",1,{7,8}}};
  const auto nb=s.mesh.n_faces()-s.mesh.n_internal_faces();
  for(auto name:{"pressure","pressure@fluent","wall_shear_magnitude","y_plus","velocity","wall_shear"}) {
    view::Field f{name,std::string_view(name).starts_with("pressure") || std::string_view(name).starts_with("wall_shear")?"Pa":std::string_view(name)=="velocity"?"m/s":"1",std::string_view(name)=="velocity" || std::string_view(name)=="wall_shear"?3u:1u,{}};
    for(Index b=0;b<nb;++b) {
      const auto c=s.mesh.owner[b+s.mesh.n_internal_faces()];
      if(f.components==3)f.values.insert(f.values.end(),{3,4,0});
      else f.values.push_back(std::string_view(name)=="pressure"?float(10+10*c):std::string_view(name)=="pressure@fluent"?float(8+10*c):5);
    }
    s.boundary.push_back(std::move(f));
  }
  view::write_snapshot(t.path("view-input.bin"),s);
  const auto copy=view::read_snapshot(t.path("view-input.bin"));
  CHECK(copy.mesh.face_nodes==s.mesh.face_nodes);CHECK(copy.cells[0].values==s.cells[0].values);CHECK(copy.reference);
  // A v1 snapshot ends before the seven optional v2 properties. Preserve
  // backwards reads rather than making old finished runs unusable.
  {std::ifstream file(t.path("view-input.bin"),std::ios::binary);std::string bytes{std::istreambuf_iterator<char>(file),{}};bytes[4]=1;bytes.resize(bytes.size()-7*4);mesh_text(t.path("legacy.bin"),bytes);}
  const auto legacy=view::read_snapshot(t.path("legacy.bin"));CHECK(legacy.cells[0].values==s.cells[0].values);CHECK_FALSE(legacy.properties.density);

  mesh_text(t.path("summary.json"),"{\"v\":1}");
  std::vector<view::Plane> planes;
  for(auto text:{"0.5,0,0,1,0,0","0,0.5,0,0,1,0","1,0,0,1,0,0","5,0,0,1,0,0","0,1,0,1,1,0"})planes.push_back(view::parse_plane(text));
  view::export_view(t.directory,copy,planes);
  auto command=quote(KUTTA_PYTHON_PATH)+" "+quote(std::filesystem::path(KUTTA_SOURCE_DIR)/"tests/io/check_view.py")+" "+quote(t.directory);
#ifdef _WIN32
  command='"'+command+'"';
#endif
  CHECK(std::system(command.c_str())==0);
  view::export_view(t.directory,copy);
  CHECK_FALSE(std::filesystem::exists(t.path("view/slice-0.bin")));
  CHECK_THROWS(view::parse_plane("0,0,0,1,0,nan"));
  CHECK_THROWS(view::parse_plane("0,0,0,0,0,0"));
  CHECK_THROWS(view::parse_plane("0,0,0,1,0,0,"));
  // Truncated, unknown-version and invalid-topology snapshots fail at the
  // file boundary rather than allocating from untrusted lengths.
  mesh_text(t.path("bad.bin"),"KUTR");CHECK_THROWS(view::read_snapshot(t.path("bad.bin")));
  {std::fstream f(t.path("view-input.bin"),std::ios::in|std::ios::out|std::ios::binary);f.seekp(4);const char version=99;f.write(&version,1);}
  CHECK_THROWS(view::read_snapshot(t.path("view-input.bin")));
}
TEST_CASE("View export triangulates a concave polygon and its nonconvex volume slice") {
  using namespace kutta::detail;
  MeshTemp t;view::Snapshot s;
  const std::vector<Vec3d> base{{0,0,0},{2,0,0},{2,1,0},{1,1,0},{1,2,0},{0,2,0}};
  std::vector<Vec3d> points=base;for(auto p:base)points.push_back({p.x,p.y,1});
  std::vector<kutta::detail::CellFaces> cells(1);cells[0].push_back({5,4,3,2,1,0});cells[0].push_back({6,7,8,9,10,11});
  for(Index i=0;i<6;++i){const auto next=(i+1)%6;cells[0].push_back({i,next,next+6,i+6});}
  s.mesh=from_cells(points,cells,{{"shell",PatchKind::Wall,"wall",1}},[](const auto&){return 0;});
  s.cells={{"pressure","Pa",1,{7}}};s.boundary={{"pressure","Pa",1,std::vector<float>(s.mesh.n_faces(),7)}};
  const auto plane=view::parse_plane("0,0,0.5,0,0,1");
  view::export_view(t.directory,s,std::span(&plane,1));
  binary::Reader r(t.path("view/slice-0.bin"));CHECK(r.scalar<std::uint32_t>()==0x5654554b);
  CHECK(r.scalar<std::uint32_t>()==1);CHECK(r.scalar<std::uint32_t>()==32);
  const auto nv=r.scalar<std::uint32_t>(),ni=r.scalar<std::uint32_t>();
  (void)r.scalar<std::uint32_t>();(void)r.scalar<std::uint32_t>();(void)r.scalar<std::uint32_t>();
  const auto p=r.array<float>(nv*3);const auto triangles=r.array<std::uint32_t>(ni);double area=0;
  for(std::size_t i=0;i<triangles.size();i+=3){Vec3d v[3];for(int j=0;j<3;++j){const auto n=triangles[i+j]*3;v[j]={p[n],p[n+1],p[n+2]};CHECK(v[j].z==doctest::Approx(.5));}area+=norm(cross(v[1]-v[0],v[2]-v[0]))*.5;}
  CHECK(area==doctest::Approx(3));
}
TEST_CASE("Zone edge polylines contain each polygon edge exactly once and no triangulation diagonal") {
  using namespace kutta::detail;
  const std::vector<Vec3d> base{{0,0,0},{2,0,0},{2,1,0},{1,1,0},{1,2,0},{0,2,0}};
  std::vector<Vec3d> points=base;for(auto p:base)points.push_back({p.x,p.y,1});
  std::vector<kutta::detail::CellFaces> cells(1);cells[0].push_back({5,4,3,2,1,0});cells[0].push_back({6,7,8,9,10,11});
  for(Index i=0;i<6;++i){const auto next=(i+1)%6;cells[0].push_back({i,next,next+6,i+6});}
  view::Snapshot prism;prism.mesh=from_cells(points,cells,{{"shell",PatchKind::Wall,"wall",1}},[](const auto&){return 0;});
  view::Snapshot box;box.mesh=make_box_hex(3,2,2,{0,0,0},{3,2,2});
  for(auto* s:{&prism,&box}) {
    MeshTemp t;view::export_view(t.directory,*s);
    for(std::size_t z=0;z<s->mesh.patches.size();++z) {
      const auto& patch=s->mesh.patches[z];CAPTURE(z);
      binary::Reader surface(t.path("view/zone-"+std::to_string(z)+".surface.bin")),edges(t.path("view/zone-"+std::to_string(z)+".edges.bin"));
      std::uint32_t sh[8],eh[8];for(auto& n:sh)n=surface.scalar<std::uint32_t>();for(auto& n:eh)n=edges.scalar<std::uint32_t>();
      CHECK(eh[0]==0x4c54554b);CHECK(eh[3]==sh[3]);CHECK(eh[5]==0);
      const auto positions=surface.array<float>(sh[3]*3);CHECK(edges.array<float>(eh[3]*3)==positions);
      const auto indices=edges.array<std::uint32_t>(eh[4]);const auto offsets=edges.array<std::uint32_t>(eh[6]+1);
      // Polygon edges by position; zone vertex order is the surface's.
      const auto key=[&](std::uint32_t a,std::uint32_t b){std::array<float,6> x{positions[a*3],positions[a*3+1],positions[a*3+2],positions[b*3],positions[b*3+1],positions[b*3+2]};
        std::array<float,3> l{x[0],x[1],x[2]},r{x[3],x[4],x[5]};return l<r?std::pair(l,r):std::pair(r,l);};
      std::map<std::array<float,3>,std::uint32_t> id;for(std::uint32_t v=0;v<sh[3];++v)id[{positions[v*3],positions[v*3+1],positions[v*3+2]}]=v;
      std::set<std::pair<std::array<float,3>,std::array<float,3>>> polygon,chained;std::size_t chained_count=0;
      for(Index f=patch.start;f<patch.start+patch.size;++f){const auto ns=nodes(s->mesh,f);for(std::size_t i=0;i<ns.size();++i){const auto a=s->mesh.points[ns[i]],b=s->mesh.points[ns[(i+1)%ns.size()]];
        polygon.insert(key(id.at({float(a.x),float(a.y),float(a.z)}),id.at({float(b.x),float(b.y),float(b.z)})));}}
      for(std::uint32_t l=0;l<eh[6];++l)for(auto i=offsets[l];i+1<offsets[l+1];++i){chained.insert(key(indices[i],indices[i+1]));++chained_count;}
      CHECK(chained==polygon);CHECK(chained_count==polygon.size());CHECK(indices.size()==polygon.size()+eh[6]);
    }
  }
}

namespace {
struct SliceRead {double area=0;std::uint64_t skipped=0;std::size_t triangles=0;bool on_plane=true;std::set<float> values;};
// Slice n of a published view: signed area along the plane normal, whether
// every vertex lies on the plane, and the distinct "cell" values it shows.
SliceRead read_slice(const std::filesystem::path& dir,std::size_t n) {
  using namespace kutta::detail;
  const auto manifest=std::get<Json::Object>(parse_json(mesh_text(dir/"view/manifest.json")).value);
  const auto& d=std::get<Json::Object>(std::get<Json::Array>(manifest.at("slices").value)[n].value);
  const auto vec=[&](const char* key){const auto& a=std::get<Json::Array>(d.at(key).value);return Vec3d{std::get<double>(a[0].value),std::get<double>(a[1].value),std::get<double>(a[2].value)};};
  const auto origin=vec("origin"),normal=vec("normal");SliceRead out;
  out.skipped=std::uint64_t(std::get<double>(d.at("skipped_cells").value));
  binary::Reader r(dir/("view/slice-"+std::to_string(n)+".bin"));
  for(int i=0;i<3;++i)(void)r.scalar<std::uint32_t>();
  const auto nv=r.scalar<std::uint32_t>(),ni=r.scalar<std::uint32_t>();for(int i=0;i<3;++i)(void)r.scalar<std::uint32_t>();
  const auto p=r.array<float>(nv*3);const auto t=r.array<std::uint32_t>(ni);out.triangles=ni/3;
  for(std::size_t i=0;i<t.size();i+=3) {
    Vec3d v[3];for(int j=0;j<3;++j){const auto k=t[i+j]*3;v[j]={p[k],p[k+1],p[k+2]};out.on_plane &= std::abs(dot(v[j]-origin,normal))<1e-6;}
    out.area+=dot(cross(v[1]-v[0],v[2]-v[0]),normal)*.5;
  }
  const auto& fields=std::get<Json::Array>(d.at("fields").value);
  for(const auto& f:fields)if(std::get<Json::Object>(f.value).at("name")==Json{"cell"}) {
    const auto values=r.array<float>(nv);out.values.insert(values.begin(),values.end());
  }
  return out;
}
std::string mesher_cube() {
  using namespace kutta::detail;
  const auto m=make_box_hex(1,1,1,{.13,.17,.1},{.73,.77,.7});
  std::ostringstream out;out<<"solid body\n";
  for(Index f=0;f<m.n_faces();++f){const auto ns=nodes(m,f);for(int k=1;k<3;++k){out<<"facet normal 0 0 0\nouter loop\n";for(auto n:{ns[0],ns[k],ns[k+1]}){const auto q=m.points[n];out<<"vertex "<<q.x<<' '<<q.y<<' '<<q.z<<'\n';}out<<"endloop\nendfacet\n";}}
  return out.str()+"endsolid body\n";
}
}
TEST_CASE("Slices of a snapped mesher mesh close through hanging faces, vertices and near-plane nodes") {
  using namespace kutta::detail;
  // The box the mesher's CLI test meshes, with prism layers (polyhedra where
  // they terminate) and a refinement box, so the core has 2:1 hanging faces.
  MeshTemp t;mesh_text(t.path("body.stl"),mesher_cube());const std::array stl{t.path("body.stl")};
  auto settings=mesher::parse_settings(R"({"v":1,"base_size":0.4,"surface":{"min_size":0.1,"max_size":0.1},"tunnel":{"upstream":0.3,"downstream":0.3,"side":0.3,"top":0.3},"layers":{"count":2,"first_height":0.005,"growth":1.2},"max_cells":200000})");
  settings.boxes.push_back({"wake",{{.7,.1,.1},{1.2,.8,.7}},.05});
  auto generated=mesher::generate(mesher::read_stl(stl,settings),settings);
  view::Snapshot s;s.mesh=std::move(generated.mesh);
  REQUIRE(std::get<double>(generated.report.at("max_face_level_difference").value)>=1);
  // Each cell shows its own index, so overlapping or duplicated cells and
  // the coincident-face owner rule are visible in the values.
  view::Field cell{"cell","1",1,{}};for(Index c=0;c<s.mesh.n_cells;++c)cell.values.push_back(float(c));s.cells={cell};
  Vec3d lo=s.mesh.points.front(),hi=lo;
  for(auto q:s.mesh.points){lo={std::min(lo.x,q.x),std::min(lo.y,q.y),std::min(lo.z,q.z)};hi={std::max(hi.x,q.x),std::max(hi.y,q.y),std::max(hi.z,q.z)};}
  const auto eps=view::slice_tolerance(s.mesh);
  // The octree level with the most nodes inside the body's height: the plane
  // passes exactly through its vertices and along its faces and edges.
  std::map<double,int> levels;for(auto q:s.mesh.points)if(q.z>.15 && q.z<.65)++levels[q.z];
  const auto level=std::max_element(levels.begin(),levels.end(),[](auto a,auto b){return a.second<b.second;})->first;
  CAPTURE(level);REQUIRE(levels[level]>50);
  const double tunnel=(hi.x-lo.x)*(hi.y-lo.y),expected=tunnel-.36;
  const auto slice_at=[&](std::initializer_list<std::string> planes){
    std::vector<view::Plane> p;for(const auto& text:planes)p.push_back(view::parse_plane(text));view::export_view(t.directory,s,p);};
  const auto z=[](double v){std::ostringstream o;o.precision(17);o<<"0,0,"<<v<<",0,0,1";return o.str();};
  slice_at({z(level),z(level+1e-4)});
  for(std::size_t n=0;n<2;++n){const auto r=read_slice(t.directory,n);CAPTURE(n);CHECK(r.skipped==0);CHECK(r.on_plane);CHECK(r.area==doctest::Approx(expected).epsilon(1e-6));}
  // Snapping leaves nodes scattered across the tolerance: some within eps of
  // the plane, some just beyond, on both sides. This failed every export
  // ("slice intersection is not a closed polygon") and, where a cell had
  // some faces within eps, left holes.
  std::uint32_t seed=1;
  for(auto& q:s.mesh.points)if(q.z==level){seed=seed*1664525u+1013904223u;q.z+=eps*(double(seed>>8)/double(1u<<24)*6-3);}
  slice_at({z(level),z(level+1e-4),z(level-1e-4),"0.43,0.47,"+std::to_string(level)+",0,0,-1"});
  for(std::size_t n=0;n<4;++n){const auto r=read_slice(t.directory,n);CAPTURE(n);CHECK(r.skipped==0);CHECK(r.on_plane);CHECK(std::abs(r.area)==doctest::Approx(expected).epsilon(1e-6));}
  // A plane through lattice vertices along cell diagonals, compared with the
  // same plane moved off every vertex.
  const auto diagonal=[&](double shift){std::ostringstream o;o.precision(17);o<<.43+shift<<",0.47,0.4,1,1,0";return o.str();};
  slice_at({diagonal(0),diagonal(1e-7)});
  const auto a=read_slice(t.directory,0),b=read_slice(t.directory,1);
  CHECK(a.skipped==0);CHECK(a.on_plane);CHECK(a.area>0);CHECK(a.area==doctest::Approx(b.area).epsilon(1e-6));
}
TEST_CASE("A slice on an interior face emits it once with the owner's values, and a cell whose faces do not close is skipped and counted") {
  using namespace kutta::detail;
  MeshTemp t;view::Snapshot s;s.mesh=make_box_hex(2,1,1,{0,0,0},{2,1,1});
  s.cells={{"cell","1",1,{0,1}}};
  for(auto text:{"1,0,0,1,0,0","1,0,0,-1,0,0"}) {
    const auto plane=view::parse_plane(text);view::export_view(t.directory,s,std::span(&plane,1));
    const auto r=read_slice(t.directory,0);CHECK(r.skipped==0);CHECK(std::abs(r.area)==doctest::Approx(1));
    CHECK(r.values==std::set<float>{float(s.mesh.owner[0])});
  }
  // Give one boundary face of cell 1 its own copy of a node: two of the
  // cell's edges are then not shared by two of its faces, so its faces no
  // longer bound a closed surface. The export still succeeds, without it.
  auto open=s;const auto plane=view::parse_plane("1.5,0,0,1,0,0");
  for(Index f=open.mesh.n_internal_faces();f<open.mesh.n_faces();++f) {
    const auto ns=nodes(open.mesh,f);
    if(open.mesh.owner[f]!=1 || std::all_of(ns.begin(),ns.end(),[&](Index n){return open.mesh.points[n].x==open.mesh.points[ns[0]].x;}))continue;
    open.mesh.points.push_back(open.mesh.points[ns[2]]);open.mesh.face_nodes[std::size_t(open.mesh.face_offsets[f]+2)]=Index(open.mesh.points.size()-1);break;
  }
  view::export_view(t.directory,open,std::span(&plane,1));
  const auto r=read_slice(t.directory,0);CHECK(r.skipped==1);CHECK(r.triangles==0);
}
