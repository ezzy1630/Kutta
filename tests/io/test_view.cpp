#include "kutta/io/view.hpp"
#include "kutta/io/binary.hpp"
#include "kutta/mesh/detail.hpp"
#include "test_support.hpp"
#include <doctest/doctest.h>
#include <cstdlib>
#include <fstream>
#include <limits>

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
