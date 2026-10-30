#include "export.hpp"
#include "protocol.hpp"
#include "kutta/post/post.hpp"
#include "kutta/core/json.hpp"
#include "kutta/core/path.hpp"
#include "kutta/core/sha256.hpp"
#include "kutta/mesh/detail.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <iostream>
#include <limits>

namespace kutta::app {
namespace {
using namespace detail;
using view::Field;
constexpr float missing=std::numeric_limits<float>::quiet_NaN();
Field vector_field(Device& d,std::string name,std::string unit,const DeviceVec3& v) {
  const auto x=d.download(v.x),y=d.download(v.y),z=d.download(v.z);Field f{std::move(name),std::move(unit),3,{}};f.values.reserve(x.size()*3);
  for(std::size_t i=0;i<x.size();++i)f.values.insert(f.values.end(),{x[i],y[i],z[i]});return f;
}
void add_derived(std::vector<Field>& fields,std::optional<double> pressure,double dynamic) {
  // Derive before interpolation: magnitude is averaged from face magnitudes,
  // not the magnitude of a vector that may cancel around a corner.
  const auto count=fields.size();
  for(std::size_t i=0;i<count;++i) {
    const auto& f=fields[i];const auto suffix=f.name.ends_with("@fluent")?"@fluent":"";
    if(f.components==3) {
      Field mag{(f.name.starts_with("velocity")?std::string("velocity_magnitude"):std::string("wall_shear_magnitude"))+suffix,f.unit,1,{}};mag.values.reserve(f.values.size()/3);
      for(std::size_t j=0;j<f.values.size();j+=3)mag.values.push_back(float(std::hypot(double(f.values[j]),double(f.values[j+1]),double(f.values[j+2]))));
      fields.push_back(std::move(mag));
    } else if(f.name.starts_with("pressure") && pressure && dynamic>0 && std::isfinite(dynamic)) {
      Field cp{std::string("pressure_coefficient")+suffix,"1",1,{}};cp.values.reserve(f.values.size());
      for(float v:f.values)cp.values.push_back(float((v-*pressure)/dynamic));fields.push_back(std::move(cp));
    }
  }
}
}
view::Fields collect_view_fields(Device& dev,const fv::SteadySolver& solver,
                     const fluent::FluentCase& c,const fluent::FluentData* reference,const setup::FluentSetup& setup,std::string_view provenance,bool volume) {
  view::Fields s;s.reference=reference!=nullptr;s.provenance=provenance;
  s.properties={setup.settings.density,setup.settings.viscosity,c.reference.density,c.reference.velocity,c.reference.pressure,c.reference.area,c.reference.length};
  const auto& m=solver.mesh();const auto nb=m.n_faces()-m.n_internal_faces();
  // Turbulence fields and y+ exist only for the active closure.
  const bool sst=setup.settings.sst.has_value(),epsilon=setup.settings.realizable_ke.has_value();
  auto cell_velocity=vector_field(dev,"velocity","m/s",solver.velocity());
  if(volume) {
    s.cells.push_back(std::move(cell_velocity));
    s.cells.push_back({"pressure","Pa",1,dev.download(solver.pressure())});
    if(sst || epsilon) s.cells.push_back({"k","m^2/s^2",1,dev.download(solver.kinetic_energy())});
    if(sst) s.cells.push_back({"omega","1/s",1,dev.download(solver.specific_dissipation())});
    if(epsilon) s.cells.push_back({"epsilon","m^2/s^3",1,dev.download(solver.dissipation())});
    if(sst || epsilon) {
      auto ratio=dev.download(sst?solver.turbulence().mu_t:solver.realizable_turbulence().mu_t);for(auto& v:ratio)v=float(v/setup.settings.viscosity);
      s.cells.push_back({"viscosity_ratio","1",1,std::move(ratio)});
    }
  }
  s.boundary.push_back({"pressure","Pa",1,dev.download(solver.boundary_pressure())});
  s.boundary.push_back(vector_field(dev,"wall_shear","Pa",solver.wall_shear()));
  // k-epsilon wall functions report y*, Fluent's wall-function coordinate.
  if(sst) s.boundary.push_back({"y_plus","1",1,dev.download(solver.turbulence().wall_y_plus)});
  if(epsilon) s.boundary.push_back({"y_star","1",1,dev.download(solver.realizable_turbulence().wall_y_star)});
  Field velocity{"velocity","m/s",3,std::vector<float>(std::size_t(nb)*3)};
  for(const auto& p:m.patches) {
    auto bc=std::find_if(setup.boundaries.begin(),setup.boundaries.end(),[&](const auto& b){return b.patch==p.name;});
    KUTTA_REQUIRE(bc!=setup.boundaries.end(),"missing export boundary condition");
    for(Index f=p.start;f<p.start+p.size;++f) {
      const auto cell=std::size_t(m.owner[f])*3;const auto& u=volume?s.cells.front().values:cell_velocity.values;
      Vec3d v{u[cell],u[cell+1],u[cell+2]};
      if(const auto* inlet=std::get_if<fv::VelocityInlet>(&bc->condition))v=std::get<Vec3d>(inlet->velocity);
      else if(std::holds_alternative<fv::Symmetry>(bc->condition)){const auto n=solver.geometry().face_area[f]/norm(solver.geometry().face_area[f]);v=v-n*dot(v,n);}
      else if(const auto* wall=std::get_if<fv::Wall>(&bc->condition))v=wall->velocity;
      const auto i=std::size_t(f-m.n_internal_faces())*3;velocity.values[i]=float(v.x);velocity.values[i+1]=float(v.y);velocity.values[i+2]=float(v.z);
    }
  }
  s.boundary.push_back(std::move(velocity));
  if(volume) {
    const auto mass_flux=dev.download(solver.mass_flux());
    s.boundary.push_back({"mass_flux","kg/s",1,std::vector<float>(mass_flux.begin()+m.n_internal_faces(),mass_flux.end())});
  }
  if(reference) {
    // Cell-zone vectors retain original Fluent-ID order through RCM; their
    // stored indices are updated by the mesh permutation.
    std::vector<std::uint64_t> ids(m.n_cells);
    for(const auto& z:m.cell_zones) {
      auto range=std::find_if(c.mesh.cell_zones.begin(),c.mesh.cell_zones.end(),[&](const auto& r){return r.zone_id==z.source_zone_id;});
      KUTTA_REQUIRE(range!=c.mesh.cell_zones.end() && range->size()==z.cells.size(),"reference cell zone mismatch");
      for(std::size_t i=0;i<z.cells.size();++i)ids[z.cells[i]]=range->first+i;
    }
    const auto sample=[&](const auto& fields,std::string_view key,std::uint64_t id,int component=0) {
      auto it=fields.find(std::string(key));if(it==fields.end())return missing;
      const auto values=it->second.at(id);return values.size()>std::size_t(component)?float(values[component]):missing;
    };
    if(volume) {
      for(const auto& [name,source,unit]:std::vector<std::array<std::string,3>>{
          {"pressure","SV_P","Pa"},{"k","SV_K","m^2/s^2"},{"omega","SV_O","1/s"},{"epsilon","SV_D","m^2/s^3"},{"viscosity_ratio","SV_MU_T","1"}}) {
        if(!reference->cell_fields.contains(source))continue;
        Field f{name+"@fluent",unit,1,{}};f.values.reserve(ids.size());
        for(auto id:ids){const auto v=sample(reference->cell_fields,source,id);f.values.push_back(name=="viscosity_ratio"?float(v/setup.settings.viscosity):v);}s.cells.push_back(std::move(f));
      }
      Field u{"velocity@fluent","m/s",3,{}};
      for(auto id:ids)for(auto key:{"SV_U","SV_V","SV_W"})u.values.push_back(sample(reference->cell_fields,key,id));s.cells.push_back(std::move(u));
    }
    Field p{"pressure@fluent","Pa",1,std::vector<float>(nb,missing)},shear{"wall_shear@fluent","Pa",3,std::vector<float>(std::size_t(nb)*3,missing)},yp{"y_plus@fluent","1",1,std::vector<float>(nb,missing)},uv{"velocity@fluent","m/s",3,std::vector<float>(std::size_t(nb)*3,missing)};
    for(const auto& patch:m.patches) {
      auto range=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& r){return r.zone_id==patch.source_zone_id;});
      KUTTA_REQUIRE(range!=c.mesh.face_zones.end() && range->size()==std::uint64_t(patch.size),"reference face zone mismatch");
      for(Index f=patch.start;f<patch.start+patch.size;++f) {
        const auto old=solver.permutation().faces_new_to_old[f];KUTTA_REQUIRE(old>=patch.start && old<patch.start+patch.size,"reference boundary permutation mismatch");
        const auto id=range->first+std::uint64_t(old-patch.start);const auto b=f-m.n_internal_faces();
        p.values[b]=sample(reference->face_fields,"SV_P",id);
        if(std::isnan(p.values[b]))p.values[b]=sample(reference->cell_fields,"SV_P",ids[m.owner[f]]);
        yp.values[b]=sample(reference->face_fields,"SV_WALL_YPLUS",id);
        for(int axis=0;axis<3;++axis) {
          // Fluent stores integrated stress on the fluid (N/face). Normalize
          // by the physical face area before surface averaging or deltas.
          shear.values[std::size_t(b)*3+axis]=float(sample(reference->face_fields,"SV_WALL_SHEAR",id,axis)/norm(solver.geometry().face_area[f]));
          const char* keys[]={"SV_U","SV_V","SV_W"};float v=sample(reference->face_fields,keys[axis],id);
          if(std::isnan(v))v=sample(reference->cell_fields,keys[axis],ids[m.owner[f]]);
          uv.values[std::size_t(b)*3+axis]=v;
        }
      }
    }
    for(auto* f:{&p,&shear,&yp,&uv})s.boundary.push_back(std::move(*f));
  }
  double dynamic=0;
  if(c.reference.density && c.reference.velocity && *c.reference.density>0 && *c.reference.velocity>0)dynamic=.5**c.reference.density**c.reference.velocity**c.reference.velocity;
  add_derived(s.cells,c.reference.pressure,dynamic);add_derived(s.boundary,c.reference.pressure,dynamic);
  return s;
}
void save_view_input(const std::filesystem::path& path,Device& dev,const fv::SteadySolver& solver,
                     const fluent::FluentCase& c,const fluent::FluentData* reference,const setup::FluentSetup& setup,std::string_view provenance) {
  view::write_snapshot(path,solver.mesh(),collect_view_fields(dev,solver,c,reference,setup,provenance,true));
}
int export_view(int argc,char** argv) {
  if(argc<3){std::cerr<<"Usage: kutta export-view DIR [--slice ox,oy,oz,nx,ny,nz]... [--post spec.json | --streamlines spec.json]\n";return 2;}
  std::vector<view::Plane> planes;std::optional<post::Request> request;
  try {for(int i=3;i<argc;i+=2){
    KUTTA_REQUIRE(i+1<argc,"missing export option value");const std::string_view key=argv[i];
    if(key=="--slice")planes.push_back(view::parse_plane(argv[i+1]));
    else if(key=="--post" || key=="--streamlines"){KUTTA_REQUIRE(!request,"only one post-processing specification is allowed");request=post::read_request(argv[i+1],key=="--streamlines");}
    else throw Error("unknown export option: "+std::string(key));
  }}
  catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
  try {
    const auto start=std::chrono::steady_clock::now();const std::filesystem::path directory=argv[2];
    KUTTA_REQUIRE(std::filesystem::exists(directory/"summary.json"),"run has no completion marker summary.json");
    std::string cache;
    if(request) {
      // The build identity is part of every cache key, so a rebuilt CLI never
      // reuses work computed by different code.
      const auto stats=post::export_view(directory,planes,*request,std::string(KUTTA_VERSION)+"+"+KUTTA_COMMIT);
      cache=" (post cache: reused "+std::to_string(stats.reused)+", computed "+std::to_string(stats.computed)+")";
    } else view::export_view(directory,view::read_snapshot(directory/"view-input.bin"),planes);
    std::cout<<"Wrote "<<(directory/"view").string()<<" in "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" s"<<cache<<"\n";return 0;
  }catch(const std::exception& e){std::cerr<<"kutta export-view: "<<e.what()<<'\n';return 1;}
}
int export_geometry(int argc,char** argv) {
  bool json=false;std::filesystem::path input,out;
  for(int i=2;i<argc;++i) {
    const std::string_view arg=argv[i];
    if(arg=="--json" && !json)json=true;
    else if(arg=="--out" && out.empty() && i+1<argc)out=argv[++i];
    else if(!arg.starts_with("--") && input.empty())input=argv[i];
    else input.clear(),out.clear(),i=argc;
  }
  if(input.empty() || out.empty()) {
    const auto message="Usage: kutta export-geometry CASE --out DIR [--json]";
    if(json)error_event(message);
    std::cerr<<message<<'\n';return 2;
  }
  try {
    const auto start=std::chrono::steady_clock::now();
    // Hash concurrently with the mesh read; both are bounded by file I/O.
    auto hash=std::async(std::launch::async,[&]{return sha256_file(input);});
    const auto mesh=read_mesh(input);
    std::filesystem::create_directories(out);
    const auto stage=out/(".geometry-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::uint64_t triangles=0;const auto target=out/"geometry",old=out/(".geometry-old-"+stage.filename().string().substr(10));
    try {
      const auto written=view::write_surface_view(stage,mesh,view::Fields{},true,{{"geometry_only",true},{"case",Json::Object{{"sha256",hash.get()}}}});
      for(const auto& z:std::get<Json::Array>(written.at("zones").value))triangles+=std::uint64_t(std::get<double>(std::get<Json::Object>(z.value).at("indices").value))/3;
      if(std::filesystem::exists(target))std::filesystem::rename(target,old);
      try {std::filesystem::rename(stage,target);}catch(...){if(std::filesystem::exists(old))std::filesystem::rename(old,target);throw;}
      std::filesystem::remove_all(old);
    }catch(...){std::error_code ec;std::filesystem::remove_all(stage,ec);throw;}
    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    const auto manifest=utf8_path(std::filesystem::absolute(target/"manifest.json"));
    if(json)std::cout<<dump_json(Json::Object{{"v",1},{"event","geometry"},{"manifest",manifest},{"zones",double(mesh.patches.size())},{"triangles",double(triangles)},{"seconds",seconds}})<<std::endl;
    else std::cout<<"Wrote "<<manifest<<" ("<<mesh.patches.size()<<" zones, "<<triangles<<" triangles) in "<<seconds<<" s\n";
    return 0;
  }catch(const std::exception& e){if(json)error_event(e.what());std::cerr<<"kutta export-geometry: "<<e.what()<<'\n';return 1;}
}
} // namespace kutta::app
