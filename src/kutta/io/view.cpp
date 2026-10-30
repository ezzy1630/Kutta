#include "view.hpp"
#include "view_geometry.hpp"
#include "binary.hpp"
#include "kutta/core/json.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>

namespace kutta::view {
namespace {
using namespace detail;
constexpr float missing=std::numeric_limits<float>::quiet_NaN();
Vec3d polygon_area(std::span<const Vec3d> p) {
  Vec3d area;for(std::size_t i=1;i+1<p.size();++i)area=area+cross(p[i]-p[0],p[i+1]-p[0])*.5;return area;
}
bool visible(std::string_view name,bool wall) {
  const auto at=name.find('@');name=name.substr(0,at);
  if(name=="mass_flux")return false;
  if(name.starts_with("wall_shear") || name.starts_with("skin_friction") || name=="y_plus" || name=="y_star")return wall;
  if(name=="velocity" || name=="velocity_magnitude")return !wall;
  return true;
}
Geometry surface(const PolyMesh& m,const Fields& s,const Patch& patch,std::vector<Index>& local) {
  Geometry g;std::vector<Index> vertices;std::vector<double> area(patch.size);
  for(Index f=patch.start;f<patch.start+patch.size;++f) {
    std::vector<Vec3d> points;std::vector<std::uint32_t> ids;
    for(auto n:nodes(m,f)) {
      if(local[n]<0){local[n]=index(vertices.size());vertices.push_back(n);}
      ids.push_back(std::uint32_t(local[n]));points.push_back(m.points[n]);
    }
    area[f-patch.start]=norm(polygon_area(points));
    for(auto i:triangulate(points))g.indices.push_back(ids[i]);
  }
  for(auto n:vertices){const auto p=m.points[n];g.positions.insert(g.positions.end(),{float(p.x),float(p.y),float(p.z)});}
  for(const auto& f:s.boundary)if(visible(f.name,patch.kind==PatchKind::Wall)) {
    Field out{f.name,f.unit,f.components,{}};out.values.resize(vertices.size()*f.components);
    std::vector<double> sums(out.values.size()),weights(vertices.size());std::vector<bool> absent(vertices.size());
    for(Index face=patch.start;face<patch.start+patch.size;++face) {
      const auto b=std::size_t(face-m.n_internal_faces())*f.components;bool valid=true;
      for(std::uint32_t c=0;c<f.components;++c)valid &= std::isfinite(f.values[b+c]);
      for(auto n:nodes(m,face)) {
        const auto v=local[n];absent[v]=absent[v] || !valid;
        if(valid){const auto a=area[face-patch.start];weights[v]+=a;for(std::uint32_t c=0;c<f.components;++c)sums[std::size_t(v)*f.components+c]+=a*f.values[b+c];}
      }
    }
    for(std::size_t v=0;v<vertices.size();++v)for(std::uint32_t c=0;c<f.components;++c)out.values[v*f.components+c]=absent[v] || weights[v]==0?missing:float(sums[v*f.components+c]/weights[v]);
    if(std::any_of(out.values.begin(),out.values.end(),[](float v){return std::isfinite(v);}))g.fields.push_back(std::move(out));
  }
  for(auto n:vertices)local[n]=-1;
  append_deltas(g.fields);return g;
}
void append_polygon(Geometry& g,std::vector<Vec3d> ring,const Plane& plane,Index cell,std::span<const Field> fields) {
  if(ring.size()<3 || norm(polygon_area(ring))==0)return;
  if(dot(polygon_area(ring),plane.normal)<0)std::reverse(ring.begin(),ring.end());
  const auto base=std::uint32_t(g.positions.size()/3);
  KUTTA_REQUIRE(std::uint64_t(base)+ring.size()<=UINT32_MAX,"view exceeds uint32 vertices");
  for(auto i:triangulate(ring))g.indices.push_back(base+i);
  for(auto p:ring)g.positions.insert(g.positions.end(),{float(p.x),float(p.y),float(p.z)});
  for(std::size_t i=0;i<fields.size();++i){const auto& f=fields[i];for(std::size_t v=0;v<ring.size();++v)for(std::uint32_t c=0;c<f.components;++c)g.fields[i].values.push_back(f.values[std::size_t(cell)*f.components+c]);}
}
Geometry slice(const Snapshot& s,const Plane& plane,const kutta::CellFaces& cf,double eps) {
  const auto& m=s.mesh;Geometry g;
  for(const auto& f:s.cells)g.fields.push_back({f.name,f.unit,f.components,{}});
  std::vector<double> distance;distance.reserve(m.points.size());for(auto p:m.points)distance.push_back(dot(p-plane.origin,plane.normal));
  std::vector<bool> candidate(m.n_cells);
  for(Index f=0;f<m.n_faces();++f) {
    double lo=std::numeric_limits<double>::infinity(),hi=-lo;
    for(auto n:nodes(m,f)){lo=std::min(lo,distance[n]);hi=std::max(hi,distance[n]);}
    if(lo<=eps && hi>=-eps){candidate[m.owner[f]]=true;if(f<m.n_internal_faces())candidate[m.neighbour[f]]=true;}
  }
  for(Index cell=0;cell<m.n_cells;++cell)if(candidate[cell]) {
    bool coplanar=false;
    // A plane coincident with an existing face is emitted once, using its
    // owner cell's values. Neighbour copies would z-fight in the viewer.
    for(auto i=cf.offsets[cell];i<cf.offsets[cell+1];++i) {
      const auto f=cf.faces[i];const auto ns=nodes(m,f);
      if(std::all_of(ns.begin(),ns.end(),[&](Index n){return std::abs(distance[n])<=eps;})) {
        coplanar=true;if(m.owner[f]==cell){std::vector<Vec3d> p;for(auto n:ns)p.push_back(m.points[n]);append_polygon(g,std::move(p),plane,cell,s.cells);}
      }
    }
    if(coplanar)continue;
    bool positive=false,negative=false;
    for(auto i=cf.offsets[cell];i<cf.offsets[cell+1];++i)for(auto n:nodes(m,cf.faces[i])){positive |= distance[n]>eps;negative |= distance[n]<-eps;}
    if(!positive || !negative)continue; // tangent point/edge, zero section area
    std::vector<Vec3d> points;std::vector<std::array<std::uint32_t,2>> edges;
    const auto point=[&](Vec3d p){for(std::size_t i=0;i<points.size();++i)if(norm(points[i]-p)<=eps)return std::uint32_t(i);points.push_back(p);return std::uint32_t(points.size()-1);};
    for(auto i=cf.offsets[cell];i<cf.offsets[cell+1];++i) {
      const auto ns=nodes(m,cf.faces[i]);std::vector<Vec3d> p;for(auto n:ns)p.push_back(m.points[n]);
      const auto triangles=triangulate(p);
      for(std::size_t t=0;t<triangles.size();t+=3) {
        std::vector<std::uint32_t> hits;
        const auto hit=[&](Vec3d v){const auto id=point(v);if(std::find(hits.begin(),hits.end(),id)==hits.end())hits.push_back(id);};
        for(int j=0;j<3;++j) {
          const auto a=triangles[t+j],b=triangles[t+(j+1)%3];const auto da=distance[ns[a]],db=distance[ns[b]];
          if(std::abs(da)<=eps)hit(p[a]);
          if((da>eps && db<-eps) || (da<-eps && db>eps))hit(p[a]+(p[b]-p[a])*(da/(da-db)));
        }
        if(hits.size()==2){std::array edge{std::min(hits[0],hits[1]),std::max(hits[0],hits[1])};if(std::find(edges.begin(),edges.end(),edge)==edges.end())edges.push_back(edge);}
      }
    }
    std::vector<bool> used(edges.size());
    for(std::size_t first=0;first<edges.size();++first)if(!used[first]) {
      std::vector<Vec3d> ring;auto start=edges[first][0],current=start;std::size_t edge=first;
      do {
        ring.push_back(points[current]);used[edge]=true;current=edges[edge][0]==current?edges[edge][1]:edges[edge][0];
        if(current==start)break;
        edge=0;for(;edge<edges.size();++edge)if(!used[edge] && (edges[edge][0]==current || edges[edge][1]==current))break;
        KUTTA_REQUIRE(edge<edges.size(),"slice intersection is not a closed polygon");
      }while(ring.size()<=edges.size());
      KUTTA_REQUIRE(current==start,"slice polygon did not close");append_polygon(g,std::move(ring),plane,cell,s.cells);
    }
  }
  std::erase_if(g.fields,[](const auto& f){return !f.values.empty() && std::none_of(f.values.begin(),f.values.end(),[](float v){return std::isfinite(v);});});
  append_deltas(g.fields);return g;
}
struct Range {double lo=std::numeric_limits<double>::infinity(),hi=-lo;void add(double n){if(std::isfinite(n)){lo=std::min(lo,n);hi=std::max(hi,n);}}};
// Each mesh edge once, chained greedily into edge-disjoint polylines over the
// zone's own vertices (identical order), so drawing it shows polygon edges
// without triangulation diagonals in about one index per edge.
Geometry edges(const PolyMesh& m,const Patch& patch,std::vector<Index>& local) {
  Geometry g;std::vector<Index> vertices;std::vector<std::uint64_t> keys;
  for(Index f=patch.start;f<patch.start+patch.size;++f) {
    const auto ns=nodes(m,f);
    for(auto n:ns)if(local[n]<0){local[n]=index(vertices.size());vertices.push_back(n);}
    for(std::size_t i=0;i<ns.size();++i) {
      const auto a=std::uint64_t(local[ns[i]]),b=std::uint64_t(local[ns[(i+1)%ns.size()]]);
      if(a!=b)keys.push_back(std::min(a,b)<<32|std::max(a,b));
    }
  }
  for(auto n:vertices){const auto p=m.points[n];g.positions.insert(g.positions.end(),{float(p.x),float(p.y),float(p.z)});local[n]=-1;}
  std::sort(keys.begin(),keys.end());keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
  const auto nv=vertices.size();std::vector<std::uint32_t> offsets(nv+1),adjacent(keys.size()*2),cursor;
  for(auto k:keys){++offsets[(k>>32)+1];++offsets[(k&0xffffffffu)+1];}
  for(std::size_t v=0;v<nv;++v)offsets[v+1]+=offsets[v];
  cursor.assign(offsets.begin(),offsets.end()-1);
  for(std::uint32_t e=0;e<keys.size();++e){adjacent[cursor[keys[e]>>32]++]=e;adjacent[cursor[keys[e]&0xffffffffu]++]=e;}
  cursor.assign(offsets.begin(),offsets.end()-1);
  std::vector<bool> used(keys.size());g.line_offsets.push_back(0);
  const auto walk=[&](std::uint32_t v) {
    for(bool first=true;;first=false) {
      while(cursor[v]<offsets[v+1] && used[adjacent[cursor[v]]])++cursor[v];
      if(cursor[v]==offsets[v+1])break;
      const auto e=adjacent[cursor[v]++];used[e]=true;
      if(first)g.indices.push_back(v);
      v=std::uint32_t(keys[e]>>32)==v?std::uint32_t(keys[e]&0xffffffffu):std::uint32_t(keys[e]>>32);g.indices.push_back(v);
    }
    if(g.indices.size()!=g.line_offsets.back())g.line_offsets.push_back(std::uint32_t(g.indices.size()));
  };
  // Start at odd-degree vertices first: open chains end there, so fewer lines.
  for(std::uint32_t v=0;v<nv;++v)if((offsets[v+1]-offsets[v])%2)walk(v);
  for(std::uint32_t v=0;v<nv;++v)while(cursor[v]<offsets[v+1]){const auto before=g.indices.size();walk(v);if(g.indices.size()==before)break;}
  return g;
}
Json bounds(std::span<const Vec3d> points) {
  Vec3d lo=points.front(),hi=lo;for(auto p:points){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
  return Json::Array{Json::Array{lo.x,lo.y,lo.z},Json::Array{hi.x,hi.y,hi.z}};
}
}
Json::Object manifest_header(bool reference,const std::string& provenance,std::span<const Vec3d> points) {
  KUTTA_REQUIRE(!points.empty(),"view mesh has no points");
  return {{"v",1},{"format",Json::Object{{"name","kutta-view"},{"version",1},{"endian","little"},{"header_bytes",32}}},{"reference",reference},{"provenance",provenance},{"bounds",bounds(points)}};
}
Json::Object write_zone(const std::filesystem::path& dir,std::size_t ordinal,const PolyMesh& mesh,const Fields& fields,std::vector<Index>& local,bool with_edges) {
  const auto& patch=mesh.patches[ordinal];const auto name="zone-"+std::to_string(ordinal);
  auto zone=write_geometry(dir/(name+".surface.bin"),surface(mesh,fields,patch,local)).descriptor;
  zone.emplace("id",patch.source_zone_id);zone.emplace("name",patch.name);zone.emplace("type",patch.source_type);zone.emplace("faces",double(patch.size));
  if(with_edges)zone.emplace("edges",write_geometry(dir/(name+".edges.bin"),edges(mesh,patch,local)).descriptor);
  return zone;
}
double slice_tolerance(const PolyMesh& mesh) {
  Vec3d lo=mesh.points.front(),hi=lo;
  for(auto p:mesh.points){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
  // Include absolute-coordinate roundoff without discarding small meshes.
  return std::max(norm(hi-lo)*1e-10,std::max(norm(lo),norm(hi))*8*std::numeric_limits<double>::epsilon());
}
Json::Object write_slice(const std::filesystem::path& file,const Snapshot& s,const Plane& p,const kutta::CellFaces& cf,double eps) {
  auto descriptor=write_geometry(file,slice(s,p,cf,eps)).descriptor;
  descriptor.emplace("origin",Json::Array{p.origin.x,p.origin.y,p.origin.z});descriptor.emplace("normal",Json::Array{p.normal.x,p.normal.y,p.normal.z});
  return descriptor;
}
void write_manifest(const std::filesystem::path& path,Json::Object manifest) {
  // Ranges combine every geometry descriptor's own field ranges, so cached
  // and freshly computed parts yield the same catalog.
  std::map<std::string,std::pair<Json::Object,Range>> catalog;
  for(auto key:{"zones","slices","streamlines","isosurfaces","probes"}) {
    const auto it=manifest.find(key);if(it==manifest.end())continue;
    for(const auto& entry:std::get<Json::Array>(it->second.value))for(const auto& f:std::get<Json::Array>(std::get<Json::Object>(entry.value).at("fields").value)) {
      const auto& o=std::get<Json::Object>(f.value);auto& c=catalog[std::get<std::string>(o.at("name").value)];
      c.first={{"name",o.at("name")},{"unit",o.at("unit")},{"components",o.at("components")},{"location","vertex"}};
      for(auto bound:{"min","max"})if(const auto* n=std::get_if<double>(&o.at(bound).value))c.second.add(*n);
    }
  }
  Json::Array fields;
  for(auto& [name,c]:catalog){(void)name;c.first.emplace("min",std::isfinite(c.second.lo)?Json(c.second.lo):Json(nullptr));c.first.emplace("max",std::isfinite(c.second.hi)?Json(c.second.hi):Json(nullptr));fields.emplace_back(std::move(c.first));}
  manifest.insert_or_assign("fields",std::move(fields));
  std::ofstream out(path,std::ios::binary);out.exceptions(std::ios::failbit|std::ios::badbit);out<<dump_json(manifest)<<'\n';out.close();
}
Json::Object write_surface_view(const std::filesystem::path& directory,const PolyMesh& mesh,const Fields& fields,bool with_edges,Json::Object extra) {
  KUTTA_REQUIRE(!mesh.points.empty(),"surface view mesh has no points");
  for(const auto& field:fields.boundary)KUTTA_REQUIRE((field.components==1 || field.components==3) && field.values.size()==std::size_t(mesh.n_faces()-mesh.n_internal_faces())*field.components,"invalid surface field shape: "+field.name);
  KUTTA_REQUIRE(std::filesystem::create_directory(directory),"surface staging directory already exists");
  auto manifest=manifest_header(fields.reference,fields.provenance,mesh.points);
  std::vector<Index> local(mesh.points.size(),-1);Json::Array zones;
  for(std::size_t i=0;i<mesh.patches.size();++i)zones.emplace_back(write_zone(directory,i,mesh,fields,local,with_edges));
  manifest.emplace("zones",std::move(zones));manifest.emplace("slices",Json::Array{});
  for(auto& [key,value]:extra)manifest.insert_or_assign(key,std::move(value));
  write_manifest(directory/"manifest.json",manifest);
  return manifest;
}
std::filesystem::path staging(const std::filesystem::path& directory) {
  for(;;) {
    const auto stage=directory/(".view-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    if(std::filesystem::create_directory(stage))return stage;
  }
}
void publish(const std::filesystem::path& stage,const std::filesystem::path& directory) {
  const auto target=directory/"view",old=directory/(".view-old-"+stage.filename().string().substr(6));
  if(std::filesystem::exists(target))std::filesystem::rename(target,old);
  try {std::filesystem::rename(stage,target);}catch(...){if(std::filesystem::exists(old))std::filesystem::rename(old,target);throw;}
  std::filesystem::remove_all(old);
}
void export_view(const std::filesystem::path& directory,const Snapshot& s,std::span<const Plane> planes) {
  KUTTA_REQUIRE(!s.mesh.points.empty(),"view snapshot has no points");
  std::filesystem::create_directories(directory);const auto stage=staging(directory);
  try {
    auto manifest=manifest_header(s.reference,s.provenance,s.mesh.points);
    std::vector<Index> local(s.mesh.points.size(),-1);Json::Array zones,slices;
    for(std::size_t i=0;i<s.mesh.patches.size();++i)zones.emplace_back(write_zone(stage,i,s.mesh,s,local,true));
    if(!planes.empty()) {
      const auto cf=build_cell_faces(s.mesh);const auto eps=slice_tolerance(s.mesh);
      for(std::size_t i=0;i<planes.size();++i)slices.emplace_back(write_slice(stage/("slice-"+std::to_string(i)+".bin"),s,planes[i],cf,eps));
    }
    manifest.emplace("zones",std::move(zones));manifest.emplace("slices",std::move(slices));
    write_manifest(stage/"manifest.json",std::move(manifest));
    publish(stage,directory);
  }catch(...){std::error_code ec;std::filesystem::remove_all(stage,ec);throw;}
}
} // namespace kutta::view
