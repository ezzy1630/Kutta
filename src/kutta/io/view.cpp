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
#include <optional>

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
// Triangulates before touching g, so a ring that cannot be triangulated
// leaves the geometry unchanged and returns false.
bool append_polygon(Geometry& g,std::vector<Vec3d> ring,const Plane& plane,Index cell,std::span<const Field> fields) {
  if(ring.size()<3 || norm(polygon_area(ring))==0)return true;
  if(dot(polygon_area(ring),plane.normal)<0)std::reverse(ring.begin(),ring.end());
  const auto base=std::uint32_t(g.positions.size()/3);
  KUTTA_REQUIRE(std::uint64_t(base)+ring.size()<=UINT32_MAX,"view exceeds uint32 vertices");
  std::vector<std::uint32_t> triangles;
  try {triangles=triangulate(ring);}catch(const Error&){return false;}
  for(auto i:triangles)g.indices.push_back(base+i);
  for(auto p:ring)g.positions.insert(g.positions.end(),{float(p.x),float(p.y),float(p.z)});
  for(std::size_t i=0;i<fields.size();++i){const auto& f=fields[i];for(std::size_t v=0;v<ring.size();++v)for(std::uint32_t c=0;c<f.components;++c)g.fields[i].values.push_back(f.values[std::size_t(cell)*f.components+c]);}
  return true;
}
// Section of one cell by the plane as closed rings, or nullopt when the cell's
// faces do not bound a closed surface (an edge not shared by exactly two of its
// faces). Each node is strictly above the plane (distance > eps) or not; a
// node within eps counts as below. A face of more than three nodes is cut as
// the triangle fan about its node average, the solver's face geometry, so the
// cut follows warped and concave faces (2:1 hanging faces whose mid-edge node
// is not collinear) instead of chording across them, and both cells sharing a
// face cut the same surface. Crossings lie only where a fan triangle's edge
// joins an above and a below point and are identified by that edge (mesh edge
// by node pair, spoke by face and node), so on a closed cell surface every
// crossing has exactly two segments and the rings always close, whatever the
// plane does at vertices, on edges or near the faces themselves.
std::optional<std::vector<std::vector<Vec3d>>> section(const PolyMesh& m,const kutta::CellFaces& cf,Index cell,std::span<const double> distance,double eps) {
  std::vector<std::uint64_t> keys;std::vector<Vec3d> points;std::vector<std::array<std::uint32_t,2>> link;
  const auto crossing=[&](std::uint64_t key,Vec3d pa,double da,Vec3d pb,double db) {
    for(std::size_t i=0;i<keys.size();++i)if(keys[i]==key)return std::uint32_t(i);
    const auto t=std::clamp(da/(da-db),0.,1.);
    keys.push_back(key);points.push_back(pa+(pb-pa)*t);link.push_back({UINT32_MAX,UINT32_MAX});
    return std::uint32_t(points.size()-1);
  };
  const auto edge=[&](Index a,Index b){if(a>b)std::swap(a,b);return crossing(std::uint64_t(a)<<32|std::uint64_t(b),m.points[a],distance[a],m.points[b],distance[b]);};
  bool closed=true;
  const auto connect=[&](std::uint32_t a,std::uint32_t b) {
    for(auto [from,to]:{std::pair{a,b},std::pair{b,a}}) {
      auto& l=link[from];if(l[0]==UINT32_MAX)l[0]=to;else if(l[1]==UINT32_MAX)l[1]=to;else closed=false;
    }
  };
  for(auto i=cf.offsets[cell];i<cf.offsets[cell+1];++i) {
    const auto f=cf.faces[i];const auto ns=nodes(m,f);
    bool any_above=false,any_below=false;for(auto n:ns){any_above |= distance[n]>eps;any_below |= distance[n]<=eps;}
    if(!any_above || !any_below)continue;
    if(ns.size()==3) {
      std::uint32_t hit[2];int count=0;
      for(std::size_t j=0;j<3;++j){const auto a=ns[j],b=ns[(j+1)%3];if((distance[a]>eps)!=(distance[b]>eps))hit[count++]=edge(a,b);}
      connect(hit[0],hit[1]);continue;
    }
    // The fan centre's distance is the mean of its nodes', exactly as the
    // centre is their mean, so it never leaves their range.
    Vec3d centre;double dc=0;for(auto n:ns){centre=centre+m.points[n];dc+=distance[n];}
    centre=centre*(1./double(ns.size()));dc/=double(ns.size());const bool above_c=dc>eps;
    const auto spoke=[&](Index n){return crossing(std::uint64_t(1)<<63|std::uint64_t(f)<<32|std::uint64_t(n),centre,dc,m.points[n],distance[n]);};
    for(std::size_t j=0;j<ns.size();++j) {
      const auto a=ns[j],b=ns[(j+1)%ns.size()];const bool above_a=distance[a]>eps,above_b=distance[b]>eps;
      if(above_a==above_b && above_a==above_c)continue;
      std::uint32_t hit[2];int count=0;
      if(above_a!=above_b)hit[count++]=edge(a,b);
      if(above_b!=above_c)hit[count++]=spoke(b);
      if(above_c!=above_a)hit[count++]=spoke(a);
      connect(hit[0],hit[1]);
    }
  }
  if(!closed)return std::nullopt;
  for(const auto& l:link)if(l[1]==UINT32_MAX)return std::nullopt;
  std::vector<std::vector<Vec3d>> rings;std::vector<bool> used(points.size());
  for(std::uint32_t start=0;start<points.size();++start)if(!used[start]) {
    std::vector<Vec3d> ring;std::uint32_t previous=UINT32_MAX,current=start;
    do {
      used[current]=true;
      // Crossings at t=0 or 1 on edges sharing a node coincide; keep one.
      if(ring.empty() || norm(ring.back()-points[current])>eps)ring.push_back(points[current]);
      const auto next=link[current][0]==previous?link[current][1]:link[current][0];
      previous=current;current=next;
    }while(current!=start && !used[current]);
    if(current!=start)return std::nullopt;
    while(ring.size()>1 && norm(ring.back()-ring.front())<=eps)ring.pop_back();
    // A ring through the same point twice (the fan centre of a face lying
    // almost in the plane, between two above-plane parts of that face)
    // splits there into simple rings.
    std::vector<std::vector<Vec3d>> pending{std::move(ring)};
    while(!pending.empty()) {
      auto r=std::move(pending.back());pending.pop_back();bool split=false;
      for(std::size_t a=0;a<r.size() && !split;++a)for(std::size_t b=a+1;b<r.size();++b)if(norm(r[a]-r[b])<=eps) {
        std::vector<Vec3d> inner(r.begin()+std::ptrdiff_t(a),r.begin()+std::ptrdiff_t(b)),outer(r.begin()+std::ptrdiff_t(b),r.end());
        outer.insert(outer.end(),r.begin(),r.begin()+std::ptrdiff_t(a));
        pending.push_back(std::move(inner));pending.push_back(std::move(outer));split=true;break;
      }
      if(split)continue;
      // Spoke crossings on a planar face lie on the chord between its edge
      // crossings; drop such collinear points so slices stay as small as
      // before (this also removes zero-width spikes).
      for(std::size_t j=0;r.size()>=3 && j<r.size();) {
        const auto& p0=r[(j+r.size()-1)%r.size()];const auto& p1=r[j];const auto& p2=r[(j+1)%r.size()];
        if(norm(cross(p1-p0,p2-p1))<=1e-9*norm(p1-p0)*norm(p2-p1)){r.erase(r.begin()+std::ptrdiff_t(j));j=0;}else ++j;
      }
      if(r.size()>=3)rings.push_back(std::move(r));
    }
  }
  return rings;
}
template<class Append>
std::uint64_t slice_polygons(const PolyMesh& m,const Plane& plane,const kutta::CellFaces& cf,double eps,
                             Append append) {
  std::uint64_t skipped=0;
  std::vector<double> distance;distance.reserve(m.points.size());for(auto p:m.points)distance.push_back(dot(p-plane.origin,plane.normal));
  std::vector<bool> candidate(m.n_cells);
  for(Index f=0;f<m.n_faces();++f) {
    double lo=std::numeric_limits<double>::infinity(),hi=-lo;
    for(auto n:nodes(m,f)){lo=std::min(lo,distance[n]);hi=std::max(hi,distance[n]);}
    if(lo>eps || hi<-eps) continue;
    candidate[m.owner[f]]=true;
    if(f<m.n_internal_faces())candidate[m.neighbour[f]]=true;
  }
  std::vector<Index> coplanar;
  for(Index cell=0;cell<m.n_cells;++cell)if(candidate[cell]) {
    // Faces lying in the plane (every node within eps).
    coplanar.clear();
    for(auto i=cf.offsets[cell];i<cf.offsets[cell+1];++i) {
      const auto ns=nodes(m,cf.faces[i]);
      if(std::all_of(ns.begin(),ns.end(),[&](Index n){return std::abs(distance[n])<=eps;}))coplanar.push_back(cf.faces[i]);
    }
    const auto face_polygon=[&](Index f){std::vector<Vec3d> p;for(auto n:nodes(m,f))p.push_back(m.points[n]);return p;};
    bool above=false,below=false,ok=true;
    for(auto i=cf.offsets[cell];i<cf.offsets[cell+1];++i)for(auto n:nodes(m,cf.faces[i])){above |= distance[n]>eps;below |= distance[n]<=eps;}
    if(!above) {
      // Wholly on or below the plane: only a boundary face lying in the
      // plane has no cell above it to emit it.
      for(auto f:coplanar)if(f>=m.n_internal_faces())ok &= append(face_polygon(f),cell);
      skipped+=!ok;continue;
    }
    if(!below)continue;
    const auto rings=section(m,cf,cell,distance,eps);
    if(!rings){++skipped;continue;}
    std::vector<const std::vector<Vec3d>*> kept;double section_area=0;
    for(const auto& ring:*rings) {
      // A cell touching the plane only along a vertex or an edge (within eps)
      // yields a ring of zero width: no section area, nothing to draw.
      double diameter=0;for(const auto& p:ring)diameter=std::max(diameter,norm(p-ring.front()));
      const auto area=norm(polygon_area(ring));
      if(area>eps*diameter){kept.push_back(&ring);section_area+=area;}
    }
    // When the section is exactly this cell's faces lying in the plane, emit
    // those faces with their owner's values, once: a plane coincident with an
    // existing face keeps its owner's values, as before, and neighbour copies
    // would z-fight in the viewer. Otherwise (a face only partly within eps,
    // as snapped meshes have) the cell's own section is the consistent cut.
    double face_area=0;for(auto f:coplanar)face_area+=norm(polygon_area(face_polygon(f)));
    if(!coplanar.empty() && std::abs(face_area-section_area)<=1e-6*std::max(face_area,section_area))
      for(auto f:coplanar)ok &= append(face_polygon(f),m.owner[f]);
    else for(const auto* ring:kept)ok &= append(*ring,cell);
    skipped+=!ok;
  }
  return skipped;
}
Geometry slice(const Snapshot& s,const Plane& plane,const kutta::CellFaces& cf,double eps,std::uint64_t& skipped) {
  Geometry g;
  for(const auto& f:s.cells)g.fields.push_back({f.name,f.unit,f.components,{}});
  skipped=slice_polygons(s.mesh,plane,cf,eps,[&](const auto& ring,Index cell) {
    return append_polygon(g,ring,plane,cell,s.cells);
  });
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
  std::uint64_t skipped=0;auto descriptor=write_geometry(file,slice(s,p,cf,eps,skipped)).descriptor;
  // Cells whose faces do not close, or whose section cannot be triangulated,
  // are left out instead of failing the whole export (v1 additive).
  descriptor.emplace("skipped_cells",double(skipped));
  descriptor.emplace("origin",Json::Array{p.origin.x,p.origin.y,p.origin.z});descriptor.emplace("normal",Json::Array{p.normal.x,p.normal.y,p.normal.z});
  return descriptor;
}
Json::Object mesh_gap_section(const PolyMesh& mesh,Vec3d start,Vec3d end) {
  const auto gap=end-start;
  const auto width=norm(gap);
  KUTTA_REQUIRE(width>0 && std::isfinite(width),"mesh section requires distinct endpoints");
  auto vertical=gap/width;
  if(vertical.z<0)vertical=vertical*(-1);
  auto normal=cross(vertical,Vec3d{1,0,0});
  if(norm(normal)<0.1)normal=cross(vertical,Vec3d{0,1,0});
  normal=normal/norm(normal);
  const auto horizontal=cross(normal,vertical);
  const Plane plane{start+gap*.5,normal};
  const double radius=width*4;
  Json::Array polygons;
  std::size_t points=0;
  const auto skipped=slice_polygons(mesh,plane,build_cell_faces(mesh),slice_tolerance(mesh),
    [&](const auto& ring,Index) {
      const auto project=[&](Vec3d point) {
        const auto local=point-plane.origin;
        return std::pair{dot(local,horizontal),dot(local,vertical)};
      };
      auto [min_x,min_y]=project(ring.front());
      auto max_x=min_x,max_y=min_y;
      for(const auto p:ring) {
        const auto [x,y]=project(p);
        min_x=std::min(min_x,x);max_x=std::max(max_x,x);
        min_y=std::min(min_y,y);max_y=std::max(max_y,y);
      }
      if(min_x>radius || max_x<-radius || min_y>radius || max_y<-radius) return true;
      points+=ring.size();
      KUTTA_REQUIRE(points<=80'000,"section contains too many vertices to display");
      Json::Array polygon;
      for(const auto p:ring) {
        const auto [x,y]=project(p);
        polygon.emplace_back(Json::Array{x,y});
      }
      polygons.emplace_back(std::move(polygon));
      return true;
    });
  return {{"v",1},{"event","mesh-section"},{"radius_m",radius},{"width_m",width},
          {"polygons",std::move(polygons)},{"skipped_cells",double(skipped)}};
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
