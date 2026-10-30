#include "view.hpp"
#include "binary.hpp"
#include "kutta/core/json.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>

namespace kutta::view {
namespace {
using namespace detail;
constexpr float missing=std::numeric_limits<float>::quiet_NaN();
struct Geometry {
  std::vector<float> positions;
  std::vector<std::uint32_t> indices;
  std::vector<Field> fields;
};
Vec3d polygon_area(std::span<const Vec3d> p) {
  Vec3d area;for(std::size_t i=1;i+1<p.size();++i)area=area+cross(p[i]-p[0],p[i+1]-p[0])*.5;return area;
}
// Ear clipping in the dominant projection preserves concave Fluent polygons.
// Convex faces take a linear fan path. Warped faces use the same projection;
// triangles retain original 3D vertices and outward winding.
std::vector<std::uint32_t> triangulate(std::span<const Vec3d> p) {
  KUTTA_REQUIRE(p.size()>=3,"view polygon has fewer than three vertices");
  const auto area=polygon_area(p);KUTTA_REQUIRE(norm(area)>0,"view polygon has zero area");
  int axis=0;if(std::abs(area.y)>std::abs(area.x))axis=1;if(std::abs(area.z)>(axis==0?std::abs(area.x):std::abs(area.y)))axis=2;
  const double direction=(axis==0?area.x:axis==1?area.y:area.z)>0?1:-1;
  const auto turn=[&](Vec3d a,Vec3d b,Vec3d c){const auto v=cross(b-a,c-a);return direction*(axis==0?v.x:axis==1?v.y:v.z);};
  double scale=0;for(auto v:p)scale=std::max(scale,norm(v-p[0]));const double eps=scale*scale*1e-13;
  bool convex=true;for(std::size_t i=0;i<p.size();++i)convex &= turn(p[i],p[(i+1)%p.size()],p[(i+2)%p.size()])>eps;
  std::vector<std::uint32_t> out;out.reserve((p.size()-2)*3);
  if(convex){for(std::uint32_t i=1;i+1<p.size();++i)out.insert(out.end(),{0,i,i+1});return out;}
  std::vector<std::uint32_t> ring(p.size());std::iota(ring.begin(),ring.end(),0);
  // Slice triangles can introduce collinear points along face diagonals.
  for(std::size_t i=0;ring.size()>3 && i<ring.size();) {
    const auto a=ring[(i+ring.size()-1)%ring.size()],b=ring[i],c=ring[(i+1)%ring.size()];
    if(std::abs(turn(p[a],p[b],p[c]))<=eps){ring.erase(ring.begin()+std::ptrdiff_t(i));i=0;}else ++i;
  }
  while(ring.size()>3) {
    bool found=false;
    for(std::size_t i=0;i<ring.size();++i) {
      const auto a=ring[(i+ring.size()-1)%ring.size()],b=ring[i],c=ring[(i+1)%ring.size()];
      if(turn(p[a],p[b],p[c])<=eps)continue;
      bool inside=false;for(auto j:ring)if(j!=a && j!=b && j!=c && turn(p[a],p[b],p[j])>=-eps && turn(p[b],p[c],p[j])>=-eps && turn(p[c],p[a],p[j])>=-eps){inside=true;break;}
      if(inside)continue;
      out.insert(out.end(),{a,b,c});ring.erase(ring.begin()+std::ptrdiff_t(i));found=true;break;
    }
    KUTTA_REQUIRE(found,"cannot triangulate self-intersecting or degenerate view polygon");
  }
  out.insert(out.end(),ring.begin(),ring.end());return out;
}
void add_deltas(std::vector<Field>& fields) {
  const auto count=fields.size();
  for(std::size_t i=0;i<count;++i) {
    const auto& f=fields[i];if(!f.name.ends_with("@fluent"))continue;
    const auto name=f.name.substr(0,f.name.size()-7);
    auto base=std::find_if(fields.begin(),fields.begin()+std::ptrdiff_t(count),[&](const auto& b){return b.name==name;});
    if(base==fields.begin()+std::ptrdiff_t(count))continue;
    Field d{name+"@delta",f.unit,f.components,{}};d.values.reserve(f.values.size());
    for(std::size_t j=0;j<f.values.size();++j)d.values.push_back(base->values[j]-f.values[j]);fields.push_back(std::move(d));
  }
}
bool visible(std::string_view name,bool wall) {
  const auto at=name.find('@');name=name.substr(0,at);
  return name=="pressure" || name=="pressure_coefficient" ||
    (wall?(name=="wall_shear" || name=="wall_shear_magnitude" || name=="y_plus"):(name=="velocity" || name=="velocity_magnitude"));
}
Geometry surface(const Snapshot& s,const Patch& patch,std::vector<Index>& local) {
  const auto& m=s.mesh;Geometry g;std::vector<Index> vertices;std::vector<double> area(patch.size);
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
  add_deltas(g.fields);return g;
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
  add_deltas(g.fields);return g;
}
struct Range {double lo=std::numeric_limits<double>::infinity(),hi=-lo;void add(double n){if(std::isfinite(n)){lo=std::min(lo,n);hi=std::max(hi,n);}}};
Range range(const Field& f) {
  Range out;for(std::size_t i=0;i<f.values.size();i+=f.components) {
    const auto n=f.components==1?double(f.values[i]):std::hypot(double(f.values[i]),double(f.values[i+1]),double(f.values[i+2]));out.add(n);
  }
  return out;
}
void range_json(std::ostream& out,Range r) {out<<",\"min\":";json_number(out,r.lo);out<<",\"max\":";json_number(out,r.hi);}
void bounds_json(std::ostream& out,std::span<const Vec3d> points) {
  Vec3d lo=points.front(),hi=lo;for(auto p:points){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
  out<<"[["<<lo.x<<','<<lo.y<<','<<lo.z<<"],["<<hi.x<<','<<hi.y<<','<<hi.z<<"]]";
}
void payload(const std::filesystem::path& path,const Geometry& g,std::ostream& manifest,std::map<std::string,std::pair<Field,Range>>& catalog) {
  KUTTA_REQUIRE(g.positions.size()/3<=UINT32_MAX && g.indices.size()<=UINT32_MAX && g.fields.size()<=UINT32_MAX,"view exceeds uint32 capacity");
  binary::Writer w(path);w.scalar(std::uint32_t(0x5654554b));w.scalar(std::uint32_t(1));w.scalar(std::uint32_t(32));
  w.scalar(std::uint32_t(g.positions.size()/3));w.scalar(std::uint32_t(g.indices.size()));w.scalar(std::uint32_t(g.fields.size()));w.scalar(std::uint32_t(0));w.scalar(std::uint32_t(0));
  for(float p:g.positions)KUTTA_REQUIRE(std::isfinite(p),"position cannot be represented as float32");
  w.array<float>(g.positions);w.array<std::uint32_t>(g.indices);
  manifest<<",\"file\":"<<json_string(path.filename().string())<<",\"vertices\":"<<g.positions.size()/3<<",\"indices\":"<<g.indices.size()<<",\"fields\":[";
  std::uint64_t offset=32+g.positions.size()*4+g.indices.size()*4;bool first=true;
  for(const auto& f:g.fields) {
    KUTTA_REQUIRE(f.values.size()==g.positions.size()/3*f.components,"view vertex field shape mismatch");
    w.array<float>(f.values);auto r=range(f);
    if(!first)manifest<<',';first=false;
    manifest<<"{\"name\":"<<json_string(f.name)<<",\"unit\":"<<json_string(f.unit)<<",\"components\":"<<f.components<<",\"location\":\"vertex\",\"offset_bytes\":"<<offset;range_json(manifest,r);manifest<<'}';
    offset+=f.values.size()*4;
    auto& entry=catalog[f.name];entry.first.name=f.name;entry.first.unit=f.unit;entry.first.components=f.components;entry.second.add(r.lo);entry.second.add(r.hi);
  }
  manifest<<"]}";w.close();
}
}
void export_view(const std::filesystem::path& directory,const Snapshot& s,std::span<const Plane> planes) {
  KUTTA_REQUIRE(!s.mesh.points.empty(),"view snapshot has no points");
  const auto token=std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
  const auto stage=directory/(".view-"+token),target=directory/"view",old=directory/(".view-old-"+token);
  std::filesystem::create_directories(stage);
  try {
    std::ofstream manifest(stage/"manifest.json");manifest.exceptions(std::ios::failbit|std::ios::badbit);manifest.imbue(std::locale::classic());manifest<<std::setprecision(17);
    manifest<<"{\"v\":1,\"format\":{\"name\":\"kutta-view\",\"version\":1,\"endian\":\"little\",\"header_bytes\":32},\"reference\":"<<(s.reference?"true":"false")<<",\"provenance\":"<<json_string(s.provenance)<<",\"bounds\":";bounds_json(manifest,s.mesh.points);
    std::map<std::string,std::pair<Field,Range>> catalog;std::vector<Index> local(s.mesh.points.size(),-1);
    manifest<<",\"zones\":[";bool first=true;
    for(std::size_t i=0;i<s.mesh.patches.size();++i) {
      const auto& p=s.mesh.patches[i];auto g=surface(s,p,local);
      if(!first)manifest<<',';first=false;manifest<<"{\"id\":"<<p.source_zone_id<<",\"name\":"<<json_string(p.name)<<",\"type\":"<<json_string(p.source_type)<<",\"faces\":"<<p.size;
      payload(stage/("zone-"+std::to_string(i)+".surface.bin"),g,manifest,catalog);
    }
    manifest<<"],\"slices\":[";
    if(!planes.empty()) {
      const auto cf=build_cell_faces(s.mesh);Vec3d lo=s.mesh.points.front(),hi=lo;
      for(auto p:s.mesh.points){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
      // Include absolute-coordinate roundoff without discarding small meshes.
      const double eps=std::max(norm(hi-lo)*1e-10,std::max(norm(lo),norm(hi))*8*std::numeric_limits<double>::epsilon());
      for(std::size_t i=0;i<planes.size();++i) {
        auto g=slice(s,planes[i],cf,eps);const auto& p=planes[i];if(i)manifest<<',';
        manifest<<"{\"origin\":["<<p.origin.x<<','<<p.origin.y<<','<<p.origin.z<<"],\"normal\":["<<p.normal.x<<','<<p.normal.y<<','<<p.normal.z<<']';
        payload(stage/("slice-"+std::to_string(i)+".bin"),g,manifest,catalog);
      }
    }
    manifest<<"],\"fields\":[";first=true;
    for(const auto& [name,entry]:catalog){if(!first)manifest<<',';first=false;manifest<<"{\"name\":"<<json_string(name)<<",\"unit\":"<<json_string(entry.first.unit)<<",\"components\":"<<entry.first.components<<",\"location\":\"vertex\"";range_json(manifest,entry.second);manifest<<'}';}
    manifest<<"]}\n";manifest.close();
    if(std::filesystem::exists(target))std::filesystem::rename(target,old);
    try {std::filesystem::rename(stage,target);}catch(...){if(std::filesystem::exists(old))std::filesystem::rename(old,target);throw;}
    std::filesystem::remove_all(old);
  }catch(...){std::error_code ec;std::filesystem::remove_all(stage,ec);throw;}
}
} // namespace kutta::view
