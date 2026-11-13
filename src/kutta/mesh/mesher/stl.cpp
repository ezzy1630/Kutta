#include "external_aero.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/core/path.hpp"
#include <bit>
#include <fstream>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>

namespace kutta::mesher {
using namespace detail;
namespace {
struct Raw { Vec3d a,b,c; Index zone; };
double coord(Vec3d p,int a) { return a==0?p.x:a==1?p.y:p.z; }
std::uint32_t u32(const unsigned char* b) { return std::uint32_t(b[0]) | (std::uint32_t(b[1])<<8) | (std::uint32_t(b[2])<<16) | (std::uint32_t(b[3])<<24); }
Vec3d point(const unsigned char* p,double scale) { return {std::bit_cast<float>(u32(p))*scale,std::bit_cast<float>(u32(p+4))*scale,std::bit_cast<float>(u32(p+8))*scale}; }
Box empty_box() { const auto inf=std::numeric_limits<double>::infinity(); return {{inf,inf,inf},{-inf,-inf,-inf}}; }
void extend(Box& b,Vec3d p) { b.min={std::min(b.min.x,p.x),std::min(b.min.y,p.y),std::min(b.min.z,p.z)}; b.max={std::max(b.max.x,p.x),std::max(b.max.y,p.y),std::max(b.max.z,p.z)}; }
std::string clean(std::string s) {
  const auto first=s.find_first_not_of(" \t\r\n"); if(first==std::string::npos) return {};
  return s.substr(first,s.find_last_not_of(" \t\r\n")-first+1);
}
}
Surface read_stl(std::span<const std::filesystem::path> paths,const Settings& settings) {
  KUTTA_REQUIRE(!paths.empty(),"STL: at least one input is required");
  Surface out; std::vector<Raw> raw; std::set<std::string> names;
  const auto scale=settings.units=="mm"?0.001:1.0;
  auto zone=[&](std::string name) {
    KUTTA_REQUIRE(!name.empty(),"STL: empty file/solid name");
    const auto base=name; int i=2; while(!names.insert(name).second) name=base+"-"+std::to_string(i++);
    out.zones.push_back(name); return index(out.zones.size()-1);
  };
  for(const auto& path:paths) {
    std::ifstream in(path,std::ios::binary); KUTTA_REQUIRE(in,"STL: cannot open "+path.string());
    KUTTA_REQUIRE(std::filesystem::file_size(path)<=std::uintmax_t(1024)*1024*1024,"STL: input exceeds 1 GiB");
    std::vector<unsigned char> bytes((std::istreambuf_iterator<char>(in)),{});
    KUTTA_REQUIRE(bytes.size()<=std::size_t(1024)*1024*1024,"STL: input exceeds 1 GiB");
    if(bytes.size()>=84 && 84ULL+50ULL*u32(bytes.data()+80)==bytes.size()) {
      const auto z=zone(utf8_path(path.stem()));
      for(std::size_t p=84;p<bytes.size();p+=50) raw.push_back({point(bytes.data()+p+12,scale),point(bytes.data()+p+24,scale),point(bytes.data()+p+36,scale),z});
    } else {
      std::istringstream text(std::string(bytes.begin(),bytes.end())); text.imbue(std::locale::classic());
      std::string token; bool open=false; Index z=-1;
      auto expect=[&](std::string_view wanted) { KUTTA_REQUIRE(bool(text>>token)&&token==wanted,"STL "+path.filename().string()+": expected "+std::string(wanted)); };
      auto real=[&]() { std::string v; KUTTA_REQUIRE(bool(text>>v),"STL: truncated number"); const auto d=parse_double(v); KUTTA_REQUIRE(d && std::isfinite(*d),"STL: invalid number "+v); return *d; };
      while(text>>token) {
        if(token=="solid") { KUTTA_REQUIRE(!open,"STL: nested solid"); std::string name; std::getline(text,name); name=clean(name); z=zone(name.empty()?utf8_path(path.stem()):name); open=true; }
        else if(token=="endsolid") { KUTTA_REQUIRE(open,"STL: endsolid without solid"); std::getline(text,token); open=false; }
        else {
          KUTTA_REQUIRE(open && token=="facet","STL: expected solid or facet");
          expect("normal"); (void)real(); (void)real(); (void)real(); expect("outer"); expect("loop");
          std::array<Vec3d,3> p;
          for(auto& v:p) { expect("vertex"); v={real()*scale,real()*scale,real()*scale}; }
          expect("endloop"); expect("endfacet"); raw.push_back({p[0],p[1],p[2],z});
        }
      }
      KUTTA_REQUIRE(!open,"STL: missing endsolid");
    }
  }
  KUTTA_REQUIRE(!raw.empty(),"STL: no triangles");
  out.bounds=empty_box(); for(const auto& t:raw) for(auto p:{t.a,t.b,t.c}) { KUTTA_REQUIRE(finite(p),"STL: nonfinite coordinate"); extend(out.bounds,p); }
  const auto diagonal=norm(out.bounds.max-out.bounds.min);
  KUTTA_REQUIRE(diagonal>0 && diagonal<1e7,"STL: invalid body bounds; check units");
  const double tolerance=diagonal*1e-7;
  using Key=std::array<std::int64_t,3>;
  std::map<Key,std::vector<Index>> bins;
  std::size_t welded=0,degenerate=0,duplicates=0;
  auto vertex=[&](Vec3d p) {
    const auto r=(p-out.bounds.min)/tolerance;
    Key k{std::int64_t(std::floor(r.x)),std::int64_t(std::floor(r.y)),std::int64_t(std::floor(r.z))};
    for(int a=-1;a<=1;++a) for(int b=-1;b<=1;++b) for(int c=-1;c<=1;++c) {
      const auto it=bins.find({k[0]+a,k[1]+b,k[2]+c}); if(it==bins.end()) continue;
      for(auto v:it->second) if(norm(out.points[v]-p)<=tolerance) { ++welded; return v; }
    }
    const auto v=index(out.points.size()); out.points.push_back(p); bins[k].push_back(v); return v;
  };
  std::set<std::array<Index,4>> unique;
  for(const auto& t:raw) {
    Triangle tri{{vertex(t.a),vertex(t.b),vertex(t.c)},t.zone,0,{}};
    const auto& n=tri.nodes;
    if(n[0]==n[1] || n[1]==n[2] || n[2]==n[0] || norm(cross(t.b-t.a,t.c-t.a))<=tolerance*tolerance) { ++degenerate; continue; }
    auto sorted=n; std::sort(sorted.begin(),sorted.end());
    if(!unique.insert({t.zone,sorted[0],sorted[1],sorted[2]}).second) { ++duplicates; continue; }
    out.triangles.push_back(tri);
  }
  KUTTA_REQUIRE(!out.triangles.empty(),"STL: all triangles degenerate");
  using Edge=std::array<Index,3>; // zone, lower node, upper node
  std::map<Edge,std::vector<std::pair<Index,bool>>> edges;
  for(Index t=0;t<index(out.triangles.size());++t) for(int e=0;e<3;++e) {
    const auto& tri=out.triangles[t]; auto a=tri.nodes[e],b=tri.nodes[(e+1)%3];
    edges[{tri.zone,std::min(a,b),std::max(a,b)}].push_back({t,a<b});
  }
  std::vector<std::vector<std::pair<Index,bool>>> adj(out.triangles.size());
  std::size_t boundary=0,nonmanifold=0;
  for(const auto& [edge,uses]:edges) {
    (void)edge;
    if(uses.size()==1) ++boundary;
    else if(uses.size()!=2) ++nonmanifold;
    if(uses.size()==2) { const auto a=uses[0],b=uses[1]; adj[a.first].push_back({b.first,a.second==b.second}); adj[b.first].push_back({a.first,a.second==b.second}); }
  }
  std::vector<int> flip(out.triangles.size(),-1); Index components=0; std::size_t reversed=0,conflicts=0;
  for(Index start=0;start<index(out.triangles.size());++start) if(flip[start]<0) {
    std::vector<Index> members{start}; flip[start]=0;
    for(std::size_t q=0;q<members.size();++q) for(const auto [b,change]:adj[members[q]]) {
      const auto want=flip[members[q]]^int(change);
      if(flip[b]<0) { flip[b]=want; members.push_back(b); } else if(flip[b]!=want) ++conflicts;
    }
    double volume=0;
    for(auto i:members) { const auto& t=out.triangles[i]; const auto a=out.points[t.nodes[0]]-out.bounds.min,b=out.points[t.nodes[1]]-out.bounds.min,c=out.points[t.nodes[2]]-out.bounds.min; volume+=dot(a,cross(b,c))*(flip[i]?-1:1); }
    for(auto i:members) { auto& t=out.triangles[i]; if(bool(flip[i])^(volume<0)) { std::swap(t.nodes[1],t.nodes[2]); ++reversed; } t.component=components; }
    ++components;
  }
  const auto feature_cos=std::cos(settings.curvature_degrees*std::acos(-1.)/180);
  auto normal=[&](Index i) { const auto& n=out.triangles[i].nodes; const auto v=cross(out.points[n[1]]-out.points[n[0]],out.points[n[2]]-out.points[n[0]]); return v/norm(v); };
  for(const auto& [edge,uses]:edges) {
    (void)edge;
    const bool feature=uses.size()!=2 || dot(normal(uses[0].first),normal(uses[1].first))<feature_cos;
    if(feature) for(const auto& [t,dir]:uses) {
      (void)dir; auto& tri=out.triangles[t];
      for(int e=0;e<3;++e) if(std::min(tri.nodes[e],tri.nodes[(e+1)%3])==edge[1] && std::max(tri.nodes[e],tri.nodes[(e+1)%3])==edge[2]) tri.feature_edges[e]=true;
    }
  }
  out.repairs={{"input_triangles",double(raw.size())},{"triangles",double(out.triangles.size())},{"welded_vertex_occurrences",double(welded)},
    {"weld_tolerance_m",tolerance},{"removed_degenerate_triangles",double(degenerate)},{"removed_duplicate_triangles",double(duplicates)},
    {"reoriented_triangles",double(reversed)},{"open_edges",double(boundary)},{"nonmanifold_edges",double(nonmanifold)},
    {"orientation_conflicts",double(conflicts/2)},{"components",double(components)},
    {"hole_filling",false},{"inside_method","union of connected shells; majority of three ray parities (no flood fill)"}};
  return out;
}

namespace {
double box_distance(Box b,Vec3d p) {
  const Vec3d d{std::max({b.min.x-p.x,0.,p.x-b.max.x}),std::max({b.min.y-p.y,0.,p.y-b.max.y}),std::max({b.min.z-p.z,0.,p.z-b.max.z})}; return dot(d,d);
}
// Closest point regions of a triangle, including its edges and vertices.
double triangle_distance(Vec3d p,Vec3d a,Vec3d b,Vec3d c) {
  const auto ab=b-a,ac=c-a,ap=p-a; const auto d1=dot(ab,ap),d2=dot(ac,ap);
  if(d1<=0 && d2<=0) return dot(ap,ap);
  const auto bp=p-b; const auto d3=dot(ab,bp),d4=dot(ac,bp);
  if(d3>=0 && d4<=d3) return dot(bp,bp);
  const auto vc=d1*d4-d3*d2;
  if(vc<=0 && d1>=0 && d3<=0) { const auto d=ap-ab*(d1/(d1-d3)); return dot(d,d); }
  const auto cp=p-c; const auto d5=dot(ab,cp),d6=dot(ac,cp);
  if(d6>=0 && d5<=d6) return dot(cp,cp);
  const auto vb=d5*d2-d1*d6;
  if(vb<=0 && d2>=0 && d6<=0) { const auto d=ap-ac*(d2/(d2-d6)); return dot(d,d); }
  const auto va=d3*d6-d5*d4;
  if(va<=0 && d4-d3>=0 && d5-d6>=0) { const auto d=bp-(c-b)*((d4-d3)/(d4-d3+d5-d6)); return dot(d,d); }
  const auto n=cross(ab,ac); const auto d=dot(ap,n); return d*d/dot(n,n);
}
bool ray_box(Box b,Vec3d p,Vec3d d) {
  double lo=0,hi=std::numeric_limits<double>::infinity();
  for(int a=0;a<3;++a) { auto t0=(coord(b.min,a)-coord(p,a))/coord(d,a),t1=(coord(b.max,a)-coord(p,a))/coord(d,a); if(t0>t1) std::swap(t0,t1); lo=std::max(lo,t0); hi=std::min(hi,t1); }
  return hi>=lo;
}
}
struct SurfaceQuery::Impl {
  const Surface& s;
  struct Node { Box b; int first=0,count=0,left=-1,right=-1; };
  std::vector<Node> nodes;
  std::vector<Index> order;
  explicit Impl(const Surface& surface):s(surface) { order.resize(s.triangles.size()); std::iota(order.begin(),order.end(),0); build(0,int(order.size())); }
  Box bounds(Index t) const { Box b=empty_box(); for(auto n:s.triangles[t].nodes) extend(b,s.points[n]); return b; }
  int build(int first,int count) {
    const auto id=int(nodes.size()); nodes.emplace_back(); Box b=empty_box();
    for(int i=first;i<first+count;++i) { const auto tb=bounds(order[i]); extend(b,tb.min); extend(b,tb.max); }
    nodes[id]={b,first,count,-1,-1};
    if(count>8) {
      const auto d=b.max-b.min; const int a=d.x>d.y?(d.x>d.z?0:2):(d.y>d.z?1:2); const auto mid=first+count/2;
      std::nth_element(order.begin()+first,order.begin()+mid,order.begin()+first+count,[&](Index x,Index y) { const auto bx=bounds(x),by=bounds(y); return coord(bx.min+bx.max,a)<coord(by.min+by.max,a); });
      const auto left=build(first,count/2),right=build(mid,count-count/2); nodes[id].left=left; nodes[id].right=right;
    }
    return id;
  }
  void nearest(int id,Vec3d p,double& distance,Index& tri,Index excluded=-1) const {
    const auto& n=nodes[id]; if(box_distance(n.b,p)>distance) return;
    if(n.left>=0) {
      const bool left=box_distance(nodes[n.left].b,p)<box_distance(nodes[n.right].b,p);
      nearest(left?n.left:n.right,p,distance,tri,excluded); nearest(left?n.right:n.left,p,distance,tri,excluded); return;
    }
    for(int i=n.first;i<n.first+n.count;++i) { const auto t=order[i]; const auto& tr=s.triangles[t]; if(tr.component==excluded) continue; const auto& v=tr.nodes; const auto d=triangle_distance(p,s.points[v[0]],s.points[v[1]],s.points[v[2]]); if(d<distance) { distance=d; tri=t; } }
  }
  void hits(int id,Vec3d p,Vec3d d,std::vector<std::pair<Index,double>>& out) const {
    const auto& n=nodes[id]; if(!ray_box(n.b,p,d)) return;
    if(n.left>=0) { hits(n.left,p,d,out); hits(n.right,p,d,out); return; }
    for(int i=n.first;i<n.first+n.count;++i) {
      const auto& tri=s.triangles[order[i]]; const auto a=s.points[tri.nodes[0]],ab=s.points[tri.nodes[1]]-a,ac=s.points[tri.nodes[2]]-a;
      const auto h=cross(d,ac); const auto det=dot(ab,h);
      if(std::abs(det)<1e-14*norm(ab)*norm(ac)) continue;
      const auto q=p-a; const auto u=dot(q,h)/det; if(u<0 || u>1) continue;
      const auto r=cross(q,ab); const auto v=dot(d,r)/det; if(v<0 || u+v>1) continue;
      const auto t=dot(ac,r)/det; if(t>1e-12*norm(s.bounds.max-s.bounds.min)) out.push_back({tri.component,t});
    }
  }
  bool inside(Vec3d p) const {
    // A union of shells cannot contain a point outside their global bounds.
    if(p.x<s.bounds.min.x || p.x>s.bounds.max.x || p.y<s.bounds.min.y || p.y>s.bounds.max.y || p.z<s.bounds.min.z || p.z>s.bounds.max.z) return false;
    int votes=0;
    constexpr Vec3d rays[]={{1,.3713906763541037,.6947465906068658},{.21913721,1,.51768143},{.61391717,.27131971,1}};
    for(auto ray:rays) {
      std::vector<std::pair<Index,double>> h; hits(0,p,ray,h); std::sort(h.begin(),h.end());
      bool in=false; std::size_t j=0;
      while(j<h.size()) {
        const auto component=h[j].first; int crossings=0; double prev=-1;
        while(j<h.size() && h[j].first==component) { if(prev<0 || std::abs(h[j].second-prev)>1e-10*std::max(1.,std::abs(prev))) { ++crossings; prev=h[j].second; } ++j; }
        in=in || crossings%2;
      }
      votes+=in;
    }
    return votes>=2;
  }
};
SurfaceQuery::SurfaceQuery(const Surface& s):impl_(std::make_unique<Impl>(s)) {}
SurfaceQuery::~SurfaceQuery()=default;
SurfaceQuery::SurfaceQuery(SurfaceQuery&&) noexcept=default;
SurfaceQuery& SurfaceQuery::operator=(SurfaceQuery&&) noexcept=default;
Sample SurfaceQuery::nearest(Vec3d p) const {
  double d=std::numeric_limits<double>::infinity(); Index tri=-1; impl_->nearest(0,p,d,tri);
  return {std::sqrt(d),tri,false};
}
Sample SurfaceQuery::sample(Vec3d p) const { auto s=nearest(p); s.inside=impl_->inside(p); return s; }
double SurfaceQuery::other_distance(Vec3d p,Index component) const {
  double d=std::numeric_limits<double>::infinity(); Index tri=-1; impl_->nearest(0,p,d,tri,component); return std::sqrt(d);
}
}  // namespace kutta::mesher
