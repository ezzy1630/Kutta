#include "external_aero.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "kutta/core/memory.hpp"
#include <chrono>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

namespace kutta::mesher {
using namespace detail;
namespace {
using Grid=std::array<int,3>;
struct GridHash {
  std::size_t operator()(const Grid& k) const { std::size_t h=0; for(auto v:k) h=(h*1000003)^std::hash<int>{}(v); return h; }
};
struct Octant { Grid lo; int size,level,child=-1; Sample sample{}; };
struct Tree {
  const Surface& surface; const Settings& settings; const SurfaceQuery& query;
  Box bounds; double unit; int root_size; Grid counts;
  std::vector<Octant> nodes;
  Index leaves=0;
  Tree(const Surface& sf,const Settings& st,const SurfaceQuery& q):surface(sf),settings(st),query(q) {
    const auto length=sf.bounds.max.x-sf.bounds.min.x;
    KUTTA_REQUIRE(length>0,"mesh: body must have positive x length");
    bounds={{sf.bounds.min.x-st.upstream*length,st.symmetry_y?0:sf.bounds.min.y-st.side*length,st.ground?0:sf.bounds.min.z-st.side*length},
      {sf.bounds.max.x+st.downstream*length,sf.bounds.max.y+st.side*length,sf.bounds.max.z+st.top*length}};
    KUTTA_REQUIRE(!st.ground || sf.bounds.min.z>=-1e-7*length,"mesh: ground is z=0 but geometry extends below it; translate the STL or disable ground");
    KUTTA_REQUIRE(!st.symmetry_y || sf.bounds.max.y>0,"mesh: y symmetry retains y>=0, but geometry is wholly outside it");
    auto minimum=st.min_size; for(const auto& b:st.boxes) minimum=std::min(minimum,b.size);
    int depth=0; while(st.base_size/std::pow(2.,depth)>minimum) ++depth;
    KUTTA_REQUIRE(depth<=16,"mesh: requested size range exceeds 16 octree levels");
    root_size=1<<(depth+1); unit=st.base_size/root_size; // extra bit for leaf centres
    const auto extent=(bounds.max-bounds.min)/st.base_size;
    KUTTA_REQUIRE(extent.x>0 && extent.y>0 && extent.z>0 && extent.x<100000 && extent.y<100000 && extent.z<100000,"mesh: invalid tunnel dimensions");
    counts={int(std::ceil(extent.x)),int(std::ceil(extent.y)),int(std::ceil(extent.z))};
    for(int a=0;a<3;++a) KUTTA_REQUIRE(std::int64_t(counts[a])*root_size<std::numeric_limits<int>::max()-1,"mesh: requested tunnel/size range exceeds integer lattice capacity");
    const auto roots=std::int64_t(counts[0])*counts[1]*counts[2];
    KUTTA_REQUIRE(roots<=st.max_cells,"mesh/max_cells: tunnel base grid alone exceeds the cell budget; increase base_size");
    bounds.max=bounds.min+Vec3d{double(counts[0]),double(counts[1]),double(counts[2])}*st.base_size;
    nodes.reserve(std::size_t(roots));
    for(int k=0;k<counts[2];++k) for(int j=0;j<counts[1];++j) for(int i=0;i<counts[0];++i) nodes.push_back({{i*root_size,j*root_size,k*root_size},root_size,0,-1,{}});
    leaves=Index(roots);
  }
  Vec3d position(Grid g) const { return bounds.min+Vec3d{double(g[0]),double(g[1]),double(g[2])}*unit; }
  Vec3d center(const Octant& n) const { return position({n.lo[0]+n.size/2,n.lo[1]+n.size/2,n.lo[2]+n.size/2}); }
  int locate(Grid p) const {
    for(int a=0;a<3;++a) if(p[a]<0 || p[a]>=counts[a]*root_size) return -1;
    int id=(p[2]/root_size*counts[1]+p[1]/root_size)*counts[0]+p[0]/root_size;
    while(nodes[id].child>=0) { const auto& n=nodes[id]; int c=0; for(int a=0;a<3;++a) if(p[a]>=n.lo[a]+n.size/2) c|=1<<a; id=n.child+c; }
    return id;
  }
  void split(int id) {
    const auto n=nodes[id]; KUTTA_REQUIRE(n.size>2,"mesh: balance exceeded minimum octree size");
    KUTTA_REQUIRE(leaves<=settings.max_cells-7,"mesh/max_cells: octree budget exceeded; increase sizes or max_cells");
    const auto first=int(nodes.size()); nodes[id].child=first;
    for(int i=0;i<8;++i) { auto p=n.lo; for(int a=0;a<3;++a) p[a]+=(i>>a&1)*n.size/2; nodes.push_back({p,n.size/2,n.level+1,-1,{}}); }
    leaves+=7;
  }
  void refine() {
    for(std::size_t i=0;i<nodes.size();++i) {
      auto& n=nodes[i]; if(n.child>=0) continue;
      const auto p=center(n); n.sample=query.sample(p);
      const auto h=n.size*unit; double target=settings.base_size;
      if(n.sample.distance<1.75*h) {
        target=settings.max_size;
        const auto& nearest=surface.triangles[n.sample.triangle];
        for(int e=0;e<3;++e) if(nearest.feature_edges[e]) {
          const auto a=surface.points[nearest.nodes[e]],b=surface.points[nearest.nodes[(e+1)%3]],edge=b-a;
          const auto t=std::clamp(dot(p-a,edge)/dot(edge,edge),0.,1.);
          if(norm(p-(a+edge*t))<1.75*h) target=settings.min_size;
        }
        const auto other=query.other_distance(p,surface.triangles[n.sample.triangle].component);
        // A conservative local slot estimate; distance sum overestimates gaps
        // away from the medial axis, hence the 1.5-cell proximity band above.
        if(other<settings.proximity_cells*settings.max_size)
          target=std::min(target,std::max(settings.min_size,(other+n.sample.distance)/settings.proximity_cells));
      }
      for(const auto& box:settings.boxes) {
        const auto a=position(n.lo),b=a+Vec3d{h,h,h};
        if(a.x<box.bounds.max.x && b.x>box.bounds.min.x && a.y<box.bounds.max.y && b.y>box.bounds.min.y && a.z<box.bounds.max.z && b.z>box.bounds.min.z) target=std::min(target,box.size);
      }
      // Limit adjacent target sizes by growing at most one cell size away
      // from the surface, before the exact face balance pass.
      target=std::min(target,std::max(settings.max_size,n.sample.distance*.6));
      if(h>target*(1+1e-12) && n.size>2) split(int(i));
    }
    // Each fine leaf asks every face neighbour for no more than a 2:1 jump.
    // Sampling at its face centre always finds a coarser covering neighbour.
    bool changed=true;
    while(changed) {
      changed=false;
      for(std::size_t i=0;i<nodes.size();++i) if(nodes[i].child<0) {
        const auto n=nodes[i];
        for(int x=-1;x<=1;++x) for(int y=-1;y<=1;++y) for(int z=-1;z<=1;++z) {
          if(x==0 && y==0 && z==0) continue;
          const Grid direction{x,y,z}; auto p=n.lo;
          for(int a=0;a<3;++a) p[a]+=direction[a]<0?-1:direction[a]>0?n.size:n.size/2;
          const auto j=locate(p);
          if(j>=0 && nodes[j].size>2*n.size) { split(j); changed=true; }
        }
      }
    }
    // New balance leaves have no stored surface sample yet.
    for(auto& n:nodes) if(n.child<0) n.sample=query.sample(center(n));
  }
  int max_level_jump() const {
    int jump=0;
    for(const auto& n:nodes) if(n.child<0) for(int x=-1;x<=1;++x) for(int y=-1;y<=1;++y) for(int z=-1;z<=1;++z) {
      if(x==0 && y==0 && z==0) continue;
      const Grid direction{x,y,z}; auto p=n.lo;
      for(int a=0;a<3;++a) p[a]+=direction[a]<0?-1:direction[a]>0?n.size:n.size/2;
      const auto j=locate(p); if(j>=0) jump=std::max(jump,std::abs(n.level-nodes[j].level));
    }
    return jump;
  }
  std::vector<std::array<Grid,4>> face_quads(const Octant& n,int axis,int side) const {
    auto mid=n.lo; for(int a=0;a<3;++a) mid[a]+=n.size/2;
    mid[axis]=n.lo[axis]+(side?n.size:-1); const auto neighbour=locate(mid);
    // Looking at four quadrants also catches a fine neighbour beside a coarse
    // neighbour at the sampled centre (a balanced face may be partly split).
    bool divide=neighbour>=0 && nodes[neighbour].size<n.size;
    const int u=(axis+1)%3,v=(axis+2)%3;
    for(int i=0;i<2;++i) for(int j=0;j<2;++j) { auto p=mid; p[u]=n.lo[u]+(2*i+1)*n.size/4; p[v]=n.lo[v]+(2*j+1)*n.size/4; const auto other=locate(p); divide=divide || (other>=0 && nodes[other].size<n.size); }
    const int count=divide?2:1,size=n.size/count;
    std::vector<std::array<Grid,4>> out;
    for(int j=0;j<count;++j) for(int i=0;i<count;++i) {
      auto a=n.lo; a[axis]+=side*n.size; a[u]+=i*size; a[v]+=j*size;
      auto b=a,c=a,d=a; b[u]+=size; c[u]+=size; c[v]+=size; d[v]+=size;
      out.push_back({a,b,c,d});
    }
    return out;
  }
};
struct FaceKey { std::array<Index,8> n{-1,-1,-1,-1,-1,-1,-1,-1}; bool operator==(const FaceKey&) const=default; };
struct FaceHash { std::size_t operator()(const FaceKey& k) const { std::size_t h=0; for(auto v:k.n) h=h*1000003+std::hash<Index>{}(v); return h; } };
struct Builder {
  const Tree& tree; const Surface& surface; const SurfaceQuery& query; const Settings& settings;
  std::vector<Vec3d> points;
  std::vector<double> distance;
  std::unordered_map<Grid,Index,GridHash> vertices;
  std::unordered_map<std::uint64_t,Index> intersections;
  std::vector<Face> faces;
  std::unordered_map<FaceKey,Index,FaceHash> pending;
  Index cells=0,hexes=0,cut_cells=0;
  const Progress& progress;
  explicit Builder(const Tree& t,const Progress& p):tree(t),surface(t.surface),query(t.query),settings(t.settings),progress(p) {}
  Index vertex(Grid g) {
    const auto it=vertices.find(g); if(it!=vertices.end()) return it->second;
    const auto p=tree.position(g); const auto s=query.sample(p); const auto id=index(points.size());
    double local_size=tree.settings.base_size;
    for(int x=-1;x<=0;++x) for(int y=-1;y<=0;++y) for(int z=-1;z<=0;++z) {
      const auto leaf=tree.locate({g[0]+x,g[1]+y,g[2]+z});
      if(leaf>=0) local_size=std::min(local_size,tree.nodes[leaf].size*tree.unit);
    }
    // Snap the isosurface to a nearby existing vertex to prevent arbitrarily
    // small cut slivers. The shared sample limits displacement to 1% of the
    // finest incident leaf edge and is identical on both sides of a face.
    points.push_back(p); const auto d=s.distance<local_size*.01?0:s.distance*(s.inside?-1:1); distance.push_back(d); vertices.emplace(g,id); return id;
  }
  Index crossing(Index a,Index b) {
    if(distance[a]==0) return a; if(distance[b]==0) return b;
    if(a>b) std::swap(a,b);
    const auto key=(std::uint64_t(a)<<32)|std::uint32_t(b);
    const auto it=intersections.find(key); if(it!=intersections.end()) return it->second;
    const auto t=distance[a]/(distance[a]-distance[b]); const auto id=index(points.size());
    points.push_back(points[a]+(points[b]-points[a])*t); distance.push_back(0); intersections.emplace(key,id); return id;
  }
  void add_cell(std::vector<std::vector<Index>> polygons,bool hex=false,bool cut=false) {
    std::vector<Index> all; for(const auto& f:polygons) all.insert(all.end(),f.begin(),f.end());
    std::sort(all.begin(),all.end()); all.erase(std::unique(all.begin(),all.end()),all.end());
    if(all.size()<4) return;
    Vec3d c{}; const auto origin=points[all[0]]; for(auto v:all) c=c+(points[v]-origin); c=origin+c/double(all.size());
    KUTTA_REQUIRE(cells<settings.max_cells,"mesh/max_cells: cut-cell budget exceeded; increase sizes or max_cells");
    for(auto& poly:polygons) {
      poly.erase(std::unique(poly.begin(),poly.end()),poly.end()); if(poly.size()>1 && poly.front()==poly.back()) poly.pop_back();
      if(poly.size()<3) continue;
      Vec3d area{},fc{}; for(std::size_t j=0;j<poly.size();++j) { area=area+cross(points[poly[j]]-c,points[poly[(j+1)%poly.size()]]-c); fc=fc+(points[poly[j]]-c); }
      KUTTA_REQUIRE(norm(area)>0,"mesh: degenerate clipped face; refine or translate geometry away from grid coincidence");
      if(dot(area,fc)<0) std::reverse(poly.begin(),poly.end());
      FaceKey key; KUTTA_REQUIRE(poly.size()<=8,"mesh: unexpected construction polygon"); std::copy(poly.begin(),poly.end(),key.n.begin()); std::sort(key.n.begin(),key.n.end());
      const auto it=pending.find(key);
      if(it==pending.end()) { const auto id=index(faces.size()); faces.push_back({std::move(poly),cells,-1,-1}); pending.emplace(key,id); }
      else { auto& f=faces[it->second]; KUTTA_REQUIRE(f.neighbour<0 && f.owner!=cells,"mesh: nonmanifold face"); f.neighbour=cells; pending.erase(it); }
    }
    ++cells; hexes+=hex; cut_cells+=cut;
  }
  using Poly = std::vector<std::vector<Index>>;
  Poly tetra(std::array<Index,4> t) {
    const int positive=int(std::count_if(t.begin(),t.end(),[&](Index v) { return distance[v]>0; }));
    if(!positive) return {};
    constexpr int sides[4][3]={{0,1,2},{0,3,1},{1,3,2},{2,3,0}};
    std::vector<std::vector<Index>> polygons;
    for(const auto& side:sides) {
      std::vector<Index> f;
      for(int e=0;e<3;++e) { const auto a=t[side[e]],b=t[side[(e+1)%3]]; if(distance[a]>=0) f.push_back(a); if((distance[a]<0 && distance[b]>0)||(distance[a]>0 && distance[b]<0)) f.push_back(crossing(a,b)); }
      f.erase(std::unique(f.begin(),f.end()),f.end()); if(f.size()>1 && f.front()==f.back()) f.pop_back(); if(f.size()>=3) polygons.push_back(std::move(f));
    }
    std::vector<Index> cap;
    if(positive<4) {
      for(int a=0;a<4;++a) {
        if(distance[t[a]]==0) cap.push_back(t[a]);
        for(int b=a+1;b<4;++b) if(distance[t[a]]*distance[t[b]]<0) cap.push_back(crossing(t[a],t[b]));
      }
      std::sort(cap.begin(),cap.end()); cap.erase(std::unique(cap.begin(),cap.end()),cap.end());
      if(cap.size()>=3) {
        const auto o=points[cap[0]],u=points[cap[1]]-o;
        Vec3d normal{}; for(std::size_t i=2;i<cap.size();++i) if(norm(normal)==0) normal=cross(u,points[cap[i]]-o);
        if(norm(normal)>0) {
          const auto v=cross(normal,u); Vec3d center{}; for(auto n:cap) center=center+(points[n]-o); center=o+center/double(cap.size());
          std::sort(cap.begin(),cap.end(),[&](Index a,Index b) { const auto x=points[a]-center,y=points[b]-center; return std::atan2(dot(x,v)/norm(v),dot(x,u)/norm(u))<std::atan2(dot(y,v)/norm(v),dot(y,u)/norm(u)); });
          // When a whole original tetra face is on the surface, it is already
          // in the clipped original faces, so do not add it a second time.
          auto sorted=cap; std::sort(sorted.begin(),sorted.end()); bool duplicate=false;
          for(auto f:polygons) { std::sort(f.begin(),f.end()); duplicate=duplicate || f==sorted; }
          if(!duplicate) polygons.push_back(cap);
        }
      }
    }
    return polygons;
  }
  // Merge only when the union remains convex. This removes the artificial
  // tetra subdivision at planar walls and keeps re-entrant solid corners split.
  // No volume is discarded and no face is moved by agglomeration.
  void agglomerate(std::vector<Poly>& polys,double h) const {
    auto vertices_of=[](const Poly& p) {
      std::vector<Index> ids; for(const auto& f:p) ids.insert(ids.end(),f.begin(),f.end());
      std::sort(ids.begin(),ids.end()); ids.erase(std::unique(ids.begin(),ids.end()),ids.end()); return ids;
    };
    for(auto& p:polys) {
      const auto ids=vertices_of(p); Vec3d center{}; for(auto v:ids) center=center+points[v]; center=center/double(ids.size());
      for(auto& f:p) { Vec3d area{},fc{}; for(std::size_t i=0;i<f.size();++i) { area=area+cross(points[f[i]]-center,points[f[(i+1)%f.size()]]-center); fc=fc+points[f[i]]-center; } if(dot(area,fc)<0) std::reverse(f.begin(),f.end()); }
    }
    bool changed=true;
    while(changed) {
      changed=false;
      for(std::size_t i=0;i<polys.size();++i) if(!polys[i].empty()) for(std::size_t j=i+1;j<polys.size();++j) if(!polys[j].empty()) {
        std::vector<std::pair<std::size_t,std::size_t>> common;
        for(std::size_t a=0;a<polys[i].size();++a) for(std::size_t b=0;b<polys[j].size();++b) {
          auto x=polys[i][a],y=polys[j][b]; std::sort(x.begin(),x.end()); std::sort(y.begin(),y.end()); if(x==y) common.push_back({a,b});
        }
        if(common.empty()) continue;
        Poly merged;
        for(std::size_t a=0;a<polys[i].size();++a) if(std::none_of(common.begin(),common.end(),[&](auto p) { return p.first==a; })) merged.push_back(polys[i][a]);
        for(std::size_t b=0;b<polys[j].size();++b) if(std::none_of(common.begin(),common.end(),[&](auto p) { return p.second==b; })) merged.push_back(polys[j][b]);
        const auto ids=vertices_of(merged); bool convex=true;
        for(const auto& f:merged) {
          const auto origin=points[f[0]]; Vec3d area{}; for(std::size_t v=1;v+1<f.size();++v) area=area+cross(points[f[v]]-origin,points[f[v+1]]-origin);
          for(auto v:ids) if(dot(area,points[v]-origin)>1e-10*h*norm(area)) { convex=false; break; }
          if(!convex) break;
        }
        if(convex) { polys[i]=std::move(merged); polys[j].clear(); changed=true; }
      }
    }
  }
  PolyMesh build() {
    std::unordered_set<Grid,GridHash> corners;
    for(const auto& n:tree.nodes) if(n.child<0) for(int c=0;c<8;++c) {
      auto g=n.lo; for(int a=0;a<3;++a) g[a]+=(c>>a&1)*n.size; corners.insert(g);
    }
    std::vector<std::int8_t> kind(tree.nodes.size(),-1);
    // Decide shared face representation before publishing either side. Whole
    // fluid neighbours use one planar polygon; only cut interfaces need fans.
    for(std::size_t id=0;id<tree.nodes.size();++id) {
      const auto& n=tree.nodes[id]; if(n.child>=0) continue;
      if(n.sample.distance>n.size*tree.unit*.87) { if(!n.sample.inside) kind[id]=0; continue; }
      bool negative=false,positive=false;
      auto sample=[&](Grid g) { const auto d=distance[vertex(g)]; negative=negative || d<0; positive=positive || d>0; };
      sample({n.lo[0]+n.size/2,n.lo[1]+n.size/2,n.lo[2]+n.size/2});
      for(int a=0;a<3;++a) for(int side=0;side<2;++side) for(const auto& g:tree.face_quads(n,a,side)) {
        Grid mid{}; for(int axis=0;axis<3;++axis) mid[axis]=(g[0][axis]+g[2][axis])/2; sample(mid);
        for(int j=0;j<4;++j) { sample(g[j]); Grid edge{}; for(int axis=0;axis<3;++axis) edge[axis]=(g[j][axis]+g[(j+1)%4][axis])/2; if(corners.contains(edge)) sample(edge); }
      }
      if(positive) kind[id]=negative?1:0;
    }
    auto last_progress=std::chrono::steady_clock::now();
    for(std::size_t id=0;id<tree.nodes.size();++id) if(kind[id]>=0) {
      if(id%1024==0 && progress && std::chrono::steady_clock::now()-last_progress>std::chrono::seconds(1)) {
        progress("trim",.5+.24*double(id)/tree.nodes.size(),"Clipping and merging convex surface cells"); last_progress=std::chrono::steady_clock::now();
      }
      const auto& n=tree.nodes[id]; const bool cut=kind[id]==1;
      std::vector<std::array<Index,3>> triangles; Poly polygons; int quad_count=0;
      for(int a=0;a<3;++a) for(int side=0;side<2;++side) for(const auto& g:tree.face_quads(n,a,side)) {
        ++quad_count; Grid mid{}; for(int axis=0;axis<3;++axis) mid[axis]=(g[0][axis]+g[2][axis])/2;
        std::vector<Index> rim;
        for(int j=0;j<4;++j) {
          rim.push_back(vertex(g[j])); Grid edge{};
          for(int axis=0;axis<3;++axis) edge[axis]=(g[j][axis]+g[(j+1)%4][axis])/2;
          if(corners.contains(edge)) rim.push_back(vertex(edge));
        }
        auto across=mid; across[a]+=side?0:-1; const auto neighbour=tree.locate(across);
        const bool fan=cut || (neighbour>=0 && kind[neighbour]==1);
        if(fan) {
          const auto face_center=vertex(mid);
          for(std::size_t j=0;j<rim.size();++j) {
            const std::array<Index,3> t{face_center,rim[j],rim[(j+1)%rim.size()]};
            if(cut) triangles.push_back(t); else polygons.push_back({t[0],t[1],t[2]});
          }
        } else polygons.push_back(std::move(rim));
      }
      if(cut) {
        const auto center=vertex({n.lo[0]+n.size/2,n.lo[1]+n.size/2,n.lo[2]+n.size/2});
        std::vector<Poly> polys;
        for(const auto& t:triangles) { auto p=tetra({center,t[0],t[1],t[2]}); if(!p.empty()) polys.push_back(std::move(p)); }
        agglomerate(polys,n.size*tree.unit);
        for(auto& p:polys) if(!p.empty()) add_cell(std::move(p),false,true);
      } else add_cell(std::move(polygons),quad_count==6,false);
    }
    KUTTA_REQUIRE(cells>0,"mesh: no fluid cells remain");
    std::vector<Patch> patches{{"inlet",PatchKind::VelocityInlet,"velocity-inlet",3},{"outlet",PatchKind::PressureOutlet,"pressure-outlet",4},
      {tree.settings.symmetry_y?"symmetry":"side-low",PatchKind::Symmetry,"symmetry",5},{"side-high",PatchKind::Symmetry,"symmetry",6},
      {tree.settings.ground?"ground":"bottom",tree.settings.ground?PatchKind::Wall:PatchKind::Symmetry,tree.settings.ground?"wall":"symmetry",7},{"top",PatchKind::Symmetry,"symmetry",8}};
    for(const auto& name:surface.zones) {
      KUTTA_REQUIRE(std::none_of(patches.begin(),patches.end(),[&](const Patch& p) { return p.name==name; }),"STL zone name collides with tunnel boundary: "+name);
      patches.push_back({name,PatchKind::Wall,"wall",int(patches.size()+3)});
    }
    for(const auto& [key,id]:pending) {
      (void)key; auto& f=faces[id]; Vec3d c{}; for(auto n:f.nodes) c=c+points[n]; c=c/double(f.nodes.size());
      f.patch=-1;
      for(int a=0;a<3;++a) for(int side=0;side<2;++side) {
        const auto value=[&](Vec3d p) { return a==0?p.x:a==1?p.y:p.z; };
        const auto plane=value(side?tree.bounds.max:tree.bounds.min);
        if(std::all_of(f.nodes.begin(),f.nodes.end(),[&](Index v) { return std::abs(value(points[v])-plane)<tree.unit*1e-6; })) f.patch=2*a+side;
      }
      if(f.patch<0) {
        KUTTA_REQUIRE(std::all_of(f.nodes.begin(),f.nodes.end(),[&](Index v) { return distance[v]==0; }),"mesh: unmatched internal face is not on the shared isosurface");
        const auto s=query.sample(c);
        KUTTA_REQUIRE(s.distance<tree.settings.base_size,"mesh: unpaired internal face away from surface");
        f.patch=6+surface.triangles[s.triangle].zone;
      }
    }
    // Remove empty zones (for example a wheel wholly on the discarded side of
    // a symmetry model); keep their absence explicit in the report.
    std::vector<Index> counts(patches.size()); for(const auto& f:faces) if(f.neighbour<0) ++counts[f.patch];
    std::vector<Index> map(patches.size(),-1); std::vector<Patch> used;
    for(std::size_t i=0;i<patches.size();++i) if(counts[i]) { map[i]=index(used.size()); used.push_back(patches[i]); }
    for(auto& f:faces) if(f.neighbour<0) f.patch=map[f.patch];
    // Compact unused centres/inside vertices before export.
    std::vector<Index> remap(points.size(),-1); std::vector<Vec3d> used_points;
    for(auto& f:faces) for(auto& n:f.nodes) { if(remap[n]<0) { remap[n]=index(used_points.size()); used_points.push_back(points[n]); } n=remap[n]; }
    return pack(std::move(used_points),cells,std::move(faces),std::move(used));
  }
};
Json vec(Vec3d p) { return Json::Array{p.x,p.y,p.z}; }
}
Result generate(const Surface& surface,const Settings& settings,double speed,const Progress& progress) {
  const auto start=std::chrono::steady_clock::now();
  const auto emit=[&](std::string_view stage,double f,std::string_view m) { if(progress) progress(stage,f,m); };
  KUTTA_REQUIRE(settings.layer_count==0,"mesh/layers/count: prism layers are not available in v1; set count to 0 for a wall-function mesh");
  KUTTA_REQUIRE(std::isfinite(speed) && speed>0,"mesh: speed must be positive and finite");
  SurfaceQuery query(surface);
  emit("octree",.15,"Refining surface, slots and volume boxes");
  Tree tree(surface,settings,query); tree.refine();
  const auto level_jump=tree.max_level_jump(); KUTTA_REQUIRE(level_jump<=1,"mesh: octree 2:1 balance invariant failed");
  emit("trim",.5,"Clipping a conforming shared signed-distance field");
  Builder builder(tree,progress); auto mesh=builder.build();
  emit("layers",.75,"No prism layers requested; wall-function sizing is reported");
  emit("quality",.8,"Checking closure, positive volume and convexity");
  validate(mesh); const auto geometry=compute_geometry(mesh); const auto quality=orthogonal_quality(mesh,geometry);
  const auto worst=std::min_element(quality.begin(),quality.end());
  const auto min_quality=*worst; const auto worst_cell=Index(worst-quality.begin());
  const auto poor_cells=std::count_if(quality.begin(),quality.end(),[](double q) { return q<=.05; });
  double max_skew=0,max_ratio=1; Index nonconvex=0;
  // Every face plane must bound every vertex of its cell. This detects a
  // nonconvex cut cell instead of treating a positive signed volume as proof.
  const auto cell_faces=build_cell_faces(mesh);
  for(Index c=0;c<mesh.n_cells;++c) {
    bool convex=true; std::vector<Index> vertices;
    for(auto i=cell_faces.offsets[c];i<cell_faces.offsets[c+1];++i) { const auto ns=nodes(mesh,cell_faces.faces[i]); vertices.insert(vertices.end(),ns.begin(),ns.end()); }
    std::sort(vertices.begin(),vertices.end()); vertices.erase(std::unique(vertices.begin(),vertices.end()),vertices.end());
    for(auto i=cell_faces.offsets[c];i<cell_faces.offsets[c+1];++i) {
      const auto f=cell_faces.faces[i]; const auto normal=geometry.face_area[f]*cell_faces.signs[i];
      const auto tolerance=1e-8*std::cbrt(geometry.cell_volume[c])*norm(normal);
      for(auto v:vertices)
        if(dot(normal,mesh.points[v]-geometry.face_centroid[f])>tolerance) convex=false;
    }
    nonconvex+=!convex;
  }
  for(Index f=0;f<mesh.n_internal_faces();++f) {
    const auto o=mesh.owner[f],n=mesh.neighbour[f]; const auto link=geometry.cell_centroid[n]-geometry.cell_centroid[o];
    const auto t=dot(geometry.face_centroid[f]-geometry.cell_centroid[o],geometry.face_area[f])/dot(link,geometry.face_area[f]);
    const auto skew=norm(geometry.face_centroid[f]-(geometry.cell_centroid[o]+link*t))/std::max(norm(link),1e-30); max_skew=std::max(max_skew,skew);
    max_ratio=std::max(max_ratio,std::max(geometry.cell_volume[o]/geometry.cell_volume[n],geometry.cell_volume[n]/geometry.cell_volume[o]));
  }
  KUTTA_REQUIRE(nonconvex==0,"mesh: nonconvex cells detected; nothing written");
  KUTTA_REQUIRE(min_quality>1e-8 && std::isfinite(max_skew) && std::isfinite(max_ratio),"mesh: cut-cell quality below hard floor (1e-8); refine or translate geometry");
  Json::Array zones,missing;
  for(const auto& p:mesh.patches) zones.push_back(Json::Object{{"name",p.name},{"faces",double(p.size)},{"type",p.source_type}});
  for(const auto& name:surface.zones) if(std::none_of(mesh.patches.begin(),mesh.patches.end(),[&](const Patch& p) { return p.name==name; })) missing.push_back(name);
  double min_wall=std::numeric_limits<double>::infinity(),max_wall=0,sum_wall=0; Index wall_faces=0;
  for(const auto& p:mesh.patches) if(p.kind==PatchKind::Wall) for(Index f=p.start;f<p.start+p.size;++f) {
    const auto d=std::abs(dot(geometry.face_centroid[f]-geometry.cell_centroid[mesh.owner[f]],geometry.face_area[f]))/norm(geometry.face_area[f]);
    min_wall=std::min(min_wall,d); max_wall=std::max(max_wall,d); sum_wall+=d; ++wall_faces;
  }
  const double nu=1.7894e-5/1.225,length=surface.bounds.max.x-surface.bounds.min.x;
  const double cf=.026/std::pow(speed*length/nu,1./7),ustar=speed*std::sqrt(cf/2);
  const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  const auto usage=memory::usage();
  Json::Object report{{"v",1},{"cells",double(mesh.n_cells)},{"faces",double(mesh.n_faces())},{"nodes",double(mesh.n_points())},
    {"cell_types",Json::Object{{"hex",double(builder.hexes)},{"poly",double(mesh.n_cells-builder.hexes)},{"prism",0}}},
    {"cut_cells",double(builder.cut_cells)},{"octree_leaves",double(tree.leaves)},{"max_face_level_difference",level_jump},
    {"bounds",Json::Array{vec(tree.bounds.min),vec(tree.bounds.max)}},{"body_bounds",Json::Array{vec(surface.bounds.min),vec(surface.bounds.max)}},
    {"zones",std::move(zones)},{"missing_surface_zones",std::move(missing)},{"repairs",surface.repairs},
    {"quality",Json::Object{{"cells_at_or_below_0_05",double(poor_cells)},{"worst_cell",Json::Object{{"index",double(worst_cell)},{"centroid",vec(geometry.cell_centroid[worst_cell])},{"volume",geometry.cell_volume[worst_cell]}}},{"orthogonal_min",min_quality},{"orthogonal_mean",std::accumulate(quality.begin(),quality.end(),0.)/quality.size()},
      {"skewness_max",max_skew},{"adjacent_volume_ratio_max",max_ratio},{"nonconvex_cells",double(nonconvex)},{"negative_cells",0},
      {"volume_min",*std::min_element(geometry.cell_volume.begin(),geometry.cell_volume.end())},{"hard_orthogonal_floor",1e-8}}},
    {"layers",Json::Object{{"requested",settings.layer_count},{"generated",0},{"status","wall-function mesh; prism layers unavailable"}}},
    {"estimated_first_cell_y_plus",wall_faces?Json(Json::Object{{"speed_m_s",speed},{"min",min_wall*ustar/nu},{"mean",sum_wall/wall_faces*ustar/nu},{"max",max_wall*ustar/nu},{"method","flat-plate Cf=0.026/Re_L^(1/7); centroid wall distance; estimate, not solved y+"}}):Json(nullptr)},
    {"wall_seconds",elapsed},{"peak_memory_bytes",usage.peak?Json(double(usage.peak)):Json(nullptr)},
    {"geometry_method","shared piecewise-linear signed distance; cut tetrahedra near STL, whole octree cells elsewhere"},
    {"isosurface_snap_fraction",.01},{"geometry_resolution_m",settings.min_size},{"warnings",Json::Array{"Sub-cell features and leaks larger than the weld tolerance are not repaired. Resolve slots with surface/min_size; inspect missing_surface_zones.","No prism layers: drag and downforce require independent grid and wall-resolution validation.","Hex counts include octree cubes with triangulated planar faces; Fluent stores these as polyhedra."}}};
  return {std::move(mesh),std::move(report)};
}
}  // namespace kutta::mesher
