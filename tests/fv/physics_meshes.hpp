#pragma once
#include <cmath>
#include <map>
#include <numbers>
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::test {
// Deterministic clipped Voronoi prisms: actual non-Cartesian polyhedra,
// including non-orthogonal faces, with the same rectangular boundary as hex.
inline PolyMesh physics_poly_box(int nx,int ny) {
  using namespace detail;
  std::vector<Vec3d> seeds;
  for(int j=0;j<ny;++j) for(int i=0;i<nx;++i)
    seeds.push_back({(i+.5+.13*std::sin(3*i+7*j))/nx,(j+.5+.13*std::cos(5*i+2*j))/ny,0});
  std::vector<Vec3d> points;std::map<std::array<long long,3>,Index> ids;
  auto vertex=[&](Vec3d p) {
    const std::array<long long,3> key{std::llround(p.x*1e10),std::llround(p.y*1e10),std::llround(p.z*1e10)};
    const auto [it,added]=ids.emplace(key,Index(points.size()));if(added) points.push_back(p);return it->second;
  };
  std::vector<CellFaces> cells;
  for(auto a:seeds) {
    std::vector<Vec3d> poly{{0,0,0},{1,0,0},{1,1,0},{0,1,0}};
    for(auto b:seeds) {
      auto normal=b-a;if(norm(normal)<1e-12) continue;
      const auto mid=(a+b)*.5;std::vector<Vec3d> clipped;
      for(std::size_t k=0;k<poly.size();++k) {
        const auto p=poly[k],q=poly[(k+1)%poly.size()];const auto dp=dot(p-mid,normal),dq=dot(q-mid,normal);
        if(dp<=1e-14) clipped.push_back(p);
        if((dp<0 && dq>0) || (dp>0 && dq<0)) clipped.push_back(p+(q-p)*(dp/(dp-dq)));
      }
      poly=std::move(clipped);
    }
    CellFaces faces;std::vector<Index> bottom,top;
    for(auto p:poly) {bottom.push_back(vertex(p));p.z=.1;top.push_back(vertex(p));}
    faces.push_back(bottom);faces.push_back(top);
    for(std::size_t k=0;k<poly.size();++k) {const auto j=(k+1)%poly.size();faces.push_back({bottom[k],bottom[j],top[j],top[k]});}
    cells.push_back(std::move(faces));
  }
  std::vector<Patch> patches{{"xmin"},{"xmax"},{"ymin"},{"ymax"},{"zmin"},{"zmax"}};
  auto mesh=from_cells(std::move(points),cells,std::move(patches),[&](const std::vector<Index>& ns) {
    // from_cells owns the moved points; recover coordinates through the key map.
    for(int side=0;side<6;++side) {
      bool all=true;for(auto n:ns) {
        auto it=std::find_if(ids.begin(),ids.end(),[&](const auto& kv){return kv.second==n;});
        const auto value=it->first[side/2];const auto target=side%2?(side/2==2?1000000000LL:10000000000LL):0LL;
        all=all && std::abs(value-target)<2;
      }
      if(all) return Index(side);
    }
    throw Error("Voronoi verification mesh has an unmatched internal face");
  });
  validate(mesh);return mesh;
}
// Each poly ring cell spans two sectors, with both arc faces split to remain
// conformal with the adjacent hex ring. Poly cells have eight boundary faces.
inline PolyMesh physics_annulus(int nr,int nt,bool poly) {
  using namespace detail;
  std::vector<Vec3d> points;
  for(int z=0;z<2;++z) for(int j=0;j<=nr;++j) for(int k=0;k<nt;++k) {
    const auto r=1.+double(j)/nr,t=2*std::numbers::pi*k/nt;
    points.push_back({r*std::cos(t),r*std::sin(t),.1*z});
  }
  auto id=[&](int z,int r,int k) {return Index((z*(nr+1)+r)*nt+(k%nt));};
  std::vector<CellFaces> cells;
  for(int r=0;r<nr;++r) {
    const auto step=poly && r%2?2:1;
    for(int k=0;k<nt;k+=step) {
      std::vector<std::pair<int,int>> ring;
      for(int j=0;j<=step;++j) ring.emplace_back(r,k+j);
      for(int j=step;j>=0;--j) ring.emplace_back(r+1,k+j);
      CellFaces faces;std::vector<Index> bottom,top;
      for(auto [a,b]:ring) {bottom.push_back(id(0,a,b));top.push_back(id(1,a,b));}
      faces.push_back(bottom);faces.push_back(top);
      for(std::size_t j=0;j<ring.size();++j) {
        auto q=(j+1)%ring.size();faces.push_back({bottom[j],bottom[q],top[q],top[j]});
      }
      cells.push_back(std::move(faces));
    }
  }
  // Keep a copy for boundary classification; no field lives on the host in a solve.
  const auto coordinates=points;
  auto mesh=from_cells(std::move(points),cells,{{"inner"},{"outer"},{"zmin"},{"zmax"}},
    [&](const std::vector<Index>& ns) {
      bool bottom=true,top=true,inner=true;
      for(auto n:ns) {auto p=coordinates[n];bottom=bottom && p.z==0;top=top && p.z==.1;inner=inner && std::hypot(p.x,p.y)<1.000001;}
      return Index(bottom?2:top?3:inner?0:1);
    });
  validate(mesh);return mesh;
}
inline void physics_zone(PolyMesh& mesh,int id=7) {
  CellZone z{"physics","fluid",id,{}};for(Index c=0;c<mesh.n_cells;++c) z.cells.push_back(c);mesh.cell_zones={std::move(z)};
}
} // namespace kutta::test
