// Fluent's Green-Gauss cell-based gradient (operators.hpp, green_gauss_cell).
#include <cmath>
#include <iostream>
#include <map>

#include "../backend_helpers.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

using namespace kutta;
using namespace kutta::detail;

namespace {
double field(Vec3d p) {return std::sin(1.7*p.x)+.6*p.y*p.y-.3*p.x*p.z+.2*p.z;}
struct Values {std::vector<float> cell,boundary;};
Values sample(const PolyMesh& mesh,const MeshGeometry& g) {
  Values v;
  for(auto p:g.cell_centroid) v.cell.push_back(float(field(p)));
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) v.boundary.push_back(float(field(g.face_centroid[f])));
  return v;
}
// Theory Guide 2024 R2 Eq. 23-30 with Eq. 23-31, evaluated literally: the
// face value is the mean of the two cell values, or the boundary value.
std::vector<Vec3d> reference(const PolyMesh& mesh,const MeshGeometry& g,const Values& v) {
  std::vector<Vec3d> out(mesh.n_cells);
  for(Index f=0;f<mesh.n_faces();++f) {
    const auto a=g.face_area[f];
    if(f<mesh.n_internal_faces()) {
      const double face=.5*(double(v.cell[mesh.owner[f]])+double(v.cell[mesh.neighbour[f]]));
      out[mesh.owner[f]]=out[mesh.owner[f]]+a*face;out[mesh.neighbour[f]]=out[mesh.neighbour[f]]-a*face;
    } else out[mesh.owner[f]]=out[mesh.owner[f]]+a*double(v.boundary[f-mesh.n_internal_faces()]);
  }
  for(Index c=0;c<mesh.n_cells;++c) out[c]=out[c]/g.cell_volume[c];
  return out;
}
std::array<std::vector<float>,3> gradient(Device& dev,const PolyMesh& mesh,const MeshGeometry& g,const Values& v) {
  const auto m=upload_mesh(dev,mesh,g);const auto s=fv::make_stencil(dev,mesh,g);
  const auto fit=fv::green_gauss_cell(dev,mesh,g);
  auto cell=dev.upload<float>(v.cell),boundary=dev.upload<float>(v.boundary);auto out=fv::vector_field(dev,mesh.n_cells);
  fv::gradient(dev,m,s,fit,cell,boundary,out);
  return {dev.download(out.x),dev.download(out.y),dev.download(out.z)};
}
void matches_definition(const PolyMesh& mesh) {
  const auto g=compute_geometry(mesh);const auto v=sample(mesh,g);const auto expected=reference(mesh,g,v);
  // Float rounding of rows and values against per-cell gradient magnitudes.
  double scale=0;for(const auto& e:expected) scale=std::max(scale,norm(e));
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto x=gradient(*dev,mesh,g,v);double worst=0;
    for(Index c=0;c<mesh.n_cells;++c) worst=std::max(worst,norm(Vec3d{x[0][c],x[1][c],x[2][c]}-expected[c]));
    std::cout<<"FV_GREEN_GAUSS backend="<<backend_name(dev->backend())<<" cells="<<mesh.n_cells<<" max_error="<<worst<<" scale="<<scale<<'\n';
    CHECK(worst<2e-5*scale);
  }
}
// Polygons listed counter-clockwise; each directed edge has its cell on the
// left, the convention extrude_2d expects. All boundary edges form patch 0.
PolyMesh planar(std::vector<Vec3d> points,const std::vector<std::vector<Index>>& cells) {
  PlanarMesh p;p.points=std::move(points);p.n_cells=Index(cells.size());
  std::map<std::pair<Index,Index>,std::size_t> seen;
  for(Index c=0;c<Index(cells.size());++c) for(std::size_t k=0;k<cells[c].size();++k) {
    const auto a=cells[c][k],b=cells[c][(k+1)%cells[c].size()];
    if(auto it=seen.find({b,a});it!=seen.end()) {p.neighbour[it->second]=c;p.edge_patch[it->second]=-1;continue;}
    seen.emplace(std::pair{a,b},p.edges.size());
    p.edges.push_back({a,b});p.owner.push_back(c);p.neighbour.push_back(-1);p.edge_patch.push_back(0);
  }
  p.patches.push_back({"outer",PatchKind::Wall,"wall"});
  p.cell_zones.push_back({"fluid","fluid",1,{}});for(Index c=0;c<p.n_cells;++c) p.cell_zones[0].cells.push_back(c);
  return extrude_2d(p,.37);
}
PolyMesh mixed_planar() {
  // Quadrilaterals and triangles on a distorted grid, as in the cylinder mesh.
  std::vector<Vec3d> pts;
  for(int j=0;j<4;++j) for(int i=0;i<4;++i) pts.push_back({i+.23*std::sin(1.3*i+2.1*j),j+.19*std::cos(.7*i+1.9*j),0});
  const auto id=[](int i,int j){return Index(4*j+i);};
  std::vector<std::vector<Index>> cells;
  for(int j=0;j<3;++j) for(int i=0;i<3;++i) {
    if((i+j)%2) cells.push_back({id(i,j),id(i+1,j),id(i+1,j+1),id(i,j+1)});
    else {cells.push_back({id(i,j),id(i+1,j),id(i+1,j+1)});cells.push_back({id(i,j),id(i+1,j+1),id(i,j+1)});}
  }
  return planar(pts,cells);
}
}

TEST_CASE("FV Green-Gauss cell gradient matches Fluent's face-average definition") {
  matches_definition(make_box_hex(4,3,3,{0,0,0},{1,.8,1.3}));
  matches_definition(make_box_tet(3,3,2,{0,0,0},{1,1,1}));
  matches_definition(mixed_planar());
}

TEST_CASE("FV Green-Gauss cell gradient is linear exact on a uniform hex box") {
  // Face centroids bisect the centroid links, so the face mean is exact.
  const auto mesh=make_box_hex(5,4,3,{0,0,0},{1,1,1});const auto g=compute_geometry(mesh);
  Values v;
  const auto linear=[](Vec3d p){return 1.3+.7*p.x-.4*p.y+.2*p.z;};
  for(auto p:g.cell_centroid) v.cell.push_back(float(linear(p)));
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) v.boundary.push_back(float(linear(g.face_centroid[f])));
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto x=gradient(*dev,mesh,g,v);
    for(Index c=0;c<mesh.n_cells;++c) {
      CHECK(std::abs(x[0][c]-.7)<2e-5);CHECK(std::abs(x[1][c]+.4)<2e-5);CHECK(std::abs(x[2][c]-.2)<2e-5);
    }
  }
}

TEST_CASE("FV Green-Gauss on a one-cell extrusion is Fluent's planar gradient") {
  // A planar Fluent case has no z faces. Symmetry caps carry the cell value
  // so they add nothing, and the side faces' area per volume is the edge's
  // length per cell area: the 3D operator is the 2D one, for any thickness.
  const auto mesh=mixed_planar();const auto g=compute_geometry(mesh);
  Values v=sample(mesh,g);
  for(const auto& p:mesh.patches) if(p.name=="front" || p.name=="back")
    for(Index f=p.start;f<p.start+p.size;++f) v.boundary[f-mesh.n_internal_faces()]=v.cell[mesh.owner[f]];
  // 2D evaluation on the side faces only: edge normal times length.
  std::vector<Vec3d> planar(mesh.n_cells);std::vector<double> area(mesh.n_cells);
  for(Index f=0;f<mesh.n_faces();++f) {
    const auto a=g.face_area[f];
    if(std::abs(a.z)>1e-12*norm(a)) continue;
    const auto edge=Vec3d{a.x,a.y,0}/.37;
    if(f<mesh.n_internal_faces()) {
      const double face=.5*(double(v.cell[mesh.owner[f]])+double(v.cell[mesh.neighbour[f]]));
      planar[mesh.owner[f]]=planar[mesh.owner[f]]+edge*face;planar[mesh.neighbour[f]]=planar[mesh.neighbour[f]]-edge*face;
    } else planar[mesh.owner[f]]=planar[mesh.owner[f]]+edge*double(v.boundary[f-mesh.n_internal_faces()]);
  }
  for(Index c=0;c<mesh.n_cells;++c) area[c]=g.cell_volume[c]/.37;
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto x=gradient(*dev,mesh,g,v);
    for(Index c=0;c<mesh.n_cells;++c) {
      const auto expected=planar[c]/area[c];
      // Side-face area vectors carry z components at double rounding only.
      CHECK(std::abs(x[2][c])<1e-12*(1+norm(expected)));
      CHECK(std::abs(x[0][c]-expected.x)<2e-5*(1+norm(expected)));
      CHECK(std::abs(x[1][c]-expected.y)<2e-5*(1+norm(expected)));
    }
  }
}
