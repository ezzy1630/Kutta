#include <iostream>
#include <limits>
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"
#include "../backend_helpers.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;

namespace {
void linear_exactness(const PolyMesh& mesh) {
  auto g=compute_geometry(mesh);
  const auto field=[](Vec3d p) { return float(1.3+0.7*p.x-0.4*p.y+0.2*p.z); };
  std::vector<float> cell,boundary;
  for(auto p:g.cell_centroid) cell.push_back(field(p));
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) boundary.push_back(field(g.face_centroid[f]));
  std::vector<float> baseline;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g); auto s=fv::make_stencil(*dev,mesh,g);
    auto dc=dev->upload<float>(cell), db=dev->upload<float>(boundary);
    auto grad=fv::vector_field(*dev,mesh.n_cells); auto psi=dev->alloc<float>(mesh.n_cells);
    fv::gradient(*dev,m,s,dc,db,grad);
    auto x=dev->download(grad.x),y=dev->download(grad.y),z=dev->download(grad.z);
    double max_error=0;
    for(Index c=0;c<mesh.n_cells;++c) max_error=std::max({max_error,std::abs(x[c]-.7),std::abs(y[c]+.4),std::abs(z[c]-.2)});
    std::cout<<"FV_LSQ backend="<<backend_name(dev->backend())<<" cells="<<mesh.n_cells<<" max_error="<<max_error<<'\n';
    const auto cf=build_cell_faces(mesh);
    const auto inv=dev->download(s.lsq_inverse);
    using namespace kutta::detail;
    const auto exact=[](Vec3d p) {return 1.3+.7*p.x-.4*p.y+.2*p.z;};
    for(Index c=0;c<mesh.n_cells;++c) {
      double bound[3]{}, reference[3]{};
      for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
        const auto f=cf.faces[j];
        const auto n=f<mesh.n_internal_faces()?(mesh.owner[f]==c?mesh.neighbour[f]:mesh.owner[f]):-1;
        const auto point=n<0?g.face_centroid[f]:g.cell_centroid[n];
        const auto d=point-g.cell_centroid[c];
        const double value=n<0?boundary[f-mesh.n_internal_faces()]:cell[n];
        const double rounding=std::abs(value-exact(point))+std::abs(cell[c]-exact(g.cell_centroid[c]));
        const double weight[3]{(inv[c]*d.x+inv[c+mesh.n_cells]*d.y+inv[c+2*mesh.n_cells]*d.z)/dot(d,d),
          (inv[c+mesh.n_cells]*d.x+inv[c+3*mesh.n_cells]*d.y+inv[c+4*mesh.n_cells]*d.z)/dot(d,d),
          (inv[c+2*mesh.n_cells]*d.x+inv[c+4*mesh.n_cells]*d.y+inv[c+5*mesh.n_cells]*d.z)/dot(d,d)};
        for(int k=0;k<3;++k) {bound[k]+=std::abs(weight[k])*rounding;reference[k]+=weight[k]*(value-cell[c]);}
      }
      const double actual[3]{x[c],y[c],z[c]}, expected[3]{.7,-.4,.2};
      for(int k=0;k<3;++k) {
        CHECK(std::abs(actual[k]-reference[k])<2e-6);
        CHECK(std::abs(actual[k]-expected[k])<bound[k]+2e-6);
      }
    }
    if(baseline.empty()) baseline=x;
    else for(Index c=0;c<mesh.n_cells;++c) CHECK(x[c]==doctest::Approx(baseline[c]).epsilon(1e-6));
    fv::limit_gradient(*dev,m,s,dc,db,grad,psi);
    auto limiter=dev->download(psi);
    auto rx=dev->download(s.face_delta.x),ry=dev->download(s.face_delta.y),rz=dev->download(s.face_delta.z);
    auto adj=dev->download(s.adjacent),offsets=dev->download(m.cell_face_offsets);
    for(Index c=0;c<mesh.n_cells;++c) {
      float lo=cell[c],hi=cell[c];
      for(Index j=offsets[c];j<offsets[c+1];++j) {
        float v=adj[j]<0?boundary[-adj[j]-1]:cell[adj[j]];lo=std::min(lo,v);hi=std::max(hi,v);
      }
      for(Index j=offsets[c];j<offsets[c+1];++j) {
        float reconstructed=cell[c]+limiter[c]*(x[c]*rx[j]+y[c]*ry[j]+z[c]*rz[j]);
        CHECK(reconstructed>=lo-2e-6f); CHECK(reconstructed<=hi+2e-6f);
      }
    }
  }
}
}
TEST_CASE("FV least squares is linear exact and limiter bounds face reconstructions") {
  linear_exactness(make_box_hex(5,4,3,{0,0,0},{1,1,1}));
  linear_exactness(make_box_tet(4,3,2,{0,0,0},{1,1,1}));
}
TEST_CASE("FV least squares is linear exact on real Fluent polyhedra") {
  const std::filesystem::path path="/Users/ezzy/.Others/Kutta/validation/fluent_examples/box/poly.msh.h5";
#if KUTTA_HAS_HDF5
  if(std::filesystem::exists(path)) linear_exactness(read_fluent_cff(path));
  else MESSAGE("SKIP: external Fluent poly box absent");
#else
  MESSAGE("SKIP: HDF5 is unavailable");
#endif
}
