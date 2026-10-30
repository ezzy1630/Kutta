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
  const std::filesystem::path path=KUTTA_SOURCE_DIR "/validation/fluent_examples/box/poly.msh.h5";
#if KUTTA_HAS_HDF5
  if(std::filesystem::exists(path)) linear_exactness(read_fluent_cff(path));
  else MESSAGE("SKIP: external Fluent poly box absent");
#else
  MESSAGE("SKIP: HDF5 is unavailable");
#endif
}

TEST_CASE("FV fused velocity gradients equal scalar gradient kernels") {
  auto mesh=make_box_tet(3,3,3,{0,0,0},{1,1,1});auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);
    auto q=fv::vector_field(*dev,m.n_cells),bc=fv::vector_field(*dev,m.n_boundary_faces());
    std::vector<float> cell,face;
    for(auto p:g.cell_centroid) cell.push_back(float(p.x*p.x-.3*p.y+std::sin(p.z)));
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) {
      auto p=g.face_centroid[f];face.push_back(float(p.x*p.x-.3*p.y+std::sin(p.z)));
    }
    dev->upload<float>(q.x,cell);dev->upload<float>(bc.x,face);
    dev->fill(q.y,2);dev->fill(bc.y,2);dev->fill(q.z,-3);dev->fill(bc.z,-3);
    auto u=fv::vector_field(*dev,m.n_cells),v=fv::vector_field(*dev,m.n_cells),w=fv::vector_field(*dev,m.n_cells);
    auto reference=fv::vector_field(*dev,m.n_cells);
    fv::vector_gradient(*dev,m,s,q,bc,u,v,w);
    fv::gradient(*dev,m,s,q.x,bc.x,reference);
    CHECK(dev->download(u.x)==dev->download(reference.x));
    CHECK(dev->download(u.y)==dev->download(reference.y));
    CHECK(dev->download(u.z)==dev->download(reference.z));
    CHECK(dev->max_abs(v.x,m.n_cells)==0);CHECK(dev->max_abs(v.y,m.n_cells)==0);CHECK(dev->max_abs(v.z,m.n_cells)==0);
    CHECK(dev->max_abs(w.x,m.n_cells)==0);CHECK(dev->max_abs(w.y,m.n_cells)==0);CHECK(dev->max_abs(w.z,m.n_cells)==0);
  }
}

TEST_CASE("FV extrapolated boundaries cannot preserve seeded null-space gradients") {
  using namespace kutta::detail;
  auto host=make_box_tet(4,4,4,{0,0,0},{1,1,1});
  const auto geom=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,host,geom);auto s=fv::make_stencil(*dev,host,geom);
    auto cell=dev->alloc<float>(m.n_cells);dev->fill(cell,3);
    auto bc=dev->alloc<float>(m.n_boundary_faces());
    auto grad=fv::vector_field(*dev,m.n_cells,1);
    std::vector<float> boundary(m.n_boundary_faces());
    // Reproduce the former feedback loop on a constant field.
    for(int step=0;step<200;++step) {
      const auto x=dev->download(grad.x),y=dev->download(grad.y),z=dev->download(grad.z);
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) {
        const auto c=host.owner[f];const auto d=geom.face_centroid[f]-geom.cell_centroid[c];
        boundary[f-m.n_internal_faces]=float(3+x[c]*d.x+y[c]*d.y+z[c]*d.z);
      }
      dev->upload<float>(bc,boundary);fv::gradient(*dev,m,s,cell,bc,grad);
    }
    CHECK(dev->max_abs(grad.x,m.n_cells)>.01);
    const std::vector<Index> independent(m.n_boundary_faces(),0);
    auto lsq=fv::boundary_least_squares(*dev,host,geom,independent);
    fv::gradient(*dev,m,s,lsq,cell,bc,grad);
    CHECK(dev->max_abs(grad.x,m.n_cells)==0);
    CHECK(dev->max_abs(grad.y,m.n_cells)==0);
    CHECK(dev->max_abs(grad.z,m.n_cells)==0);
  }
}
