#include "../fluent_examples.hpp"
#include <iostream>
#include <limits>
#include "kernels.hpp"
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
  const auto path=fluent_examples()/"box/poly.msh.h5";
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
    const std::vector<fv::BoundaryObservation> independent(m.n_boundary_faces(),fv::BoundaryObservation::Omitted);
    auto lsq=fv::boundary_least_squares(*dev,host,geom,independent);
    fv::gradient(*dev,m,s,lsq,cell,bc,grad);
    CHECK(dev->max_abs(grad.x,m.n_cells)==0);
    CHECK(dev->max_abs(grad.y,m.n_cells)==0);
    CHECK(dev->max_abs(grad.z,m.n_cells)==0);
  }
}

TEST_CASE("FV projected symmetry observations preserve tangential affine gradients") {
  using namespace kutta::detail;
  auto mesh=make_box_tet(3,3,3,{0,0,0},{1,1,1});const auto geom=compute_geometry(mesh);
  const auto value=[](Vec3d p) {return float(1+.3*p.x-.2*p.y);};
  std::vector<float> cell,boundary;
  for(auto x:geom.cell_centroid) cell.push_back(value(x));
  std::vector<fv::BoundaryObservation> observations;
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) {
    const bool symmetry=std::abs(geom.face_area[f].z)>0;
    observations.push_back(symmetry?fv::BoundaryObservation::NormalFootValue:fv::BoundaryObservation::FaceValue);
    boundary.push_back(symmetry?cell[mesh.owner[f]]:value(geom.face_centroid[f]));
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,geom);auto stencil=fv::make_stencil(*dev,mesh,geom);
    auto lsq=fv::boundary_least_squares(*dev,mesh,geom,observations);
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary);
    auto grad=fv::vector_field(*dev,m.n_cells);
    fv::gradient(*dev,m,stencil,lsq,q,bc,grad);
    const auto x=dev->download(grad.x),y=dev->download(grad.y),z=dev->download(grad.z);
    for(Index c=0;c<m.n_cells;++c) {
      CHECK(x[c]==doctest::Approx(.3).epsilon(1e-5));
      CHECK(y[c]==doctest::Approx(-.2).epsilon(1e-5));
      CHECK(std::abs(z[c])<1e-5);
    }
  }
}

TEST_CASE("FV jump derivative satisfies independent face-row normal equations") {
  using namespace kutta::detail;
  auto mesh=make_box_tet(3,2,2,{0,0,0},{1,1,1});const auto geom=compute_geometry(mesh);
  const auto cf=build_cell_faces(mesh);
  const auto value=[](Vec3d p) {return p.x*p.x+.3*p.y*p.z;};
  std::vector<float> cell,boundary;std::vector<Vec3d> base;
  for(auto p:geom.cell_centroid) {cell.push_back(float(value(p)));base.push_back({2*p.x,.3*p.z,.3*p.y});}
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) boundary.push_back(float(value(geom.face_centroid[f])));
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,geom);auto stencil=fv::make_stencil(*dev,mesh,geom);
    auto fit=fv::face_jump_least_squares(*dev,mesh,geom);
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary);
    auto g=fv::vector_field(*dev,mesh.n_cells),out=fv::vector_field(*dev,mesh.n_cells);
    std::vector<float> gx,gy,gz;
    for(auto v:base) {gx.push_back(float(v.x));gy.push_back(float(v.y));gz.push_back(float(v.z));}
    dev->upload<float>(g.x,gx);dev->upload<float>(g.y,gy);dev->upload<float>(g.z,gz);
    fv::jump_corrected_gradient(*dev,m,stencil,fit,q,bc,g,out);
    auto x=dev->download(out.x),y=dev->download(out.y),z=dev->download(out.z);
    double nonzero=0;
    for(Index c=0;c<mesh.n_cells;++c) {
      Vec3d lhs{},rhs{};const Vec3d correction{x[c]-gx[c],y[c]-gy[c],z[c]-gz[c]};
      for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
        const auto f=cf.faces[j];const auto r=geom.face_centroid[f]-geom.cell_centroid[c];
        const Vec3d gc{gx[c],gy[c],gz[c]};const double own=cell[c]+dot(gc,r);
        double next=boundary[0];
        if(f<mesh.n_internal_faces()) {
          const auto n=mesh.owner[f]==c?mesh.neighbour[f]:mesh.owner[f];
          next=cell[n]+dot(Vec3d{gx[n],gy[n],gz[n]},geom.face_centroid[f]-geom.cell_centroid[n]);
        } else next=boundary[f-mesh.n_internal_faces()];
        lhs=lhs+r*(dot(r,correction)/dot(r,r));rhs=rhs+r*(.5*(next-own)/dot(r,r));
      }
      nonzero+=norm(correction);
      CHECK(norm(lhs-rhs)<3e-6);
    }
    CHECK(nonzero>.01);
  }
}

TEST_CASE("fused limiter and scaling preserve staged reconstruction") {
  const auto host=make_box_tet(4,3,2,{0,0,0},{1,.5,.25});
  const auto geom=compute_geometry(host);
  std::vector<float> cell,boundary;
  const auto value=[](Vec3d p) {return float(std::sin(11*p.x)+.7*p.y*p.z);};
  for(auto p:geom.cell_centroid) cell.push_back(value(p));
  for(Index f=host.n_internal_faces();f<host.n_faces();++f) boundary.push_back(value(geom.face_centroid[f]));
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,host,geom);auto s=fv::make_stencil(*dev,host,geom);
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary);
    auto g=fv::vector_field(*dev,m.n_cells),staged=fv::vector_field(*dev,m.n_cells),fused=fv::vector_field(*dev,m.n_cells);
    auto psi=dev->alloc<float>(m.n_cells),fused_psi=dev->alloc<float>(m.n_cells);
    fv::gradient(*dev,m,s,q,bc,g);
    fv::limit_gradient(*dev,m,s,q,bc,g,psi);
    kern::fv_scale_gradient(*dev,m.n_cells,g.x,g.y,g.z,psi,staged.x,staged.y,staged.z);
    kern::fv_limiter_scale(*dev,m.n_cells,m.cell_face_offsets,s.adjacent,
      s.face_delta.x,s.face_delta.y,s.face_delta.z,q,bc,g.x,g.y,g.z,fused_psi,fused.x,fused.y,fused.z,-1.f);
    CHECK(dev->download(psi)==dev->download(fused_psi));
    CHECK(dev->download(staged.x)==dev->download(fused.x));
    CHECK(dev->download(staged.y)==dev->download(fused.y));
    CHECK(dev->download(staged.z)==dev->download(fused.z));
  }
}

// SolverSettings::Limiter::Differentiable. Venkatakrishnan's function never
// exceeds delta/d, so with a zero threshold it still bounds face values; a
// threshold far above the field's variation leaves the gradient unlimited.
TEST_CASE("differentiable limiter bounds at zero threshold and passes small variations") {
  const auto host=make_box_tet(4,3,2,{0,0,0},{1,.5,.25});
  const auto geom=compute_geometry(host);
  std::vector<float> cell,boundary;
  const auto value=[](Vec3d p) {return float(std::sin(11*p.x)+.7*p.y*p.z);};
  for(auto p:geom.cell_centroid) cell.push_back(value(p));
  for(Index f=host.n_internal_faces();f<host.n_faces();++f) boundary.push_back(value(geom.face_centroid[f]));
  const auto [lo_it,hi_it]=std::minmax_element(cell.begin(),cell.end());
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,host,geom);auto s=fv::make_stencil(*dev,host,geom);
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary);
    auto g=fv::vector_field(*dev,m.n_cells),fused=fv::vector_field(*dev,m.n_cells);
    auto psi=dev->alloc<float>(m.n_cells),fused_psi=dev->alloc<float>(m.n_cells);
    std::vector<Buffer<float>> scratch;for(int k=0;k<4;++k) scratch.push_back(dev->alloc<float>(m.n_cells));
    auto shifted=dev->upload<float>([&]{auto v=cell;for(auto& x:v) x=3*x-2;return v;}());
    const Buffer<float>* fields[]{&q,&shifted};
    const auto range=fv::field_ranges(*dev,fields,m.n_cells,scratch);
    CHECK(range[0].range==doctest::Approx(*hi_it-*lo_it).epsilon(1e-6));CHECK(range[1].range==doctest::Approx(3*(*hi_it-*lo_it)).epsilon(1e-6));
    fv::gradient(*dev,m,s,q,bc,g);
    fv::limit_gradient(*dev,m,s,q,bc,g,psi,nullptr,nullptr,0.f);
    const auto limiter=dev->download(psi);
    const auto x=dev->download(g.x),y=dev->download(g.y),z=dev->download(g.z);
    const auto rx=dev->download(s.face_delta.x),ry=dev->download(s.face_delta.y),rz=dev->download(s.face_delta.z);
    const auto adj=dev->download(s.adjacent),offsets=dev->download(m.cell_face_offsets);
    int limited=0;
    for(Index c=0;c<m.n_cells;++c) {
      float lo=cell[c],hi=cell[c];
      for(Index j=offsets[c];j<offsets[c+1];++j) {
        const float v=adj[j]<0?boundary[-adj[j]-1]:cell[adj[j]];lo=std::min(lo,v);hi=std::max(hi,v);
      }
      for(Index j=offsets[c];j<offsets[c+1];++j) {
        const float reconstructed=cell[c]+limiter[c]*(x[c]*rx[j]+y[c]*ry[j]+z[c]*rz[j]);
        CHECK(reconstructed>=lo-2e-6f);CHECK(reconstructed<=hi+2e-6f);
      }
      limited+=limiter[c]<.99f;
    }
    CHECK(limited>0);
    // The fused kernel and the float-pair limiter with zero low words agree.
    const float eps2=std::pow(10.f*(*hi_it-*lo_it),2.f);
    fv::limit_gradient(*dev,m,s,q,bc,g,psi,nullptr,nullptr,eps2);
    kern::fv_limiter_scale(*dev,m.n_cells,m.cell_face_offsets,s.adjacent,
      s.face_delta.x,s.face_delta.y,s.face_delta.z,q,bc,g.x,g.y,g.z,fused_psi,fused.x,fused.y,fused.z,eps2);
    CHECK(dev->download(psi)==dev->download(fused_psi));
    auto low=dev->alloc<float>(m.n_cells);dev->fill(low,0);
    auto kind=dev->alloc<Index>(boundary.size());dev->fill(kind,0);
    fv::limit_gradient(*dev,m,s,q,bc,g,fused_psi,&low,&kind,eps2);
    CHECK(dev->download(psi)==dev->download(fused_psi));
    for(float v:dev->download(psi)) CHECK(v>.97f);
  }
}
