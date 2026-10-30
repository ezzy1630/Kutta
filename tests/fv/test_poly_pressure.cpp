#include <iostream>
#include <array>
#include <limits>
#include "../backend_helpers.hpp"
#include "physics_meshes.hpp"
#include "kernels.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/mesh/topology.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/linalg/csr.hpp"

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
struct ErrorNorm {
  double square=0,weight=0,maximum=0;
  void add(double e,double w=1) {
    square+=w*e*e;weight+=w;
    maximum=std::isfinite(e)?std::max(maximum,std::abs(e)):std::numeric_limits<double>::infinity();
  }
  double rms() const {return std::sqrt(square/weight);}
};
}
TEST_CASE("FV poly annulus manufactured pressure reconstruction") {
  KUTTA_FOR_EACH_BACKEND(dev) for(int nr:{8,16,32}) {
    const auto mesh=test::physics_annulus(nr,8*nr,true);const auto g=compute_geometry(mesh);
    const auto m=upload_mesh(*dev,mesh,g);const auto st=make_stencil(*dev,mesh,g);
    const auto cf=build_cell_faces(mesh);
    std::vector<bool> wall(mesh.n_cells,false);
    for(const auto& patch:mesh.patches) if(patch.name=="inner" || patch.name=="outer")
      for(Index f=patch.start;f<patch.start+patch.size;++f) wall[mesh.owner[f]]=true;
    for(bool linear:{true,false}) for(bool exact_boundary:{false,true}) {
      auto value=[&](Vec3d x) {return linear? .7*x.x-.3*x.y : .5*(x.x*x.x+x.y*x.y);};
      auto derivative=[&](Vec3d x) {return linear?Vec3d{.7,-.3,0}:Vec3d{x.x,x.y,0};};
      std::vector<float> cells,bc;
      for(auto x:g.cell_centroid) cells.push_back(float(value(x)));
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) bc.push_back(float(value(g.face_centroid[f])));
      auto p=dev->upload<float>(cells),b=dev->upload<float>(bc);auto grad=vector_field(*dev,m.n_cells);
      std::vector<BoundaryObservation> rows(m.n_boundary_faces(),exact_boundary?BoundaryObservation::FaceValue:BoundaryObservation::Omitted);
      const auto lsq=boundary_least_squares(*dev,mesh,g,rows);gradient(*dev,m,st,lsq,p,b,grad);
      const auto gx=dev->download(grad.x),gy=dev->download(grad.y),gz=dev->download(grad.z);
      ErrorNorm interior,boundary;
      for(Index c=0;c<m.n_cells;++c) {
        const auto e=Vec3d{gx[c],gy[c],gz[c]}-derivative(g.cell_centroid[c]);
        (wall[c]?boundary:interior).add(norm(e),g.cell_volume[c]);
      }
      auto fp=dev->alloc<float>(m.n_faces),jump=dev->alloc<float>(m.n_faces);
      std::vector<Index> kinds;
      for(Index f=m.n_internal_faces;f<m.n_faces;++f)
        kinds.push_back(std::abs(g.face_area[f].z)>.8*norm(g.face_area[f])?3:0);
      const auto boundary_kind=dev->upload<Index>(kinds);
      auto reconstructed_boundary=dev->alloc<float>(m.n_boundary_faces());
      kern::fv_pressure_boundary(*dev,m.n_boundary_faces(),m.owner,boundary_kind,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,p,grad.x,grad.y,grad.z,
        b,reconstructed_boundary,m.n_internal_faces);
      kern::fv_pressure_faces(*dev,m.n_faces,m.owner,m.neighbour,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
        p,grad.x,grad.y,grad.z,reconstructed_boundary,fp,jump,m.n_internal_faces);
      auto face=dev->download(fp),jumps=dev->download(jump);
      std::array<ErrorNorm,6> face_errors,jump_errors;
      double max_skew=0;int bounded=0;
      for(Index f=0;f<m.n_internal_faces;++f) {
        const auto p0=mesh.owner[f],p1=mesh.neighbour[f];const auto x=g.face_centroid[f];
        const auto a=g.face_area[f];const auto d=g.cell_centroid[p1]-g.cell_centroid[p0];
        const auto skew=x-(g.cell_centroid[p0]+g.cell_centroid[p1])*.5;
        const double ratio=norm(skew)/norm(d);max_skew=std::max(max_skew,ratio);if(ratio>1) ++bounded;
        const bool radial=std::abs(dot(a,Vec3d{x.x,x.y,0}))>.8*norm(a)*std::hypot(x.x,x.y);
        const int kind=radial?((cf.offsets[p0+1]-cf.offsets[p0])==(cf.offsets[p1+1]-cf.offsets[p1])?1:2):0;
        face_errors[kind].add(face[f]-value(x),norm(a));jump_errors[kind].add(jumps[f],norm(a));
      }
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) {
        const auto a=g.face_area[f],x=g.face_centroid[f];
        const int kind=std::abs(a.z)>.8*norm(a)?5:std::hypot(x.x,x.y)<1.5?3:4;
        face_errors[kind].add(face[f]-value(x),norm(a));
      }
      std::cout<<"POLY_PRESSURE backend="<<backend_name(dev->backend())<<" nr="<<nr<<" linear="<<linear<<" exact_bc="<<exact_boundary
        <<" grad_interior="<<interior.rms()<<" grad_wall="<<boundary.rms()<<" skew_max="<<max_skew<<" bounded="<<bounded;
      for(int k=0;k<6;++k) if(face_errors[k].weight>0) {
        std::cout<<" face"<<k<<"="<<face_errors[k].rms();
        if(jump_errors[k].weight>0) std::cout<<" jump"<<k<<"="<<jump_errors[k].rms();
      }
      std::cout<<std::endl;
      CHECK(bounded==0);
      if(linear) {
        CHECK(interior.maximum<2e-5);CHECK(boundary.maximum<2e-5);
        for(const auto& e:face_errors) if(e.weight>0) CHECK(e.maximum<2e-6);
        for(const auto& e:jump_errors) if(e.weight>0) CHECK(e.maximum<2e-6);
      }
    }
  }
}

TEST_CASE("FV poly quadratic viscous face flux and documented secondary gradient") {
  KUTTA_FOR_EACH_BACKEND(dev) for(int nr:{8,16,32}) {
    const auto mesh=test::physics_annulus(nr,8*nr,true);const auto geo=compute_geometry(mesh);
    const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
    // A divergence-free quadratic velocity: (x^2, -2xy, 0).
    for(int component=0;component<2;++component) for(bool exact_gradient:{false,true}) {
      const auto value=[&](Vec3d x) {return component==0?x.x*x.x:-2*x.x*x.y;};
      const auto derivative=[&](Vec3d x) {return component==0?Vec3d{2*x.x,0,0}:Vec3d{-2*x.y,-2*x.x,0};};
      std::vector<float> cell,bc,dx,dy;
      for(auto x:geo.cell_centroid) {cell.push_back(float(value(x)));dx.push_back(float(derivative(x).x));dy.push_back(float(derivative(x).y));}
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) bc.push_back(float(value(geo.face_centroid[f])));
      auto q=dev->upload<float>(cell),b=dev->upload<float>(bc),zero=dev->alloc<float>(m.n_cells),flux=dev->alloc<float>(m.n_faces),mu=dev->alloc<float>(m.n_faces);
      dev->fill(zero,0);dev->fill(flux,0);dev->fill(mu,1);
      auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);
      TransportEquation equation(*dev,m,st,mesh,geo);TransportSettings controls;controls.limit=false;
      if(exact_gradient) {dev->upload<float>(equation.grad.x,dx);dev->upload<float>(equation.grad.y,dy);dev->fill(equation.grad.z,0);}
      else gradient(*dev,m,st,q,b,equation.grad);
      equation.assemble_reconstructed(q,b,kind,flux,mu,zero,controls);
      const auto correction=dev->download(equation.face_correction),gx=dev->download(equation.grad.x),gy=dev->download(equation.grad.y),gz=dev->download(equation.grad.z);
      const auto weights=dev->download(m.weight),diffusion=dev->download(st.diffusion);
      std::array<ErrorNorm,5> errors;ErrorNorm published_error;
      for(Index f=0;f<m.n_faces;++f) {
        const auto p=mesh.owner[f];const auto a=geo.face_area[f],x=geo.face_centroid[f];
        const auto d=(f<m.n_internal_faces?geo.cell_centroid[mesh.neighbour[f]]:x)-geo.cell_centroid[p];
        const double alpha=dot(a,a)/dot(a,d);
        Vec3d avg{gx[p],gy[p],gz[p]};int face_type=2;
        if(f<m.n_internal_faces) {
          const auto n=mesh.neighbour[f];
          avg=avg*weights[f]+Vec3d{gx[n],gy[n],gz[n]}*(1-weights[f]);
          const bool radial=std::abs(dot(a,Vec3d{x.x,x.y,0}))>.8*norm(a)*std::hypot(x.x,x.y);
          face_type=radial?1:0;
        } else {
          if(std::abs(a.z)>.8*norm(a)) face_type=4;
          else face_type=std::hypot(x.x,x.y)<1.5?2:3;
        }
        const double scalar_jump=(f<m.n_internal_faces?cell[mesh.neighbour[f]]:bc[f-m.n_internal_faces])-double(cell[p]);
        const double actual=double(diffusion[f])*scalar_jump+correction[f];
        // Independent host-double evaluation of Fluent CM Eq. 3-3, using
        // the same stored float observations and the geometry-derived mean.
        const auto secondary=a-d*alpha;
        const double published=alpha*scalar_jump+dot(avg,secondary);
        published_error.add((actual-published)/norm(a),norm(a));
        errors[face_type].add((actual-dot(derivative(x),a))/norm(a),norm(a));
      }
      std::cout<<"POLY_VISCOUS backend="<<backend_name(dev->backend())<<" nr="<<nr<<" component="<<component<<" exact_grad="<<exact_gradient
        <<" published_Linf="<<published_error.maximum;
      for(int k=0;k<5;++k) if(errors[k].weight) std::cout<<" normal_derivative"<<k<<"="<<errors[k].rms();
      std::cout<<std::endl;
      CHECK(published_error.maximum<1e-4);
    }
  }
}

namespace {
std::vector<float> predict_faces(Device& dev,const DeviceMesh& m,const Buffer<Index>& kind,
    const Buffer<float>& q,const Buffer<float>& bc,const DeviceVec3& g,const Buffer<float>& ap,Buffer<float>& out) {
  auto spare1=dev.alloc<float>(m.n_faces),spare2=dev.alloc<float>(m.n_faces);
  kern::fv_momentum_faces(dev,m.n_faces,m.owner,m.neighbour,kind,
    m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
    q,q,q,g.x,g.y,g.z,g.x,g.y,g.z,g.x,g.y,g.z,ap,bc,bc,bc,out,spare1,spare2,m.n_internal_faces);
  return dev.download(out);
}
}
TEST_CASE("FV face projection preserves the normal stencil for quadratic fields") {
  KUTTA_FOR_EACH_BACKEND(dev) for(int variant=0;variant<5;++variant) {
    const auto mesh=variant==0?make_box_hex(4,4,4,{0,0,0},{1,1,1}):variant==1?make_box_tet(3,3,3,{0,0,0},{1,1,1}):
      test::physics_annulus(8<<(variant-2),8*(8<<(variant-2)),true);
    const auto geo=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
    const auto value=[](Vec3d x) {return x.x*x.x-2*x.x*x.y+3*x.y*x.y+.5*x.z*x.z+.3*x.x*x.z-.7*x.y*x.z+.2*x.x-.4*x.y+.6*x.z;};
    const auto derivative=[](Vec3d x) {return Vec3d{2*x.x-2*x.y+.3*x.z+.2,-2*x.x+6*x.y-.7*x.z-.4,x.z+.3*x.x-.7*x.y+.6};};
    std::vector<float> cell,boundary,diagonal,dx,dy,dz;
    for(auto x:geo.cell_centroid) {cell.push_back(float(value(x)));const auto g=derivative(x);dx.push_back(float(g.x));dy.push_back(float(g.y));dz.push_back(float(g.z));diagonal.push_back(1+.1f*(cell.size()%17));}
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(geo.face_centroid[f])));
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary),ap=dev->upload<float>(diagonal);
    const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    auto reconstruction=std::make_shared<ProjectedFaceReconstruction>(*dev,mesh,geo,rows);
    CHECK(reconstruction->fallback_faces()==0);
    if(variant==0) CHECK(reconstruction->corrected_faces()==0);
    auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);
    auto zero=dev->alloc<float>(m.n_cells),mass=dev->alloc<float>(m.n_faces),mu=dev->alloc<float>(m.n_faces);
    dev->fill(zero,0);dev->fill(mass,0);dev->fill(mu,1);
    TransportEquation equation(*dev,m,st,mesh,geo,{},reconstruction);TransportSettings controls;controls.limit=false;
    // Exact gradients isolate reconstruction from the cell-gradient operator.
    dev->upload<float>(equation.grad.x,dx);dev->upload<float>(equation.grad.y,dy);dev->upload<float>(equation.grad.z,dz);
    equation.assemble_reconstructed(q,bc,kind,mass,mu,zero,controls);
    auto face_values=dev->alloc<float>(m.n_faces);
    const auto predictor=predict_faces(*dev,m,kind,q,bc,equation.grad,ap,face_values);
    auto phi=dev->alloc<float>(m.n_cells);dev->fill(phi,0);
    reconstruction->momentum_values(q,bc,ap,face_values,&phi);
    CHECK(dev->download(face_values)==predictor);
    reconstruction->momentum_values(q,bc,ap,face_values);
    const auto face=dev->download(face_values),correction=dev->download(equation.face_correction),diffusion=dev->download(st.diffusion);
    auto grad=vector_field(*dev,m.n_faces);reconstruction->gradient(q,bc,grad);
    const auto gx=dev->download(grad.x),gy=dev->download(grad.y),gz=dev->download(grad.z);const auto mask=dev->download(reconstruction->mask());
    ErrorNorm flux_error,value_error,gradient_error;std::array<ErrorNorm,5> analytic_flux;
    for(Index f=0;f<m.n_faces;++f) {
      const auto p=mesh.owner[f];const auto area=geo.face_area[f],x=geo.face_centroid[f],n=area/norm(area);
      const double dp=dot(x-geo.cell_centroid[p],n),dn=f<m.n_internal_faces?dot(geo.cell_centroid[mesh.neighbour[f]]-x,n):0;
      const auto pf=x-n*dp,nf=x+n*dn;
      const double jump=(f<m.n_internal_faces?cell[mesh.neighbour[f]]:boundary[f-m.n_internal_faces])-double(cell[p]);
      const double flux=double(diffusion[f])*jump+correction[f];
      flux_error.add(flux/norm(area)-(value(nf)-value(pf))/(dp+dn),norm(area));
      if(f<m.n_internal_faces) value_error.add(face[f]-(dn*value(pf)+dp*value(nf))/(dp+dn),norm(area));
      if(mask[f]) gradient_error.add(norm(Vec3d{gx[f],gy[f],gz[f]}-derivative(x)),norm(area));
      int type=0;
      if(f<m.n_internal_faces) type=std::abs(dot(area,Vec3d{x.x,x.y,0}))>.8*norm(area)*std::hypot(x.x,x.y)?1:0;
      else type=std::abs(n.z)>.8?4:std::hypot(x.x,x.y)<1.5?2:3;
      analytic_flux[type].add(flux/norm(area)-dot(derivative(x),n),norm(area));
    }
    INFO("backend "<<backend_name(dev->backend())<<" variant "<<variant);
    std::cout<<"PROJECTED_QUADRATIC backend="<<backend_name(dev->backend())<<" variant="<<variant
      <<" normal_stencil_Linf="<<flux_error.maximum<<" normal_interpolation_Linf="<<value_error.maximum<<" gradient_Linf="<<gradient_error.maximum;
    for(int k=0;k<5;++k) if(analytic_flux[k].weight) std::cout<<" analytic_flux"<<k<<"="<<analytic_flux[k].rms();
    std::cout<<std::endl;
    CHECK(flux_error.maximum<2e-4);CHECK(value_error.maximum<2e-5);CHECK(gradient_error.maximum<2e-4);
  }
}

TEST_CASE("FV projected wall gradient carries the assembled viscous flux") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool poly:{false,true}) {
    const auto mesh=poly?test::physics_annulus(8,64,true):make_box_tet(3,3,3,{0,0,0},{1,1,1});
    const auto geo=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
    auto value=[](Vec3d x){return std::sin(1.3*x.x)*std::exp(.7*x.y)+x.z*x.z*x.x;};
    std::vector<float> cells,boundary;for(auto x:geo.cell_centroid) cells.push_back(float(value(x)));
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(geo.face_centroid[f])));
    auto q=dev->upload<float>(cells),bc=dev->upload<float>(boundary),c=dev->alloc<float>(m.n_faces);
    const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    ProjectedFaceReconstruction fit(*dev,mesh,geo,rows);auto raw=vector_field(*dev,m.n_faces),consistent=vector_field(*dev,m.n_faces);
    fit.gradient(q,bc,raw);fit.diffusion_correction(q,bc,c);fit.flux_gradient(q,bc,m,st,consistent);
    const auto x=dev->download(raw.x),y=dev->download(raw.y),z=dev->download(raw.z);
    const auto fx=dev->download(consistent.x),fy=dev->download(consistent.y),fz=dev->download(consistent.z);
    const auto correction=dev->download(c),diffusion=dev->download(st.diffusion);const auto mask=dev->download(fit.mask());
    ErrorNorm before,after,tangent;
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) if(mask[f]) {
      const auto a=geo.face_area[f],n=a/norm(a),g0=Vec3d{x[f],y[f],z[f]},g1=Vec3d{fx[f],fy[f],fz[f]};
      const double flux=double(diffusion[f])*(boundary[f-m.n_internal_faces]-double(cells[mesh.owner[f]]))+correction[f];
      before.add((dot(g0,a)-flux)/norm(a));after.add((dot(g1,a)-flux)/norm(a));
      const auto delta=g1-g0;tangent.add(norm(delta-n*dot(delta,n)));
    }
    std::cout<<"PROJECTED_WALL backend="<<backend_name(dev->backend())<<" poly="<<poly<<" raw_Linf="<<before.maximum
      <<" consistent_Linf="<<after.maximum<<" tangent_change_Linf="<<tangent.maximum<<std::endl;
    CHECK(after.maximum<2e-5);CHECK(tangent.maximum<2e-6);
  }
}

TEST_CASE("FV face projection is invariant under owner neighbour reversal") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool poly:{false,true}) {
    const auto original=poly?test::physics_annulus(8,64,true):make_box_tet(3,3,3,{0,0,0},{1,1,1});
    auto reversed=original;
    // Reverse cell numbering while preserving the mesh invariant owner < neighbour.
    // Every internal face consequently reverses its physical ownership.
    for(auto& c:reversed.owner) c=reversed.n_cells-1-c;
    for(auto& c:reversed.neighbour) c=reversed.n_cells-1-c;
    for(auto& zone:reversed.cell_zones) for(auto& c:zone.cells) c=reversed.n_cells-1-c;
    for(Index f=0;f<reversed.n_internal_faces();++f) {
      std::swap(reversed.owner[f],reversed.neighbour[f]);
      std::reverse(reversed.face_nodes.begin()+reversed.face_offsets[f],reversed.face_nodes.begin()+reversed.face_offsets[f+1]);
    }
    std::vector<float> corrections[2],faces[2];
    for(int side=0;side<2;++side) {
      const auto& mesh=side?reversed:original;const auto geo=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
      auto value=[](Vec3d x){return std::sin(x.x+2*x.y)+.3*x.x*x.z*x.z;};
      std::vector<float> cells,boundary,diagonal;
      for(auto x:geo.cell_centroid) {cells.push_back(float(value(x)));diagonal.push_back(1+.1f*((side?mesh.n_cells-Index(cells.size())+1:Index(cells.size()))%17));}
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(geo.face_centroid[f])));
      auto q=dev->upload<float>(cells),bc=dev->upload<float>(boundary),ap=dev->upload<float>(diagonal),c=dev->alloc<float>(m.n_faces),out=dev->alloc<float>(m.n_faces);
      auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);auto grad=vector_field(*dev,m.n_cells);gradient(*dev,m,st,q,bc,grad);
      const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
      ProjectedFaceReconstruction fit(*dev,mesh,geo,rows);fit.diffusion_correction(q,bc,c);
      predict_faces(*dev,m,kind,q,bc,grad,ap,out);fit.momentum_values(q,bc,ap,out);
      corrections[side]=dev->download(c);faces[side]=dev->download(out);
    }
    ErrorNorm correction_error,value_error;
    for(Index f=0;f<original.n_faces();++f) {
      correction_error.add(corrections[0][f]+(f<original.n_internal_faces()?1:-1)*corrections[1][f]);
      value_error.add(faces[0][f]-faces[1][f]);
    }
    INFO("backend "<<backend_name(dev->backend())<<" poly "<<poly);
    CHECK(correction_error.maximum<2e-6);CHECK(value_error.maximum<2e-6);
  }
}
