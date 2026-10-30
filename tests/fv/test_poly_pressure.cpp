#include <iostream>
#include <array>
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
  void add(double e,double w=1) {square+=w*e*e;weight+=w;maximum=std::max(maximum,std::abs(e));}
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
      kern::fv_pressure_faces(*dev,m.n_faces,m.owner,m.neighbour,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
        p,grad.x,grad.y,grad.z,b,fp,jump,m.n_internal_faces);
      auto face=dev->download(fp),jumps=dev->download(jump);
      std::array<ErrorNorm,3> face_errors,jump_errors;
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
      std::cout<<"POLY_PRESSURE backend="<<backend_name(dev->backend())<<" nr="<<nr<<" linear="<<linear<<" exact_bc="<<exact_boundary
        <<" grad_interior="<<interior.rms()<<" grad_wall="<<boundary.rms()<<" skew_max="<<max_skew<<" bounded="<<bounded;
      for(int k=0;k<3;++k) if(face_errors[k].weight>0) std::cout<<" face"<<k<<"="<<face_errors[k].rms()<<" jump"<<k<<"="<<jump_errors[k].rms();
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
