// Manufactured divergence-free velocity (x^2, -2xy, 0). Production diffusion
// kernels, independent analytic face flux, and the same hex/poly refinement.
#include <array>
#include <iomanip>
#include <iostream>
#include "../tests/fv/physics_meshes.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/mesh/topology.hpp"
using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
struct Norm {
  double square=0,weight=0,maximum=0;
  void add(double e,double a) {square+=a*e*e;weight+=a;maximum=std::max(maximum,std::abs(e));}
};
}
int main(int argc,char** argv) {
  if(argc>2 || (argc==2 && std::string_view(argv[1])!="--baseline")) {std::cerr<<"Usage: bench_quadratic_viscous [--baseline]\n";return 2;}
  const bool baseline=argc==2;
  constexpr const char* names[]{"azimuthal","radial_regular","radial_hanging","inner_wall","outer_wall","cap"};
  std::cout<<std::setprecision(10);
  for(auto backend:Device::available()) for(bool poly:{false,true}) for(int nr:{8,16,32}) {
    auto dev=Device::create(backend);const auto host=test::physics_annulus(nr,8*nr,poly);const auto geo=compute_geometry(host);
    const auto m=upload_mesh(*dev,host,geo);const auto stencil=make_stencil(*dev,host,geo);const auto cf=build_cell_faces(host);
    std::vector<int> types(m.n_faces);std::vector<Index> kinds(m.n_boundary_faces(),0);
    std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    for(Index f=0;f<m.n_faces;++f) {
      const auto a=geo.face_area[f],x=geo.face_centroid[f];
      if(f<m.n_internal_faces) {
        const bool radial=std::abs(dot(a,Vec3d{x.x,x.y,0}))>.8*norm(a)*std::hypot(x.x,x.y);
        const auto p=host.owner[f],n=host.neighbour[f];
        types[f]=!radial?0:(cf.offsets[p+1]-cf.offsets[p]==cf.offsets[n+1]-cf.offsets[n]?1:2);
      } else {
        types[f]=std::abs(a.z)>.8*norm(a)?5:std::hypot(x.x,x.y)<1.5?3:4;
        if(types[f]==5) {kinds[f-m.n_internal_faces]=2;rows[f-m.n_internal_faces]=BoundaryObservation::NormalFootValue;}
      }
    }
    const auto fit=baseline?std::shared_ptr<ProjectedFaceReconstruction>{}:std::make_shared<ProjectedFaceReconstruction>(*dev,host,geo,rows);
    TransportEquation equation(*dev,m,stencil,host,geo,{},fit);TransportSettings controls;controls.limit=false;
    auto kind=dev->upload<Index>(kinds);auto zero=dev->alloc<float>(m.n_cells),mass=dev->alloc<float>(m.n_faces),mu=dev->alloc<float>(m.n_faces);
    dev->fill(zero,0);dev->fill(mass,0);dev->fill(mu,1);const auto lsq=boundary_least_squares(*dev,host,geo,rows);
    const auto diffusion=dev->download(stencil.diffusion);
    for(int component=0;component<2;++component) for(bool exact:{false,true}) {
      auto value=[&](Vec3d x){return component==0?x.x*x.x:-2*x.x*x.y;};
      auto derivative=[&](Vec3d x){return component==0?Vec3d{2*x.x,0,0}:Vec3d{-2*x.y,-2*x.x,0};};
      std::vector<float> cells,boundary,dx,dy;
      for(auto x:geo.cell_centroid) {cells.push_back(float(value(x)));const auto g=derivative(x);dx.push_back(float(g.x));dy.push_back(float(g.y));}
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(geo.face_centroid[f])));
      auto q=dev->upload<float>(cells),bc=dev->upload<float>(boundary);
      if(exact) {dev->upload<float>(equation.grad.x,dx);dev->upload<float>(equation.grad.y,dy);dev->fill(equation.grad.z,0);}
      else gradient(*dev,m,stencil,lsq,q,bc,equation.grad);
      equation.assemble_reconstructed(q,bc,kind,mass,mu,zero,controls);
      const auto correction=dev->download(equation.face_correction);std::array<Norm,6> errors;
      for(Index f=0;f<m.n_faces;++f) {
        const auto a=geo.face_area[f];double flux=0;
        if(types[f]!=5) {
          const double jump=(f<m.n_internal_faces?cells[host.neighbour[f]]:boundary[f-m.n_internal_faces])-double(cells[host.owner[f]]);
          flux=diffusion[f]*jump+correction[f];
        }
        errors[types[f]].add((flux-dot(derivative(geo.face_centroid[f]),a))/norm(a),norm(a));
      }
      for(int type=0;type<6;++type) if(errors[type].weight)
        std::cout<<"QUADRATIC_OPERATOR backend="<<backend_name(backend)<<" poly="<<poly<<" nr="<<nr<<" component="<<component
          <<" gradient="<<(exact?"exact":"lsq")<<" operator=viscous_normal_derivative location="<<names[type]
          <<" rms="<<std::sqrt(errors[type].square/errors[type].weight)<<" maximum="<<errors[type].maximum<<'\n';
    }
  }
}
