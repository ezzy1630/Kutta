// Exact-state term isolation: this benchmark never changes SteadySolver.
#include <iostream>
#include "../tests/fv/physics_meshes.hpp"
#include "kernels.hpp"
#include "kutta/fv/transport.hpp"
using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
struct ErrorNorm {
  double square=0,weight=0;
  void add(double error,double volume) {square+=volume*error*error;weight+=volume;}
  double rms() const {return std::sqrt(square/weight);}
};
}
int main() {
  for(auto backend:Device::available()) for(int nr:{8,16,32}) {
    auto dev=Device::create(backend);
    const auto mesh=test::physics_annulus(nr,8*nr,true);const auto geo=compute_geometry(mesh);
    const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
    std::vector<float> ux,uy,bx,by,mass(m.n_faces,0);std::vector<Index> kind(m.n_boundary_faces(),1);
    auto exact=[](Vec3d x) {const double angular=-1./3+4./(3*(x.x*x.x+x.y*x.y));return Vec3d{-angular*x.y,angular*x.x,0};};
    auto pressure=[](Vec3d x) {const double r=std::hypot(x.x,x.y);return r*r/18.-8./9*std::log(r)-8./(9*r*r);};
    for(auto x:geo.cell_centroid) {auto u=exact(x);ux.push_back(float(u.x));uy.push_back(float(u.y));}
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) {auto u=exact(geo.face_centroid[f]);bx.push_back(float(u.x));by.push_back(float(u.y));}
    for(const auto& patch:mesh.patches) if(patch.name=="inner" || patch.name=="outer")
      for(Index f=patch.start;f<patch.start+patch.size;++f) kind[f-m.n_internal_faces]=0;
    for(Index f=0;f<m.n_internal_faces;++f) mass[f]=float(dot(exact(geo.face_centroid[f]),geo.face_area[f]));
    auto u=dev->upload<float>(ux),v=dev->upload<float>(uy),bu=dev->upload<float>(bx),bv=dev->upload<float>(by),flux=dev->upload<float>(mass);auto bk=dev->upload<Index>(kind);
    auto zero=dev->alloc<float>(m.n_cells),mu=dev->alloc<float>(m.n_faces),res=dev->alloc<float>(m.n_cells);dev->fill(zero,0);
    // Compare pressure-force truncation with diffusion on the identical state.
    std::vector<float> pc,pb,pgx,pgy;
    for(auto x:geo.cell_centroid) {pc.push_back(float(pressure(x)));const double a=-1./3+4./(3*(x.x*x.x+x.y*x.y));pgx.push_back(float(a*a*x.x));pgy.push_back(float(a*a*x.y));}
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) pb.push_back(float(pressure(geo.face_centroid[f])));
    auto pp=dev->upload<float>(pc),bp=dev->upload<float>(pb),pf=dev->alloc<float>(m.n_faces),pj=dev->alloc<float>(m.n_faces);
    auto pg=vector_field(*dev,m.n_cells),force=vector_field(*dev,m.n_cells);
    std::vector<BoundaryObservation> omitted(m.n_boundary_faces(),BoundaryObservation::Omitted);
    const auto pressure_lsq=boundary_least_squares(*dev,mesh,geo,omitted);
    for(int mode=0;mode<3;++mode) {
      if(mode==0) gradient(*dev,m,st,pressure_lsq,pp,bp,pg);
      else {dev->upload<float>(pg.x,pgx);dev->upload<float>(pg.y,pgy);dev->fill(pg.z,0);}
      kern::fv_pressure_faces(*dev,m.n_faces,m.owner,m.neighbour,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
        m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,pp,pg.x,pg.y,pg.z,bp,pf,pj,m.n_internal_faces);
      if(mode==2) {std::vector<float> exact_faces;for(auto x:geo.face_centroid) exact_faces.push_back(float(pressure(x)));dev->upload<float>(pf,exact_faces);}
      kern::fv_momentum_source(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
        m.area.x,m.area.y,m.area.z,pf,m.volume,zero,zero,zero,force.x,force.y,force.z);
      const auto fx=dev->download(force.x),fy=dev->download(force.y);
      ErrorNorm radial,tangential;
      for(Index c=0;c<m.n_cells;++c) {const auto x=geo.cell_centroid[c];const double r=std::hypot(x.x,x.y),vol=geo.cell_volume[c];
        radial.add(-(fx[c]*x.x+fy[c]*x.y)/r-std::pow(norm(exact(x)),2)/r,vol);
        tangential.add((-fx[c]*x.y+fy[c]*x.x)/r,vol);}
      std::cout<<"POLY_PRESSURE_FORCE backend="<<backend_name(dev->backend())<<" nr="<<nr<<" mode="<<mode
        <<" radial_defect="<<radial.rms()<<" tangential_defect="<<tangential.rms()<<std::endl;
    }

    TransportSettings settings;settings.limit=false;
    std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::Omitted);
    for(Index i=0;i<m.n_boundary_faces();++i) if(kind[i]==0) rows[i]=BoundaryObservation::FaceValue;
    const auto lsq=boundary_least_squares(*dev,mesh,geo,rows);
    TransportEquation equation(*dev,m,st,mesh,geo);
    for(int gradient_mode:{0,1,2}) for(bool convection:{false,true}) {
      std::vector<std::vector<float>> defect;
      for(int component=0;component<2;++component) {
        const auto& field=component==0?u:v;const auto& boundary=component==0?bu:bv;
        if(gradient_mode==1) {
          std::vector<float> dx,dy;
          for(auto x:geo.cell_centroid) {
            const double r2=x.x*x.x+x.y*x.y,a=-1./3+4./(3*r2),ax=-8*x.x/(3*r2*r2),ay=-8*x.y/(3*r2*r2);
            dx.push_back(float(component==0?-x.y*ax:a+x.x*ax));dy.push_back(float(component==0?-a-x.y*ay:x.x*ay));
          }
          dev->upload<float>(equation.grad.x,dx);dev->upload<float>(equation.grad.y,dy);dev->fill(equation.grad.z,0);
        } else gradient(*dev,m,st,lsq,field,boundary,equation.grad);
        if(gradient_mode==2) {
          kern::fv_pressure_faces(*dev,m.n_faces,m.owner,m.neighbour,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
            m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,field,equation.grad.x,equation.grad.y,equation.grad.z,boundary,pf,pj,m.n_internal_faces);
          kern::fv_momentum_source(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
            m.area.x,m.area.y,m.area.z,pf,m.volume,zero,zero,zero,force.x,force.y,force.z);
          kern::scale_copy(*dev,m.n_cells,equation.grad.x,force.x,-1);kern::scale_copy(*dev,m.n_cells,equation.grad.y,force.y,-1);kern::scale_copy(*dev,m.n_cells,equation.grad.z,force.z,-1);
        }
        dev->fill(mu,convection?0:1);
        if(!convection) dev->fill(flux,0);else dev->upload<float>(flux,mass);
        equation.assemble_reconstructed(field,boundary,bk,flux,mu,zero,settings);
        residual(*dev,equation.matrix,field,equation.rhs,res);defect.push_back(dev->download(res));
      }
      ErrorNorm radial,tangential;
      for(Index c=0;c<m.n_cells;++c) {
        const auto x=geo.cell_centroid[c];const double r=std::hypot(x.x,x.y),vol=geo.cell_volume[c];
        // residual() is b-Au: diffusion, or negative convection.
        const Vec3d force{defect[0][c]/vol,defect[1][c]/vol,0};
        const double centrifugal=std::pow(norm(exact(x)),2)/r;
        radial.add(dot(force,Vec3d{x.x/r,x.y/r,0})-(convection?centrifugal:0),vol);
        tangential.add(dot(force,Vec3d{-x.y/r,x.x/r,0}),vol);
      }
      std::cout<<"POLY_MOMENTUM backend="<<backend_name(dev->backend())<<" nr="<<nr<<" gradient_mode="<<gradient_mode<<" convection="<<convection
        <<" radial_defect="<<radial.rms()<<" tangential_defect="<<tangential.rms()<<std::endl;
    }
  }
}
