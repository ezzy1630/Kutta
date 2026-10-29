#include <cmath>
#include <iostream>
#include <numbers>

#include "../backend_helpers.hpp"
#include "kernels.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/io/mesh_io.hpp"

using namespace kutta;

TEST_CASE("FV MMS pressure LSQ gradient convergence") {
  const double pi=std::numbers::pi;
  const auto exact=[pi](Vec3d p) {return float(std::cos(pi*p.x)*std::cos(pi*p.y));};
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    double old=0;
    for(int n:{4,8,16}) {
      auto mesh=tet?make_box_tet(n,n,n,{0,0,0},{1,1,1}):make_box_hex(n,n,n,{0,0,0},{1,1,1});
      auto g=compute_geometry(mesh);auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);
      std::vector<float> cell,bc;
      for(auto p:g.cell_centroid) cell.push_back(exact(p));
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) bc.push_back(exact(g.face_centroid[f]));
      auto q=dev->upload<float>(cell),b=dev->upload<float>(bc);auto gradient=fv::vector_field(*dev,m.n_cells);
      fv::gradient(*dev,m,s,q,b,gradient);
      auto x=dev->download(gradient.x),y=dev->download(gradient.y),z=dev->download(gradient.z);
      auto fp=dev->alloc<float>(m.n_faces),jump=dev->alloc<float>(m.n_faces);
      kern::fv_pressure_faces(*dev,m.n_faces,m.owner,m.neighbour,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
        m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
        q,gradient.x,gradient.y,gradient.z,b,fp,jump,m.n_internal_faces);
      auto zero=dev->alloc<float>(m.n_cells);dev->fill(zero,0);
      auto pressure_force=fv::vector_field(*dev,m.n_cells);
      kern::fv_momentum_source(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
        m.area.x,m.area.y,m.area.z,fp,m.volume,zero,zero,zero,
        pressure_force.x,pressure_force.y,pressure_force.z);
      auto fx=dev->download(pressure_force.x),fy=dev->download(pressure_force.y),fz=dev->download(pressure_force.z);
      double error=0,volume=0,gauss_error=0,difference=0;
      for(Index c=0;c<m.n_cells;++c) {
        const auto p=g.cell_centroid[c];
        const double ex=-pi*std::sin(pi*p.x)*std::cos(pi*p.y),ey=-pi*std::cos(pi*p.x)*std::sin(pi*p.y);
        error+=g.cell_volume[c]*(std::pow(x[c]-ex,2)+std::pow(y[c]-ey,2)+z[c]*z[c]);
        gauss_error+=g.cell_volume[c]*(std::pow(fx[c]+ex,2)+std::pow(fy[c]+ey,2)+fz[c]*fz[c]);
        difference+=g.cell_volume[c]*(std::pow(fx[c]+x[c],2)+std::pow(fy[c]+y[c],2)+std::pow(fz[c]+z[c],2));
        volume+=g.cell_volume[c];
      }
      error=std::sqrt(error/volume);
      std::cout<<"FV_PRESSURE_LSQ backend="<<backend_name(dev->backend())<<" tet="<<tet<<" n="<<n<<" L2="<<error
        <<" momentum_gradient_L2="<<std::sqrt(gauss_error/volume)<<" momentum_minus_RC_gradient_L2="<<std::sqrt(difference/volume);
      if(old>0) {
        const double order=std::log2(old/error);
        std::cout<<" order="<<order;
        if(n==16) CHECK(order>=.9);
      }
      std::cout<<std::endl;old=error;
    }
  }
}
