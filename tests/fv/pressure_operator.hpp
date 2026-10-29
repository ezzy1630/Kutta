#pragma once

#include <cmath>
#include <numbers>
#include "kernels.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::test {
struct PressureOperatorError { double exact_face, momentum, lsq; };
// The actual LSQ -> Second Order face -> momentum-source kernels, supplied
// exact boundary pressure to isolate interior truncation from BC error.
inline PressureOperatorError pressure_operator_error(Device& dev,int n,bool tet) {
  using namespace kutta::detail;
  const double pi=std::numbers::pi;
  const auto exact=[pi](Vec3d p) {return float(std::cos(pi*p.x)*std::cos(pi*p.y));};

  auto mesh=tet?make_box_tet(n,n,n,{0,0,0},{1,1,1}):make_box_hex(n,n,1,{0,0,0},{1,1,1./n});
  auto g=compute_geometry(mesh);auto m=upload_mesh(dev,mesh,g);auto s=fv::make_stencil(dev,mesh,g);
  std::vector<float> cell,bc;
  for(auto p:g.cell_centroid) cell.push_back(exact(p));
  for(Index f=m.n_internal_faces;f<m.n_faces;++f) bc.push_back(exact(g.face_centroid[f]));
  auto q=dev.upload<float>(cell),b=dev.upload<float>(bc);auto gradient=fv::vector_field(dev,m.n_cells);
  fv::gradient(dev,m,s,q,b,gradient);
  auto x=dev.download(gradient.x),y=dev.download(gradient.y),z=dev.download(gradient.z);
  auto fp=dev.alloc<float>(m.n_faces),jump=dev.alloc<float>(m.n_faces);
  kern::fv_pressure_faces(dev,m.n_faces,m.owner,m.neighbour,
    m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
    m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
    q,gradient.x,gradient.y,gradient.z,b,fp,jump,m.n_internal_faces);
  auto zero=dev.alloc<float>(m.n_cells);dev.fill(zero,0);
  auto pressure_force=fv::vector_field(dev,m.n_cells);
  kern::fv_momentum_source(dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
    m.area.x,m.area.y,m.area.z,fp,m.volume,zero,zero,zero,
    pressure_force.x,pressure_force.y,pressure_force.z);
  auto fx=dev.download(pressure_force.x),fy=dev.download(pressure_force.y),fz=dev.download(pressure_force.z);
  // Independent double-precision midpoint-face Gauss oracle. It does not
  // use reconstructed pressure or the device geometry/gradient kernels.
  std::vector<Vec3d> exact_gauss(m.n_cells);
  for(Index face=0;face<m.n_faces;++face) {
    const auto p=g.face_centroid[face];
    const auto force=g.face_area[face]*(std::cos(pi*p.x)*std::cos(pi*p.y));
    exact_gauss[mesh.owner[face]]=exact_gauss[mesh.owner[face]]+force;
    if(face<m.n_internal_faces)
      exact_gauss[mesh.neighbour[face]]=exact_gauss[mesh.neighbour[face]]-force;
  }
  double error=0,volume=0,gauss_error=0,exact_error=0;
  for(Index c=0;c<m.n_cells;++c) {
    const auto p=g.cell_centroid[c];
    const double ex=-pi*std::sin(pi*p.x)*std::cos(pi*p.y),ey=-pi*std::cos(pi*p.x)*std::sin(pi*p.y);
    const auto eg=exact_gauss[c]*(1/g.cell_volume[c]);
    exact_error+=g.cell_volume[c]*(std::pow(eg.x-ex,2)+std::pow(eg.y-ey,2)+eg.z*eg.z);
    error+=g.cell_volume[c]*(std::pow(x[c]-ex,2)+std::pow(y[c]-ey,2)+z[c]*z[c]);
    gauss_error+=g.cell_volume[c]*(std::pow(fx[c]+ex,2)+std::pow(fy[c]+ey,2)+fz[c]*fz[c]);
    volume+=g.cell_volume[c];
  }
  error=std::sqrt(error/volume);
  return {std::sqrt(exact_error/volume),std::sqrt(gauss_error/volume),error};
}
inline double order(double coarse,double fine) {return std::log2(coarse/fine);}
}  // namespace kutta::test
