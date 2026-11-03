#include <chrono>
#include <iostream>
#include <bit>
#include "kernels.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/topology.hpp"
using namespace kutta;
int main(int argc,char** argv) {
 auto mesh=read_fluent_msh(argv[1]);renumber_cells_rcm(mesh);
 auto geo=compute_geometry(mesh);auto dev=Device::create(Backend::Metal);
 auto m=upload_mesh(*dev,mesh,geo);auto s=fv::make_stencil(*dev,mesh,geo);auto fit=fv::face_jump_least_squares(*dev,mesh,geo);
 std::vector<float> q,bc;
 for(auto p:geo.cell_centroid) q.push_back(float(std::sin(p.x)+.3*p.y*p.z));
 for(Index f=m.n_internal_faces;f<m.n_faces;++f) {auto p=geo.face_centroid[f];bc.push_back(float(std::sin(p.x)+.3*p.y*p.z));}
 auto field=dev->upload<float>(q),boundary=dev->upload<float>(bc);
 auto g=fv::vector_field(*dev,m.n_cells),out=fv::vector_field(*dev,m.n_cells);
 fv::gradient(*dev,m,s,field,boundary,g);
 auto run=[&](bool candidate) {
  if(candidate) fv::jump_corrected_gradient(*dev,m,s,fit,field,boundary,g,out);
  else kern::fv_jump_gradient(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,s.adjacent,m.cell_face_sign,
   fit.delta.x,fit.delta.y,fit.delta.z,fit.inverse,s.face_delta.x,s.face_delta.y,s.face_delta.z,
   m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
   field,boundary,g.x,g.y,g.z,out.x,out.y,out.z,m.n_cells);
 };
 run(false);auto x=dev->download(out.x),y=dev->download(out.y),z=dev->download(out.z);
 run(true);if(x!=dev->download(out.x)||y!=dev->download(out.y)||z!=dev->download(out.z)) return 2;
 for(int rep=0;rep<5;++rep) for(int order=0;order<2;++order) {
  bool candidate=(order+rep)%2;run(candidate);dev->synchronize();
  auto start=std::chrono::steady_clock::now();
  for(int i=0;i<200;++i) run(candidate);dev->synchronize();
  std::cout<<"JUMP_TRIAL rep="<<rep<<" candidate="<<candidate<<" cells="<<m.n_cells<<" ms="<<std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count()/200<<std::endl;
 }
}
