#pragma once
#include "tmr.hpp"
#include "kernels.hpp"

namespace kutta::validation::ke_channel {
using namespace fv;
inline constexpr double re_tau=587.19;
inline std::vector<std::array<double,2>> dns(const std::filesystem::path& path) {
  std::ifstream input(path);KUTTA_REQUIRE(input,"MKM DNS file missing");std::string line;
  std::vector<std::array<double,2>> result;
  while(std::getline(input,line)) {if(line.empty() || line[0]=='#') continue;std::istringstream in(line);double y,yp,u;KUTTA_REQUIRE(bool(in>>y>>yp>>u),"invalid MKM row");result.push_back({y,u});}
  KUTTA_REQUIRE(result.size()>100 && result.back()[0]==1,"incomplete MKM profile");return result;
}
struct Result {int steps=0;bool converged=false;double seconds=0,defect=0,change=0,tau=0;std::vector<std::array<double,2>> profile;};
// The fully-developed 1D reduction solves the U, k and epsilon transport
// equations on Device with dp/dx=-1, rho=1, delta=1. There is no pressure or
// convection equation in this reduction. It is distinct from a 3D car solve.
inline Result run(Device& d,int count) {
  const auto start=std::chrono::steady_clock::now();
  auto host=make_box_hex(1,count,1,{0,0,0},{1,2,1});
  for(auto& p:host.patches) p.kind=(p.name=="ymin" || p.name=="ymax")?PatchKind::Wall:PatchKind::Symmetry;
  auto geometry=compute_geometry(host);auto mesh=upload_mesh(d,host,geometry);auto stencil=make_stencil(d,host,geometry);
  const auto n=mesh.n_cells,nf=mesh.n_faces,nb=mesh.n_boundary_faces();
  auto scalar=[&](std::size_t size,float value){auto b=d.alloc<float>(size);d.fill(b,value);return b;};
  auto velocity=vector_field(d,n),wall_velocity=vector_field(d,nb),gu=vector_field(d,n),zero_vector=vector_field(d,n);
  for(auto* v:{&velocity,&wall_velocity,&gu,&zero_vector}) {d.fill(v->x,0);d.fill(v->y,0);d.fill(v->z,0);}
  auto k=scalar(n,float(1/std::sqrt(.09))),epsilon=scalar(n,1),rho=scalar(n,1),mu=scalar(n,float(1/re_tau));
  std::vector<float> initial_u(n),initial_e(n);
  for(Index c=0;c<n;++c) {double y=std::min(geometry.cell_centroid[c].y,2-geometry.cell_centroid[c].y);initial_u[c]=float(std::log(9.793*re_tau*y)/.4187);initial_e[c]=float(1/(.4187*y));}
  d.upload<float>(velocity.x,initial_u);d.upload<float>(epsilon,initial_e);
  auto face_mu=scalar(nf,float(1/re_tau)),mass=scalar(nf,0),boundary=scalar(nb,0),body=scalar(n,1);
  auto diagonal=scalar(n,0),raw=scalar(n,0),absolute=scalar(n,0),scale=scalar(n,0),fluent=scalar(n,0);
  auto kind=d.alloc<Index>(nb),wall=d.alloc<Index>(nb),zero=d.alloc<Index>(n);d.fill(kind,1);d.fill(zero,0);
  std::vector<Index> mask(nb,0),flow_kind(nb,3);std::vector<BoundaryObservation> urows(nb,BoundaryObservation::Omitted),trows(nb,BoundaryObservation::Omitted);
  for(const auto& p:host.patches) if(p.kind==PatchKind::Wall) for(Index f=p.start;f<p.start+p.size;++f) {auto b=f-mesh.n_internal_faces;mask[b]=1;flow_kind[b]=0;urows[b]=BoundaryObservation::FaceValue;}
  d.upload<Index>(wall,mask);auto flow=d.upload<Index>(flow_kind);
  auto ulsq=boundary_least_squares(d,host,geometry,urows),tlsq=boundary_least_squares(d,host,geometry,trows);
  TransportEquation momentum(d,mesh,stencil,host,geometry),ke(d,mesh,stencil,host,geometry),ee(d,mesh,stencil,host,geometry);
  RealizableKeModel model(d,mesh);TransportSettings settings;settings.relaxation=.7f;settings.convection=ConvectionScheme::FirstOrderUpwind;
  const auto update=[&] {
    gradient(d,mesh,stencil,ulsq,velocity.x,boundary,gu);
    model.update({velocity,{gu,zero_vector,zero_vector},rho,mu,k,epsilon,wall_velocity});
    const auto& f=model.fields();kern::fv_effective_viscosity(d,nf,mesh.owner,mesh.neighbour,flow,mesh.weight,f.mu_t,face_mu,mesh.n_internal_faces,float(1/re_tau));
    kern::rke_channel_wall(d,n,mesh.cell_face_offsets,mesh.cell_faces,wall,mesh.area.x,mesh.area.y,mesh.area.z,f.wall_coefficient,diagonal,mesh.n_internal_faces);
  };
  const auto assemble=[&](int equation) {
    const auto& f=model.fields();auto& eq=equation==0?momentum:(equation==1?ke:ee);auto& q=equation==0?velocity.x:(equation==1?k:epsilon);
    gradient(d,mesh,stencil,equation==0?ulsq:tlsq,q,boundary,eq.grad);
    eq.assemble_reconstructed(q,boundary,kind,mass,equation==0?face_mu:(equation==1?f.face_diffusivity_k:f.face_diffusivity_epsilon),
      equation==0?body:(equation==1?f.k_source:f.epsilon_source),settings,
      equation==0?TransportTerms{nullptr,&diagonal}:TransportTerms{equation==1?&f.k_sink:&f.epsilon_sink});
    if(equation==2) kern::fv_fixed_rows(d,n,eq.matrix.row_ptr,eq.matrix.diag,f.fixed_epsilon_mask,f.fixed_epsilon,eq.matrix.val,eq.rhs);
  };
  // Stops on convergence or after 10,000 steps, never on wall time.
  Result result;std::vector<float> previous=initial_u;
  for(int step=1;step<=10000;++step) {
    for(int equation=0;equation<3;++equation) {
      update();assemble(equation);auto& eq=equation==0?momentum:(equation==1?ke:ee);auto& q=equation==0?velocity.x:(equation==1?k:epsilon);
      eq.solve(q,{1e-5,1e-12,300});if(equation>0) kern::fv_floor(d,n,q,equation==1?1e-14f:1e-20f);
    }
    if(step%100!=0) continue;
    update();result.defect=0;
    for(int equation=0;equation<3;++equation) {
      assemble(equation);auto& eq=equation==0?momentum:(equation==1?ke:ee);auto& q=equation==0?velocity.x:(equation==1?k:epsilon);
      residual(d,eq.matrix,q,eq.rhs,raw);
      kern::fv_scalar_monitor(d,n,q,eq.diagonal,raw,equation==2?model.fields().fixed_epsilon_mask:zero,mesh.volume,absolute,scale,fluent,1,0);
      const auto numerator=d.sum(absolute,n),denominator=d.sum(scale,n);result.defect=std::max(result.defect,numerator/std::max(denominator,1e-30));
    }
    auto current=d.download(velocity.x);result.change=0;
    for(Index c=0;c<n;++c) result.change=std::max(result.change,std::abs(double(current[c]-previous[c]))/std::max(double(current[c]),1e-20));
    previous=current;result.steps=step;result.converged=result.defect<=1e-5 && result.change<=1e-4;
    std::cout<<"KE channel progress backend="<<backend_name(d.backend())<<" n="<<n<<" steps="<<step<<" defect="<<result.defect<<" change="<<result.change<<std::endl;
    if(result.converged) break;
  }
  update();const auto u=d.download(velocity.x),stress=d.download(model.fields().wall_tau);
  result.tau=0;int walls=0;for(Index b=0;b<nb;++b) if(mask[b]) {result.tau+=stress[b];++walls;}result.tau/=walls;
  for(Index c=0;c<n;++c) if(geometry.cell_centroid[c].y<=1) result.profile.push_back({geometry.cell_centroid[c].y,u[c]});
  std::sort(result.profile.begin(),result.profile.end());d.synchronize();result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return result;
}
} // namespace kutta::validation::ke_channel
