// Exact-state operator errors. All discretization calls use production kernels;
// the host only supplies analytic observations and reduces diagnostic errors.
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include "../tests/fv/physics_meshes.hpp"
#include "kernels.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/mesh/topology.hpp"

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
struct Norm {
  double square=0,weight=0,maximum=0;
  void add(double value,double w) {
    square+=w*value*value;weight+=w;maximum=std::isfinite(value)?std::max(maximum,std::abs(value)):std::numeric_limits<double>::infinity();
  }
  double rms() const {return std::sqrt(square/weight);}
};
Vec3d velocity(Vec3d x) {
  const double a=-1./3+4./(3*(x.x*x.x+x.y*x.y));return {-a*x.y,a*x.x,0};
}
Vec3d derivative(Vec3d x,int component) {
  const double r2=x.x*x.x+x.y*x.y,a=-1./3+4./(3*r2),ax=-8*x.x/(3*r2*r2),ay=-8*x.y/(3*r2*r2);
  return component==0?Vec3d{-x.y*ax,-a-x.y*ay,0}:Vec3d{a+x.x*ax,x.x*ay,0};
}
double pressure(Vec3d x) {
  const double r=std::hypot(x.x,x.y);return r*r/18.-8./9*std::log(r)-8./(9*r*r);
}
Vec3d pressure_gradient(Vec3d x) {
  const double a=-1./3+4./(3*(x.x*x.x+x.y*x.y));return {a*a*x.x,a*a*x.y,0};
}
std::vector<Vec3d> download(Device& d,const DeviceVec3& v) {
  const auto x=d.download(v.x),y=d.download(v.y),z=d.download(v.z);
  std::vector<Vec3d> out(x.size());for(std::size_t i=0;i<x.size();++i) out[i]={x[i],y[i],z[i]};return out;
}
void upload(Device& d,const std::vector<Vec3d>& values,DeviceVec3& out) {
  std::vector<float> x,y,z;
  for(auto p:values) {x.push_back(float(p.x));y.push_back(float(p.y));z.push_back(float(p.z));}
  d.upload<float>(out.x,x);d.upload<float>(out.y,y);d.upload<float>(out.z,z);
}
}
int main(int argc,char** argv) {
  const bool baseline=argc>1 && std::string_view(argv[1])=="--baseline";
  std::cout<<std::setprecision(10);
  constexpr const char* face_names[]{"azimuthal","radial_regular","radial_hanging","inner_wall","outer_wall","cap"};
  constexpr const char* cell_names[]{"interior","inner_wall","outer_wall"};
  for(auto backend:Device::available()) for(bool poly:{false,true}) for(int nr:{8,16,32}) {
    auto dev=Device::create(backend);
    const auto host=test::physics_annulus(nr,8*nr,poly);const auto geo=compute_geometry(host);
    const auto m=upload_mesh(*dev,host,geo);const auto st=make_stencil(*dev,host,geo);
    const auto cf=build_cell_faces(host);
    std::vector<int> faces(m.n_faces),cells(m.n_cells,0);
    std::vector<Index> scalar_kind(m.n_boundary_faces(),0),flow_kind(m.n_boundary_faces(),0);
    std::vector<BoundaryObservation> velocity_rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    for(Index f=0;f<m.n_faces;++f) {
      const auto x=geo.face_centroid[f],a=geo.face_area[f];
      if(f<m.n_internal_faces) {
        const bool radial=std::abs(dot(a,Vec3d{x.x,x.y,0}))>.8*norm(a)*std::hypot(x.x,x.y);
        const auto p=host.owner[f],n=host.neighbour[f];
        faces[f]=!radial?0:(cf.offsets[p+1]-cf.offsets[p]==cf.offsets[n+1]-cf.offsets[n]?1:2);
      } else {
        faces[f]=std::abs(a.z)>.8*norm(a)?5:std::hypot(x.x,x.y)<1.5?3:4;
        if(faces[f]==5) {
          scalar_kind[f-m.n_internal_faces]=2;flow_kind[f-m.n_internal_faces]=3;
          velocity_rows[f-m.n_internal_faces]=BoundaryObservation::NormalFootValue;
        } else cells[host.owner[f]]=faces[f]-2;
      }
    }
    auto integrate=[&](Index f,auto function) {
      const auto start=host.face_offsets[f],end=host.face_offsets[f+1];
      const auto a=host.points[host.face_nodes[start]];Vec3d sum{};
      constexpr double w[]{.225,.132394152788506,.132394152788506,.132394152788506,.125939180544827,.125939180544827,.125939180544827};
      constexpr double r[]{1./3,.059715871789770,.470142064105115,.470142064105115,.797426985353087,.101286507323456,.101286507323456};
      constexpr double t[]{1./3,.470142064105115,.059715871789770,.470142064105115,.101286507323456,.797426985353087,.101286507323456};
      for(Index j=start+1;j+1<end;++j) {
        const auto b=host.points[host.face_nodes[j]],c=host.points[host.face_nodes[j+1]];
        const auto normal=geo.face_area[f]/norm(geo.face_area[f]);
        const double size=.5*dot(cross(b-a,c-a),normal);if(size==0) continue;
        for(int k=0;k<7;++k) sum=sum+function(a*r[k]+b*t[k]+c*(1-r[k]-t[k]),normal)*(size*w[k]);
      }
      return sum;
    };
    std::vector<Vec3d> integrated_viscous(m.n_faces),integrated_convection(m.n_faces),integrated_pressure(m.n_faces),reference_gradient(m.n_cells);
    for(Index f=0;f<m.n_faces;++f) {
      integrated_viscous[f]=integrate(f,[](Vec3d x,Vec3d n){return Vec3d{dot(derivative(x,0),n),dot(derivative(x,1),n),0};});
      integrated_convection[f]=integrate(f,[](Vec3d x,Vec3d n){const auto u=velocity(x);return u*dot(u,n);});
      integrated_pressure[f]=integrate(f,[](Vec3d x,Vec3d n){return n*pressure(x);});
      reference_gradient[host.owner[f]]=reference_gradient[host.owner[f]]+integrated_pressure[f];
      if(f<m.n_internal_faces) reference_gradient[host.neighbour[f]]=reference_gradient[host.neighbour[f]]-integrated_pressure[f];
    }
    for(Index c=0;c<m.n_cells;++c) reference_gradient[c]=reference_gradient[c]/geo.cell_volume[c];
    auto sk=dev->upload<Index>(scalar_kind),fk=dev->upload<Index>(flow_kind);
    auto zero=dev->alloc<float>(m.n_cells),mu=dev->alloc<float>(m.n_faces),mass=dev->alloc<float>(m.n_faces);
    dev->fill(zero,0);
    auto q=vector_field(*dev,m.n_cells),bc=vector_field(*dev,m.n_boundary_faces());
    std::vector<Vec3d> qhost;std::vector<float> phost,pbhost,mhost(m.n_faces,0);
    for(auto x:geo.cell_centroid) {qhost.push_back(velocity(x));phost.push_back(float(pressure(x)));}
    upload(*dev,qhost,q);
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) pbhost.push_back(float(pressure(geo.face_centroid[f])));
    for(Index f=0;f<m.n_internal_faces;++f) mhost[f]=float(dot(velocity(geo.face_centroid[f]),geo.face_area[f]));
    auto p=dev->upload<float>(phost),pb=dev->upload<float>(pbhost),boundary_p=dev->alloc<float>(m.n_boundary_faces());
    auto pf=dev->alloc<float>(m.n_faces),jump=dev->alloc<float>(m.n_faces);
    auto pg=vector_field(*dev,m.n_cells),force=vector_field(*dev,m.n_cells);
    const auto vfit=boundary_least_squares(*dev,host,geo,velocity_rows);
    const std::vector<BoundaryObservation> omitted(m.n_boundary_faces(),BoundaryObservation::Omitted);
    const auto pfit=boundary_least_squares(*dev,host,geo,omitted);
    const auto reconstruction=baseline?std::shared_ptr<ProjectedFaceReconstruction>{}:std::make_shared<ProjectedFaceReconstruction>(*dev,host,geo,velocity_rows);
    std::cout<<"TC_RECONSTRUCTION backend="<<backend_name(backend)<<" poly="<<poly<<" nr="<<nr
      <<" corrected="<<(reconstruction?reconstruction->corrected_faces():0)<<" fallback="<<(reconstruction?reconstruction->fallback_faces():0)<<'\n';
    TransportEquation equation(*dev,m,st,host,geo,{},reconstruction);
    auto gu=vector_field(*dev,m.n_cells),gv=vector_field(*dev,m.n_cells),gw=vector_field(*dev,m.n_cells);
    auto face_velocity=vector_field(*dev,m.n_faces),face_zero=vector_field(*dev,m.n_faces);
    dev->fill(face_zero.x,0);dev->fill(face_zero.y,0);dev->fill(face_zero.z,0);
    auto rc_flux=dev->alloc<float>(m.n_faces),coefficient=dev->alloc<float>(m.n_faces);
    for(bool actual_wall:{false,true}) for(bool exact_grad:{false,true}) {
      std::vector<Vec3d> boundary;
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) {
        const auto x=geo.face_centroid[f],a=geo.face_area[f];auto u=velocity(x);
        if(actual_wall && faces[f]!=5) u=faces[f]==3?Vec3d{-x.y,x.x,0}:Vec3d{};
        u=u-a*(dot(u,a)/dot(a,a));boundary.push_back(u);
      }
      upload(*dev,boundary,bc);
      if(exact_grad) {
        std::vector<Vec3d> values[3];
        for(auto x:geo.cell_centroid) {values[0].push_back(derivative(x,0));values[1].push_back(derivative(x,1));values[2].push_back({});}
        upload(*dev,values[0],gu);upload(*dev,values[1],gv);upload(*dev,values[2],gw);
        std::vector<Vec3d> g;for(auto x:geo.cell_centroid) g.push_back(pressure_gradient(x));upload(*dev,g,pg);
      } else {
        vector_gradient(*dev,m,st,vfit,q,bc,gu,gv,gw);gradient(*dev,m,st,pfit,p,pb,pg);
      }
      auto label=[&](const char* op,const char* location,const Norm& error) {
        if(!error.weight) return;
        std::cout<<"TC_OPERATOR backend="<<backend_name(backend)<<" poly="<<poly<<" nr="<<nr
          <<" wall_data="<<(actual_wall?"physical":"analytic")<<" gradient="<<(exact_grad?"exact":"lsq")
          <<" operator="<<op<<" location="<<location<<" rms="<<error.rms()<<" maximum="<<error.maximum<<'\n';
      };
      auto cell_errors=[&](const char* op,const std::vector<Vec3d>& values,auto exact) {
        std::array<Norm,3> errors;
        for(Index c=0;c<m.n_cells;++c) errors[cells[c]].add(norm(values[c]-exact(geo.cell_centroid[c])),geo.cell_volume[c]);
        for(int k=0;k<3;++k) label(op,cell_names[k],errors[k]);
      };
      cell_errors("pressure_lsq_gradient",download(*dev,pg),pressure_gradient);
      kern::fv_pressure_boundary(*dev,m.n_boundary_faces(),m.owner,fk,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
        p,pg.x,pg.y,pg.z,pb,boundary_p,m.n_internal_faces);
      kern::fv_pressure_faces(*dev,m.n_faces,m.owner,m.neighbour,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
        m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,p,pg.x,pg.y,pg.z,boundary_p,pf,jump,m.n_internal_faces);
      const auto pressures=dev->download(pf),jumps=dev->download(jump),diffusion=dev->download(st.diffusion);
      std::array<Norm,6> pressure_errors,jump_errors,skew_errors;
      for(Index f=0;f<m.n_faces;++f) {
        const auto a=geo.face_area[f],x=geo.face_centroid[f];
        pressure_errors[faces[f]].add(pressures[f]-pressure(x),norm(a));
        if(f<m.n_internal_faces) jump_errors[faces[f]].add(jumps[f]/norm(geo.cell_centroid[host.neighbour[f]]-geo.cell_centroid[host.owner[f]]),norm(a));
        else {
          const auto r=x-geo.cell_centroid[host.owner[f]],n=a/norm(a);skew_errors[faces[f]].add(norm(r-n*dot(r,n))/std::abs(dot(r,n)),norm(a));
        }
      }
      for(int k=0;k<6;++k) {label("pressure_face",face_names[k],pressure_errors[k]);label("rc_jump_over_link",face_names[k],jump_errors[k]);label("boundary_tangent_over_normal",face_names[k],skew_errors[k]);}
      kern::fv_momentum_source(*dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,m.area.x,m.area.y,m.area.z,pf,m.volume,
        zero,zero,zero,force.x,force.y,force.z);
      auto pressure_force=download(*dev,force);for(auto& x:pressure_force) x=x*(-1);
      cell_errors("pressure_gauss_gradient",pressure_force,pressure_gradient);

      const Buffer<float>* field[]{&q.x,&q.y};const Buffer<float>* bfield[]{&bc.x,&bc.y};const DeviceVec3* gfield[]{&gu,&gv};
      std::vector<Vec3d> viscous(m.n_faces),convection(m.n_faces),viscous_exact(m.n_faces);
      TransportSettings controls;controls.limit=false;
      for(int component=0;component<2;++component) {
        dev->upload<float>(equation.grad.x,dev->download(gfield[component]->x));
        dev->upload<float>(equation.grad.y,dev->download(gfield[component]->y));
        dev->fill(equation.grad.z,0);
        for(bool conv:{false,true}) {
          dev->fill(mu,conv?0:1);if(conv) dev->upload<float>(mass,mhost);else dev->fill(mass,0);
          equation.assemble_reconstructed(*field[component],*bfield[component],sk,mass,mu,zero,controls);
          const auto correction=dev->download(equation.face_correction),values=dev->download(*field[component]),bvalues=dev->download(*bfield[component]);
          for(Index f=0;f<m.n_faces;++f) {
            const auto c=host.owner[f];double flux=0;
            if(conv) {
              const auto up=f<m.n_internal_faces && mhost[f]<0?host.neighbour[f]:c;
              flux=mhost[f]*values[up]-correction[f];
            } else if(f<m.n_internal_faces || scalar_kind[f-m.n_internal_faces]==0) {
              flux=diffusion[f]*((f<m.n_internal_faces?values[host.neighbour[f]]:bvalues[f-m.n_internal_faces])-values[c])+correction[f];
            }
            auto& target=conv?convection[f]:viscous[f];if(component==0) target.x=flux;else target.y=flux;
            const double exact_flux=dot(derivative(geo.face_centroid[f],component),geo.face_area[f]);
            if(component==0) viscous_exact[f].x=exact_flux;else viscous_exact[f].y=exact_flux;
          }
        }
      }
      std::array<Norm,6> viscous_errors,convection_errors,integrated_viscous_errors,integrated_convection_errors;
      std::vector<Vec3d> viscous_div(m.n_cells),convective_div(m.n_cells),interior_error(m.n_cells),wall_error(m.n_cells);
      for(Index f=0;f<m.n_faces;++f) {
        const auto a=geo.face_area[f];const double area=norm(a);const auto expected=velocity(geo.face_centroid[f]);
        integrated_viscous_errors[faces[f]].add(norm(viscous[f]-integrated_viscous[f])/area,area);
        integrated_convection_errors[faces[f]].add(norm(convection[f]-integrated_convection[f])/area,area);
        viscous_errors[faces[f]].add(norm(viscous[f]-viscous_exact[f])/area,area);
        convection_errors[faces[f]].add(norm(convection[f]-expected*dot(expected,a))/area,area);
        const auto p0=host.owner[f];viscous_div[p0]=viscous_div[p0]+viscous[f];convective_div[p0]=convective_div[p0]-convection[f];
        auto& error=f<m.n_internal_faces?interior_error:wall_error;error[p0]=error[p0]+viscous[f]-viscous_exact[f];
        if(f<m.n_internal_faces) {const auto n=host.neighbour[f];viscous_div[n]=viscous_div[n]-viscous[f];convective_div[n]=convective_div[n]+convection[f];error[n]=error[n]-(viscous[f]-viscous_exact[f]);}
      }
      for(Index c=0;c<m.n_cells;++c) {const double v=geo.cell_volume[c];viscous_div[c]=viscous_div[c]/v;convective_div[c]=convective_div[c]/v;interior_error[c]=interior_error[c]/v;wall_error[c]=wall_error[c]/v;}
      for(int k=0;k<6;++k) {label("viscous_integrated_normal_flux",face_names[k],integrated_viscous_errors[k]);label("convective_integrated_normal_flux",face_names[k],integrated_convection_errors[k]);label("viscous_normal_derivative",face_names[k],viscous_errors[k]);label("convective_normal_flux",face_names[k],convection_errors[k]);}
      const auto null=[](Vec3d){return Vec3d{};};
      cell_errors("viscous_divergence",viscous_div,null);cell_errors("viscous_interior_face_error_divergence",interior_error,null);cell_errors("viscous_wall_face_error_divergence",wall_error,null);
      cell_errors("negative_convective_divergence",convective_div,pressure_gradient);
      std::array<Norm,3> pressure_truncation,convection_truncation;
      for(Index c=0;c<m.n_cells;++c) {
        pressure_truncation[cells[c]].add(norm(pressure_force[c]-reference_gradient[c]),geo.cell_volume[c]);
        convection_truncation[cells[c]].add(norm(convective_div[c]-reference_gradient[c]),geo.cell_volume[c]);
      }
      for(int k=0;k<3;++k) {label("pressure_fv_truncation",cell_names[k],pressure_truncation[k]);label("convection_fv_truncation",cell_names[k],convection_truncation[k]);}


      // Use the assembled physical diagonal, including convection and diffusion,
      // for the actual Rhie-Chow flux. No pseudo-time or relaxation contribution.
      dev->fill(mu,1);dev->upload<float>(mass,mhost);
      dev->copy(equation.grad.x,gu.x);dev->copy(equation.grad.y,gu.y);dev->copy(equation.grad.z,gu.z);
      equation.assemble_reconstructed(q.x,bc.x,sk,mass,mu,zero,controls);
      kern::fv_momentum_faces(*dev,m.n_faces,m.owner,m.neighbour,fk,
        m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
        q.x,q.y,q.z,gu.x,gu.y,gu.z,gv.x,gv.y,gv.z,gw.x,gw.y,gw.z,equation.diagonal,bc.x,bc.y,bc.z,
        face_velocity.x,face_velocity.y,face_velocity.z,m.n_internal_faces);
      if(reconstruction) {
        reconstruction->momentum_values(q.x,bc.x,equation.diagonal,face_velocity.x);
        reconstruction->momentum_values(q.y,bc.y,equation.diagonal,face_velocity.y);
        reconstruction->momentum_values(q.z,bc.z,equation.diagonal,face_velocity.z);
      }
      kern::fv_rhie_chow(*dev,m.n_faces,m.owner,m.neighbour,fk,m.area.x,m.area.y,m.area.z,st.diffusion,m.volume,equation.diagonal,equation.consistent_inverse,
        face_zero.x,face_zero.y,face_zero.z,jump,rc_flux,coefficient,m.n_internal_faces,1);
      const auto predictor=download(*dev,face_velocity);const auto rc=dev->download(rc_flux);
      std::array<Norm,6> predictor_errors,rc_errors,combined_errors;
      for(Index f=0;f<m.n_internal_faces;++f) {
        const auto a=geo.face_area[f];const double area=norm(a),e=dot(predictor[f]-velocity(geo.face_centroid[f]),a)/area;
        predictor_errors[faces[f]].add(e,area);rc_errors[faces[f]].add(rc[f]/area,area);combined_errors[faces[f]].add(e+rc[f]/area,area);
      }
      for(int k=0;k<6;++k) {label("momentum_face_normal_velocity",face_names[k],predictor_errors[k]);label("rc_normal_velocity",face_names[k],rc_errors[k]);label("combined_face_normal_velocity",face_names[k],combined_errors[k]);}
      std::vector<float> coupled_mass(m.n_faces,0);std::vector<double> continuity(m.n_cells,0);
      for(Index f=0;f<m.n_internal_faces;++f) {
        coupled_mass[f]=float(dot(predictor[f],geo.face_area[f])+rc[f]);
        continuity[host.owner[f]]+=coupled_mass[f];continuity[host.neighbour[f]]-=coupled_mass[f];
      }
      std::array<Norm,3> continuity_error;
      for(Index c=0;c<m.n_cells;++c) continuity_error[cells[c]].add(continuity[c]/geo.cell_volume[c],geo.cell_volume[c]);
      for(int k=0;k<3;++k) label("continuity_divergence",cell_names[k],continuity_error[k]);
      std::vector<Vec3d> coupled_convection(m.n_faces),coupled_divergence(m.n_cells);
      dev->fill(mu,0);dev->upload<float>(mass,coupled_mass);
      for(int component=0;component<2;++component) {
        dev->copy(equation.grad.x,gfield[component]->x);dev->copy(equation.grad.y,gfield[component]->y);dev->copy(equation.grad.z,gfield[component]->z);
        equation.assemble_reconstructed(*field[component],*bfield[component],sk,mass,mu,zero,controls);
        const auto correction=dev->download(equation.face_correction),values=dev->download(*field[component]);
        for(Index f=0;f<m.n_internal_faces;++f) {
          const auto up=coupled_mass[f]<0?host.neighbour[f]:host.owner[f];const double flux=coupled_mass[f]*values[up]-correction[f];
          if(component==0) coupled_convection[f].x=flux;else coupled_convection[f].y=flux;
        }
      }
      std::array<Norm,6> coupled_face_error;std::array<Norm,3> coupled_truncation;
      for(Index f=0;f<m.n_faces;++f) {
        const double area=norm(geo.face_area[f]);coupled_face_error[faces[f]].add(norm(coupled_convection[f]-integrated_convection[f])/area,area);
        coupled_divergence[host.owner[f]]=coupled_divergence[host.owner[f]]-coupled_convection[f];
        if(f<m.n_internal_faces) coupled_divergence[host.neighbour[f]]=coupled_divergence[host.neighbour[f]]+coupled_convection[f];
      }
      for(Index c=0;c<m.n_cells;++c) coupled_truncation[cells[c]].add(norm(coupled_divergence[c]/geo.cell_volume[c]-reference_gradient[c]),geo.cell_volume[c]);
      for(int k=0;k<6;++k) label("convective_coupled_normal_flux",face_names[k],coupled_face_error[k]);
      for(int k=0;k<3;++k) label("convection_coupled_fv_truncation",cell_names[k],coupled_truncation[k]);

    }
  }
}
