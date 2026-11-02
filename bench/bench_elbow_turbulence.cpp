// Saved-state constitutive diagnostics and regional field comparisons on the
// exact Fluent mesh. CSV cell order is the original source mesh order.
#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
std::vector<double> original(Device& d,const SteadySolver& s,const Buffer<float>& b) {
  const auto x=d.download(b);std::vector<double> out(x.size());
  for(std::size_t i=0;i<x.size();++i) out[s.permutation().cells_new_to_old[i]]=x[i];
  return out;
}
void compare(const char* name,const std::vector<double>& x,const std::vector<double>& ref,
    const std::vector<bool>& wall) {
  for(int region=0;region<3;++region) {
    double n=0,d=0;std::size_t count=0;
    for(std::size_t i=0;i<x.size();++i) if(region==0 || wall[i]==(region==1)) {
      n+=(x[i]-ref[i])*(x[i]-ref[i]);d+=ref[i]*ref[i];++count;
    }
    std::cout<<"FIELD "<<name<<" region="<<(region==0?"all":region==1?"wall":"interior")
      <<" values="<<count<<" L2="<<std::sqrt(n/d)<<'\n';
  }
}
}
int main(int argc,char** argv) {
  try {
    KUTTA_REQUIRE((argc==8 || argc==9),"Usage: bench_elbow_turbulence case data out cpu|metal saved|cold iterations approximations [cell|jump]");
    const auto c=fluent::read_case(argv[1]);const auto mesh=read_mesh(argv[1]);
    const auto data=fluent::read_data(argv[2],c);const auto saved=setup::fluent_fields(c,data,mesh);
    auto problem=setup::fluent_setup(c,mesh,setup::parse_approximations(argv[7]));
    const bool ke=!saved.epsilon.empty();
    KUTTA_REQUIRE(problem.settings.sst || problem.settings.realizable_ke,"requires SST or k-epsilon");
    if(argc==9) {
      const std::string_view strain=argv[8];
      KUTTA_REQUIRE(ke && problem.settings.realizable_ke->model.variant==KeVariant::Standard,"strain override requires standard k-epsilon");
      KUTTA_REQUIRE(strain=="cell" || strain=="jump","invalid strain reconstruction");
      problem.settings.realizable_ke->standard_strain_reconstruction=strain=="jump"?KeStrainReconstruction::JumpCorrected:KeStrainReconstruction::CellGradient;
    }
    KUTTA_REQUIRE(std::string_view(argv[4])=="cpu" || std::string_view(argv[4])=="metal","invalid backend");
    KUTTA_REQUIRE(std::string_view(argv[5])=="saved" || std::string_view(argv[5])=="cold","invalid start");
    auto d=Device::create(std::string_view(argv[4])=="cpu"?Backend::Cpu:Backend::Metal);
    SteadySolver s(*d,mesh,problem.boundaries,problem.settings);
    setup::initialize(s,std::string_view(argv[5])=="saved"?setup::Initialization{saved}:setup::initial_fields(c,mesh));
    const auto iterations=std::stoi(argv[6]);KUTTA_REQUIRE(iterations>=0,"invalid iterations");
    std::filesystem::create_directories(argv[3]);const std::filesystem::path out=argv[3];
    std::cout<<std::setprecision(12)<<"CASE sha256="<<sha256_file(argv[1])<<" data="<<sha256_file(argv[2])
      <<" commit="<<KUTTA_BENCH_COMMIT<<" backend="<<argv[4]<<" start="<<argv[5]<<" iterations="<<iterations<<'\n';
    if(iterations) s.iterate(iterations);
    std::ofstream order(out/"order.csv");for(auto v:s.permutation().cells_new_to_old) order<<v<<'\n';order.close();
    const auto r=s.monitor();
    std::cout<<"RESIDUAL iteration="<<r.iteration<<" continuity_raw="<<r.continuity.raw
      <<" u="<<r.u.scaled<<" v="<<r.v.scaled<<" w="<<r.w.scaled<<" k="<<r.k.scaled
      <<" second="<<(ke?r.epsilon.scaled:r.omega.scaled)<<'\n';
    std::vector<bool> wall(mesh.n_cells,false);
    for(const auto& p:mesh.patches) if(p.kind==PatchKind::Wall)
      for(Index f=p.start;f<p.start+p.size;++f) wall[mesh.owner[f]]=true;
    std::vector<std::pair<std::string,std::vector<double>>> fields;
    auto add=[&](const char* name,const Buffer<float>& b){fields.emplace_back(name,original(*d,s,b));};
    add("u",s.velocity().x);add("v",s.velocity().y);add("w",s.velocity().z);add("p",s.pressure());
    add("k",s.kinetic_energy());add(ke?"epsilon":"omega",ke?s.dissipation():s.specific_dissipation());
    if(ke) {
      const auto& t=s.realizable_turbulence();
      add("mu_t",t.mu_t);add("strain",t.strain);add("pk",t.k_source);add("dk",t.k_sink);
      add("ps",t.epsilon_source);add("ds",t.epsilon_sink);add("gk",t.diffusivity_k);add("gs",t.diffusivity_epsilon);
      add("fixed",t.fixed_epsilon);
    } else {
      const auto& t=s.turbulence();
      add("mu_t",t.mu_t);add("strain",t.strain);add("pk",t.k_source);add("dk",t.k_sink);
      add("ps",t.omega_source);add("ds",t.omega_sink);add("gk",t.diffusivity_k);add("gs",t.diffusivity_omega);
      add("fixed",t.fixed_omega);add("f1",t.f1);add("f2",t.f2);
    }
    std::ofstream csv(out/"cells.csv");KUTTA_REQUIRE(csv,"cannot open cells.csv");csv<<std::setprecision(17);
    for(std::size_t j=0;j<fields.size();++j) csv<<(j?",":"")<<fields[j].first;
    csv<<'\n';for(Index i=0;i<mesh.n_cells;++i) {
      for(std::size_t j=0;j<fields.size();++j) csv<<(j?",":"")<<fields[j].second[i];csv<<'\n';
    }
    std::vector<double> mut(mesh.n_cells),vel(3*mesh.n_cells),refvel(vel.size());std::vector<bool> wallvec(vel.size());
    for(Index i=0;i<mesh.n_cells;++i) {
      mut[i]=data.cell_fields.at("SV_MU_T").at(std::uint64_t(i)+1)[0];
      for(int j=0;j<3;++j) {vel[3*i+j]=fields[j].second[i];refvel[3*i+j]=j==0?saved.velocity[i].x:j==1?saved.velocity[i].y:saved.velocity[i].z;wallvec[3*i+j]=wall[i];}
    }
    compare("velocity",vel,refvel,wallvec);compare("pressure",fields[3].second,saved.pressure,wall);
    compare("k",fields[4].second,saved.k,wall);compare(ke?"epsilon":"omega",fields[5].second,ke?saved.epsilon:saved.omega,wall);
    compare("mu_t",fields[6].second,mut,wall);
    const auto bp=d->download(s.boundary_pressure());const auto& m=s.mesh();const auto& g=s.geometry();
    double inlet=0,outlet=0,ia=0,oa=0,ir=0,orr=0;
    std::ofstream wc(out/"walls.csv");wc<<std::setprecision(17)<<"cell,dy,ystar,tau,fixed,production\n";
    const auto ys=d->download(ke?s.realizable_turbulence().wall_y_star:s.turbulence().wall_y_plus);
    const auto tau=d->download(ke?s.realizable_turbulence().wall_tau:s.turbulence().wall_tau);
    const auto dy=d->download(s.device_mesh().wall_normal_distance);

    for(const auto& p:m.patches) {
      const auto z=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& x){return x.zone_id==p.source_zone_id;});
      for(Index f=p.start;f<p.start+p.size;++f) {
        const auto b=f-m.n_internal_faces();const double a=norm(g.face_area[f]);
        const auto id=s.permutation().faces_new_to_old[f];
        const auto& op=*std::find_if(mesh.patches.begin(),mesh.patches.end(),[&](const auto& x){return x.name==p.name;});
        const auto ps=data.face_fields.at("SV_P").at(z->first+id-op.start)[0];
        if(p.kind==PatchKind::VelocityInlet) {inlet+=a*bp[b];ia+=a;ir+=a*ps;}
        if(p.kind==PatchKind::PressureOutlet) {outlet+=a*bp[b];oa+=a;orr+=a*ps;}
        if(p.kind==PatchKind::Wall) wc<<s.permutation().cells_new_to_old[m.owner[f]]<<','<<dy[b]
          <<','<<ys[b]<<','<<tau[b]<<','<<fields[14].second[s.permutation().cells_new_to_old[m.owner[f]]]<<','<<fields[8].second[s.permutation().cells_new_to_old[m.owner[f]]]<<'\n';
      }
    }
    std::cout<<"DROP ours="<<inlet/ia-outlet/oa<<" saved="<<ir/ia-orr/oa<<" relative="<<std::abs((inlet/ia-outlet/oa)/(ir/ia-orr/oa)-1)<<'\n';
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
