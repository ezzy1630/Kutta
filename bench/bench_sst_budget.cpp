// Frozen-flow SST budget using the production model and scalar transport kernels.
// Diagnostic only: never changes SteadySolver or its defaults. Output is restricted
// reference-derived data and must stay outside version control.
#include "kutta/core/sha256.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"
#include "kernels.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
int count(std::string_view text) {
  std::size_t used=0;const int value=std::stoi(std::string(text),&used);
  KUTTA_REQUIRE(used==text.size() && value>=0,"invalid iteration count");return value;
}
std::vector<float> saved_flux(const std::filesystem::path& path,const PolyMesh& mesh,
    const MeshPermutation& permutation,const std::string& case_hash,const std::string& data_hash) {
  std::ifstream in(path);std::string line;
  KUTTA_REQUIRE(bool(std::getline(in,line)) && line=="# case_sha256="+case_hash,"flux case hash mismatch");
  KUTTA_REQUIRE(bool(std::getline(in,line)) && line=="# data_sha256="+data_hash,"flux data hash mismatch");
  KUTTA_REQUIRE(bool(std::getline(in,line)) && line=="owner,neighbour,flux,nodes","invalid flux columns");
  std::map<std::vector<Index>,Index> identities;
  for(Index f=0;f<mesh.n_faces();++f) {
    const auto ns=detail::nodes(mesh,f);std::vector<Index> key(ns.begin(),ns.end());std::sort(key.begin(),key.end());
    KUTTA_REQUIRE(identities.emplace(std::move(key),f).second,"duplicate mesh face identity");
  }
  std::vector<float> result(mesh.n_faces());std::vector<bool> seen(mesh.n_faces());
  while(std::getline(in,line)) {
    std::istringstream row(line);std::string a,b,q,ns;
    KUTTA_REQUIRE(bool(std::getline(row,a,',')) && bool(std::getline(row,b,',')) && bool(std::getline(row,q,',')) && bool(std::getline(row,ns)),"incomplete flux row");
    const auto parse_id=[](const std::string& value) {std::size_t used=0;const int id=std::stoi(value,&used);KUTTA_REQUIRE(used==value.size(),"invalid flux cell id");return id;};
    const int owner=parse_id(a),neighbour=parse_id(b);std::size_t used=0;const double flux=std::stod(q,&used);
    KUTTA_REQUIRE(used==q.size() && std::isfinite(flux) && std::abs(flux)<=std::numeric_limits<float>::max(),"invalid flux value");
    std::istringstream nodes(ns);std::vector<Index> key;Index node;
    while(nodes>>node) key.push_back(node);
    KUTTA_REQUIRE(nodes.eof() && !key.empty(),"invalid flux node ids");std::sort(key.begin(),key.end());
    const auto it=identities.find(key);KUTTA_REQUIRE(it!=identities.end(),"unknown flux face");const auto f=it->second;
    KUTTA_REQUIRE(!seen[f],"duplicate flux face");seen[f]=true;
    const auto p=permutation.cells_new_to_old[mesh.owner[f]];
    const auto n=f<mesh.n_internal_faces()?permutation.cells_new_to_old[mesh.neighbour[f]]:-1;
    KUTTA_REQUIRE((owner==p && neighbour==n) || (n>=0 && owner==n && neighbour==p),"flux adjacency mismatch");
    result[f]=float(owner==p?flux:-flux);
  }
  KUTTA_REQUIRE(in.eof() && std::all_of(seen.begin(),seen.end(),[](bool x){return x;}),"incomplete flux coverage");return result;
}
void vector_upload(Device& d,DeviceVec3& field,const std::vector<Vec3d>& values) {
  std::vector<float> x(values.size()),y(x.size()),z(x.size());
  for(std::size_t i=0;i<x.size();++i) {x[i]=float(values[i].x);y[i]=float(values[i].y);z[i]=float(values[i].z);}
  d.upload<float>(field.x,x);d.upload<float>(field.y,y);d.upload<float>(field.z,z);
}
}
int main(int argc,char** argv) {
 try {
  KUTTA_REQUIRE(argc==8 || argc==9,"Usage: bench_sst_budget case data output cpu|metal iterations kutta|saved flux.csv|- [compact-omega|no-cross|saved-wall-omega|molecular-omega|plain-strain|source-jump|saved-limited-strain]");
  const int iterations=count(argv[5]);const std::string backend=argv[4],flux_mode=argv[6],variant=argc==9?argv[8]:"default";
  KUTTA_REQUIRE(backend=="cpu" || backend=="metal","invalid backend");
  KUTTA_REQUIRE(flux_mode=="kutta" || flux_mode=="saved","invalid flux mode");
  KUTTA_REQUIRE(variant=="default" || variant=="compact-omega" || variant=="no-cross" || variant=="saved-wall-omega" || variant=="molecular-omega" || variant=="plain-strain" || variant=="source-jump" || variant=="saved-limited-strain","invalid variant");
  const std::filesystem::path out=argv[3];
  KUTTA_REQUIRE(!std::filesystem::exists(out/"identity.txt") && !std::filesystem::exists(out/"cells.csv"),"choose a fresh diagnostic output directory");
  const auto case_hash=sha256_file(argv[1]),data_hash=sha256_file(argv[2]);
  const auto c=fluent::read_case(argv[1]);const auto input=read_mesh(argv[1]);const auto data=fluent::read_data(argv[2],c);
  KUTTA_REQUIRE(c.mesh.cell_zones.size()==1 && c.mesh.cell_zones[0].first==1 && c.mesh.cell_zones[0].last==std::uint64_t(input.n_cells),"requires one dense source cell zone");
  const auto fields=setup::fluent_fields(c,data,input);const auto problem=setup::fluent_setup(c,input,setup::parse_approximations("wfgc"));
  KUTTA_REQUIRE(problem.settings.sst && !problem.settings.sst->transition,"requires plain SST");
  KUTTA_REQUIRE(problem.settings.gradient==GradientScheme::LeastSquares,"requires LSQ gradients");
  KUTTA_REQUIRE(!problem.settings.sst->compensated_transport && problem.settings.sst->steady_tolerance==0,"probe requires ordinary scalar transport");
  auto device=Device::create(backend=="cpu"?Backend::Cpu:Backend::Metal);auto& d=*device;
  SteadySolver reference(d,input,problem.boundaries,problem.settings);setup::initialize(reference,setup::Initialization{fields});
  const auto monitor=reference.monitor();const auto& mesh=reference.mesh();const auto& geometry=reference.geometry();const auto& m=reference.device_mesh();
  const auto& perm=reference.permutation();const auto nc=m.n_cells,nb=m.n_boundary_faces(),nf=m.n_faces;
  auto stencil=make_stencil(d,mesh,geometry);
  auto k=d.alloc<float>(nc),omega=d.alloc<float>(nc),flux=d.alloc<float>(nf),rho=d.alloc<float>(nc),mu=d.alloc<float>(nc);
  d.copy(k,reference.kinetic_energy());d.copy(omega,reference.specific_dissipation());d.copy(flux,reference.mass_flux());
  if(flux_mode=="saved") {const auto f=saved_flux(argv[7],mesh,perm,case_hash,data_hash);d.upload<float>(flux,f);}
  d.fill(rho,float(problem.settings.density));d.fill(mu,float(problem.settings.viscosity));
  std::vector<Index> kinds(nb),specs(nb,0);std::vector<float> a(nb),b(nb),saved_wall(nb);
  std::vector<Vec3d> boundary_velocity(nb),wall_velocity(nb);
  std::vector<BoundaryObservation> krows(nb,BoundaryObservation::Omitted),wrows(krows),srows(krows);
  for(const auto& patch:mesh.patches) {
    const auto it=std::find_if(problem.boundaries.begin(),problem.boundaries.end(),[&](const auto& p){return p.patch==patch.name;});
    KUTTA_REQUIRE(it!=problem.boundaries.end(),"missing patch condition");const auto& bc=it->condition;
    const auto original=std::find_if(input.patches.begin(),input.patches.end(),[&](const auto& p){return p.name==patch.name;});
    const auto zone=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& z){return z.zone_id==patch.source_zone_id;});
    KUTTA_REQUIRE(original!=input.patches.end() && zone!=c.mesh.face_zones.end(),"missing source zone");
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const Index bf=f-m.n_internal_faces,cell=mesh.owner[f];
      const Vec3d u=fields.velocity[perm.cells_new_to_old[cell]];
      const auto n=geometry.face_area[f]/norm(geometry.face_area[f]);
      const std::optional<TurbulenceCondition>* turbulence=nullptr;
      if(const auto* wall=std::get_if<Wall>(&bc)) {
        KUTTA_REQUIRE(!wall->rotation && norm(wall->velocity)==0,"requires stationary walls");kinds[bf]=0;
        wrows[bf]=BoundaryObservation::FaceValue;
        const auto id=zone->first+perm.faces_new_to_old[f]-original->start;
        saved_wall[bf]=float(data.face_fields.at("SV_O").at(id)[0]);
      } else if(const auto* inlet=std::get_if<VelocityInlet>(&bc)) {
        KUTTA_REQUIRE(std::holds_alternative<Vec3d>(inlet->velocity),"requires constant inlet velocity");
        kinds[bf]=1;boundary_velocity[bf]=std::get<Vec3d>(inlet->velocity);turbulence=&inlet->turbulence;
        krows[bf]=wrows[bf]=BoundaryObservation::FaceValue;
      } else if(const auto* outlet=std::get_if<PressureOutlet>(&bc)) {
        kinds[bf]=4;boundary_velocity[bf]=u;turbulence=&outlet->turbulence;
      } else {
        KUTTA_REQUIRE(std::holds_alternative<Symmetry>(bc),"unsupported diagnostic boundary");kinds[bf]=3;
        boundary_velocity[bf]=u-n*dot(u,n);krows[bf]=wrows[bf]=BoundaryObservation::NormalFootValue;
      }
      srows[bf]=kinds[bf]==0?BoundaryObservation::Omitted:BoundaryObservation::FaceValue;
      if(turbulence) {
        KUTTA_REQUIRE(turbulence->has_value(),"missing turbulence recipe");
        if(const auto* v=std::get_if<TurbulenceValues>(&**turbulence)) {a[bf]=float(v->k);b[bf]=float(v->omega);}
        else if(const auto* v=std::get_if<IntensityLength>(&**turbulence)) {specs[bf]=1;a[bf]=float(v->intensity);b[bf]=float(v->length);}
        else if(const auto* v=std::get_if<IntensityViscosityRatio>(&**turbulence)) {specs[bf]=2;a[bf]=float(v->intensity);b[bf]=float(v->ratio);}
        else KUTTA_REQUIRE(false,"unsupported turbulence recipe");
      }
    }
  }
  const auto hk=d.download(flux);
  for(Index bf=0;bf<nb;++bf) KUTTA_REQUIRE(kinds[bf]!=4 || hk[m.n_internal_faces+bf]>=0,"diagnostic does not implement outlet backflow velocity");
  auto kind=d.upload<Index>(kinds),spec=d.upload<Index>(specs);
  auto ta=d.upload<float>(a),tb=d.upload<float>(b);
  auto fixed_face=d.upload<float>(saved_wall),bk=d.alloc<float>(nb),bw=d.alloc<float>(nb);
  auto transport_kind=d.alloc<Index>(nb);
  auto ub=vector_field(d,nb),uw=vector_field(d,nb);vector_upload(d,ub,boundary_velocity);vector_upload(d,uw,wall_velocity);
  auto kfit=boundary_least_squares(d,mesh,geometry,krows,&stencil),wfit=boundary_least_squares(d,mesh,geometry,wrows,&stencil);
  auto sfit=boundary_least_squares(d,mesh,geometry,srows,&stencil);auto jump=face_jump_least_squares(d,mesh,geometry);
  DeviceVec3 base[]{vector_field(d,nc),vector_field(d,nc),vector_field(d,nc)};
  DeviceVec3 gs[]{vector_field(d,nc),vector_field(d,nc),vector_field(d,nc)};
  const auto& velocity=reference.velocity();vector_gradient(d,m,stencil,sfit,velocity,ub,base[0],base[1],base[2]);
  const Buffer<float>* v[]{&velocity.x,&velocity.y,&velocity.z};const Buffer<float>* vb[]{&ub.x,&ub.y,&ub.z};
  for(int j=0;j<3;++j) {
    if(variant=="plain-strain") {d.copy(gs[j].x,base[j].x);d.copy(gs[j].y,base[j].y);d.copy(gs[j].z,base[j].z);}
    else jump_corrected_gradient(d,m,stencil,jump,*v[j],*vb[j],base[j],gs[j]);
  }
  if(variant=="saved-limited-strain") {
    const auto sk=d.download(reference.turbulence().strain),f2=d.download(reference.turbulence().f2);
    std::vector<float> factor(nc,1);
    for(Index j=0;j<nc;++j) {
      const auto old=perm.cells_new_to_old[j];const double mt=data.cell_fields.at("SV_MU_T").at(std::uint64_t(old)+1)[0];
      if(mt*fields.omega[old]/(problem.settings.density*fields.k[old])<1-1e-9) {
        KUTTA_REQUIRE(mt>0 && f2[j]>0 && sk[j]>0,"cannot infer saved limited strain");
        factor[j]=float(problem.settings.density*.31*fields.k[old]/(mt*f2[j]*sk[j]));
      }
    }
    for(auto& g:gs) for(auto* q:{&g.x,&g.y,&g.z}) {
      auto values=d.download(*q);for(Index j=0;j<nc;++j) values[j]*=factor[j];d.upload<float>(*q,values);
    }
  }
  auto observed=wrows;for(Index bf=0;bf<nb;++bf) if(kinds[bf]==0) observed[bf]=BoundaryObservation::Omitted;
  auto reconstruction=variant=="compact-omega"?std::shared_ptr<const WalkFaceReconstruction>{}:walk_reconstruction({},d,mesh,geometry,m,stencil,observed);
  auto matrix=std::make_shared<CsrMatrix>(upload(d,stencil.pattern));
  TransportEquation keq(d,m,stencil,mesh,geometry,matrix),weq(d,m,stencil,mesh,geometry,matrix,reconstruction);
  auto options=problem.settings.sst->model;if(variant=="molecular-omega") options.omega_production_molecular_ratio=.001f;
  SstModel model(d,m,options);const auto& t=model.fields();
  auto source_k=vector_field(d,nc),source_w=vector_field(d,nc);
  auto no_cross_source=d.alloc<float>(nc),no_cross_sink=d.alloc<float>(nc),unused_source=d.alloc<float>(nc),unused_sink=d.alloc<float>(nc),unused_production=d.alloc<float>(nc),zero_cross=d.alloc<float>(nc);d.fill(zero_cross,0);
  const auto update=[&] {
    const auto evaluate=[&] {model.update({velocity,{gs[0],gs[1],gs[2]},rho,mu,k,omega,variant=="source-jump"?source_k:keq.grad,variant=="source-jump"?source_w:weq.grad,uw});};
    evaluate();
    kern::fv_turbulence_boundary(d,nb,m.owner,kind,spec,ta,tb,flux,ub.x,ub.y,ub.z,k,omega,
      variant=="saved-wall-omega"?fixed_face:t.wall_face_omega,bk,bw,transport_kind,m.n_internal_faces,float(problem.settings.density),float(problem.settings.viscosity));
    gradient(d,m,stencil,kfit,k,bk,keq.grad);gradient(d,m,stencil,wfit,omega,bw,weq.grad);
    if(variant=="source-jump") {
      jump_corrected_gradient(d,m,stencil,jump,k,bk,keq.grad,source_k);
      jump_corrected_gradient(d,m,stencil,jump,omega,bw,weq.grad,source_w);
    }
    evaluate();
    if(variant=="no-cross") {
      // Re-evaluate sources with exactly the same F1 and viscosity; remove only D_omega.
      kern::sst_sources(d,nc,rho,mu,k,omega,t.strain,t.vorticity,t.f1,zero_cross,t.mu_t,t.fixed_omega_mask,t.wall_production,
        unused_source,unused_sink,no_cross_source,no_cross_sink,unused_production,
        options.low_re,options.gamma_set==SstGammaSet::Menter2003,options.kato_launder?1:options.vorticity_production?2:0,
        options.production_limiter,options.limit_omega_production,options.wall_production,options.production_clip,
        options.limit_cross_diffusion,options.cross_diffusion_clip,options.k_floor,options.omega_floor,options.omega_production_molecular_ratio);
    }
  };
  auto settings=[&](bool w) {TransportSettings st;const auto& o=*problem.settings.sst;st.density=float(problem.settings.density);
    st.pseudo_time_step=o.pseudo_transient?float(monitor.pseudo_time_step):0;st.limit=problem.settings.limit;
    st.convection=w?o.omega_convection:o.k_convection;st.relaxation=float(w?o.omega_relaxation:o.k_relaxation);return st;};
  auto assemble=[&](bool w) {auto& eq=w?weq:keq;eq.assemble_reconstructed(w?omega:k,w?bw:bk,transport_kind,flux,w?t.face_diffusivity_omega:t.face_diffusivity_k,
    w?(variant=="no-cross"?no_cross_source:t.omega_source):t.k_source,settings(w),{w?(variant=="no-cross"?&no_cross_sink:&t.omega_sink):&t.k_sink});
    if(w) kern::fv_fixed_rows(d,nc,eq.matrix.row_ptr,eq.matrix.diag,t.fixed_omega_mask,t.fixed_omega,eq.matrix.val,eq.rhs);};
  std::filesystem::create_directories(out);
  std::cout<<std::setprecision(12)<<"IDENTITY commit="<<KUTTA_BENCH_COMMIT<<" case="<<case_hash<<" data="<<data_hash<<" backend="<<backend<<" flux="<<flux_mode<<" variant="<<variant<<" iterations="<<iterations<<'\n';
  update();
  if(variant=="default") {
    for(auto pair:{std::pair{&t.mu_t,&reference.turbulence().mu_t},std::pair{&t.strain,&reference.turbulence().strain},std::pair{&t.f1,&reference.turbulence().f1}}) {
      const auto x=d.download(*pair.first),y=d.download(*pair.second);double e=0,r=0;
      for(std::size_t j=0;j<x.size();++j) {e+=std::pow(x[j]-y[j],2);r+=double(y[j])*y[j];}
      std::cout<<"CLOSURE_CONTROL relative_L2="<<std::sqrt(e/r)<<'\n';KUTTA_REQUIRE(std::sqrt(e/r)<2e-6,"frozen probe differs from solver closure");
    }
  }
  auto residual=d.alloc<float>(nc),zero_cell=d.alloc<float>(nc),zero_face=d.alloc<float>(nf);d.fill(zero_cell,0);d.fill(zero_face,0);
  for(int it=0;it<iterations;++it) {
    update();assemble(false);keq.solve(k,problem.settings.sst->linear);kern::fv_floor(d,nc,k,options.k_floor);
    update();assemble(true);weq.solve(omega,problem.settings.sst->linear,&t.fixed_omega_mask);kern::fv_floor(d,nc,omega,options.omega_floor);
    if((it+1)%100==0) std::cout<<"ITERATION "<<it+1<<std::endl;
  }
  update();
  std::vector<std::pair<std::string,std::vector<float>>> columns;
  auto add=[&](const std::string& name,const Buffer<float>& q){columns.emplace_back(name,d.download(q));};
  add("k",k);add("omega",omega);add("mu_t",t.mu_t);add("strain",t.strain);add("vorticity",t.vorticity);add("f1",t.f1);add("f2",t.f2);add("cdkw",t.cdkw);
  add("pk",t.k_source);add("dk",t.k_sink);add("ps",variant=="no-cross"?no_cross_source:t.omega_source);add("ds",variant=="no-cross"?no_cross_sink:t.omega_sink);add("gk",t.diffusivity_k);add("gw",t.diffusivity_omega);add("fixed",t.fixed_omega);
  add("kx",keq.grad.x);add("ky",keq.grad.y);add("kz",keq.grad.z);add("wx",weq.grad.x);add("wy",weq.grad.y);add("wz",weq.grad.z);
  if(variant=="source-jump") {add("source_kx",source_k.x);add("source_ky",source_k.y);add("source_kz",source_k.z);add("source_wx",source_w.x);add("source_wy",source_w.y);add("source_wz",source_w.z);}
  for(bool w:{false,true}) {
    auto& eq=w?weq:keq;auto& q=w?omega:k;const auto& qb=w?bw:bk;const auto& diffusion=w?t.face_diffusivity_omega:t.face_diffusivity_k;
    assemble(w);eq.physical_residual(q,residual,w?&t.fixed_omega_mask:nullptr);add(w?"rw":"rk",residual);
    const auto rr=d.download(residual);const auto mask=d.download(t.fixed_omega_mask);double raw=0;
    for(Index j=0;j<nc;++j) if(!w || !mask[j]) raw+=std::abs(rr[j]);
    std::cout<<"NATIVE_RESIDUAL "<<(w?"omega":"k")<<" raw="<<raw;
    if(iterations==0 && variant=="default" && flux_mode=="kutta") {
      const double expected=w?monitor.omega.raw:monitor.k.raw;
      std::cout<<" solver_raw="<<expected<<" relative="<<std::abs(raw/expected-1);
      KUTTA_REQUIRE(std::abs(raw/expected-1)<2e-5,"budget assembly differs from solver monitor");
    }
    std::cout<<'\n';
    auto st=settings(w);st.relaxation=1;st.pseudo_time_step=0;
    eq.assemble_reconstructed(q,qb,transport_kind,flux,zero_face,zero_cell,st);eq.physical_residual(q,residual);add(w?"cw":"ck",residual);
    eq.assemble_reconstructed(q,qb,transport_kind,zero_face,diffusion,zero_cell,st);eq.physical_residual(q,residual);add(w?"dw":"dk_transport",residual);
  }
  std::ofstream csv(out/"cells.csv");KUTTA_REQUIRE(csv,"cannot write cells");csv<<std::setprecision(17)<<"cell";
  for(const auto& column:columns) csv<<','<<column.first;csv<<'\n';
  for(Index j=0;j<nc;++j) {csv<<perm.cells_new_to_old[j];for(const auto& column:columns) csv<<','<<column.second[j];csv<<'\n';}
  KUTTA_REQUIRE(csv,"failed writing cells");
  const auto kh=d.download(k),wh=d.download(omega);
  for(int region=0;region<3;++region) for(bool w:{false,true}) {
    const auto wall=d.download(t.fixed_omega_mask);double e=0,r=0;std::size_t n=0;
    for(Index j=0;j<nc;++j) if(region==0 || bool(wall[j])==(region==1)) {const auto old=perm.cells_new_to_old[j];const double ref=w?fields.omega[old]:fields.k[old];const double value=w?wh[j]:kh[j];e+=(value-ref)*(value-ref);r+=ref*ref;++n;}
    std::cout<<"FIELD "<<(w?"omega":"k")<<" region="<<region<<" cells="<<n<<" L2="<<std::sqrt(e/r)<<'\n';
  }
  std::ofstream metadata(out/"identity.txt");KUTTA_REQUIRE(metadata,"cannot write identity");metadata<<case_hash<<'\n'<<data_hash<<'\n'<<KUTTA_BENCH_COMMIT<<'\n'<<backend<<'\n'<<flux_mode<<'\n'<<variant<<'\n'<<iterations<<'\n';KUTTA_REQUIRE(metadata,"failed writing identity");
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
