#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <numeric>

#include "../backend_helpers.hpp"
#include "../fluent_examples.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"

using namespace kutta;
using namespace kutta::detail;
using namespace kutta::fv;
namespace {
struct ErrorNorm {double l2=0,maximum=0,relative_maximum=0;};
ErrorNorm error(std::span<const double> values,std::span<const double> saved,int components=1) {
  REQUIRE(values.size()==saved.size());double square=0,denominator=0,maximum=0,reference=0;
  for(std::size_t i=0;i<values.size();i+=components) {
    double d=0,r=0;
    for(int j=0;j<components;++j) {d+=std::pow(values[i+j]-saved[i+j],2);r+=saved[i+j]*saved[i+j];}
    square+=d;denominator+=r;maximum=std::max(maximum,std::sqrt(d));reference=std::max(reference,std::sqrt(r));
  }
  REQUIRE(denominator>0);return {std::sqrt(square/denominator),maximum,maximum/reference};
}
std::vector<double> interleave(const std::vector<Vec3d>& values) {
  std::vector<double> out;for(auto v:values) {out.push_back(v.x);out.push_back(v.y);out.push_back(v.z);}return out;
}
std::vector<double> original(Device& dev,const SteadySolver& solver,const Buffer<float>& q) {
  const auto values=dev.download(q);std::vector<double> out(values.size());
  for(std::size_t c=0;c<values.size();++c) out[solver.permutation().cells_new_to_old[c]]=values[c];return out;
}
std::vector<double> original_velocity(Device& dev,const SteadySolver& solver) {
  const auto x=original(dev,solver,solver.velocity().x),y=original(dev,solver,solver.velocity().y),z=original(dev,solver,solver.velocity().z);
  std::vector<double> out;for(std::size_t c=0;c<x.size();++c) {out.push_back(x[c]);out.push_back(y[c]);out.push_back(z[c]);}return out;
}
void print_error(const char* field,ErrorNorm e) {
  std::cout<<"V6_FIELD "<<field<<" relative_L2="<<e.l2<<" max_absolute="<<e.maximum<<" relative_max="<<e.relative_maximum<<'\n';
}
void raw_flag(fluent::RawSettings& raw,const std::string& key,bool value) {
  for(auto& e:raw.alist.elements) if(!e.elements.empty() && e.elements[0].atom_text()==key) {
    e=fluent::parse_sexpr("("+key+(value?" #t)":" #f)"));return;
  }
  FAIL("missing fixture raw flag "<<key);
}
std::vector<std::uint64_t> face_ids(const fluent::FluentCase& c,const PolyMesh& mesh) {
  std::vector<std::uint64_t> ids(mesh.n_faces(),0);
  for(const auto& patch:mesh.patches) {
    const auto it=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& z){return z.zone_id==patch.source_zone_id;});
    REQUIRE(it!=c.mesh.face_zones.end());REQUIRE(it->size()==std::uint64_t(patch.size));
    for(Index j=0;j<patch.size;++j) ids[patch.start+j]=it->first+j;
  }
  return ids;
}
struct ZoneResult {double pressure=0,area=0;Vec3d pressure_force,shear_force;};
std::map<std::string,ZoneResult> saved_zones(const fluent::FluentCase& c,const fluent::FluentData& d,const PolyMesh& mesh) {
  const auto geom=compute_geometry(mesh);const auto ids=face_ids(c,mesh);std::map<std::string,ZoneResult> out;
  for(const auto& patch:mesh.patches) {
    auto& z=out[patch.name];
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto p=d.face_fields.at("SV_P").at(ids[f])[0],area=norm(geom.face_area[f]);
      z.pressure+=p*area;z.area+=area;z.pressure_force=z.pressure_force+geom.face_area[f]*p;
      if(patch.kind==PatchKind::Wall) {
        const auto shear=d.face_fields.at("SV_WALL_SHEAR").at(ids[f]);REQUIRE(shear.size()==3);
        z.shear_force=z.shear_force-Vec3d{shear[0],shear[1],shear[2]};
      }
    }
    z.pressure/=z.area;
  }
  return out;
}
}

TEST_CASE("Fluent setup converts SST inlet specifications and rejects unsupported physics") {
  const auto cold=setup::turbulence_values(IntensityLength{.05,.07*.1016},1,1.225,1.7894e-5);
  CHECK(cold.k==doctest::Approx(.00375));CHECK(cold.omega==doctest::Approx(15.72038792));
  const auto hot=setup::turbulence_values(IntensityLength{.05,.07*.0254},1.2,1.225,1.7894e-5);
  CHECK(hot.omega==doctest::Approx(75.45786201));
  const auto vr=setup::turbulence_values(IntensityViscosityRatio{.05,4},1,1.225,1.7894e-5);
  CHECK(vr.omega==doctest::Approx(1.225*.00375/(1.7894e-5*4)));
#if KUTTA_HAS_HDF5
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)) {MESSAGE("SKIP: external Fluent case absent");return;}
  auto c=fluent::read_case(path);const auto mesh=read_fluent_cff(path);
  CHECK_THROWS_WITH(setup::fluent_setup(c,mesh),doctest::Contains("recon/cell-lsf?"));
  const auto allowed=setup::parse_approximations("wfgc");
  const auto problem=setup::fluent_setup(c,mesh,allowed);
  CHECK(problem.approximations==allowed);
  CHECK(problem.settings.density==1.225);CHECK(problem.settings.viscosity==1.7894e-5);
  CHECK(std::holds_alternative<HybridInitialization>(setup::initial_fields(c,mesh)));
#endif
}

TEST_CASE("FV V6 Fluent mixing elbow fixed point and cold start" * doctest::test_suite("slow")) {
#if !KUTTA_HAS_HDF5
  MESSAGE("SKIP: V6 requires HDF5");
#else
  const auto root=fluent_examples()/"mixing_elbow";
  const auto cp=root/"mixing_elbow.cas.h5",dp=root/"mixing_elbow.dat.h5";
  if(!std::filesystem::exists(cp) || !std::filesystem::exists(dp)) {MESSAGE("SKIP: external Fluent mixing-elbow pair absent");return;}
  const auto c=fluent::read_case(cp);const auto mesh=read_fluent_cff(cp);const auto d=fluent::read_data(dp,c);
  const auto saved=setup::fluent_fields(c,d,mesh);const auto reference_velocity=interleave(saved.velocity);
  const auto zones=saved_zones(c,d,mesh);
  std::cout<<std::setprecision(10)<<"V6 Fluent="<<*c.solver_version<<" cells="<<mesh.n_cells<<" saved_iterations="<<(d.iteration_count?std::to_string(*d.iteration_count):"absent")
    <<" last_residual_iteration="<<(d.last_residual_iteration?std::to_string(*d.last_residual_iteration):"absent")<<" reports="<<c.reports.size()<<std::endl;
  for(const auto& [name,h]:d.residuals) if(!h.rows.empty()) {
    const auto row=h.rows.back();std::cout<<"V6_FLUENT_RESIDUAL "<<name<<" raw="<<row[0]<<" scale="<<row[1]<<" scaled="<<row[0]/row[1]<<'\n';
  }
  // Keep the strict original-case error visible AND failing. Diagnostic runs
  // below opt into one discrete approximation explicitly; they are
  // not a claim that the requested same-settings solve has been achieved.
  CHECK_NOTHROW(setup::fluent_setup(c,mesh));
  CHECK_NOTHROW(setup::initial_fields(c,mesh));
  auto diagnostic=c;
  const bool standard_proxy=std::getenv("KUTTA_V6_STANDARD_PROXY")!=nullptr;
  if(standard_proxy) raw_flag(diagnostic.variables,"hyb-init?",false);
  const auto allowed=setup::parse_approximations("wfgc");
  const auto setup=setup::fluent_setup(diagnostic,mesh,allowed);const auto cold=setup::initial_fields(diagnostic,mesh);
  std::cout<<"V6_DIAGNOSTIC "<<setup::approximation_notice(setup.approximations)
    <<"; cold start uses "<<(standard_proxy?"stored standard proxy":"stored hybrid initialization settings")<<std::endl;
  int limit=4000;
  if(const auto* value=std::getenv("KUTTA_V6_ITERATIONS")) {const auto parsed=parse_double(value);REQUIRE(parsed);REQUIRE(*parsed>0);REQUIRE(*parsed<=100000);REQUIRE(std::floor(*parsed)==*parsed);limit=int(*parsed);}
  int solved=0;
  KUTTA_FOR_EACH_BACKEND(dev) {
    if(const auto* selected=std::getenv("KUTTA_V6_BACKEND");selected && std::string_view(backend_name(dev->backend()))!=selected) continue;
    for(bool seeded:{true,false}) {
      ++solved;
      const auto mode=seeded?"seeded":(standard_proxy?"cold-standard-proxy":"cold-hybrid");
      auto settings=setup.settings;settings.monitor_interval=50;
      SteadySolver solver(*dev,mesh,setup.boundaries,settings);
      setup::initialize(solver,seeded?setup::Initialization{saved}:cold);
      const auto r0=solver.monitor();
      if(seeded) {
        // This fixture has one contiguous source cell zone, independently
        // checked by the data reader and setup mapping above.
        REQUIRE(c.mesh.cell_zones.size()==1);REQUIRE(c.mesh.cell_zones[0].first==1);
        std::vector<double> distance,mut;
        for(Index cell=0;cell<mesh.n_cells;++cell) {
          distance.push_back(d.cell_fields.at("SV_WALL_DIST").at(std::uint64_t(cell)+1)[0]);
          mut.push_back(d.cell_fields.at("SV_MU_T").at(std::uint64_t(cell)+1)[0]);
        }
        print_error("initial_wall_distance",error(original(*dev,solver,solver.device_mesh().wall_distance),distance));
        const auto mut_error=error(original(*dev,solver,solver.turbulence().mu_t),mut);
        print_error("initial_mu_t",mut_error);
        CHECK(mut_error.l2<.011); // round-1 derivative-gradient ratchet
        CHECK(r0.continuity.raw<8e-4); // round-1 flux ratchet; original gates below remain strict
        // Kutta's SST terms at Fluent's saved fields, for tools/forensics/k_budget.py.
        if(const auto* directory=std::getenv("KUTTA_V6_WALL_DIR")) {
          std::filesystem::create_directories(directory);
          const auto& t=solver.turbulence();
          const auto S=original(*dev,solver,t.strain),f1=original(*dev,solver,t.f1),f2=original(*dev,solver,t.f2);
          const auto mt=original(*dev,solver,t.mu_t),pk=original(*dev,solver,t.production_k);
          std::ofstream out(std::filesystem::path(directory)/(std::string(backend_name(dev->backend()))+"-saved-sst.csv"));
          REQUIRE(out);out<<std::setprecision(17)<<"strain,mu_t,f1,f2,production_k\n";
          for(std::size_t i=0;i<S.size();++i) out<<S[i]<<','<<mt[i]<<','<<f1[i]<<','<<f2[i]<<','<<pk[i]<<'\n';
        }
      }
      std::cout<<"V6_INITIAL backend="<<backend_name(dev->backend())<<" mode="<<mode<<" u="<<r0.u.scaled<<" v="<<r0.v.scaled
        <<" w="<<r0.w.scaled<<" pseudo_time_step="<<r0.pseudo_time_step<<" continuity_raw="<<r0.continuity.raw<<" k="<<r0.k.scaled<<" omega="<<r0.omega.scaled<<std::endl;
      if(seeded) {
        for(const auto& [name,value]:std::vector<std::pair<std::string,ResidualValue>>{
            {"x-velocity",r0.u},{"y-velocity",r0.v},{"z-velocity",r0.w},{"k",r0.k},{"omega",r0.omega}})
          std::cout<<"V6_SAVED_RESIDUAL backend="<<backend_name(dev->backend())<<" equation="<<name
            <<" raw="<<value.raw<<" steady_scaled="<<value.scaled<<" Fluent_scaled="<<value.fluent_scaled
            <<" normalization="<<value.fluent_normalization<<'\n';
      }
      // The absolute mass scale depends on mesh, density and throughput. Use
      // the saved run's final unscaled continuity as the comparison target,
      // alongside a tighter 1e-6 steady transport defect for every equation.
      const double continuity_target=d.residuals.at("continuity").rows.back()[0];
      const auto start=std::chrono::steady_clock::now();bool converged=false;
      for(int outer=0;outer<limit;outer+=50) {
        solver.iterate(std::min(50,limit-outer));const auto r=solver.monitor();
        const double maximum=std::max({r.u.scaled,r.v.scaled,r.w.scaled,r.k.scaled,r.omega.scaled});
        std::cout<<"V6_ITER backend="<<backend_name(dev->backend())<<" mode="<<mode<<" iteration="<<r.iteration
          <<" u="<<r.u.scaled<<" v="<<r.v.scaled<<" w="<<r.w.scaled<<" k="<<r.k.scaled<<" omega="<<r.omega.scaled
          <<" continuity="<<r.continuity.raw<<" corrected_mass="<<r.corrected_mass_imbalance<<" max_u="<<r.maximum_velocity<<std::endl;
        REQUIRE(std::isfinite(maximum));REQUIRE(r.maximum_velocity<100);
        if(maximum<1e-6 && r.continuity.raw<continuity_target) {converged=true;break;}
      }
      const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
      const auto last=solver.report();
      const double continuity_scale=last.continuity.raw/last.continuity.scaled;
      std::cout<<"V6_INITIAL_CONTINUITY backend="<<backend_name(dev->backend())<<" mode="<<mode
        <<" scaled="<<r0.continuity.raw/continuity_scale<<" first_five_scale="<<continuity_scale<<std::endl;
      std::cout<<"V6_RESULT backend="<<backend_name(dev->backend())<<" mode="<<mode<<" iterations="<<solver.report().iteration<<" seconds="<<seconds<<" converged="<<converged<<std::endl;
      const auto eu=error(original_velocity(*dev,solver),reference_velocity,3);
      const auto ep=error(original(*dev,solver,solver.pressure()),saved.pressure);
      const auto ek=error(original(*dev,solver,solver.kinetic_energy()),saved.k);
      const auto eo=error(original(*dev,solver,solver.specific_dissipation()),saved.omega);
      print_error("velocity",eu);print_error("pressure",ep);print_error("k",ek);print_error("omega",eo);
      CHECK(converged);CHECK(eu.l2<.005);
      if(seeded) CHECK(std::max({eu.l2,ep.l2,ek.l2,eo.l2})<.001);
      const auto bp=dev->download(solver.boundary_pressure());const auto& g=solver.geometry();
      double inlet=0,inlet_ref=0,inlet_area=0,outlet=0,outlet_ref=0,outlet_area=0;
      for(const auto& patch:solver.mesh().patches) {
        double p=0,area=0;
        for(Index f=patch.start;f<patch.start+patch.size;++f) {const auto a=norm(g.face_area[f]);p+=a*bp[f-solver.mesh().n_internal_faces()];area+=a;}
        std::cout<<"V6_ZONE_PRESSURE "<<patch.name<<" ours="<<p/area<<" Fluent="<<zones.at(patch.name).pressure<<'\n';
        if(patch.kind==PatchKind::VelocityInlet) {inlet+=p;inlet_ref+=zones.at(patch.name).pressure*area;inlet_area+=area;}
        if(patch.kind==PatchKind::PressureOutlet) {outlet+=p;outlet_ref+=zones.at(patch.name).pressure*area;outlet_area+=area;}
        if(patch.kind==PatchKind::Wall) {
          const std::vector<std::string> names{patch.name};const auto force=solver.forces(names);const auto& ref=zones.at(patch.name);
          for(const auto& pair:std::vector<std::pair<const char*,std::pair<Vec3d,Vec3d>>>{{"pressure",{force.pressure,ref.pressure_force}},
              {"shear",{force.viscous,ref.shear_force}},{"total",{force.force(),ref.pressure_force+ref.shear_force}}}) {
            const auto [a,b]=pair.second;
            std::cout<<"V6_FORCE "<<patch.name<<' '<<pair.first<<" ours="<<a.x<<','<<a.y<<','<<a.z<<" Fluent="<<b.x<<','<<b.y<<','<<b.z<<" relative="<<norm(a-b)/norm(b)<<'\n';
            if(std::string_view(pair.first)!="pressure") CHECK(norm(a-b)/norm(b)<.005);
          }
        }
      }
      const auto drop=inlet/inlet_area-outlet/outlet_area,reference=inlet_ref/inlet_area-outlet_ref/outlet_area;
      std::cout<<"V6_PRESSURE_DROP ours="<<drop<<" Fluent="<<reference<<" relative="<<std::abs(drop/reference-1)<<std::endl;
      CHECK(std::abs(drop/reference-1)<.005);
      // Optional converged wall and cell fields for tools/forensics/wall_shear.py,
      // which splits the wall-force gates face by face against Fluent's saved data.
      if(const auto* directory=std::getenv("KUTTA_V6_WALL_DIR")) {
        std::filesystem::create_directories(directory);
        const auto stem=std::filesystem::path(directory)/(std::string(backend_name(dev->backend()))+"-"+mode);
        const auto u=original(*dev,solver,solver.velocity().x),v=original(*dev,solver,solver.velocity().y),w=original(*dev,solver,solver.velocity().z);
        const auto p=original(*dev,solver,solver.pressure()),k=original(*dev,solver,solver.kinetic_energy());
        const auto omega=original(*dev,solver,solver.specific_dissipation());
        std::ofstream cells(stem.string()+"-cells.csv");REQUIRE(cells);cells<<std::setprecision(17)<<"u,v,w,p,k,omega\n";
        for(std::size_t i=0;i<p.size();++i) cells<<u[i]<<','<<v[i]<<','<<w[i]<<','<<p[i]<<','<<k[i]<<','<<omega[i]<<'\n';
        // Force on the wall per face, -shear*|A| (SteadySolver::wall_shear),
        // and the wall cell in Fluent's cell order.
        const auto sx=dev->download(solver.wall_shear().x),sy=dev->download(solver.wall_shear().y),sz=dev->download(solver.wall_shear().z);
        const auto& m=solver.mesh();const auto& order=solver.permutation().cells_new_to_old;
        std::ofstream walls(stem.string()+"-walls.csv");REQUIRE(walls);walls<<std::setprecision(17)<<"patch,x,y,z,fx,fy,fz,cell\n";
        for(const auto& patch:m.patches) if(patch.kind==PatchKind::Wall) for(Index f=patch.start;f<patch.start+patch.size;++f) {
          const auto b=f-m.n_internal_faces();const double a=norm(g.face_area[f]);const auto x=g.face_centroid[f];
          walls<<patch.name<<','<<x.x<<','<<x.y<<','<<x.z<<','<<-sx[b]*a<<','<<-sy[b]*a<<','<<-sz[b]*a<<','<<order[m.owner[f]]<<'\n';
        }
      }
      // Optional fixed-point comparison, after (and separate from) every
      // unchanged V6 gate above. Both cold starts reach the same 1,200 steps.
      if(const auto* directory=std::getenv("KUTTA_V6_COMPARE_DIR")) {
        solver.iterate(std::max(0,1200-solver.report().iteration));const auto r=solver.monitor();
        std::filesystem::create_directories(directory);
        auto path=std::filesystem::path(directory)/(std::string(backend_name(dev->backend()))+"-"+mode+".csv");
        std::ofstream out(path);REQUIRE(out);out<<"# "<<setup::approximation_notice(setup.approximations)<<'\n';
        out<<std::setprecision(17)<<"u,v,w,p,k,omega\n";
        const auto u=dev->download(solver.velocity().x),v=dev->download(solver.velocity().y),w=dev->download(solver.velocity().z);
        const auto p=dev->download(solver.pressure()),k=dev->download(solver.kinetic_energy()),omega=dev->download(solver.specific_dissipation());
        for(std::size_t c=0;c<p.size();++c) out<<u[c]<<','<<v[c]<<','<<w[c]<<','<<p[c]<<','<<k[c]<<','<<omega[c]<<'\n';
        std::cout<<"V6_PLATEAU backend="<<backend_name(dev->backend())<<" mode="<<mode<<" iteration="<<r.iteration
          <<" u="<<r.u.scaled<<" v="<<r.v.scaled<<" w="<<r.w.scaled<<" k="<<r.k.scaled<<" omega="<<r.omega.scaled
          <<" continuity="<<r.continuity.raw<<std::endl;
      }
    }
  }
  CHECK(solved>0);
  std::cout<<"V6_DIAGNOSTIC "<<setup::approximation_notice(setup.approximations)<<std::endl;
#endif
}
