// Fluent corpus pairs beyond V6 (validation/results/case_coverage/REPORT.md).
// Every pair skips with a message when KUTTA_FLUENT_EXAMPLES lacks it.
#include <cmath>
#include <iomanip>
#include <iostream>

#include "../backend_helpers.hpp"
#include "../fluent_examples.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"

using namespace kutta;
using namespace kutta::detail;
using namespace kutta::fv;
using setup::Approximation;
namespace {
struct Pair {
  const char* name;std::filesystem::path cas,dat;
  std::vector<Approximation> approximations;
};
const std::filesystem::path elbow_steady="other_folders/result_files/fluent-mixing_elbow_steady-state";
const std::filesystem::path nozzle_run="flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP";
Pair pptx() {return {"PPTX elbow","examples/PPTX-report-generation/elbow.cas.h5","examples/PPTX-report-generation/elbow.dat.h5",{Approximation::wfgc}};}
Pair dash() {return {"Dash elbow","examples/Dash-App-for-Post-processing/elbow.cas.h5","examples/Dash-App-for-Post-processing/elbow.dat.h5",{Approximation::wfgc}};}
Pair nozzle() {return {"nozzle",nozzle_run/"nozzle_para_named_1.cas.h5",nozzle_run/"nozzle_para_named_1.dat.h5",{Approximation::wfgc}};}
Pair elbow_param() {return {"elbow_param","mixing_elbow/elbow_param.cas.h5","mixing_elbow/elbow_param.dat.h5",{Approximation::compressibility,Approximation::ke_options}};}
Pair elbow_two() {return {"elbow-2",elbow_steady/"elbow-2.cas.h5",elbow_steady/"elbow-2-00050.dat.h5",{Approximation::compressibility,Approximation::ke_options}};}

bool present(const Pair& p) {
  const bool ok=std::filesystem::exists(fluent_examples()/p.cas) && std::filesystem::exists(fluent_examples()/p.dat);
  if(!ok) MESSAGE("SKIP: external Fluent pair absent: "<<std::string(p.name)<<" (set KUTTA_FLUENT_EXAMPLES)");
  return ok;
}
double l2(std::span<const double> a,std::span<const double> b) {
  REQUIRE(a.size()==b.size());double n=0,d=0;
  for(std::size_t i=0;i<a.size();++i) {n+=(a[i]-b[i])*(a[i]-b[i]);d+=b[i]*b[i];}
  REQUIRE(d>0);return std::sqrt(n/d);
}
std::vector<double> interleave(const std::vector<Vec3d>& values) {
  std::vector<double> out;for(auto v:values) {out.push_back(v.x);out.push_back(v.y);out.push_back(v.z);}return out;
}
std::vector<double> original(Device& dev,const SteadySolver& solver,const Buffer<float>& q) {
  const auto values=dev.download(q);std::vector<double> out(values.size());
  for(std::size_t c=0;c<values.size();++c) out[solver.permutation().cells_new_to_old[c]]=values[c];
  return out;
}
struct Comparison {double velocity=0,pressure=0,k=0,second=0,pressure_drop=0,wall_force=0;};
// Cold start from the case's own initialization, then the V6 conventions:
// unweighted relative L2 per field, area-weighted inlet minus outlet pressure,
// and total wall force against Fluent's saved SV_P and SV_WALL_SHEAR.
Comparison solve_and_compare(Device& dev,const Pair& pair,int iterations) {
  const auto cp=fluent_examples()/pair.cas,dp=fluent_examples()/pair.dat;
  const auto c=fluent::read_case(cp);const auto mesh=read_mesh(cp);const auto d=fluent::read_data(dp,c);
  const auto saved=setup::fluent_fields(c,d,mesh);const auto problem=setup::fluent_setup(c,mesh,pair.approximations);
  SteadySolver solver(dev,mesh,problem.boundaries,problem.settings);
  setup::initialize(solver,setup::initial_fields(c,mesh));
  solver.iterate(iterations);
  const auto report=solver.monitor();
  INFO(pair.name<<" continuity="<<report.continuity.fluent_scaled);
  const bool ke=!saved.epsilon.empty();
  setup::CellFields ours;
  const auto x=original(dev,solver,solver.velocity().x),y=original(dev,solver,solver.velocity().y),z=original(dev,solver,solver.velocity().z);
  for(std::size_t i=0;i<x.size();++i) ours.velocity.push_back({x[i],y[i],z[i]});
  ours.pressure=original(dev,solver,solver.pressure());ours.k=original(dev,solver,solver.kinetic_energy());
  (ke?ours.epsilon:ours.omega)=original(dev,solver,ke?solver.dissipation():solver.specific_dissipation());
  // Initialization at the final fields rebuilds the boundary reconstruction.
  if(ke) solver.initialize_ke_fields(ours.velocity,ours.pressure,ours.k,ours.epsilon);
  else solver.initialize_fields(ours.velocity,ours.pressure,ours.k,ours.omega);
  Comparison out{l2(interleave(ours.velocity),interleave(saved.velocity)),l2(ours.pressure,saved.pressure),l2(ours.k,saved.k),
    ke?l2(ours.epsilon,saved.epsilon):l2(ours.omega,saved.omega)};
  // Fluent face IDs follow the original mesh; the solver's patches are renumbered.
  const auto reference=c.reference.pressure.value_or(0);const auto original_geometry=compute_geometry(mesh);
  double inlet_saved=0,inlet_area=0,outlet_saved=0,outlet_area=0;Vec3d force_saved;
  for(const auto& patch:mesh.patches) {
    const auto zone=std::find_if(c.mesh.face_zones.begin(),c.mesh.face_zones.end(),[&](const auto& z){return z.zone_id==patch.source_zone_id;});
    REQUIRE(zone!=c.mesh.face_zones.end());REQUIRE(zone->size()==std::uint64_t(patch.size));
    for(Index j=0;j<patch.size;++j) {
      const auto f=patch.start+j;const auto id=zone->first+j;const double a=norm(original_geometry.face_area[f]),ps=d.face_fields.at("SV_P").at(id)[0];
      if(patch.kind==PatchKind::VelocityInlet) {inlet_saved+=a*ps;inlet_area+=a;}
      if(patch.kind==PatchKind::PressureOutlet) {outlet_saved+=a*ps;outlet_area+=a;}
      if(patch.kind==PatchKind::Wall) {
        const auto shear=d.face_fields.at("SV_WALL_SHEAR").at(id);REQUIRE(shear.size()==3);
        force_saved=force_saved+original_geometry.face_area[f]*(ps-reference)-Vec3d{shear[0],shear[1],shear[2]};
      }
    }
  }
  const auto bp=dev.download(solver.boundary_pressure());const auto& g=solver.geometry();const auto& m=solver.mesh();
  double inlet=0,outlet=0;Vec3d force;
  for(const auto& patch:m.patches) {
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const double a=norm(g.face_area[f])*bp[f-m.n_internal_faces()];
      if(patch.kind==PatchKind::VelocityInlet) inlet+=a;
      if(patch.kind==PatchKind::PressureOutlet) outlet+=a;
    }
    if(patch.kind==PatchKind::Wall) {const std::vector<std::string> names{patch.name};force=force+solver.forces(names,reference).force();}
  }
  REQUIRE(inlet_area>0);REQUIRE(outlet_area>0);
  const double drop=inlet/inlet_area-outlet/outlet_area,drop_saved=inlet_saved/inlet_area-outlet_saved/outlet_area;
  out.pressure_drop=std::abs(drop/drop_saved-1);out.wall_force=norm(force-force_saved)/norm(force_saved);
  std::cout<<std::setprecision(10)<<"FLUENT_COVERAGE "<<pair.name<<" backend="<<backend_name(dev.backend())<<" iterations="<<iterations
    <<" continuity="<<report.continuity.fluent_scaled<<" velocity_L2="<<out.velocity<<" pressure_L2="<<out.pressure<<" k_L2="<<out.k
    <<(ke?" epsilon_L2=":" omega_L2=")<<out.second<<" pressure_drop="<<out.pressure_drop<<" wall_force="<<out.wall_force<<std::endl;
  return out;
}
}

TEST_CASE("Fluent corpus pairs run only with their named approximations and map their numerics") {
  for(const auto& pair:{pptx(),dash(),nozzle(),elbow_param(),elbow_two()}) {
    INFO(pair.name);
    if(!present(pair)) continue;
    const auto cp=fluent_examples()/pair.cas;
    const auto c=fluent::read_case(cp);const auto mesh=read_mesh(cp);
    CHECK_THROWS_AS(setup::fluent_setup(c,mesh),setup::RejectedSettings);
    const auto problem=setup::fluent_setup(c,mesh,pair.approximations);
    // Applicable approximations are recorded, never added silently.
    CHECK(problem.approximations==pair.approximations);
    // Reports Kutta does not evaluate are announced, not silently dropped.
    const auto unevaluated=std::count_if(c.reports.begin(),c.reports.end(),[](const auto& r){
      return r.kind==fluent::ReportDefinition::Kind::other || r.kind==fluent::ReportDefinition::Kind::moment;});
    CHECK(std::count_if(problem.convergence_notes.begin(),problem.convergence_notes.end(),[](const auto& n){return n.starts_with("report ");})==unevaluated);
    const auto& s=problem.settings;
    if(c.viscous.model==fluent::ViscousModel::ke_standard) {
      REQUIRE(s.realizable_ke);CHECK(s.realizable_ke->model.variant==KeVariant::Standard);
      CHECK(s.realizable_ke->model.wall==KeWallTreatment::Standard);
      CHECK(s.pressure_interpolation==PressureInterpolation::Standard);
      CHECK(s.pressure_correction==PressureCorrection::Simple);
      CHECK(!s.pseudo_transient);CHECK(s.convection==ConvectionScheme::FirstOrderUpwind);
    } else {
      REQUIRE(s.sst);
      CHECK(s.face_velocity==FaceVelocityWeights::Distance);
      CHECK(s.pseudo_transient==(std::string(pair.name)=="nozzle"));
    }
  }
}

// Fluent saved elbow_param after 89 iterations (scaled continuity 3.9e-4) and
// the nozzle after 10. The other pairs' saved fields are unconverged or belong
// to another material, so they are reported, not gated. Bounds are about 1.5x
// the measured CPU cold-start values in the report.
TEST_CASE("Fluent standard k-epsilon elbow and SST nozzle match converged Fluent fields" * doctest::test_suite("slow")) {
  KUTTA_FOR_EACH_BACKEND(dev) {
    if(const auto pair=elbow_param();present(pair)) {
      const auto r=solve_and_compare(*dev,pair,300);
      CHECK(r.velocity<.015);CHECK(r.pressure<.1);CHECK(r.k<.14);CHECK(r.second<.27);
      CHECK(r.pressure_drop<.015);CHECK(r.wall_force<.006);
    }
    if(const auto pair=nozzle();present(pair)) {
      const auto r=solve_and_compare(*dev,pair,300);
      CHECK(r.velocity<.025);CHECK(r.pressure<.009);CHECK(r.k<.08);CHECK(r.second<.085);
      CHECK(r.pressure_drop<.008);CHECK(r.wall_force<.01);
    }
  }
}
