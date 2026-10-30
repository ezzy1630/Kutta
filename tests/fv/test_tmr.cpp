#include "../backend_helpers.hpp"
#include "../../validation/cases/tmr.hpp"

using namespace kutta;
using namespace kutta::validation::tmr;
namespace {
const std::filesystem::path root=std::filesystem::path(KUTTA_SOURCE_DIR)/"validation/data/nasa_tmr";
double finest(const std::vector<Zone>& data,std::size_t zone=0) {
  REQUIRE(zone<data.size());REQUIRE_FALSE(data[zone].rows.empty());REQUIRE(data[zone].rows.front().size()==4);
  return data[zone].rows.front()[3];
}
double relative(double value,double reference) {return std::abs(value/reference-1);}
}
TEST_CASE("TMR reference reader preserves independent code zones and interpolation") {
  const auto cf=read_reference(root/"2DZP_flatplate/sst/cf_convergence_sstv.dat");
  REQUIRE(cf.size()==2);CHECK(cf[0].name.find("CFL3D")!=std::string::npos);
  CHECK(finest(cf)==doctest::Approx(.00269085355));CHECK(cf[1].rows.size()==5);
  CHECK(interpolate({{0,2},{1,4},{2,3}},.5)==3);
  CHECK(interpolate({{0,2},{1,4},{2,3}},2)==3);
  CHECK_THROWS(interpolate({{0,2},{1,4}},-1));
  const auto cp=read_reference(root/"2DBC_bump/sst2003/cp_bump_sst2003.dat");
  REQUIRE(cp.size()==2);CHECK(std::isfinite(reference_station(cp.front(),.75)));
}
TEST_CASE("TMR uniform bump startup stays bounded on skew wall cells") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    Options options;options.max_iterations=24;
    const auto result=run(*dev,root,Case::Bump,0,options);
    CHECK_MESSAGE(result.failure.empty(),result.failure);
    CHECK(result.iterations==24);
    CHECK(std::isfinite(result.residual.maximum_velocity));
    CHECK(result.residual.maximum_velocity<2);
    CHECK(std::isfinite(result.values.drag));
  }
}
TEST_CASE("V5 TMR flat plate SST accuracy and three-grid trend" * doctest::test_suite("slow")) {
  const double ref_cf=finest(read_reference(root/"2DZP_flatplate/sst/cf_convergence_sstv.dat"));
  const double ref_cd=finest(read_reference(root/"2DZP_flatplate/sst/drag_convergence_sstv.dat"));
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::array<double,3> error{};
    for(int level=0;level<3;++level) {
      INFO("level: "<<level);const auto r=run(*dev,root,Case::FlatPlate,level);
      CHECK_MESSAGE(r.failure.empty(),r.failure);CHECK(r.converged);CHECK(r.iterations<=3000);CHECK(std::isfinite(r.values.drag));
      if(!r.failure.empty()) {error[level]=std::numeric_limits<double>::infinity();if(level<2) continue;}
      error[level]=std::max(relative(r.values.cf[0],ref_cf),relative(r.values.drag,ref_cd));
      if(level==2) {CHECK(relative(r.values.cf[0],ref_cf)<=.05);CHECK(relative(r.values.drag,ref_cd)<=.07);}
      std::cout<<"V5 plate backend="<<backend_name(dev->backend())<<" cells="<<r.cells<<" iterations="<<r.iterations<<" seconds="<<r.solve_seconds<<" setup_s="<<r.setup_seconds<<" converged="<<r.converged<<" budget="<<r.stopped_for_budget<<" change="<<r.load_change<<" error="<<error[level]<<" cf="<<r.values.cf[0]<<" drag="<<r.values.drag<<" continuity="<<r.residual.continuity.scaled<<" continuity_raw="<<r.residual.continuity.raw<<" corrected_mass="<<r.residual.corrected_mass_imbalance<<" ru="<<r.residual.u.scaled<<" rv="<<r.residual.v.scaled<<" rk="<<r.residual.k.scaled<<" rw="<<r.residual.omega.scaled<<std::endl;
    }
    CHECK((std::all_of(error.begin(),error.end(),[](double value){return std::isfinite(value);}) && error[2]<error[0]));
  }
}
TEST_CASE("V5 TMR bump SST station accuracy and three-grid trend" * doctest::test_suite("slow")) {
  const auto cf=read_reference(root/"2DBC_bump/sst2003/cf_convergence_sst2003.dat");
  const auto cp=read_reference(root/"2DBC_bump/sst2003/cp_bump_sst2003.dat");
  // Convergence zones are peak, front, aft, each CFL3D followed by FUN3D.
  const std::array<double,3> ref_cf{finest(cf,2),finest(cf,0),finest(cf,4)};
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::array<double,3> error{};
    for(int level=0;level<3;++level) {
      INFO("level: "<<level);const auto r=run(*dev,root,Case::Bump,level);
      std::cout<<"V5 bump measurement backend="<<backend_name(dev->backend())<<" cells="<<r.cells<<" iterations="<<r.iterations<<" seconds="<<r.solve_seconds<<" setup_s="<<r.setup_seconds<<" converged="<<r.converged<<" budget="<<r.stopped_for_budget<<" change="<<r.load_change<<" drag="<<r.values.drag;
      for(auto value:r.values.cf) std::cout<<" cf="<<value;
      for(auto value:r.values.cp) std::cout<<" cp="<<value;
      std::cout<<" continuity="<<r.residual.continuity.scaled<<" failure="<<r.failure<<std::endl;
      CHECK_MESSAGE(r.failure.empty(),r.failure);CHECK(r.converged);CHECK(r.iterations<=3000);CHECK(std::isfinite(r.values.drag));
      if(!r.failure.empty()) {error[level]=std::numeric_limits<double>::infinity();if(level<2) continue;}
      for(int station=0;station<3;++station) {
        const double cf_error=relative(r.values.cf[station],ref_cf[station]);
        const double cp_error=std::abs(r.values.cp[station]-reference_station(cp.front(),bump_stations[station]));
        error[level]=std::max({error[level],cf_error,cp_error/.03});
        if(level==2) {INFO("station: "<<bump_stations[station]);CHECK(cf_error<=.1);CHECK(cp_error<=.03);}
      }
      std::cout<<"V5 bump backend="<<backend_name(dev->backend())<<" cells="<<r.cells<<" iterations="<<r.iterations<<" seconds="<<r.solve_seconds<<" setup_s="<<r.setup_seconds<<" converged="<<r.converged<<" budget="<<r.stopped_for_budget<<" change="<<r.load_change<<" error="<<error[level]<<std::endl;
    }
    CHECK((std::all_of(error.begin(),error.end(),[](double value){return std::isfinite(value);}) && error[2]<error[0]));
  }
}
