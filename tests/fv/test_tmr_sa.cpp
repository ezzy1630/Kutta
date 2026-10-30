#include "../backend_helpers.hpp"
#include "../../validation/cases/tmr_sa.hpp"
using namespace kutta;
using namespace kutta::validation::tmr;
namespace {
const std::filesystem::path root=std::filesystem::path(KUTTA_SOURCE_DIR)/"validation/data/nasa_tmr";
double reference(const std::vector<Zone>& zones,int zone=0,int column=3) {return zones.at(zone).rows.at(0).at(column);}
double relative(double a,double b) {return std::abs(a/b-1);}
}
TEST_CASE("SA TMR pinned reference values and variant identity") {
  CHECK(reference(read_reference(root/"2DZP_flatplate/sa/cf_convergence.dat"))==doctest::Approx(.00270562153));
  CHECK(reference(read_reference(root/"2DBC_bump/sa/force_convergence.dat"),0,4)==doctest::Approx(.003572382));
}
TEST_CASE("V7 SA TMR flat plate and bump accuracy convergence and three-grid trend" * doctest::test_suite("slow")) {
  for(auto c:{Case::FlatPlate,Case::Bump}) {
    const bool plate=c==Case::FlatPlate;
    const auto cf=read_reference(root/(plate?"2DZP_flatplate/sa/cf_convergence.dat":"2DBC_bump/sa/cf_convergence_sa.dat"));
    const auto drag=read_reference(root/(plate?"2DZP_flatplate/sa/drag_convergence.dat":"2DBC_bump/sa/force_convergence.dat"));
    const auto cp=plate?std::vector<Zone>{}:read_reference(root/"2DBC_bump/sa/cp_bump.dat");
    KUTTA_FOR_EACH_BACKEND(dev) {
      std::array<double,3> errors{};
      for(int level=0;level<3;++level) {
        INFO("case "<<(plate?"plate":"bump")<<" level "<<level);
        const auto r=sa::run(*dev,root,c,level);const auto& v=r.values;
        CHECK_MESSAGE(r.failure.empty(),r.failure);CHECK(r.converged);CHECK(r.iterations<=3000);CHECK(std::isfinite(v.drag));
        const double cd_error=relative(v.drag,reference(drag,0,plate?3:4));errors[level]=cd_error;
        if(level==2) CHECK(cd_error<=(plate?.07:.1));
        for(int j=0;j<(plate?1:3);++j) {
          const double cf_error=relative(v.cf[j],reference(cf,plate?0:std::array<int,3>{2,0,4}[j]));
          errors[level]=std::max(errors[level],cf_error);
          if(level==2) CHECK(cf_error<=(plate?.05:.1));
          if(!plate) {
            const double cp_error=std::abs(v.cp[j]-reference_station(cp.front(),bump_stations[j]));
            errors[level]=std::max(errors[level],cp_error/.03);if(level==2) CHECK(cp_error<=.03);
          }
        }
        if(!r.failure.empty()) errors[level]=std::numeric_limits<double>::infinity();
        std::cout<<std::setprecision(12)<<"V7 SA "<<(plate?"plate":"bump")<<" backend="<<backend_name(dev->backend())<<" level="<<level<<" cells="<<r.cells<<" steps="<<r.iterations<<" seconds="<<r.solve_seconds<<" converged="<<r.converged<<" budget="<<r.stopped_for_budget<<" load_change="<<r.load_change<<" error="<<errors[level]<<" Cd="<<v.drag;
        for(auto x:v.cf) std::cout<<" Cf="<<x;for(auto x:v.cp) std::cout<<" Cp="<<x;
        std::cout<<" continuity="<<r.residual.continuity.scaled<<" u="<<r.residual.u.scaled<<" v="<<r.residual.v.scaled<<" nu_tilde="<<r.residual.nu_tilde.scaled<<" max_yplus="<<v.max_yplus<<" failure="<<r.failure<<std::endl;
      }
      CHECK((std::all_of(errors.begin(),errors.end(),[](double e){return std::isfinite(e);}) && errors[2]<errors[0]));
    }
  }
}
