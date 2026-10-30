#include "../backend_helpers.hpp"
#include "../../validation/cases/ke_channel.hpp"
using namespace kutta;
using namespace kutta::validation;
TEST_CASE("V8 realizable KE scalable channel MKM DNS profile and momentum balance" * doctest::test_suite("slow")) {
  const auto dns=ke_channel::dns(std::filesystem::path(KUTTA_SOURCE_DIR)/"validation/data/mkm_channel/chan590.means");
  KUTTA_FOR_EACH_BACKEND(dev) {
    for(int n:{20,40,80}) {
      INFO("channel cells "<<n);const auto result=ke_channel::run(*dev,n);
      CHECK(result.converged);CHECK(std::abs(result.tau-1)<=.02);
      double error=0;
      for(const auto& point:result.profile) if(point[0]*ke_channel::re_tau>=30) {
        const double reference=tmr::interpolate(dns,point[0]);error=std::max(error,std::abs(point[1]/reference-1));
      }
      CHECK(error<=.1);
      std::cout<<std::setprecision(12)<<"V8 KE channel backend="<<backend_name(dev->backend())<<" cells="<<n<<" steps="<<result.steps<<" seconds="<<result.seconds<<" converged="<<result.converged<<" defect="<<result.defect<<" change="<<result.change<<" tau="<<result.tau<<" max_U_error="<<error<<std::endl;
    }
  }
}
