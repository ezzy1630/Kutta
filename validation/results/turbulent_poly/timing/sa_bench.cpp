// Scratch driver: one V7 SA TMR bump level, as tests/fv/test_tmr_sa.cpp runs it.
#include <iomanip>
#include <iostream>
#include "validation/cases/tmr_sa.hpp"
using namespace kutta;using namespace kutta::validation::tmr;
int main(int argc,char** argv) {
  auto dev=Device::create(std::string_view(argv[2])=="cpu"?Backend::Cpu:Backend::Metal);
  const auto r=sa::run(*dev,argv[1],Case::Bump,argv[3][0]-'0');
  std::cout<<std::setprecision(12)<<"SA_RESULT steps="<<r.iterations<<" seconds="<<r.solve_seconds<<" budget="<<r.stopped_for_budget;
  for(auto q:r.values.cf) std::cout<<" Cf="<<q;std::cout<<'\n';
}
