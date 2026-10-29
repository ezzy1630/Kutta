#include <iostream>
#include "../backend_helpers.hpp"
#include "pressure_operator.hpp"

using namespace kutta;
TEST_CASE("FV MMS pressure gradient operator convergence") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    test::PressureOperatorError old{};
    for(int n:{4,8,16}) {
      const auto e=test::pressure_operator_error(*dev,n,tet);
      std::cout<<"FV_PRESSURE_OPERATOR backend="<<backend_name(dev->backend())<<" tet="<<tet<<" n="<<n
        <<" exact_face_L2="<<e.exact_face<<" momentum_L2="<<e.momentum<<" LSQ_L2="<<e.lsq;
      if(old.exact_face>0) {
        const double exact_order=test::order(old.exact_face,e.exact_face);
        const double momentum_order=test::order(old.momentum,e.momentum);
        std::cout<<" exact_face_order="<<exact_order<<" momentum_order="<<momentum_order
          <<" LSQ_order="<<test::order(old.lsq,e.lsq);
        if(n==16) {
          CHECK(test::order(old.lsq,e.lsq)>=.9);
          CHECK(exact_order>=(tet?.9:1.9));
          if(tet) {CHECK(exact_order<=1.1);CHECK(momentum_order>=.9);}
        }
      }
      std::cout<<std::endl;old=e;
    }
  }
}
