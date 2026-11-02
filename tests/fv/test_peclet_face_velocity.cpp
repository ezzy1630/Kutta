#include "../backend_helpers.hpp"
#include "kernels.hpp"
#include "kutta/fv/operators.hpp"

using namespace kutta;

TEST_CASE("Peclet side face velocity preserves the viscous limit and affine fields") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto ints=[&](std::initializer_list<Index> x) {return dev->upload<Index>(std::vector<Index>(x));};
    const auto floats=[&](std::initializer_list<float> x) {return dev->upload<float>(std::vector<float>(x));};
    auto offsets=ints({0,2,4}),faces=ints({0,1,0,2});
    auto signs=floats({1,1,-1,1}),flux=floats({3,-2,4});
    auto outflow=dev->alloc<float>(2);
    kern::fv_convective_outflow(*dev,2,offsets,faces,signs,flux,outflow);
    CHECK(dev->download(outflow)==std::vector<float>{3,4});

    auto owner=ints({0}),neighbour=ints({1});
    auto px=floats({.5f}),nx=floats({-.5f}),zero=floats({0});
    auto q=floats({1,3}),gx=floats({2,2}),gz=floats({0,0});
    auto weights=floats({2,6}),diagonal=floats({2,6});
    auto face=floats({2,99}); // the boundary sentinel must stay untouched
    const auto blend=[&] {
      kern::fv_peclet_side_velocity(*dev,1,owner,neighbour,px,zero,zero,nx,zero,zero,
        q,gx,gz,gz,weights,diagonal,outflow,face);
    };
    blend();
    CHECK(dev->download(face)==std::vector<float>{2,99}); // affine, unequal aP

    dev->upload<float>(gx,std::vector<float>{4,0}); // side=3, point=2
    dev->fill(outflow,0);blend();
    CHECK(dev->download(face)==std::vector<float>{2,99});
    dev->upload<float>(outflow,std::vector<float>{1,1});blend();
    CHECK(dev->download(face)==std::vector<float>{2.25f,99}); // beta=1/4
    dev->upload<float>(face,std::vector<float>{2,99});
    dev->upload<float>(outflow,std::vector<float>{2,6});blend();
    CHECK(dev->download(face)==std::vector<float>{3,99}); // pure advection
  }
}
