#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include "../backend_helpers.hpp"
#include "physics_meshes.hpp"
#include "kutta/fv/face_projection.hpp"
#include "kutta/mesh/topology.hpp"
#include "kernels.hpp"

using namespace kutta;
using namespace kutta::fv;

TEST_CASE("FV paired projected predictor bounds the total field instead of its increment") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto offsets=dev->upload<Index>(std::vector<Index>{0,3}),points=dev->upload<Index>(std::vector<Index>{0,1,2});
    auto owner=dev->upload<Index>(std::vector<Index>{0}),neighbor=dev->upload<Index>(std::vector<Index>{1}),mask=dev->upload<Index>(std::vector<Index>{1});
    auto weights=dev->upload<float>(std::vector<float>{0,0,4}),link=dev->upload<float>(std::vector<float>{0,0,0});
    auto diagonal=dev->alloc<float>(3),limiter=dev->alloc<float>(3),out=dev->alloc<float>(1);
    const float epsilon=std::ldexp(1.f,-20);
    auto high=dev->upload<float>(std::vector<float>{1,2,1.5f}),low=dev->upload<float>(std::vector<float>{0,0,epsilon});
    dev->fill(diagonal,1);
    for(float theta:{1.f,.25f,0.f}) {
      dev->fill(limiter,theta);dev->fill(out,0);
      kern::fv_projected_momentum_increment(*dev,1,offsets,points,owner,neighbor,mask,weights,link,
        high,low,high,low,diagonal,limiter,out);
      // The high reconstruction is 1.5, strictly inside [1,2]. The total
      // remains inside, so its increment is 4*epsilon, although that exceeds
      // the increment field's own range [0,epsilon].
      CHECK(dev->download(out)[0]==theta*4*epsilon);
    }
    dev->upload<float>(high,std::vector<float>{1,2,2.5f});dev->fill(limiter,1);dev->fill(out,0);
    kern::fv_projected_momentum_increment(*dev,1,offsets,points,owner,neighbor,mask,weights,link,
      high,low,high,low,diagonal,limiter,out);
    // An actively clipped high reconstruction follows the same cloud point
    // as the total upper bound, including that point's low word.
    CHECK(dev->download(out)[0]==epsilon);
    dev->upload<Index>(points,std::vector<Index>{0,1,-1});
    auto boundary_high=dev->upload<float>(std::vector<float>{1.5f}),boundary_low=dev->upload<float>(std::vector<float>{epsilon});
    dev->upload<float>(high,std::vector<float>{1,2,1.5f});dev->fill(out,0);
    kern::fv_projected_momentum_increment(*dev,1,offsets,points,owner,neighbor,mask,weights,link,
      high,low,boundary_high,boundary_low,diagonal,limiter,out);
    CHECK(dev->download(out)[0]==4*epsilon);
  }
}

TEST_CASE("FV paired walk predictor matches reconstruction of a representable total field") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto host=test::physics_annulus(8,32,true);const auto geometry=compute_geometry(host);
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    std::vector<float> high,low,total,bhigh,blow,btotal;
    const auto values=[](Vec3d x) {
      const float h=float(1+.15*x.x+.2*x.y),l=float(std::ldexp(1.,-12)*std::sin(17*x.x-13*x.y));
      return std::array<float,3>{h,l,h+l};
    };
    for(auto x:geometry.cell_centroid) {auto v=values(x);high.push_back(v[0]);low.push_back(v[1]);total.push_back(v[2]);}
    for(Index f=mesh.n_internal_faces;f<mesh.n_faces;++f) {auto v=values(geometry.face_centroid[f]);bhigh.push_back(v[0]);blow.push_back(v[1]);btotal.push_back(v[2]);}
    auto h=dev->upload<float>(high),l=dev->upload<float>(low),q=dev->upload<float>(total);
    auto bh=dev->upload<float>(bhigh),bl=dev->upload<float>(blow),bq=dev->upload<float>(btotal);
    auto diagonal=dev->alloc<float>(mesh.n_cells),limiter=dev->alloc<float>(mesh.n_cells);
    auto hf=dev->alloc<float>(mesh.n_faces),lf=dev->alloc<float>(mesh.n_faces),qf=dev->alloc<float>(mesh.n_faces);
    dev->fill(diagonal,1);dev->fill(limiter,1);dev->fill(hf,0);dev->fill(lf,0);dev->fill(qf,0);
    const std::vector<BoundaryObservation> observed(mesh.n_boundary_faces(),BoundaryObservation::FaceValue);
    WalkFaceReconstruction fit(*dev,host,geometry,mesh,stencil,observed);REQUIRE(fit.corrected_faces()>0);
    fit.momentum_values(h,bh,diagonal,hf,&limiter);
    fit.momentum_increment(h,l,bh,bl,diagonal,lf,&limiter);
    fit.momentum_values(q,bq,diagonal,qf,&limiter);
    auto a=dev->download(hf),b=dev->download(lf),c=dev->download(qf);const auto mask=dev->download(fit.mask());
    double maximum=0;
    for(Index f=0;f<mesh.n_internal_faces;++f) if(mask[f]) maximum=std::max(maximum,std::abs(double(a[f])+double(b[f])-double(c[f])));
    // Full-field float rounding and the independent moment accumulation set
    // this bound; the independently clipped low field misses O(2^-12).
    INFO("backend "<<backend_name(dev->backend())<<" maximum paired predictor difference "<<maximum);
    CHECK(maximum<3e-6);
  }
}
