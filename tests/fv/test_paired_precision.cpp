#include "../backend_helpers.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kernels.hpp"

#include <cmath>
#include <vector>

using namespace kutta;
using namespace kutta::fv;

TEST_CASE("FV compensated momentum RHS retains the steady defect under strong pseudo time") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1,1,1});
  const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);
    auto stencil=make_stencil(*dev,host,geometry);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);
    equation.enable_compensated_assembly();
    auto q=dev->alloc<float>(1),source=dev->alloc<float>(1),zero=dev->alloc<float>(1);
    auto reaction=dev->alloc<float>(1),defect=dev->alloc<float>(1);
    auto boundary=dev->alloc<float>(6),flux=dev->alloc<float>(6),mu=dev->alloc<float>(6);
    auto kind=dev->alloc<Index>(6);
    dev->fill(source,.125f);dev->fill(reaction,1);dev->fill(zero,0);
    dev->fill(boundary,0);dev->fill(flux,0);dev->fill(mu,0);dev->fill(kind,1);
    for(float value:{1.f,.1f,.001f}) {
      dev->fill(q,value);
      TransportSettings settings;settings.pseudo_time_step=1e-8f;
      equation.assemble(q,boundary,kind,flux,mu,source,settings,{nullptr,&reaction});
      kern::fv_compensated_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
        equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,zero,zero,mesh.volume,
        equation.rhs,equation.rhs_low(),defect);
      // d(q)/dt = .125 - q. The steady defect is independent of the
      // pseudo-time diagonal, even when that diagonal is eight orders larger.
      CHECK(std::abs(double(dev->download(defect)[0])-(.125-double(value)))<1e-8);
    }
  }
}

TEST_CASE("FV paired momentum preserves constant diffusion across disparate face coefficients") {
  auto host=make_box_hex(4,1,1,{0,0,0},{1,1,1});
  const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);
    auto stencil=make_stencil(*dev,host,geometry);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);
    auto q=dev->alloc<float>(mesh.n_cells),zero=dev->alloc<float>(mesh.n_cells),defect=dev->alloc<float>(mesh.n_cells);
    auto boundary=dev->alloc<float>(mesh.n_boundary_faces()),flux=dev->alloc<float>(mesh.n_faces);
    auto kind=dev->alloc<Index>(mesh.n_boundary_faces());
    std::vector<float> diffusivity(mesh.n_faces,0);
    REQUIRE(mesh.n_internal_faces==3);
    diffusivity[0]=1;diffusivity[1]=1e-8f;diffusivity[2]=1e-12f;
    auto mu=dev->upload<float>(diffusivity);
    dev->fill(q,1);dev->fill(zero,0);dev->fill(boundary,1);dev->fill(flux,0);dev->fill(kind,1);
    TransportSettings settings;settings.pseudo_time_step=1;settings.relaxation=.7f;
    equation.assemble(q,boundary,kind,flux,mu,zero,settings);
    const auto linear_operator=dev->download(equation.matrix.val);
    equation.enable_compensated_assembly();
    equation.assemble(q,boundary,kind,flux,mu,zero,settings);
    CHECK(dev->download(equation.matrix.val)==linear_operator);
    kern::fv_compensated_momentum_residual(*dev,mesh.n_cells,equation.matrix.row_ptr,equation.matrix.col,
      equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,zero,zero,mesh.volume,
      equation.rhs,equation.rhs_low(),defect);
    // No gradients, source or boundary flux: the physical steady defect is zero.
    // Keeping the linear operator positive must not add a reaction to it.
    CHECK(dev->max_abs(defect,mesh.n_cells)<1e-12);
  }
}

TEST_CASE("FV pressure jump retains curvature below one pressure ulp") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto owner=dev->upload<Index>(std::vector<Index>{0});
    auto neighbour=dev->upload<Index>(std::vector<Index>{1});
    auto px=dev->upload<float>(std::vector<float>{.25f}),nx=dev->upload<float>(std::vector<float>{-.75f});
    auto zero_face=dev->alloc<float>(1),zero_cell=dev->alloc<float>(2);
    dev->fill(zero_face,0);dev->fill(zero_cell,0);
    const float slope=std::ldexp(1.f,-20),curvature=std::ldexp(1.f,-30);
    auto pressure=dev->alloc<float>(2);
    auto gradient=dev->upload<float>(std::vector<float>{slope+curvature,slope-curvature});
    auto boundary=dev->alloc<float>(0),face=dev->alloc<float>(1),jump=dev->alloc<float>(1);
    for(float offset:{0.f,1.f}) {
      dev->upload<float>(pressure,std::vector<float>{offset,offset+slope});
      kern::fv_pressure_faces(*dev,1,owner,neighbour,px,zero_face,zero_face,nx,zero_face,zero_face,
        pressure,gradient,zero_cell,zero_cell,boundary,face,jump,1);
      kern::fv_compensated_pressure_jump(*dev,1,owner,neighbour,px,zero_face,zero_face,nx,zero_face,zero_face,
        pressure,gradient,zero_cell,zero_cell,boundary,jump,1);
      // Adding a pressure gauge offset must preserve the curvature jump.
      CHECK(std::abs(double(dev->download(jump)[0])+std::ldexp(1.,-31))<1e-14);
    }
  }
}

TEST_CASE("FV paired momentum boundary preserves a constant advected diffused field") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1,1,1});const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);equation.enable_compensated_assembly();
    auto q=dev->alloc<float>(1),zero=dev->alloc<float>(1),defect=dev->alloc<float>(1);
    auto boundary=dev->alloc<float>(6);auto kind=dev->alloc<Index>(6);
    std::vector<float> flux_values(6,0),mu_values(6,0);
    flux_values[0]=-std::ldexp(1.f,-25);flux_values[1]=-flux_values[0];mu_values[0]=.5f;
    auto flux=dev->upload<float>(flux_values),mu=dev->upload<float>(mu_values);
    dev->fill(zero,0);dev->fill(kind,0);
    for(float constant:{1.f,.1f}) {
      dev->fill(q,constant);dev->fill(boundary,constant);
      TransportSettings settings;settings.pseudo_time_step=1;settings.relaxation=.7f;
      equation.assemble(q,boundary,kind,flux,mu,zero,settings);
      kern::fv_compensated_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
        equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,zero,zero,mesh.volume,
        equation.rhs,equation.rhs_low(),defect);
      // Constant scalar, zero source and balanced mass flux: neither diffusion
      // nor advection can produce a physical defect, however different D and F.
      CHECK(dev->max_abs(defect,1)<1e-12);
    }
  }
}

TEST_CASE("FV compensated source retains cancellation with an integrated RHS") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1-std::ldexp(1.,-13),1,1});const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);equation.enable_compensated_assembly();
    auto q=dev->alloc<float>(1),source=dev->alloc<float>(1),extra=dev->alloc<float>(1),defect=dev->alloc<float>(1);
    auto boundary=dev->alloc<float>(6),flux=dev->alloc<float>(6),mu=dev->alloc<float>(6);auto kind=dev->alloc<Index>(6);
    const float value=1+std::ldexp(1.f,-13);
    dev->fill(q,0);dev->fill(source,value);dev->fill(extra,-1);
    dev->fill(boundary,0);dev->fill(flux,0);dev->fill(mu,0);dev->fill(kind,1);
    TransportSettings settings;settings.pseudo_time_step=1;
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{nullptr,nullptr,&extra});
    kern::fv_compensated_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
      equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,q,q,mesh.volume,
      equation.rhs,equation.rhs_low(),defect);
    const double expected=double(value)*double(dev->download(mesh.volume)[0])-1;
    REQUIRE(std::abs(expected)>1e-10);
    CHECK(std::abs(double(dev->download(defect)[0])-expected)<1e-13);
    dev->fill(q,1);
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{&source});
    auto zero=dev->alloc<float>(1);dev->fill(zero,0);
    kern::fv_compensated_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
      equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,zero,zero,mesh.volume,
      equation.rhs,equation.rhs_low(),defect);
    // Equal production and reaction at q=1 cancel even when their products
    // with cell volume do not fit in a single float.
    CHECK(dev->max_abs(defect,1)<1e-12);
    equation.enable_compensated_solve();settings.pseudo_time_step=0;
    auto reaction=dev->alloc<float>(1);dev->fill(reaction,1);dev->fill(q,0);
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{&reaction,nullptr,&extra});
    const auto result=equation.solve(q,{1e-4,1e-20,20});
    CHECK(result.converged);
    CHECK(std::abs(double(dev->download(q)[0])-expected/double(dev->download(mesh.volume)[0]))<1e-14);
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{&reaction,nullptr,&extra});
    equation.physical_residual(q,defect);
    CHECK(dev->max_abs(defect,1)<1e-13);
  }
}

TEST_CASE("FV scalar steady stopping avoids churn and resumes for a changed source") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1,1,1});const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);
    equation.set_steady_tolerance(1e-6);
    auto q=dev->alloc<float>(1),source=dev->alloc<float>(1),sink=dev->alloc<float>(1);
    auto boundary=dev->alloc<float>(6),flux=dev->alloc<float>(6),mu=dev->alloc<float>(6);
    auto kind=dev->alloc<Index>(6);
    dev->fill(q,1);dev->fill(sink,1);dev->fill(source,1+std::ldexp(1.f,-22));
    dev->fill(boundary,0);dev->fill(flux,0);dev->fill(mu,0);dev->fill(kind,1);
    TransportSettings settings;
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{&sink});
    const auto settled=equation.solve(q,{1e-3,1e-12,20});
    CHECK(settled.iterations==0);CHECK(settled.converged);CHECK(dev->download(q)[0]==1);
    dev->fill(source,1.01f);
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{&sink});
    const auto resumed=equation.solve(q,{1e-3,1e-12,20});
    CHECK(resumed.iterations>0);CHECK(resumed.converged);
    CHECK(dev->download(q)[0]==doctest::Approx(1.01f).epsilon(1e-6));
  }
}

TEST_CASE("FV paired flux normalization retains advection of the low scalar state") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1,1,1});const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    const float flux_increment=std::ldexp(1.f,-10),scalar_increment=std::ldexp(1.f,-20);
    std::vector<float> high_values(6,0),low_values(6,0);
    high_values[0]=-1;high_values[1]=1;low_values[0]=-flux_increment;low_values[1]=flux_increment;
    auto flux=dev->upload<float>(high_values),low_flux=dev->upload<float>(low_values),previous=dev->alloc<float>(6);
    dev->fill(previous,0);
    kern::fv_paired_flux_offset(*dev,6,mesh.area.x,mesh.area.y,mesh.area.z,previous,flux,low_flux,1e-7f);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);equation.enable_compensated_assembly();
    auto q=dev->alloc<float>(1),q_low=dev->alloc<float>(1),zero=dev->alloc<float>(1),defect=dev->alloc<float>(1);
    auto boundary=dev->alloc<float>(6),mu=dev->alloc<float>(6);auto kind=dev->alloc<Index>(6);
    dev->fill(q,1);dev->fill(q_low,scalar_increment);dev->fill(zero,0);
    dev->fill(boundary,1);dev->fill(mu,0);dev->fill(kind,0);
    TransportSettings settings;settings.convection=ConvectionScheme::FirstOrderUpwind;
    equation.assemble(q,boundary,kind,flux,mu,zero,settings);
    kern::fv_compensated_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
      equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,q_low,zero,mesh.volume,
      equation.rhs,equation.rhs_low(),defect);
    // High scalar and inlet data are both one, so low-flux transport of that
    // high state cancels. The remaining defect is -F_total * scalar_increment,
    // including the flux_increment * scalar_increment cross term.
    const double expected=-(1+double(flux_increment))*double(scalar_increment);
    CHECK(std::abs(double(dev->download(defect)[0])-expected)<1e-13);
    CHECK(dev->max_abs(low_flux,6)<1e-7);
  }
}

TEST_CASE("FV paired physical reaction of low state survives strong pseudo time") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1,1,1});const auto geometry=compute_geometry(host);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    TransportEquation equation(*dev,mesh,stencil,host,geometry);equation.enable_compensated_assembly();
    auto q=dev->alloc<float>(1),low=dev->alloc<float>(1),zero=dev->alloc<float>(1);
    auto source=dev->alloc<float>(1),sink=dev->alloc<float>(1),defect=dev->alloc<float>(1);
    auto boundary=dev->alloc<float>(6),flux=dev->alloc<float>(6),mu=dev->alloc<float>(6);auto kind=dev->alloc<Index>(6);
    dev->fill(q,.125f);dev->fill(source,.125f);dev->fill(sink,1);dev->fill(zero,0);
    dev->fill(boundary,0);dev->fill(flux,0);dev->fill(mu,0);dev->fill(kind,1);
    TransportSettings settings;settings.pseudo_time_step=1e-8f;settings.relaxation=.7f;
    equation.assemble(q,boundary,kind,flux,mu,source,settings,{&sink});
    for(float correction:{-std::ldexp(1.f,-19),std::ldexp(1.f,-19)}) {
      dev->fill(low,correction);
      kern::fv_compensated_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
        equation.matrix.val,equation.matrix.diag,equation.diagonal,equation.diagonal_low(),q,low,zero,mesh.volume,
        equation.rhs,equation.rhs_low(),defect);
      // Physical source .125 - (hi+lo) is -lo. A large pseudo diagonal
      // cannot erase the reaction by subtracting two rounded O(dt^-1 * lo) products.
      CHECK(std::abs(double(dev->download(defect)[0])+double(correction))<1e-13);
    }
  }
}

TEST_CASE("FV compensated pressure force recovers Green Gauss cancellation") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto offsets=dev->upload<Index>(std::vector<Index>{0,4}),faces=dev->upload<Index>(std::vector<Index>{0,1,2,3});
    auto sign=dev->upload<float>(std::vector<float>{1,1,1,1});
    const float small=std::ldexp(1.f,-25);
    auto area=dev->upload<float>(std::vector<float>{1,small,-1,-small});
    auto zero=dev->alloc<float>(4),pressure=dev->alloc<float>(4),volume=dev->alloc<float>(1);
    auto x=dev->alloc<float>(1),y=dev->alloc<float>(1),z=dev->alloc<float>(1);
    dev->fill(zero,0);dev->fill(pressure,1);dev->fill(volume,1);
    kern::fv_compensated_pressure_force(*dev,1,offsets,faces,sign,area,zero,zero,volume,pressure,zero,x,y,z);
    // The closed surface has zero force at uniform pressure. Its rounded
    // high reduction loses the second face and leaves -small; carry it back.
    CHECK(std::abs(double(dev->download(x)[0])-double(small))<1e-14);
    CHECK(dev->download(y)[0]==0);CHECK(dev->download(z)[0]==0);
  }
}

TEST_CASE("FV paired normal backflow carries only the low flux speed") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto owner=dev->upload<Index>(std::vector<Index>{0}),kind=dev->upload<Index>(std::vector<Index>{4});
    auto area=dev->alloc<float>(1),zero=dev->alloc<float>(1),flux=dev->alloc<float>(1),low=dev->alloc<float>(1);
    auto x=dev->alloc<float>(1),y=dev->alloc<float>(1),z=dev->alloc<float>(1);auto scalar_kind=dev->alloc<Index>(1);
    dev->fill(area,2);dev->fill(zero,0);dev->fill(flux,-1);dev->fill(low,std::ldexp(1.f,-20));
    kern::fv_paired_velocity_boundary(*dev,1,owner,kind,area,zero,zero,zero,zero,zero,flux,low,x,y,z,scalar_kind,0,1);
    CHECK(dev->download(x)[0]==std::ldexp(1.f,-21));
    CHECK(dev->download(y)[0]==0);CHECK(dev->download(z)[0]==0);CHECK(dev->download(scalar_kind)[0]==0);
  }
}

TEST_CASE("FV paired pressure force faces retain curvature below the common pressure ulp") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto owner=dev->upload<Index>(std::vector<Index>{0}),neighbor=dev->upload<Index>(std::vector<Index>{1});
    auto kind=dev->upload<Index>(std::vector<Index>{1}),poor=dev->upload<Index>(std::vector<Index>{0,0});
    auto px=dev->upload<float>(std::vector<float>{.25f}),nx=dev->upload<float>(std::vector<float>{-.75f});
    auto face_zero=dev->alloc<float>(1),cell_zero=dev->alloc<float>(2),boundary=dev->upload<float>(std::vector<float>{1.25f});
    auto diagonal=dev->alloc<float>(2),face=dev->alloc<float>(1),face_low=dev->alloc<float>(1);
    dev->fill(face_zero,0);dev->fill(cell_zero,0);dev->fill(diagonal,1);
    const float slope=std::ldexp(1.f,-20),curvature=std::ldexp(1.f,-30);
    auto pressure=dev->upload<float>(std::vector<float>{1,1+slope});
    auto gradient=dev->upload<float>(std::vector<float>{slope+curvature,slope-curvature});
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,cell_zero,gradient,cell_zero,cell_zero,cell_zero,cell_zero,cell_zero,boundary,boundary,
      diagonal,face,face_low,1,0);
    const double expected=1+double(slope)*.25+double(curvature)*.5;
    CHECK(std::abs(double(dev->download(face)[0])+dev->download(face_low)[0]-expected)<1e-14);
    CHECK(dev->download(face_low)[0]==curvature*.5f);
    // Repartition the same physical pressure without changing its face force.
    auto low=dev->upload<float>(std::vector<float>{-slope,-slope});
    dev->upload<float>(pressure,std::vector<float>{1+slope,1+2*slope});
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,low,gradient,cell_zero,cell_zero,cell_zero,cell_zero,cell_zero,boundary,boundary,
      diagonal,face,face_low,1,0);
    CHECK(std::abs(double(dev->download(face)[0])+dev->download(face_low)[0]-expected)<1e-14);
    dev->fill(pressure,1);dev->fill(low,0);
    dev->upload<float>(diagonal,std::vector<float>{1000,.001f});
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,low,gradient,cell_zero,cell_zero,cell_zero,cell_zero,cell_zero,boundary,boundary,
      diagonal,face,face_low,1,1);
    CHECK(double(dev->download(face)[0])+dev->download(face_low)[0]==1);
    auto boundary_low=dev->upload<float>(std::vector<float>{curvature});
    dev->fill(kind,2);
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,low,gradient,cell_zero,cell_zero,cell_zero,cell_zero,cell_zero,boundary,boundary_low,
      diagonal,face,face_low,0,0);
    CHECK(double(dev->download(face)[0])+dev->download(face_low)[0]==1.25+double(curvature));
    dev->fill(kind,1);
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,low,gradient,cell_zero,cell_zero,cell_zero,cell_zero,cell_zero,boundary,boundary_low,
      diagonal,face,face_low,0,0);
    CHECK(std::abs(double(dev->download(face)[0])+dev->download(face_low)[0]-(1+.25*(double(slope)+curvature)))<1e-14);
    dev->upload<Index>(poor,std::vector<Index>{1,0});
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,low,gradient,cell_zero,cell_zero,cell_zero,cell_zero,cell_zero,boundary,boundary_low,
      diagonal,face,face_low,0,0);
    CHECK(double(dev->download(face)[0])+dev->download(face_low)[0]==1);
    dev->fill(poor,0);dev->upload<float>(pressure,std::vector<float>{1,1+slope});dev->fill(gradient,slope);
    const float tiny_curvature=std::ldexp(1.f,-45);
    auto low_gradient=dev->upload<float>(std::vector<float>{tiny_curvature,-tiny_curvature});
    kern::fv_paired_pressure_faces(*dev,1,owner,neighbor,kind,poor,px,face_zero,face_zero,nx,face_zero,face_zero,
      pressure,low,gradient,cell_zero,cell_zero,low_gradient,cell_zero,cell_zero,boundary,boundary_low,
      diagonal,face,face_low,1,0);
    CHECK(dev->download(face_low)[0]==.5f*tiny_curvature);
  }
}

TEST_CASE("FV compensated flux retains the physical upwind sign through hysteresis") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    const float high=std::ldexp(1.f,-20),value=-std::ldexp(1.f,-40);
    auto previous=dev->upload<float>(std::vector<float>{high}),flux=dev->upload<float>(std::vector<float>{value});
    auto low=dev->alloc<float>(1),one=dev->alloc<float>(1),zero=dev->alloc<float>(1);
    dev->fill(low,0);dev->fill(one,1);dev->fill(zero,0);
    kern::fv_paired_flux_offset(*dev,1,one,zero,zero,previous,flux,low,std::ldexp(1.f,-16));
    CHECK(dev->download(flux)[0]==value);
    CHECK(double(dev->download(flux)[0])+dev->download(low)[0]==value);
    auto owner=dev->upload<Index>(std::vector<Index>{0}),neighbor=dev->upload<Index>(std::vector<Index>{1});
    auto correction=dev->upload<float>(std::vector<float>{-2*high,0});dev->fill(flux,high);dev->fill(low,0);
    kern::fv_compensated_correct_flux(*dev,1,owner,neighbor,one,zero,zero,correction,one,flux,low,1,std::ldexp(1.f,-16));
    CHECK(dev->download(flux)[0]==-high);
  }
}

TEST_CASE("FV paired pressure gradient carries cancellation and is invariant to rebasing") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    const float small=std::ldexp(1.f,-25),curvature=std::ldexp(1.f,-40),shift=std::ldexp(1.f,-20);
    auto offsets=dev->upload<Index>(std::vector<Index>{0,4});
    auto adjacent=dev->upload<Index>(std::vector<Index>{-1,-2,-3,-4});
    auto zero=dev->alloc<float>(4),one=dev->alloc<float>(4),cell=dev->alloc<float>(1);
    dev->fill(zero,0);dev->fill(one,1);dev->fill(cell,0);
    auto inverse=dev->upload<float>(std::vector<float>{1,0,0,1,0,1});
    auto boundary=dev->upload<float>(std::vector<float>{1,small,-1,-small});
    auto low=dev->upload<float>(std::vector<float>{curvature,0,0,0});
    auto high=vector_field(*dev,1),carry=vector_field(*dev,1);
    for(int partition=0;partition<2;++partition) {
      kern::fv_paired_gradient(*dev,1,offsets,adjacent,zero,zero,zero,one,zero,zero,inverse,
        cell,cell,boundary,low,high.x,high.y,high.z,carry.x,carry.y,carry.z,1);
      CHECK(double(dev->download(high.x)[0])+dev->download(carry.x)[0]==double(curvature));
      CHECK(dev->download(high.y)[0]==0);CHECK(dev->download(high.z)[0]==0);
      dev->upload<float>(boundary,std::vector<float>{1+shift,small+shift,-1+shift,-small+shift});
      dev->upload<float>(low,std::vector<float>{curvature-shift,-shift,-shift,-shift});
    }
  }
}

TEST_CASE("FV paired Neumann pressure trace cancels its Rhie Chow jump") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto owner=dev->upload<Index>(std::vector<Index>{0}),neighbor=dev->upload<Index>(std::vector<Index>{-1});
    auto kind=dev->upload<Index>(std::vector<Index>{1});
    auto p=dev->upload<float>(std::vector<float>{.125f}),low=dev->upload<float>(std::vector<float>{std::ldexp(1.f,-30)});
    auto g=dev->upload<float>(std::vector<float>{std::ldexp(1.f,-26)}),gl=dev->upload<float>(std::vector<float>{std::ldexp(1.f,-40)});
    auto r=dev->upload<float>(std::vector<float>{.5f}),zero=dev->alloc<float>(1);
    auto bc=dev->alloc<float>(1),bl=dev->alloc<float>(1),jh=dev->alloc<float>(1),jl=dev->alloc<float>(1);
    dev->fill(zero,0);
    kern::fv_pressure_boundary(*dev,1,owner,kind,r,zero,zero,p,g,zero,zero,zero,bc,0);
    kern::fv_paired_pressure_boundary_low(*dev,1,owner,kind,r,zero,zero,p,low,g,zero,zero,gl,zero,zero,bc,bl,0);
    kern::fv_compensated_pressure_jump(*dev,1,owner,neighbor,r,zero,zero,zero,zero,zero,p,g,zero,zero,bc,jh,0);
    kern::fv_compensated_pressure_jump(*dev,1,owner,neighbor,r,zero,zero,zero,zero,zero,low,gl,zero,zero,bl,jl,0);
    CHECK(std::abs(double(dev->download(jh)[0])+dev->download(jl)[0])<1e-16);
  }
}

TEST_CASE("FV compensated Rhie Chow retains velocity dot products below a cancelling high flux") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto owner=dev->upload<Index>(std::vector<Index>{0}),neighbor=dev->upload<Index>(std::vector<Index>{1});
    auto kind=dev->upload<Index>(std::vector<Index>{1});
    auto one=dev->alloc<float>(2),zero=dev->alloc<float>(2),small=dev->alloc<float>(1),negative=dev->alloc<float>(1);
    dev->fill(one,1);dev->fill(zero,0);dev->fill(small,std::ldexp(1.f,-25));dev->fill(negative,-1);
    auto high=dev->alloc<float>(1),low=dev->alloc<float>(1);
    kern::fv_compensated_rhie_chow(*dev,1,owner,neighbor,kind,one,one,one,zero,one,one,
      one,small,negative,zero,zero,zero,zero,high,low,1,1);
    CHECK(double(dev->download(high)[0])+dev->download(low)[0]==std::ldexp(1.,-25));
    // The physical pressure curvature lives entirely in the carried state.
    const float slope=std::ldexp(1.f,-20),curve=std::ldexp(1.f,-30);
    auto p=dev->upload<float>(std::vector<float>{1,1+slope}),pl=dev->upload<float>(std::vector<float>{0,-slope});
    auto g=dev->upload<float>(std::vector<float>{slope+curve,slope-curve});
    auto gl=dev->upload<float>(std::vector<float>{-slope,-slope});
    auto px=dev->upload<float>(std::vector<float>{.25f}),nx=dev->upload<float>(std::vector<float>{-.75f});
    kern::fv_paired_pressure_jump(*dev,1,owner,neighbor,px,zero,zero,nx,zero,zero,p,pl,g,zero,zero,gl,zero,zero,
      zero,zero,high,1);
    CHECK(std::abs(double(dev->download(high)[0])+std::ldexp(1.,-31))<1e-14);
  }
}

TEST_CASE("FV compensated SIMPLEC pressure update preserves a local correction") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    const float high=std::ldexp(1.f,-20),low=-high+std::ldexp(1.f,-43),delta=std::ldexp(1.f,-60);
    auto p=dev->upload<float>(std::vector<float>{high}),pl=dev->upload<float>(std::vector<float>{low});
    auto correction=dev->upload<float>(std::vector<float>{2*delta}),one=dev->alloc<float>(1),zero=dev->alloc<float>(1);
    auto u=vector_field(*dev,1),ul=vector_field(*dev,1);dev->fill(one,1);dev->fill(zero,0);
    kern::fv_compensated_correct_velocity(*dev,1,one,one,one,zero,zero,zero,correction,
      u.x,u.y,u.z,p,ul.x,ul.y,ul.z,pl,.5f,1e-6f,16*std::numeric_limits<float>::epsilon());
    CHECK(double(dev->download(p)[0])+dev->download(pl)[0]==std::ldexp(1.,-43)+double(delta));
  }
}

TEST_CASE("FV constitutive derivatives retain a subulp physical velocity variation") {
  auto host=make_box_hex(1,1,1,{0,0,0},{1,1,1});const auto geometry=compute_geometry(host);
  const float slope=std::ldexp(1.f,-27),shift=std::ldexp(1.f,-20);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=upload_mesh(*dev,host,geometry);auto stencil=make_stencil(*dev,host,geometry);
    std::vector<BoundaryObservation> observations(6,BoundaryObservation::FaceValue);
    auto fit=boundary_least_squares(*dev,host,geometry,observations,&stencil);
    auto jump=face_jump_least_squares(*dev,host,geometry);
    auto q=dev->alloc<float>(1),low=dev->alloc<float>(1),bc=dev->alloc<float>(6),bl=dev->alloc<float>(6);
    auto physical=dev->alloc<float>(1),physical_bc=dev->alloc<float>(6);
    auto g=vector_field(*dev,1),gl=vector_field(*dev,1),out=vector_field(*dev,1);
    for(int partition=0;partition<2;++partition) {
      const float offset=partition?shift:0;dev->fill(q,1+offset);dev->fill(low,-offset);dev->fill(bc,1+offset);
      std::vector<float> boundary_low;
      for(const auto& face:geometry.face_centroid) boundary_low.push_back(slope*float(face.x-.5)-offset);
      dev->upload<float>(bl,boundary_low);
      kern::fv_paired_value(*dev,1,q,low,physical);kern::fv_paired_value(*dev,6,bc,bl,physical_bc);
      gradient(*dev,mesh,stencil,fit,physical,physical_bc,out);
      CHECK(dev->download(out.x)[0]==0); // The old order erases the physical derivative.
      paired_gradient(*dev,mesh,stencil,fit,q,low,bc,bl,g,gl);
      CHECK(double(dev->download(g.x)[0])+dev->download(gl.x)[0]==double(slope));
      paired_jump_corrected_gradient(*dev,mesh,stencil,jump,q,low,bc,bl,g,gl,out);
      CHECK(dev->download(out.x)[0]==slope);
      CHECK(dev->download(out.y)[0]==0);CHECK(dev->download(out.z)[0]==0);
    }
  }
}
