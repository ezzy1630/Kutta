#include "../backend_helpers.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kernels.hpp"

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
      kern::fv_paired_momentum_residual(*dev,1,equation.matrix.row_ptr,equation.matrix.col,
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
    kern::fv_paired_momentum_residual(*dev,mesh.n_cells,equation.matrix.row_ptr,equation.matrix.col,
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
