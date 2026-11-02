#include "kutta/fv/transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {

TransportEquation::TransportEquation(Device& dev,const DeviceMesh& m,const Stencil& s,
                                     const PolyMesh&,const MeshGeometry&,
                                     std::shared_ptr<CsrMatrix> shared_matrix,
                                     std::shared_ptr<const FaceReconstruction> diffusion_reconstruction,
                                     std::shared_ptr<KrylovWorkspace> workspace)
    :matrix_storage_(shared_matrix?std::move(shared_matrix):std::make_shared<CsrMatrix>(upload(dev,s.pattern))),
     matrix(*matrix_storage_),rhs(dev.alloc<float>(m.n_cells)),
     diagonal(dev.alloc<float>(m.n_cells)),consistent_inverse(dev.alloc<float>(m.n_cells)),
     grad(vector_field(dev,m.n_cells)),limiter(dev.alloc<float>(m.n_cells)),
     face_correction(dev.alloc<float>(m.n_faces)),dev_(dev),mesh_(m),stencil_(s) {
  KUTTA_REQUIRE(matrix.n_rows==m.n_cells && matrix.nnz==s.pattern.nnz(),"FV shared matrix pattern dimensions mismatch");
  // Jacobi/BiCGStab setup needs only a valid pattern and diagonal. Values are
  // refreshed from the assembled device matrix before every solve.
  CsrMatrixHost initial{s.pattern,std::vector<double>(s.pattern.nnz(),0)};
  for(auto j:s.pattern.diag) initial.val[j]=1;
  SolverOptions options;options.workspace=std::move(workspace);
  solver_=make_solver(transport_solver_kind,dev,initial,options);
  zero_=dev.alloc<float>(m.n_cells);dev.fill(zero_,0);
  set_reconstruction(std::move(diffusion_reconstruction));
}

void TransportEquation::set_reconstruction(std::shared_ptr<const FaceReconstruction> reconstruction) {
  diffusion_reconstruction_=std::move(reconstruction);
  if(diffusion_reconstruction_ && diffusion_correction_.size()!=std::size_t(mesh_.n_faces))
    diffusion_correction_=dev_.alloc<float>(mesh_.n_faces);
}

void TransportEquation::enable_compensated_assembly() {
  rhs_low_=dev_.alloc<float>(mesh_.n_cells);dev_.fill(rhs_low_,0);
  diagonal_low_=dev_.alloc<float>(mesh_.n_cells);dev_.fill(diagonal_low_,0);
}

void TransportEquation::enable_compensated_solve() {
  enable_compensated_assembly();
  scalar_delta_=dev_.alloc<float>(mesh_.n_cells);
  scalar_residual_=dev_.alloc<float>(mesh_.n_cells);
}

void TransportEquation::physical_residual(const Buffer<float>& field,Buffer<float>& result,const Buffer<Index>* excluded) {
  if(rhs_low_.empty()) {residual(dev_,matrix,field,rhs,result);return;}
  kern::fv_compensated_momentum_residual(dev_,mesh_.n_cells,matrix.row_ptr,matrix.col,matrix.val,matrix.diag,
    diagonal,diagonal_low_,field,zero_,zero_,mesh_.volume,rhs,rhs_low_,result);
  if(excluded) kern::fv_scalar_set_fixed(dev_,mesh_.n_cells,*excluded,zero_,result);
}

void TransportEquation::set_steady_tolerance(double tolerance) {
  KUTTA_REQUIRE(std::isfinite(tolerance) && tolerance>=0,"FV invalid scalar steady tolerance");
  if(tolerance==0) {steady_stopping_.reset();return;}
  steady_stopping_.emplace(SteadyStopping{tolerance,dev_.alloc<float>(mesh_.n_cells),
    dev_.alloc<float>(mesh_.n_cells),dev_.alloc<float>(mesh_.n_cells)});
}

void TransportEquation::assemble(const Buffer<float>& field,const Buffer<float>& bc,
    const Buffer<Index>& kind,const Buffer<float>& flux,const Buffer<float>& mu,
    const Buffer<float>& source,const TransportSettings& settings,TransportTerms terms) {
  gradient(dev_,mesh_,stencil_,field,bc,grad);
  assemble_reconstructed(field,bc,kind,flux,mu,source,settings,terms);
}

void TransportEquation::assemble_reconstructed(const Buffer<float>& field,const Buffer<float>& bc,
    const Buffer<Index>& kind,const Buffer<float>& flux,const Buffer<float>& mu,
    const Buffer<float>& source,const TransportSettings& settings,TransportTerms terms) {
  const auto& m=mesh_; const auto& s=stencil_; auto& dev=dev_;
  KUTTA_REQUIRE(field.size()==std::size_t(m.n_cells) && source.size()==field.size() &&
    bc.size()==std::size_t(m.n_boundary_faces()) && kind.size()==bc.size() &&
    flux.size()==std::size_t(m.n_faces) && mu.size()==flux.size(),"FV transport field size mismatch");
  KUTTA_REQUIRE(settings.density>0 && settings.pseudo_time_step>=0 && settings.relaxation>0 &&
    settings.relaxation<=1,"FV invalid transport relaxation");
  for(const auto* term:{terms.sink,terms.diagonal,terms.rhs})
    KUTTA_REQUIRE(!term || term->size()==field.size(),"FV transport term size mismatch");
  if(settings.limit) limit_gradient(dev,m,s,field,bc,grad,limiter,settings.limiter_low,&kind,
    settings.limiter_epsilon2); else dev.fill(limiter,1);
  // Keep the legacy dispatch (including its shader signature) for models
  // that have not enabled the projected reconstruction.
  if(diffusion_reconstruction_) {
    diffusion_reconstruction_->diffusion_correction(field,bc,diffusion_correction_);
    kern::fv_projected_transport_faces(dev,m.n_faces,m.owner,m.neighbour,m.weight,
      m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
      m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
      s.nonorthogonal.x,s.nonorthogonal.y,s.nonorthogonal.z,
      grad.x,grad.y,grad.z,limiter,flux,mu,kind,
      diffusion_reconstruction_->mask(),diffusion_correction_,face_correction,m.n_internal_faces,
      settings.convection==ConvectionScheme::SecondOrderUpwind?1:0);
  } else {
    kern::fv_transport_faces(dev,m.n_faces,m.owner,m.neighbour,m.weight,
      m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
      m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
      s.nonorthogonal.x,s.nonorthogonal.y,s.nonorthogonal.z,
      grad.x,grad.y,grad.z,limiter,flux,mu,kind,face_correction,m.n_internal_faces,
      settings.convection==ConvectionScheme::SecondOrderUpwind?1:0);
  }
  kern::fv_transport_assemble(dev,m.n_cells,m.cell_face_offsets,m.cell_faces,m.cell_face_sign,
    s.adjacent,s.matrix_entry,matrix.row_ptr,matrix.diag,s.diffusion,flux,mu,kind,bc,face_correction,
    source,m.volume,field,terms.sink?*terms.sink:zero_,terms.diagonal?*terms.diagonal:zero_,
    terms.rhs?*terms.rhs:zero_,matrix.val,rhs,rhs_low_.empty()?zero_:rhs_low_,diagonal,
    diagonal_low_.empty()?zero_:diagonal_low_,consistent_inverse,m.n_internal_faces,
    settings.density,settings.pseudo_time_step,settings.relaxation,rhs_low_.empty()?0:1);
}

SolveStats TransportEquation::solve(Buffer<float>& field,const SolveControl& requested,const Buffer<Index>* excluded) {
  auto control=requested;
  if(!scalar_delta_.empty()) physical_residual(field,scalar_residual_,excluded);
  if(steady_stopping_) {
    KUTTA_REQUIRE(!excluded || excluded->size()==field.size(),"FV scalar exclusion size mismatch");
    auto& work=*steady_stopping_;
    if(scalar_delta_.empty()) residual(dev_,matrix,field,rhs,work.residual);
    else dev_.copy(work.residual,scalar_residual_);
    kern::fv_scalar_steady_measure(dev_,mesh_.n_cells,field,diagonal,work.residual,
      excluded?*excluded:matrix.diag,work.absolute,work.scale,excluded?1:0);
    const auto values=dev_.reduce_many({Device::HostReduction::sum(work.absolute,mesh_.n_cells),
      Device::HostReduction::sum(work.scale,mesh_.n_cells),
      Device::HostReduction::dot(work.residual,work.residual,mesh_.n_cells)});
    if(values[0]<=work.tolerance*values[1] && std::isfinite(values[2]))
      control.abs_tol=std::max(control.abs_tol,std::nextafter(std::sqrt(values[2]),std::numeric_limits<double>::infinity()));
  }
  solver_->refresh(dev_,matrix);
  if(scalar_delta_.empty()) return solver_->solve(dev_,matrix,field,rhs,control);
  dev_.fill(scalar_delta_,0);
  auto result=solver_->solve(dev_,matrix,scalar_delta_,scalar_residual_,control);
  kern::axpbyz(dev_,mesh_.n_cells,field,field,scalar_delta_,1,1);
  return result;
}
}  // namespace kutta::fv
