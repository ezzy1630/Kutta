#pragma once
#include <map>
#include "kutta/fv/steady_solver.hpp"

namespace kutta::fv {
// Owns SA-only storage and assembly; shared flow-driver hooks stay small.
class SaTransport {
 public:
  SaTransport(Device&,const DeviceMesh&,const Stencil&,const PolyMesh&,const MeshGeometry&,
    const std::map<std::string,BoundaryCondition>&,const SolverSettings&,std::shared_ptr<CsrMatrix>,
    std::shared_ptr<const WalkFaceReconstruction> shared_faces={});
  void initialize(std::span<const double>,const MeshPermutation&);
  void update(const DeviceVec3&,const DeviceVec3&,const DeviceVec3&,const DeviceVec3&,const DeviceVec3&,const Buffer<Index>&,
    const Buffer<float>& mass,Buffer<float>& face_viscosity);
  SolveStats solve(float time_step);
  ResidualValue monitor(float time_step);
  const SaFields& fields() const {return model_.fields();}
  const Buffer<float>& modified_viscosity() const {return q_;}
  const SaWallFields& wall_fields() const {return wall_fields_;}
 private:
  void assemble(float time_step);
  Device& dev_;const DeviceMesh& mesh_;const Stencil& stencil_;SaSettings settings_;
  float rho_,mu_;bool limit_;
  SaModel model_;TransportEquation equation_;BoundaryLeastSquares lsq_;
  Buffer<float> q_,rho_cell_,mu_cell_,strain_,rotation_,boundary_,values_,length_,residual_,absolute_,scale_,fluent_scale_;
  Buffer<Index> spec_,kind_,zero_;
  Buffer<float> wall_mask_;SaWallFields wall_fields_;
  const Buffer<float>* mass_=nullptr;
};
} // namespace kutta::fv
