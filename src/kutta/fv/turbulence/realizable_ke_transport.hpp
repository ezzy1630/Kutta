#pragma once
#include <map>
#include "kutta/fv/steady_solver.hpp"
namespace kutta::fv {
class RealizableKeTransport {
 public:
  RealizableKeTransport(Device&,const DeviceMesh&,const Stencil&,const PolyMesh&,const MeshGeometry&,
    const std::map<std::string,BoundaryCondition>&,const SolverSettings&,std::shared_ptr<CsrMatrix>);
  void initialize(std::span<const double>,std::span<const double>,const MeshPermutation&);
  void update(const DeviceVec3&,const DeviceVec3&,SstVelocityGradients,const Buffer<Index>&,
    const Buffer<float>& mass,Buffer<float>& face_viscosity);
  SolveStats solve(bool epsilon,float time_step);
  std::array<ResidualValue,2> monitor(float time_step);
  const RealizableKeFields& fields() const {return model_.fields();}
  const Buffer<float>& kinetic() const {return k_;}
  const Buffer<float>& epsilon() const {return epsilon_;}
 private:
  void assemble(bool epsilon,float time_step);
  Device& d_;const DeviceMesh& mesh_;const Stencil& stencil_;RealizableKeSettings settings_;float rho_,mu_;bool limit_;
  RealizableKeModel model_;TransportEquation k_equation_,e_equation_;BoundaryLeastSquares k_lsq_,e_lsq_;
  Buffer<float> k_,epsilon_,rho_cell_,mu_cell_,bk_,be_,a_,b_,residual_,absolute_,scale_,fluent_scale_;
  Buffer<Index> spec_,kind_,zero_;const Buffer<float>* mass_=nullptr;
};
} // namespace kutta::fv
