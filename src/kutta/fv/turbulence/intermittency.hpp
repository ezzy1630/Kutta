#pragma once
#include <map>
#include "kutta/fv/steady_solver.hpp"

namespace kutta::fv {
// Gamma equation of Fluent's Intermittency Transition Model (TG 4.8.2); its
// k-equation coupling is SstOptions::intermittency. Owns gamma and its
// transport; the flow driver calls update() after each SstModel::update().
//
//   IntermittencyTransport t(dev, mesh, stencil, host, geom, bcs, settings, matrix, k_lsq);
//   sst.update({..., &t.gamma()});
//   t.update(sst.fields(), rho, mu, k, omega, velocity_gradient, boundary_kind, mass);
//   t.solve(time_step);  // after k and omega
class IntermittencyTransport {
 public:
  // k_lsq: the SST k fit. Gamma has k's boundary observations (fixed on
  // inflow faces, zero flux on walls), so it shares that factorization.
  IntermittencyTransport(Device&,const DeviceMesh&,const Stencil&,const PolyMesh&,const MeshGeometry&,
    const std::map<std::string,BoundaryCondition>&,const SolverSettings&,std::shared_ptr<CsrMatrix>,
    const BoundaryLeastSquares& k_lsq,std::shared_ptr<KrylovWorkspace> workspace={});
  void initialize(std::span<const double> gamma,const MeshPermutation&);
  void fill(double gamma);
  void update(const SstFields&,const Buffer<float>& rho,const Buffer<float>& mu,
    const Buffer<float>& k,const Buffer<float>& omega,const DeviceVec3 (&velocity_gradient)[3],
    const Buffer<Index>& boundary_kind,const Buffer<float>& mass);
  SolveStats solve(float time_step);
  ResidualValue monitor(float time_step);
  const Buffer<float>& gamma() const {return gamma_;}
  const IntermittencyFields& fields() const {return fields_;}
 private:
  void assemble(float time_step);
  Device& dev_;const DeviceMesh& mesh_;const Stencil& stencil_;const BoundaryLeastSquares& lsq_;
  IntermittencySettings settings_;float rho_,k_floor_,omega_floor_;bool limit_;
  TransportEquation equation_;
  IntermittencyFields fields_;
  Buffer<float> gamma_,boundary_,inflow_,residual_,absolute_,scale_,fluent_scale_;
  Buffer<Index> kind_,wall_,zero_;
  const Buffer<float>* mass_=nullptr;
};
} // namespace kutta::fv
