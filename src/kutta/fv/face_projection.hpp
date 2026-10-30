#pragma once
#include "kutta/fv/operators.hpp"

namespace kutta::fv {
// Reconstruct tangential shifts to the face-normal line, preserving the compact
// two-point normal stencil. Geometry-only quadratic weights; device reductions
// supply diffusion, momentum interpolation and flux-consistent wall gradients.
class ProjectedFaceReconstruction {
 public:
  ProjectedFaceReconstruction(Device&,const PolyMesh&,const MeshGeometry&,
                        std::span<const BoundaryObservation>,
                        std::span<const Index> nonsmooth_cells={});
  void gradient(const Buffer<float>& field,const Buffer<float>& boundary,DeviceVec3& out) const;
  void flux_gradient(const Buffer<float>& field,const Buffer<float>& boundary,const DeviceMesh&,const Stencil&,DeviceVec3& out) const;
  void diffusion_correction(const Buffer<float>& field,const Buffer<float>& boundary,Buffer<float>& out) const;
  void momentum_values(const Buffer<float>& field,const Buffer<float>& boundary,const Buffer<float>& diagonal,
                       Buffer<float>& out,const Buffer<float>* cell_limiter=nullptr) const;
  const Buffer<Index>& mask() const {return mask_;}
  std::size_t corrected_faces() const {return corrected_;}
  std::size_t fallback_faces() const {return fallback_;}
 private:
  Device& dev_;
  Index faces_,internal_;
  std::size_t corrected_=0,fallback_=0;
  Buffer<Index> offsets_,points_,owner_,mask_;
  DeviceVec3 weights_;
  Buffer<float> diffusion_weights_;
  Buffer<float> value_weights_;
  Buffer<float> link_weights_;
  Buffer<Index> neighbour_;
};
} // namespace kutta::fv
