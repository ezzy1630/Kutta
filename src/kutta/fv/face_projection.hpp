#pragma once
#include <memory>
#include <vector>

#include "kutta/fv/operators.hpp"

namespace kutta::fv {
// Face functionals that shift both adjacent cell states tangentially onto the
// face-normal line, preserving the compact two-point normal stencil.
class FaceReconstruction {
 public:
  virtual ~FaceReconstruction()=default;
  // Deferred diffusion correction per unit diffusivity; zero where mask()==0.
  virtual void diffusion_correction(const Buffer<float>& field,const Buffer<float>& boundary,Buffer<float>& out) const=0;
  // Replaces internal masked face values of the momentum-weighted predictor.
  virtual void momentum_values(const Buffer<float>& field,const Buffer<float>& boundary,const Buffer<float>& diagonal,
                               Buffer<float>& out,const Buffer<float>* cell_limiter=nullptr) const=0;
  virtual const Buffer<Index>& mask() const=0;
  virtual std::size_t corrected_faces() const=0;
};

// Reconstruct tangential shifts to the face-normal line, preserving the compact
// two-point normal stencil. Geometry-only quadratic weights; device reductions
// supply diffusion, momentum interpolation and flux-consistent wall gradients.
class ProjectedFaceReconstruction final : public FaceReconstruction {
 public:
  ProjectedFaceReconstruction(Device&,const PolyMesh&,const MeshGeometry&,
                        std::span<const BoundaryObservation>,
                        std::span<const Index> nonsmooth_cells={});
  void gradient(const Buffer<float>& field,const Buffer<float>& boundary,DeviceVec3& out) const;
  void flux_gradient(const Buffer<float>& field,const Buffer<float>& boundary,const DeviceMesh&,const Stencil&,DeviceVec3& out) const;
  void diffusion_correction(const Buffer<float>& field,const Buffer<float>& boundary,Buffer<float>& out) const override;
  void momentum_values(const Buffer<float>& field,const Buffer<float>& boundary,const Buffer<float>& diagonal,
                       Buffer<float>& out,const Buffer<float>* cell_limiter=nullptr) const override;
  const Buffer<Index>& mask() const override {return mask_;}
  std::size_t corrected_faces() const override {return corrected_;}
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

// The same normal-line functionals, fitted with walk-count weights: a sample
// is weighted by its number of walks of length <= 3 (self steps allowed) from
// the owner plus from the neighbour. Those weights make every fit moment a sum
// of per-cell moments, gathered over face neighbours one hop at a time, so the
// storage is ten floats per functional per corrected face and forty per cell,
// rather than one weight per functional for every cloud sample of every face.
class WalkFaceReconstruction final : public FaceReconstruction {
 public:
  static constexpr int hops=3;
  // The device mesh and stencil must outlive this object.
  WalkFaceReconstruction(Device&,const PolyMesh&,const MeshGeometry&,const DeviceMesh&,const Stencil&,
                         std::span<const BoundaryObservation>,std::span<const Index> nonsmooth_cells={});
  void diffusion_correction(const Buffer<float>& field,const Buffer<float>& boundary,Buffer<float>& out) const override;
  void momentum_values(const Buffer<float>& field,const Buffer<float>& boundary,const Buffer<float>& diagonal,
                       Buffer<float>& out,const Buffer<float>* cell_limiter=nullptr) const override;
  const Buffer<Index>& mask() const override {return mask_;}
  std::size_t corrected_faces() const override {return corrected_;}
  std::size_t fallback_faces() const {return fallback_;}
  // True when built from exactly these prescribed boundary values.
  bool observes(std::span<const BoundaryObservation>) const;
 private:
  // Leaves the hop-3 moments in moment_[1], lo_[1] and hi_[1].
  void moments(const Buffer<float>& field,const Buffer<float>& boundary) const;
  Device& dev_;
  const DeviceMesh& mesh_;
  const Stencil& stencil_;
  Index slots_=0;
  std::size_t corrected_=0,fallback_=0;
  std::vector<Index> observed_host_;
  Buffer<Index> observed_,slot_,mask_;
  Buffer<float> geometry_[hops+1],coefficients_;
  mutable Buffer<float> moment_[2],lo_[2],hi_[2];
};

// Reuse shared when it observes exactly these boundary values; otherwise build
// a reconstruction for them. Null when no face needs a correction.
std::shared_ptr<const WalkFaceReconstruction> walk_reconstruction(std::shared_ptr<const WalkFaceReconstruction> shared,
    Device&,const PolyMesh&,const MeshGeometry&,const DeviceMesh&,const Stencil&,std::span<const BoundaryObservation>);
} // namespace kutta::fv
