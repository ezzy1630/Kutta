#pragma once
// Construction mesh shared by the castellate, snap and layer stages. Faces use
// detail::Face (outward from owner); boundary faces carry a patch index:
// [0,6) are the tunnel planes (2*axis+side), 6+zone are STL wall zones.

#include "external_aero.hpp"
#include "kutta/mesh/detail.hpp"
#include <thread>

namespace kutta::mesher::internal {
using detail::Face;

enum CellKind : std::uint8_t { Hex = 0, Prism = 1 };
constexpr int first_wall = 6;

struct Work {
  std::vector<Vec3d> points;
  std::vector<Face> faces;
  Index cells = 0;
  std::vector<CellKind> kind;           // per cell
  std::vector<double> size;             // per cell: octree edge length (layer cells: their wall cell's)
  std::vector<std::uint8_t> planes;     // per point: bit 2*axis+side when on that tunnel plane
  Box bounds;                           // tunnel
};

// Cell -> faces CSR; sign +1 when the cell owns the face.
struct Adjacency {
  std::vector<Index> offsets, faces;
  std::vector<std::int8_t> signs;
};
Adjacency cell_faces(const Work&);
// Point -> faces CSR.
struct PointFaces { std::vector<Index> offsets, faces; };
PointFaces point_faces(const Work&);

struct FaceGeom { Vec3d centroid, area; };
// The same polygon convention as kutta::face_geometry; zero area is returned,
// not thrown, so callers can reject a trial position.
FaceGeom polygon(const std::vector<Vec3d>& points, const std::vector<Index>& nodes);

struct Geometry {
  std::vector<FaceGeom> face;
  std::vector<Vec3d> centroid;
  std::vector<double> volume, quality;  // quality: Fluent orthogonal quality
};
// Full geometry and orthogonal quality (kutta::compute_geometry conventions).
Geometry geometry(const Work&, const Adjacency&);
// Recomputes faces, then cells, then quality of `cells` and their face
// neighbours. `cells` must contain every cell of every face in `faces`.
void update(const Work&, const Adjacency&, Geometry&, const std::vector<Index>& faces, const std::vector<Index>& cells);
// Neighbours of the cells (sharing a face), excluding boundary.
std::vector<Index> with_neighbours(const Work&, const Adjacency&, const std::vector<Index>& cells);

// Moves each source point by target*scale and decays the motion over a few
// point layers into the volume (Laplacian smoothing, tunnel planes slide).
// Wherever a cell falls below min_quality, nearby source scales halve and,
// after a few rounds, drop to zero, so the result is never worse than the
// starting mesh. Scales are updated in place. Cells still bad with every
// nearby scale at zero (an invalid starting mesh) are returned as `stuck`.
struct DisplaceResult { int rounds = 0; Index first_bad = 0; std::vector<Index> stuck; };
DisplaceResult displace(Work&, const Adjacency&, const PointFaces&, const std::vector<Index>& sources,
                        const std::vector<Vec3d>& target, std::vector<double>& scale, double min_quality, const Progress&, std::string_view stage);
// displace for repeated targets on one mesh (layer collapse). The first run
// is a displace call. A later run recomputes only near sources whose target
// changed: their scales, and those of nearby sources, restart at 1, the
// smoothing is updated by the change alone (exact, as it is linear) and only
// cells on moved points are rechecked. Scales elsewhere are kept. The Work,
// adjacency and point faces must outlive it, and the Work's topology must
// not change between runs.
class Displacement {
 public:
  Displacement(Work&, const Adjacency&, const PointFaces&, std::vector<Index> sources, double min_quality);
  DisplaceResult run(const std::vector<Vec3d>& target, const Progress&, std::string_view stage);
  std::vector<double> scale;  // per point; sources only
 private:
  void apply_all();
  std::vector<Index> apply_changed();
  void refresh(const std::vector<Index>& points, std::vector<Index>& checked);
  void grow(std::vector<Index> front, std::vector<Index>& reset);
  DisplaceResult retreat(std::vector<Index>& check, const std::vector<std::uint32_t>* allowed, std::vector<Index>& outside, const Progress&, std::string_view stage);
  Work& w_;
  const Adjacency& adjacency_;
  const PointFaces& pf_;
  std::vector<Index> sources_, band_, moving_faces_, moving_cells_;
  double min_quality_;
  std::vector<int> depth_;
  std::vector<Vec3d> origin_, displacement_, value_, target_;
  std::vector<std::uint32_t> mark_, face_mark_, cell_mark_, reset_mark_;
  std::vector<Index> slot_;
  std::uint32_t stamp_ = 0, run_ = 0;
  Geometry geometry_;
  bool first_ = true, whole_ = true;
};
// Calls fn(q) for each point q joined to p by a face edge (repeats possible).
template <class F>
void point_neighbours(const Work& w, const PointFaces& pf, Index p, F&& fn) {
  for (auto i = pf.offsets[p]; i < pf.offsets[p + 1]; ++i) {
    const auto& ns = w.faces[pf.faces[i]].nodes;
    const auto k = std::size_t(std::find(ns.begin(), ns.end(), p) - ns.begin());
    fn(ns[(k + 1) % ns.size()]);
    fn(ns[(k + ns.size() - 1) % ns.size()]);
  }
}
// Moves wall points onto the STL with implicit feature-edge capture and
// smooths the displacement into the volume, retreating where a cell would
// fall below `min_quality`.
Json::Object snap(Work&, const Surface&, const SurfaceQuery&, const Settings&, const Progress&);
// Inserts prism layers on the selected STL walls by pushing the snapped mesh
// away from the wall; collapses layers locally where quality would fail.
Json::Object add_layers(Work&, const Surface&, const Settings&, const Progress&);

constexpr double snap_min_quality = 0.15;
// Near layered walls, octree cells are at most this many last-layer
// thicknesses (a smooth layer-to-core transition; 50:1 diverged in testing).
constexpr double layer_transition = 4;
// Default body-relative refinement boxes (body_boxes), as multiples of the
// surface max_size.
constexpr double near_box_size = 2, wake_box_size = 4;
// Measured Metal device memory per cell for the default steady solve
// (SIMPLEC, SST) and for Coupled, on this mesher's car meshes
// (validation/results/mesher); not measured on CUDA.
constexpr double simplec_bytes_per_cell = 2100, coupled_bytes_per_cell = 2600;
constexpr double memory_warning_bytes = 10e9;
constexpr double layer_min_quality = 0.15;
// Faces smaller than this fraction of volume^(2/3) count as degenerate.
constexpr double degenerate_face = 1e-4;

Vec3d closest_on_triangle(Vec3d p, Vec3d a, Vec3d b, Vec3d c);

template <class F>
void parallel_for(Index n, F&& f) {
  const auto threads = Index(std::max(1u, std::min(16u, std::thread::hardware_concurrency())));
  if (n < 4096 || threads == 1) { for (Index i = 0; i < n; ++i) f(i); return; }
  std::vector<std::thread> pool;
  const Index chunk = (n + threads - 1) / threads;
  for (Index t = 0; t < threads; ++t)
    pool.emplace_back([&, t] { for (Index i = t * chunk; i < std::min(n, (t + 1) * chunk); ++i) f(i); });
  for (auto& t : pool) t.join();
}
}  // namespace kutta::mesher::internal
