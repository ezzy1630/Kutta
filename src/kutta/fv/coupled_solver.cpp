#include "kutta/fv/coupled_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>

#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {
using namespace detail;
namespace {
// Symmetric weighted cell graph without self loops.
struct Graph {
  std::vector<Index> row_ptr{0}, col;
  std::vector<double> weight;
  Index rows() const { return Index(row_ptr.size()) - 1; }
};

// Greedy strongest-unassigned pairing with index tie-breaking (deterministic).
// A row only pairs through edges of at least threshold * its strongest edge.
std::vector<Index> pair_rows(const Graph& g, double threshold, Index& count) {
  std::vector<Index> map(g.rows(), -1);
  count = 0;
  for (Index i = 0; i < g.rows(); ++i) {
    if (map[i] != -1) continue;
    double strongest = 0;
    for (Index j = g.row_ptr[i]; j < g.row_ptr[i + 1]; ++j) strongest = std::max(strongest, g.weight[j]);
    Index partner = -1;
    double best = 0;
    for (Index j = g.row_ptr[i]; j < g.row_ptr[i + 1]; ++j)
      if (map[g.col[j]] == -1 && g.weight[j] > best && g.weight[j] >= threshold * strongest) {
        partner = g.col[j];
        best = g.weight[j];
      }
    map[i] = count;
    if (partner != -1) map[partner] = count;
    ++count;
  }
  return map;
}

// Member lists by counting sort: no per-row containers on large meshes.
void members_of(const std::vector<Index>& map, Index count, std::vector<Index>& offsets,
                std::vector<Index>& members) {
  offsets.assign(std::size_t(count) + 1, 0);
  for (auto g : map) ++offsets[g + 1];
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
  members.resize(map.size());
  auto cursor = offsets;
  for (Index i = 0; i < Index(map.size()); ++i) members[cursor[map[i]]++] = i;
}

Graph contract(const Graph& g, const std::vector<Index>& map, Index count) {
  std::vector<Index> offsets, members;
  members_of(map, count, offsets, members);
  Graph out;
  out.col.reserve(g.col.size());
  out.weight.reserve(g.col.size());
  std::vector<std::pair<Index, double>> row;
  for (Index c = 0; c < count; ++c) {
    row.clear();
    for (Index m = offsets[c]; m < offsets[c + 1]; ++m) {
      const auto i = members[m];
      for (Index j = g.row_ptr[i]; j < g.row_ptr[i + 1]; ++j)
        if (map[g.col[j]] != c) row.push_back({map[g.col[j]], g.weight[j]});
    }
    std::sort(row.begin(), row.end());
    for (std::size_t k = 0; k < row.size();) {
      auto end = k;
      double sum = 0;
      while (end < row.size() && row[end].first == row[k].first) sum += row[end++].second;
      out.col.push_back(row[k].first);
      out.weight.push_back(sum);
      k = end;
    }
    out.row_ptr.push_back(Index(out.col.size()));
  }
  return out;
}

struct Transfer {
  Buffer<Index> aggregate, member_offsets, members;
};
struct Level {
  Index n = 0, nnz = 0;
  Buffer<Index> rows, cols, diag;
  Buffer<float> a, c, g, h, d;  // g, h: 3*nnz, d: 4n (see fv_coupled.kern.h)
  Buffer<float> x, b, r, tmp;
};
}  // namespace

struct CoupledSolver::Impl {
  Device& dev;
  const DeviceMesh& mesh;
  const Stencil& stencil;
  const Buffer<Index>& boundary_kind;
  const BoundaryLeastSquares& lsq;
  CoupledOptions options;
  Index n = 0, reference = -1;
  // Fine level (matrix-free).
  Buffer<float> cell, coefficient, fine_r, fine_tmp;  // cell: [du|dv|dw|cd|G|H|Q|cs], 20n
  std::vector<Transfer> transfers;  // transfers[k]: level k -> k+1
  std::vector<Level> levels;        // coarse levels 1..L, stored at [0..L-1]
  // Krylov workspace over 4n.
  Buffer<float> x, b, r, shadow, p, q, s, t, z, zs, scalars, gauge;
  Buffer<Index> status;
  // Per-solve state shared by the operator kernels.
  const Buffer<float>* mass = nullptr;
  const Buffer<float>* viscosity = nullptr;
  const Buffer<float>* weight = nullptr;
  float wc = 1;
  std::size_t owned_bytes = 0;

  template <class T>
  Buffer<T> alloc(std::size_t count) {
    owned_bytes += count * sizeof(T);
    return dev.alloc<T>(count);
  }
  template <class T>
  Buffer<T> upload(const std::vector<T>& values) {
    owned_bytes += values.size() * sizeof(T);
    return dev.upload<T>(values);
  }

  Impl(Device& device, const DeviceMesh& m, const Stencil& st, const PolyMesh& host,
       const MeshGeometry& geom, const Buffer<Index>& kind, const BoundaryLeastSquares& pressure_lsq,
       Index ref, std::span<const Index> low_order, CoupledOptions opt)
      : dev(device), mesh(m), stencil(st), boundary_kind(kind), lsq(pressure_lsq), options(opt),
        n(m.n_cells), reference(ref) {
    KUTTA_REQUIRE(options.pre_sweeps >= 1 && options.post_sweeps >= 0 && options.coarsest_sweeps >= 1 &&
                      options.smoother_damping > 0 && options.smoother_damping <= 1 && options.courant > 0 &&
                      options.coarsest_rows >= 1 && options.aggregation_passes >= 1 &&
                      options.aggregation_passes <= 8,
                  "coupled: invalid options");
    const auto kinds = dev.download(kind);
    const auto nf = m.n_faces, ni = m.n_internal_faces;
    // Pressure-force sensitivities (fv_coupled.kern.h), unscaled: G takes
    // half of each internal face and the whole face where the boundary
    // pressure is extrapolated from the cell. Q uses the residual's own
    // extrapolation vectors, none in low-order cells. Incidences follow the device
    // cell_face_offsets/cell_faces order of the LSQ deltas. The flux block H
    // follows the momentum face weights and is refreshed every solve.
    const auto offsets = dev.download(m.cell_face_offsets), faces = dev.download(m.cell_faces);
    const auto signs = dev.download(m.cell_face_sign);
    // Internal rows are per entry (possibly the stencil's); boundary rows are by boundary face.
    const auto& rows = lsq.rows();
    const auto dx = dev.download(rows.x), dy = dev.download(rows.y), dz = dev.download(rows.z);
    const auto bx = dev.download(lsq.boundary_delta.x), by = dev.download(lsq.boundary_delta.y),
               bz = dev.download(lsq.boundary_delta.z);
    const auto inverse = dev.download(lsq.inverse);
    const auto download3 = [&](const DeviceVec3& v) {
      return std::array{dev.download(v.x), dev.download(v.y), dev.download(v.z)};
    };
    const auto owner_r = download3(m.owner_to_face), neighbour_r = download3(m.neighbour_to_face);
    std::vector<float> packed(20 * std::size_t(n), 0);
    const auto put = [&](int slot, Index c, double value) {packed[std::size_t(slot) * n + c] = float(value);};
    for (Index c = 0; c < n; ++c) {
      Vec3d g, delta;
      std::array<double, 9> t{};  // sum S (x_f - x_c)^T, row major
      for (Index j = offsets[c]; j < offsets[c + 1]; ++j) {
        const auto f = faces[j];
        const auto area = geom.face_area[f] * double(signs[j]);
        const auto& side = f < ni && host.owner[f] != c ? neighbour_r : owner_r;
        const bool extrapolated = low_order.empty() || low_order[c] == 0;
        const Vec3d r = extrapolated ? Vec3d{side[0][f], side[1][f], side[2][f]} : Vec3d{};
        const bool prescribed = f >= ni && (kinds[f - ni] == 2 || kinds[f - ni] == 4);
        // An unprescribed boundary row observes the cell's own pressure, so
        // its difference, and its share of d g_c / d p_c, is zero.
        if (f < ni) delta = delta + Vec3d{dx[j], dy[j], dz[j]};
        else if (prescribed) delta = delta + Vec3d{bx[f - ni], by[f - ni], bz[f - ni]};
        double weight = 0;
        if (f < ni) {g = g + area * .5; weight = .5;}
        else if (!prescribed) {g = g + area; weight = 1;}
        const double s[]{area.x, area.y, area.z}, q[]{r.x, r.y, r.z};
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) t[3 * a + b] += weight * s[a] * q[b];
      }
      const double inv[9]{inverse[c], inverse[c + n], inverse[c + 2 * n], inverse[c + n], inverse[c + 3 * n],
                          inverse[c + 4 * n], inverse[c + 2 * n], inverse[c + 4 * n], inverse[c + 5 * n]};
      std::array<double, 9> q{};
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b)
        for (int k = 0; k < 3; ++k) q[3 * a + b] += t[3 * a + k] * inv[3 * k + b];
      // d g_c / d p_c = -M^-1 sum_j d_j, so the own-cell term enters G_cc.
      g = g - Vec3d{q[0] * delta.x + q[1] * delta.y + q[2] * delta.z, q[3] * delta.x + q[4] * delta.y + q[5] * delta.z,
                    q[6] * delta.x + q[7] * delta.y + q[8] * delta.z};
      put(4, c, g.x); put(5, c, g.y); put(6, c, g.z);
      for (int k = 0; k < 9; ++k) put(10 + k, c, q[k]);
    }
    cell = upload(packed);
    coefficient = alloc<float>(nf);
    fine_r = alloc<float>(4 * std::size_t(n)); fine_tmp = alloc<float>(4 * std::size_t(n));
    for (auto* v : {&x, &b, &r, &shadow, &p, &q, &s, &t, &z, &zs}) *v = alloc<float>(4 * std::size_t(n));
    scalars = alloc<float>(24); gauge = alloc<float>(2); status = alloc<Index>(3);
    dev.fill(x, 0);

    // Aggregates follow the geometric pressure Laplacian |S|^2/|S.d|, the
    // same strength the segregated pressure hierarchy starts from.
    Graph graph;
    {
      // The stencil's sorted cell-adjacency CSR, without its diagonal.
      const auto& pattern = st.pattern;
      std::vector<double> weight(pattern.col.size(), 0);
      const auto entry = [&](Index row, Index col) {
        const auto first = pattern.col.begin() + pattern.row_ptr[row];
        return std::lower_bound(first, pattern.col.begin() + pattern.row_ptr[row + 1], col) - pattern.col.begin();
      };
      for (Index f = 0; f < ni; ++f) {
        const auto o = host.owner[f], k = host.neighbour[f];
        const auto area = geom.face_area[f];
        const auto d = geom.cell_centroid[k] - geom.cell_centroid[o];
        const double w = dot(area, area) / std::max(std::abs(dot(area, d)), 1e-300);
        weight[entry(o, k)] += w;
        weight[entry(k, o)] += w;
      }
      graph.col.reserve(pattern.col.size());
      graph.weight.reserve(pattern.col.size());
      for (Index i = 0; i < n; ++i) {
        for (Index j = pattern.row_ptr[i]; j < pattern.row_ptr[i + 1]; ++j)
          if (j != pattern.diag[i]) {graph.col.push_back(pattern.col[j]); graph.weight.push_back(weight[j]);}
        graph.row_ptr.push_back(Index(graph.col.size()));
      }
    }
    Index rows_now = n;
    while (rows_now > options.coarsest_rows && levels.size() < 24) {
      std::vector<Index> map(rows_now);
      std::iota(map.begin(), map.end(), 0);
      Graph coarse = graph;
      Index count = rows_now;
      for (int pass = 0; pass < options.aggregation_passes && count > options.coarsest_rows; ++pass) {
        Index next_count = 0;
        const auto next = pair_rows(coarse, .25, next_count);
        coarse = contract(coarse, next, next_count);
        for (auto& v : map) v = next[v];
        count = next_count;
      }
      if (count >= .85 * rows_now) break;  // disconnected or stalled graph
      Transfer transfer;
      std::vector<Index> offsets, members;
      members_of(map, count, offsets, members);
      transfer.aggregate = upload(map);
      transfer.member_offsets = upload(offsets);
      transfer.members = upload(members);
      transfers.push_back(std::move(transfer));
      // Coarse CSR pattern: the contracted graph plus the diagonal.
      std::vector<Index> row_ptr{0}, cols, diag(count);
      for (Index i = 0; i < count; ++i) {
        bool placed = false;
        for (Index j = coarse.row_ptr[i]; j <= coarse.row_ptr[i + 1]; ++j) {
          const bool end = j == coarse.row_ptr[i + 1];
          if (!placed && (end || coarse.col[j] > i)) {diag[i] = Index(cols.size()); cols.push_back(i); placed = true;}
          if (!end) cols.push_back(coarse.col[j]);
        }
        row_ptr.push_back(Index(cols.size()));
      }
      Level level;
      level.n = count;
      level.nnz = Index(cols.size());
      level.rows = upload(row_ptr); level.cols = upload(cols); level.diag = upload(diag);
      level.a = alloc<float>(cols.size()); level.c = alloc<float>(cols.size());
      level.g = alloc<float>(3 * cols.size()); level.h = alloc<float>(3 * cols.size());
      level.d = alloc<float>(4 * std::size_t(count));
      for (auto* v : {&level.x, &level.b, &level.r, &level.tmp}) *v = alloc<float>(4 * std::size_t(count));
      levels.push_back(std::move(level));
      graph = std::move(coarse);
      rows_now = count;
    }
  }

  // ---- operators -------------------------------------------------------
  void fine_spmv(const Buffer<float>& in, Buffer<float>& out) {
    const auto& m = mesh;
    kern::coupled_fine_spmv(dev, n, m.cell_face_offsets, m.cell_faces, m.cell_face_sign, stencil.adjacent,
        m.area.x, m.area.y, m.area.z, stencil.diffusion, *viscosity, *mass, coefficient, *weight,
        lsq.rows().x, lsq.rows().y, lsq.rows().z, cell, in, out, n, wc);
  }
  void fine_residual(const Buffer<float>& in, const Buffer<float>& rhs, Buffer<float>& out) {
    const auto& m = mesh;
    kern::coupled_fine_residual(dev, n, m.cell_face_offsets, m.cell_faces, m.cell_face_sign, stencil.adjacent,
        m.area.x, m.area.y, m.area.z, stencil.diffusion, *viscosity, *mass, coefficient, *weight,
        lsq.rows().x, lsq.rows().y, lsq.rows().z, cell, in, rhs, out, n, wc);
  }
  void fine_smooth(Buffer<float>& out, const Buffer<float>& rhs, int sweeps, bool first) {
    const auto& m = mesh;
    for (int k = 0; k < sweeps; ++k) {
      kern::coupled_fine_smooth(dev, n, m.cell_face_offsets, m.cell_faces, m.cell_face_sign, stencil.adjacent,
          m.area.x, m.area.y, m.area.z, stencil.diffusion, *viscosity, *mass, coefficient, *weight,
          lsq.rows().x, lsq.rows().y, lsq.rows().z, cell, out, rhs, fine_tmp, n, wc, options.smoother_damping,
          first && k == 0 ? 1 : 0);
      std::swap(out, fine_tmp);
    }
  }
  void coarse_smooth(Level& l, int sweeps, bool first) {
    for (int k = 0; k < sweeps; ++k) {
      kern::coupled_coarse_smooth(dev, l.n, l.rows, l.cols, l.diag, l.a, l.c, l.g, l.h, l.d,
          l.x, l.b, l.tmp, l.n, l.nnz, options.smoother_damping, first && k == 0 ? 1 : 0);
      std::swap(l.x, l.tmp);
    }
  }
  // Coarse level k (levels[k]) receives its RHS in levels[k].b.
  void cycle(std::size_t k) {
    auto& l = levels[k];
    if (k + 1 == levels.size()) {
      coarse_smooth(l, options.coarsest_sweeps, true);
      return;
    }
    auto& c = levels[k + 1];
    coarse_smooth(l, options.pre_sweeps, true);
    kern::coupled_coarse_residual(dev, l.n, l.rows, l.cols, l.diag, l.a, l.c, l.g, l.h, l.d,
        l.x, l.b, l.r, l.n, l.nnz);
    kern::coupled_restrict(dev, c.n, transfers[k + 1].member_offsets, transfers[k + 1].members, l.r, c.b, l.n, c.n);
    cycle(k + 1);
    kern::coupled_prolong(dev, l.n, transfers[k + 1].aggregate, c.x, l.x, l.n, c.n);
    coarse_smooth(l, options.post_sweeps, false);
  }
  // z = V-cycle(rhs), a fixed linear preconditioner.
  void precondition(const Buffer<float>& rhs, Buffer<float>& out) {
    fine_smooth(out, rhs, options.pre_sweeps, true);
    if (!levels.empty()) {
      auto& c = levels[0];
      fine_residual(out, rhs, fine_r);
      kern::coupled_restrict(dev, c.n, transfers[0].member_offsets, transfers[0].members, fine_r, c.b, n, c.n);
      cycle(0);
      kern::coupled_prolong(dev, n, transfers[0].aggregate, c.x, out, n, c.n);
    }
    fine_smooth(out, rhs, options.post_sweeps, false);
  }
  void galerkin() {
    if (levels.empty()) return;
    const auto& m = mesh;
    auto& c = levels[0];
    const auto& t = transfers[0];
    kern::coupled_galerkin_fine(dev, c.n, t.member_offsets, t.members, t.aggregate,
        m.cell_face_offsets, m.cell_faces, m.cell_face_sign, stencil.adjacent, m.area.x, m.area.y, m.area.z,
        stencil.diffusion, *viscosity, *mass, coefficient, *weight, lsq.rows().x, lsq.rows().y, lsq.rows().z, cell,
        c.rows, c.cols, c.a, c.c, c.g, c.h, c.d, n, c.n, c.nnz, wc);
    for (std::size_t k = 1; k < levels.size(); ++k) {
      auto& f = levels[k - 1];
      auto& l = levels[k];
      const auto& tr = transfers[k];
      kern::coupled_galerkin_coarse(dev, l.n, tr.member_offsets, tr.members, tr.aggregate,
          f.rows, f.cols, f.diag, f.a, f.c, f.g, f.h, f.d,
          l.rows, l.cols, l.diag, l.a, l.c, l.g, l.h, l.d, f.n, f.nnz, l.n, l.nnz);
    }
    // Smoother Schur bounds; tmp is free until the first sweep.
    for (auto& l : levels) {
      kern::coupled_coarse_gsum(dev, l.n, l.rows, l.g, l.tmp, l.n, l.nnz);
      kern::coupled_coarse_schur(dev, l.n, l.rows, l.cols, l.diag, l.h, l.tmp, l.d, l.n, l.nnz);
    }
  }

  double norm(const Buffer<float>& v) { return std::sqrt(dev.dot(v, v, v.size())); }

  // Right-preconditioned BiCGStab on the 4n block vector, reusing the device
  // scalar kernels of the scalar solvers (layout documented in scalars.kern.h).
  SolveStats bicgstab(const SolveControl& control) {
    const auto n4 = 4 * std::size_t(n);
    SolveStats stats;
    dev.fill(x, 0);
    dev.copy(r, b);
    stats.initial_residual = stats.final_residual = norm(r);
    KUTTA_REQUIRE(std::isfinite(stats.initial_residual), "coupled: nonfinite residual");
    const double target = std::max(control.abs_tol, control.rel_tol * stats.initial_residual);
    if (stats.initial_residual <= target) {stats.converged = true; return stats;}
    std::array<float, 24> values{};
    const auto pair = [&](int offset, double value) {
      values[offset] = float(value);
      values[offset + 1] = std::isfinite(values[offset]) ? float(value - double(values[offset])) : 0;
    };
    pair(20, target * target);
    pair(22, 1e6 * stats.initial_residual * stats.initial_residual);
    dev.upload<float>(scalars, values);
    const auto restart = [&] {
      kern::solver_reset(dev, 1, scalars, status, stats.iterations);
      dev.copy(shadow, r);
      dev.fill(p, 0);
      dev.fill(q, 0);
    };
    restart();
    int restarts = 0;
    while (stats.iterations < control.max_iter) {
      dev.dot_into(shadow, r, n4, scalars, 0);
      kern::solver_bicg_beta(dev, 1, scalars, status);
      kern::solver_bicg_direction(dev, n4, r, q, p, scalars, status);
      precondition(p, z);
      fine_spmv(z, q);
      dev.dot_into(shadow, q, n4, scalars, 2);
      kern::solver_bicg_alpha(dev, 1, scalars, status);
      kern::solver_bicg_difference(dev, n4, r, q, s, scalars, status);
      dev.dot_into(s, s, n4, scalars, 16);
      kern::solver_bicg_s(dev, 1, scalars, status);
      precondition(s, zs);
      fine_spmv(zs, t);
      dev.dot2_into(t, t, t, s, n4, scalars, 2);
      kern::solver_bicg_omega(dev, 1, scalars, status);
      kern::solver_bicg_step(dev, n4, x, r, z, zs, s, t, scalars, status);
      dev.dot_into(r, r, n4, scalars, 18);
      kern::solver_check(dev, 1, scalars, status, 0);
      const auto state = dev.download(status);
      stats.iterations = state[1];
      if (state[0] == 0) continue;
      // Verify a candidate against the true residual b - J x.
      fine_residual(x, b, r);
      stats.final_residual = norm(r);
      stats.converged = std::isfinite(stats.final_residual) && stats.final_residual <= target;
      if (stats.converged || state[0] == 2 || !std::isfinite(stats.final_residual) ||
          stats.final_residual > 1e3 * stats.initial_residual || ++restarts > 4)
        return bounded(stats, target);
      restart();
    }
    fine_residual(x, b, r);
    stats.final_residual = norm(r);
    stats.converged = stats.final_residual <= target;
    return bounded(stats, target);
  }
  // A Krylov breakdown must not reach the flow field: if the iterate grew
  // the residual, take one multigrid cycle instead (a bounded step).
  SolveStats bounded(SolveStats stats, double target) {
    if (std::isfinite(stats.final_residual) && stats.final_residual <= stats.initial_residual) return stats;
    precondition(b, x);
    fine_residual(x, b, r);
    stats.final_residual = norm(r);
    KUTTA_REQUIRE(std::isfinite(stats.final_residual), "coupled: nonfinite update");
    stats.converged = stats.final_residual <= target;
    return stats;
  }
};

CoupledSolver::CoupledSolver(Device& dev, const DeviceMesh& mesh, const Stencil& stencil,
                             const PolyMesh& host, const MeshGeometry& geometry,
                             const Buffer<Index>& boundary_kind, const BoundaryLeastSquares& pressure_lsq,
                             Index reference, std::span<const Index> low_order_pressure, CoupledOptions options)
    : impl_(std::make_unique<Impl>(dev, mesh, stencil, host, geometry, boundary_kind, pressure_lsq,
                                   reference, low_order_pressure, options)) {}
CoupledSolver::~CoupledSolver() = default;

void CoupledSolver::capture_diagonal(int component, const CsrMatrix& momentum, const Buffer<float>& physical) {
  KUTTA_REQUIRE(component >= 0 && component < 3, "coupled: invalid momentum component");
  auto& s = *impl_;
  kern::coupled_copy_diagonal(s.dev, s.n, momentum.row_ptr, momentum.diag, momentum.val, physical, s.cell,
                              component * s.n, float(1 / s.options.courant));
}

SolveStats CoupledSolver::solve(const CoupledResiduals& in, const SolveControl& control) {
  auto& s = *impl_;
  const auto& m = s.mesh;
  KUTTA_REQUIRE(in.density > 0 && in.velocity_scale > 0 && std::isfinite(in.velocity_scale),
                "coupled: invalid scales");
  s.mass = &in.mass;
  s.viscosity = &in.viscosity;
  s.weight = &in.face_weights;
  // Continuity rows in momentum units: volume flux times rho*U.
  const double wc = in.density * in.velocity_scale;
  s.wc = float(wc);
  kern::coupled_face_coefficient(s.dev, m.n_faces, m.owner, m.neighbour, s.boundary_kind,
      s.stencil.diffusion, m.volume, in.physical_diagonal, s.coefficient, m.n_internal_faces, float(wc * wc));
  kern::coupled_velocity_block(s.dev, s.n, m.cell_face_offsets, m.cell_faces, m.cell_face_sign, s.stencil.adjacent,
      s.boundary_kind, m.area.x, m.area.y, m.area.z, in.face_weights, s.cell, s.n, m.n_internal_faces);
  kern::coupled_rhs(s.dev, s.n, m.cell_face_offsets, m.cell_faces, m.cell_face_sign, s.coefficient, in.mass,
      in.net_outflow ? *in.net_outflow : in.mass, in.u, in.v, in.w, s.cell, s.b, s.n, s.reference,
      float(wc / in.density), in.net_outflow ? 1 : 0);
  // Smoother Schur bound; fine_tmp is free until the first sweep.
  kern::coupled_fine_gsum(s.dev, s.n, m.cell_face_offsets, s.stencil.adjacent, m.cell_faces, m.cell_face_sign,
      m.area.x, m.area.y, m.area.z, s.lsq.rows().x, s.lsq.rows().y, s.lsq.rows().z, s.cell, s.fine_tmp, s.n, float(wc));
  kern::coupled_fine_schur(s.dev, s.n, m.cell_face_offsets, s.stencil.adjacent, m.cell_faces,
      m.area.x, m.area.y, m.area.z, in.face_weights, s.fine_tmp, s.cell, s.n, float(wc));
  s.galerkin();
  return s.bicgstab(control);
}

void CoupledSolver::update(DeviceVec3& velocity, Buffer<float>& pressure, double alpha_u, double alpha_p) {
  auto& s = *impl_;
  kern::coupled_update(s.dev, s.n, s.x, velocity.x, velocity.y, velocity.z, pressure, s.n,
                       float(alpha_u), float(alpha_p * s.wc));
  if (s.reference >= 0) {
    kern::coupled_pick(s.dev, 1, pressure, s.gauge, s.reference);
    kern::coupled_gauge(s.dev, s.n, s.gauge, pressure);
  }
}

void CoupledSolver::update(DeviceVec3& velocity, Buffer<float>& pressure, const CoupledPairs& pairs,
                           double alpha_u, double alpha_p) {
  auto& s = *impl_;
  auto& lo = pairs.velocity_lo;
  kern::coupled_paired_update(s.dev, s.n, s.x, velocity.x, velocity.y, velocity.z, pressure, lo.x, lo.y, lo.z,
                              pairs.pressure_lo, s.n, float(alpha_u), float(alpha_p * s.wc), pairs.velocity_slack,
                              pairs.pressure_slack);
  if (s.reference >= 0) {
    kern::coupled_paired_pick(s.dev, 1, pressure, pairs.pressure_lo, s.gauge, s.reference);
    kern::coupled_paired_gauge(s.dev, s.n, s.gauge, pressure, pairs.pressure_lo, pairs.pressure_slack);
  }
}

std::size_t CoupledSolver::bytes() const { return impl_->owned_bytes; }

std::string CoupledSolver::describe() const {
  std::ostringstream out;
  out << "coupled block BiCGStab, aggregation AMG V-cycle, " << impl_->levels.size() + 1 << " levels (";
  out << impl_->n;
  for (const auto& l : impl_->levels) out << '/' << l.n;
  out << " rows), block Jacobi smoother";
  return out.str();
}

}  // namespace kutta::fv
