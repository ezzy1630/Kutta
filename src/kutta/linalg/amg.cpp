#include "kutta/linalg/amg.hpp"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <utility>

#include "kutta/linalg/detail.hpp"

namespace kutta {
namespace {
using linalg_detail::norm;

struct Coarsening {
  CsrMatrixHost matrix;
  std::vector<Index> aggregate, member_offsets, members, entry_offsets, entries;
};

// Greedy strongest-unassigned matching, with deterministic index tie-breaking.
// Repeating on the pair graph forms larger aggregates. Strong normal
// coupling in a boundary layer is contracted before weak tangential edges.
std::vector<Index> pairs(const CsrMatrixHost& a, double threshold, bool exclude_isolated) {
  const auto& p = a.pattern;
  // CSR transpose indices allow strength to see incoming as well as outgoing
  // edges, including structurally one-sided upwind operators. No atomics.
  std::vector<Index> offsets(std::size_t(p.n_rows) + 1), sources(p.nnz());
  for (auto col : p.col) ++offsets[col + 1];
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
  auto cursor = offsets;
  for (Index i = 0; i < p.n_rows; ++i)
    for (Index j = p.row_ptr[i]; j < p.row_ptr[i + 1]; ++j)
      sources[cursor[p.col[j]]++] = i;
  // Transposed rows are sorted because source rows were visited in order.
  auto neighbors = [&](Index i, auto visit) {
    Index out = p.row_ptr[i], incoming = offsets[i];
    while (out < p.row_ptr[i + 1] || incoming < offsets[i + 1]) {
      const Index oc = out < p.row_ptr[i + 1] ? p.col[out] : p.n_rows;
      const Index ic = incoming < offsets[i + 1] ? sources[incoming] : p.n_rows;
      const Index col = std::min(oc, ic);
      double weight = 0;
      if (oc == col) weight = std::abs(a.val[out++]);
      if (ic == col) {
        const auto first = p.col.begin() + p.row_ptr[col];
        const auto last = p.col.begin() + p.row_ptr[col + 1];
        const auto entry = std::lower_bound(first, last, i) - p.col.begin();
        weight = std::max(weight, std::abs(a.val[entry]));
        ++incoming;
      }
      if (col != i && weight != 0) visit(col, weight);
    }
  };
  std::vector<Index> map(p.n_rows, -1);
  Index next = 0;
  for (Index i = 0; i < p.n_rows; ++i) {
    if (map[i] != -1) continue;
    double strongest = 0;
    neighbors(i, [&](Index, double weight) { strongest = std::max(strongest, weight); });
    if (exclude_isolated && strongest == 0) continue;
    Index partner = -1;
    double weight = 0;
    neighbors(i, [&](Index col, double w) {
      if (map[col] == -1 && w > weight && w >= threshold * strongest) {
        partner = col;
        weight = w;
      }
    });
    map[i] = next;
    if (partner != -1) map[partner] = next;
    ++next;
  }
  return map;
}

Coarsening coarsen(const CsrMatrixHost& fine, std::vector<Index> aggregate, bool entry_map) {
  Coarsening c;
  c.aggregate = std::move(aggregate);
  const Index nc = c.aggregate.empty() ? 0 : *std::max_element(c.aggregate.begin(), c.aggregate.end()) + 1;
  auto& cp = c.matrix.pattern;
  cp.n_rows = nc;
  cp.row_ptr.assign(std::size_t(nc) + 1, 0);
  cp.diag.resize(nc);
  c.member_offsets.assign(std::size_t(nc) + 1, 0);
  for (auto g : c.aggregate) if (g >= 0) ++c.member_offsets[g + 1];
  std::partial_sum(c.member_offsets.begin(), c.member_offsets.end(), c.member_offsets.begin());
  c.members.resize(c.member_offsets.back());
  auto cursor = c.member_offsets;
  for (Index i = 0; i < fine.pattern.n_rows; ++i)
    if (c.aggregate[i] >= 0) c.members[cursor[c.aggregate[i]]++] = i;

  // Sort only one aggregate row at a time, not a global nnz-sized tuple array.
  // Source indices are tie-breakers, fixing the sum order across all backends.
  std::vector<std::pair<Index, Index>> row;
  if (entry_map) {
    c.entries.reserve(fine.pattern.nnz());
    c.entry_offsets.push_back(0);
  }
  for (Index g = 0; g < nc; ++g) {
    row.clear();
    for (Index m = c.member_offsets[g]; m < c.member_offsets[g + 1]; ++m) {
      Index i = c.members[m];
      for (Index j = fine.pattern.row_ptr[i]; j < fine.pattern.row_ptr[i + 1]; ++j)
        if (c.aggregate[fine.pattern.col[j]] >= 0)
          row.emplace_back(c.aggregate[fine.pattern.col[j]], j);
    }
    std::sort(row.begin(), row.end());
    Index previous = -1;
    for (auto [col, source] : row) {
      if (col != previous) {
        if (entry_map && previous != -1) c.entry_offsets.push_back(static_cast<Index>(c.entries.size()));
        if (col == g) cp.diag[g] = static_cast<Index>(cp.col.size());
        cp.col.push_back(col);
        c.matrix.val.push_back(0);
        previous = col;
      }
      if (entry_map) c.entries.push_back(source);
      c.matrix.val.back() += fine.val[source];
    }
    if (entry_map) c.entry_offsets.push_back(static_cast<Index>(c.entries.size()));
    cp.row_ptr[g + 1] = static_cast<Index>(cp.col.size());
  }
  return c;
}

Coarsening aggregate(const CsrMatrixHost& a, const AmgOptions& options) {
  auto map = pairs(a, options.strength_threshold, true);
  for (int pass = 1; pass < options.aggregation_passes; ++pass) {
    // Pairing only needs the intermediate operator; entry maps are built once
    // for the final transfer. Release it before rebuilding from the fine grid.
    auto intermediate = coarsen(a, std::move(map), false);
    auto next = pairs(intermediate.matrix, options.strength_threshold, false);
    map = std::move(intermediate.aggregate);
    for (auto& g : map) if (g >= 0) g = next[g];
  }
  auto c = coarsen(a, std::move(map), true);
  return c;
}

struct Transfer {
  Buffer<Index> aggregate, member_offsets, members, entry_offsets, entries;
};
struct Level {
  CsrMatrix matrix;  // fine matrix comes from the caller; levels > 0 own CSR
  Buffer<float> inverse, l1, invalid, x, b, r, p, q, z, t;
  Transfer transfer;
  Level(Device& dev, std::size_t n)
      : inverse(dev.alloc<float>(n)), l1(dev.alloc<float>(n)), invalid(dev.alloc<float>(n)),
        x(dev.alloc<float>(n)), b(dev.alloc<float>(n)), r(dev.alloc<float>(n)),
        p(dev.alloc<float>(n)), q(dev.alloc<float>(n)), z(dev.alloc<float>(n)),
        t(dev.alloc<float>(n)) {}
};

class Aggregation final : public linalg_detail::Preconditioner {
 public:
  Aggregation(SolverKind kind, Device& dev, const CsrMatrixHost& a, AmgOptions options)
      : options_(options), spd_(kind == SolverKind::AmgPcg),
        nullspace_(linalg_detail::constant_nullspace(a)) {
    levels_.emplace_back(dev, a.pattern.n_rows);
    const CsrMatrixHost* fine = &a;
    CsrMatrixHost owned;
    std::size_t total_nnz = a.pattern.nnz();
    while (fine->pattern.n_rows > options_.coarse_rows &&
           levels_.size() < std::size_t(options_.max_levels)) {
      auto c = aggregate(*fine, options_);
      // Unconnected/diagonal graphs cannot coarsen. Smoothing remains valid.
      const auto active = c.members.size();
      if (active == 0 || c.matrix.pattern.n_rows >= 0.85 * active) break;
      auto& t = levels_.back().transfer;
      t.aggregate = dev.upload<Index>(c.aggregate);
      t.member_offsets = dev.upload<Index>(c.member_offsets);
      t.members = dev.upload<Index>(c.members);
      t.entry_offsets = dev.upload<Index>(c.entry_offsets);
      t.entries = dev.upload<Index>(c.entries);
      total_nnz += c.matrix.pattern.nnz();
      levels_.emplace_back(dev, c.matrix.pattern.n_rows);
      levels_.back().matrix = upload(dev, c.matrix.pattern);
      owned = std::move(c.matrix);
      fine = &owned;
    }
    complexity_ = a.pattern.nnz() ? double(total_nnz) / a.pattern.nnz() : 0;
  }

  void refresh(Device& dev, const CsrMatrix& a) override {
    const CsrMatrix* matrix = &a;
    for (std::size_t i = 0; i < levels_.size(); ++i) {
      auto& l = levels_[i];
      kern::krylov_diagonal(dev, matrix->n_rows, matrix->row_ptr, matrix->diag, matrix->val,
                           l.inverse, l.l1, l.invalid);
      // The last Neumann operator may be the scalar zero matrix. Its
      // projected RHS is zero and its smoothing correction is zero.
      const double allowed = nullspace_ && matrix->n_rows == 1 ? 1 : 0;
      KUTTA_REQUIRE(dev.max_abs(l.invalid, matrix->n_rows) <= allowed,
                    "AMG: invalid values or zero diagonal");
      if (i + 1 < levels_.size()) {
        auto& coarse = levels_[i + 1].matrix;
        auto& t = l.transfer;
        kern::amg_galerkin(dev, coarse.nnz, t.entry_offsets, t.entries, matrix->val, coarse.val);
        matrix = &coarse;
      }
    }
  }

  void apply(Device& dev, const CsrMatrix& a, const Buffer<float>& r,
             Buffer<float>& z) override {
    cycle(dev, 0, a, r, z);
  }
  bool variable() const override { return options_.cycle == AmgCycle::K; }
  std::string describe() const override {
    std::ostringstream out;
    out << "aggregation AMG " << (options_.cycle == AmgCycle::K ? 'K' : 'V')
        << "-cycle, " << levels_.size() << " levels, op complexity "
        << std::fixed << std::setprecision(3) << complexity_;
    return out.str();
  }

 private:
  void smooth(Device& dev, Level& l, const CsrMatrix& a, const Buffer<float>& b,
              Buffer<float>& x, int sweeps) {
    for (int sweep = 0; sweep < sweeps; ++sweep) {
      residual(dev, a, x, b, l.r);
      kern::krylov_smooth(dev, a.n_rows, l.r, l.l1, x, 2.0f / 3.0f);
    }
  }
  void cycle(Device& dev, std::size_t index, const CsrMatrix& a,
             const Buffer<float>& b, Buffer<float>& x) {
    auto& l = levels_[index];
    dev.fill(x, 0.0f);
    if (index + 1 == levels_.size()) {
      // Entire coarsest solve stays on the device. A fixed sweep count
      // avoids a CPU download/synchronization and retains V-cycle linearity.
      smooth(dev, l, a, b, x, options_.coarse_sweeps);
    } else {
      smooth(dev, l, a, b, x, options_.pre_sweeps);
      residual(dev, a, x, b, l.r);
      auto& c = levels_[index + 1];
      kern::amg_restrict(dev, c.b.size(), l.transfer.member_offsets,
                         l.transfer.members, l.r, c.b);
      if (nullspace_) linalg_detail::project(dev, c.b);
      if (options_.cycle == AmgCycle::K && index + 2 < levels_.size())
        inner(dev, index + 1);
      else cycle(dev, index + 1, c.matrix, c.b, c.x);
      kern::amg_prolong(dev, a.n_rows, l.transfer.aggregate, c.x, x);
      smooth(dev, l, a, b, x, options_.post_sweeps);
    }
    if (nullspace_) linalg_detail::project(dev, x);
  }

  // Two recursive preconditioned Krylov steps on a coarse equation. Its
  // workspaces belong to that level and are distinct from cycle's residual.
  // SPD uses A-conjugate directions, nonsymmetric uses GCR orthogonalization.
  void inner(Device& dev, std::size_t index) {
    auto& l = levels_[index];
    const auto& a = l.matrix;
    const auto n = a.n_rows;
    cycle(dev, index, a, l.b, l.z);
    dev.copy(l.p, l.z);
    spmv(dev, a, l.p, l.q);
    double denom = spd_ ? dev.dot(l.p, l.q, n) : dev.dot(l.q, l.q, n);
    double numer = spd_ ? dev.dot(l.b, l.p, n) : dev.dot(l.b, l.q, n);
    if (!(denom > 0) || !linalg_detail::coefficient(numer / denom)) {
      dev.copy(l.x, l.z);
      return;
    }
    double alpha = numer / denom;
    kern::scale_copy(dev, n, l.x, l.p, static_cast<float>(alpha));
    kern::krylov_difference(dev, n, l.b, l.q, l.t, static_cast<float>(alpha));
    // This adaptive exit makes the K-cycle nonlinear; the outer PCG uses
    // flexible CG rather than assuming a stationary preconditioner.
    if (norm(dev, l.t) <= 0.25 * norm(dev, l.b)) return;
    cycle(dev, index, a, l.t, l.z);
    double beta;
    if (spd_) beta = -dev.dot(l.z, l.q, n) / denom;
    else {
      spmv(dev, a, l.z, l.r);
      beta = -dev.dot(l.r, l.q, n) / denom;
    }
    if (!linalg_detail::coefficient(beta)) return;
    kern::krylov_direction(dev, n, l.z, l.p, static_cast<float>(beta));
    spmv(dev, a, l.p, l.q);
    denom = spd_ ? dev.dot(l.p, l.q, n) : dev.dot(l.q, l.q, n);
    numer = spd_ ? dev.dot(l.t, l.p, n) : dev.dot(l.t, l.q, n);
    if (denom > 0 && linalg_detail::coefficient(numer / denom))
      kern::axpy(dev, n, l.x, l.p, static_cast<float>(numer / denom));
  }

  AmgOptions options_;
  bool spd_, nullspace_;
  double complexity_ = 0;
  std::vector<Level> levels_;
};
}  // namespace

std::unique_ptr<LinearSolver> make_amg_solver(SolverKind kind, Device& dev,
    const CsrMatrixHost& a, const AmgOptions& options) {
  KUTTA_REQUIRE(kind == SolverKind::AmgPcg || kind == SolverKind::AmgBicgstab,
                "AMG: expected an AMG solver kind");
  KUTTA_REQUIRE(options.pre_sweeps > 0 && options.post_sweeps > 0 &&
                    options.coarse_sweeps > 0 && options.coarse_rows >= 1 &&
                    options.max_levels >= 1 && options.aggregation_passes >= 1 &&
                    options.aggregation_passes <= 8 &&
                    std::isfinite(options.strength_threshold) &&
                    options.strength_threshold >= 0 && options.strength_threshold <= 1,
                "AMG: invalid options");
  KUTTA_REQUIRE(kind != SolverKind::AmgPcg || options.pre_sweeps == options.post_sweeps,
                "AMG-PCG requires equal pre/post sweeps");
  // Validate the public host CSR boundary before traversing any arrays.
  auto checked = upload(dev, a);
  return linalg_detail::krylov_solver(kind, dev, a,
      std::make_unique<Aggregation>(kind, dev, a, options), checked);
}

std::unique_ptr<LinearSolver> make_solver(SolverKind kind, Device& dev, const CsrMatrixHost& a) {
  switch (kind) {
    case SolverKind::AmgPcg:
    case SolverKind::AmgBicgstab: return make_amg_solver(kind, dev, a, AmgOptions{});
    case SolverKind::Jacobi:
    case SolverKind::PcgJacobi:
    case SolverKind::BicgstabJacobi: {
      auto initial = upload(dev, a);
      return linalg_detail::krylov_solver(kind, dev, a, nullptr, initial);
    }
  }
  throw Error("unknown SolverKind");
}
}  // namespace kutta
