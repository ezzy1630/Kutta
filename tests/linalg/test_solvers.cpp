#include <iostream>
#include <limits>
#include <map>
#include "kernels.hpp"

#include "../backend_helpers.hpp"
#include "kutta/linalg/amg.hpp"
#include "matrices.hpp"
#include "instrumented_device.hpp"

using namespace kutta;
using namespace kutta::test_linalg;

namespace {
void check_solution(Device& dev, const CsrMatrixHost& a, SolverKind kind,
                    double tolerance = 2e-6, bool neumann = false, AmgCycle cycle = AmgCycle::V) {
  auto matrix = upload(dev, a);
  auto exact = random_solution(a.pattern.n_rows);
  auto b = rhs(a, exact);
  if (neumann) {
    // Integer-valued compatible forcing avoids introducing a mean by float conversion.
    std::fill(b.begin(), b.end(), 0.0f); b.front() = 1; b.back() = -1;
  }
  auto reference = kind == SolverKind::BicgstabJacobi || kind == SolverKind::AmgBicgstab ?
      reference_dense(a, b) : reference_cg(a, b, neumann);
  AmgOptions options; options.cycle = cycle;
  options.nullspace = neumann ? NullSpace::Constant : NullSpace::None;
  auto solver = kind == SolverKind::AmgPcg || kind == SolverKind::AmgBicgstab ?
      make_amg_solver(kind, dev, a, options) :
      make_solver(kind, dev, a, {.nullspace = options.nullspace});
  auto x = dev.alloc<float>(a.pattern.n_rows);
  dev.fill(x, 0.0f);
  auto db = dev.upload<float>(b);
  auto stats = solver->solve(dev, matrix, x, db, {tolerance, 1e-12, 2000});
  INFO(solver->describe() << ", iterations " << stats.iterations << ", rel "
       << stats.final_residual / stats.initial_residual);
  CHECK(stats.converged);
  auto got = dev.download(x);
  double error = 0, scale = 0;
  for (std::size_t i = 0; i < got.size(); ++i) {
    error += std::pow(got[i] - reference[i], 2); scale += reference[i] * reference[i];
  }
  CHECK(std::sqrt(error / scale) < 2e-4);
  auto r = dev.alloc<float>(x.size());
  residual(dev, matrix, x, db, r);
  CHECK(stats.final_residual == doctest::Approx(std::sqrt(dev.dot(r, r, r.size()))).epsilon(1e-6));
}
}

TEST_CASE("CSR validates its boundary and kernels match host double arithmetic") {
  auto a = poisson(5, 4, 3, 2, true);
  auto exact = random_solution(a.pattern.n_rows);
  std::vector<float> baseline;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto matrix = upload(*dev, a);
    auto x = dev->upload<float>(exact);
    auto y = dev->alloc<float>(exact.size());
    spmv(*dev, matrix, x, y);
    auto got = dev->download(y);
    auto expect = rhs(a, exact);
    for (std::size_t i = 0; i < got.size(); ++i)
      CHECK(got[i] == doctest::Approx(expect[i]).epsilon(2e-6));
    if (baseline.empty()) baseline = got;
    else for (std::size_t i = 0; i < got.size(); ++i)
      CHECK(got[i] == doctest::Approx(baseline[i]).epsilon(2e-6));
    auto b = dev->upload<float>(expect);
    residual(*dev, matrix, x, b, y);
    CHECK(dev->max_abs(y, y.size()) < 1e-5);
    CHECK_THROWS_AS(spmv(*dev, matrix, x, x), Error);
    auto bad = a; bad.pattern.diag[0] = -1;
    CHECK_THROWS_AS(upload(*dev, bad), Error);
    bad = a; bad.pattern.col[0] = a.pattern.n_rows;
    CHECK_THROWS_AS(upload(*dev, bad), Error);
    bad = a; bad.val[0] = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(upload(*dev, bad), Error);
    CHECK_THROWS_AS(make_solver(SolverKind::AmgPcg, *dev, bad), Error);
  }
}

TEST_CASE("Jacobi PCG and AMG solve SPD problems against double CG") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    for (auto kind : {SolverKind::Jacobi, SolverKind::PcgJacobi, SolverKind::AmgPcg}) {
      check_solution(*dev, poisson(8, 7), kind);
      check_solution(*dev, poisson(6, 5, 4), kind);
      check_solution(*dev, poisson(8, 7, 1, 1, true), kind);
    }
    check_solution(*dev, poisson(12, 9, 5, 1000), SolverKind::AmgPcg);
    check_solution(*dev, poisson(12, 9, 5, 10000), SolverKind::AmgPcg);
    check_solution(*dev, poisson(8, 7, 1, 1000), SolverKind::PcgJacobi);
  }
}

TEST_CASE("BiCGStab solves upwind convection diffusion against dense double reference") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    for (auto kind : {SolverKind::BicgstabJacobi, SolverKind::AmgBicgstab})
      check_solution(*dev, poisson(20, 18, 1, 1, false, false, 3), kind);
    check_solution(*dev, poisson(20, 18, 1, 1, false, false, 3),
                   SolverKind::AmgBicgstab, 2e-6, false, AmgCycle::K);
  }
}

TEST_CASE("constant Neumann nullspace is projected and incompatible RHS is rejected") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = poisson(8, 7, 1, 1, false, true);
    for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg}) {
      check_solution(*dev, a, kind, 2e-6, true);
      auto solver = make_solver(kind, *dev, a, {.nullspace = NullSpace::Constant});
      CHECK(solver->describe().find("nullspace") != std::string::npos);
      auto matrix = upload(*dev, a);
      auto x = dev->alloc<float>(a.pattern.n_rows), b = dev->alloc<float>(a.pattern.n_rows);
      dev->fill(x, 4.0f); dev->fill(b, 1.0f);
      CHECK_THROWS_AS(solver->solve(*dev, matrix, x, b, {}), Error);
    }
  }
}

TEST_CASE("controls zero residual nonzero initial guess and breakdown are truthful") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = poisson(5, 4);
    auto matrix = upload(*dev, a);
    auto solver = make_solver(SolverKind::PcgJacobi, *dev, a);
    auto x = dev->alloc<float>(a.pattern.n_rows), b = dev->alloc<float>(a.pattern.n_rows);
    dev->fill(x, 0); dev->fill(b, 0);
    CHECK(solver->solve(*dev, matrix, x, b, {}).iterations == 0);
    CHECK(solver->solve(*dev, matrix, x, b, {}).converged);
    dev->fill(b, 1);
    auto stats = solver->solve(*dev, matrix, x, b, {0, 0, 0});
    CHECK(stats.iterations == 0); CHECK_FALSE(stats.converged);
    CHECK(stats.final_residual == doctest::Approx(std::sqrt(double(a.pattern.n_rows))));
    CHECK(solver->solve(*dev, matrix, x, b, {0, 10, 0}).converged);
    CHECK_THROWS_AS(solver->solve(*dev, matrix, x, b, {-1, 0, 10}), Error);
    dev->fill(x, 0.5f);
    stats = solver->solve(*dev, matrix, x, b, {1e-5, 1e-10, 100});
    CHECK(stats.converged);
    // Negative definite input violates PCG's SPD contract without NaN propagation.
    for (auto& v : a.val) v = -v;
    matrix = upload(*dev, a); solver->refresh(*dev, matrix); dev->fill(x, 0);
    stats = solver->solve(*dev, matrix, x, b, {});
    CHECK_FALSE(stats.converged); CHECK(std::isfinite(stats.final_residual));
    dev->fill(matrix.val, 0);
    CHECK_THROWS_AS(solver->refresh(*dev, matrix), Error);
  }
}

TEST_CASE("AMG refresh agrees with fresh setup and V and K cycles are selectable") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = poisson(14, 13, 8);
    auto matrix = upload(*dev, a);
    for (auto cycle : {AmgCycle::V, AmgCycle::K}) {
      AmgOptions options; options.cycle = cycle;
      auto solver = make_amg_solver(SolverKind::AmgPcg, *dev, a, options);
      auto scaled = a;
      for (auto& v : scaled.val) v = static_cast<float>(v * 1.37);
      // Preserve relative edge strengths (and hence P), but change the reaction
      // per cell. A stale coarse operator cannot cancel this in Krylov scaling.
      for (Index i = 0; i < a.pattern.n_rows; ++i)
        scaled.val[a.pattern.diag[i]] = static_cast<float>(
            scaled.val[a.pattern.diag[i]] + 0.031 * (1 + hash(i) % 17));
      auto new_matrix = upload(*dev, scaled);
      // Refresh values only; keep the uploaded original pattern.
      dev->copy(matrix.val, new_matrix.val);
      solver->refresh(*dev, matrix);
      auto fresh = make_amg_solver(SolverKind::AmgPcg, *dev, scaled, options);
      auto b = dev->upload<float>(rhs(scaled, random_solution(a.pattern.n_rows)));
      auto x = dev->alloc<float>(a.pattern.n_rows), y = dev->alloc<float>(a.pattern.n_rows);
      dev->fill(x, 0); dev->fill(y, 0);
      auto s = solver->solve(*dev, matrix, x, b, {1e-6, 1e-12, 400});
      auto t = fresh->solve(*dev, matrix, y, b, {1e-6, 1e-12, 400});
      CHECK(s.converged); CHECK(t.converged); CHECK(s.iterations == t.iterations);
      CHECK(s.final_residual == t.final_residual);
      CHECK(dev->download(x) == dev->download(y));
      dev->copy(matrix.val, upload(*dev, a).val);
    }
  }
}

TEST_CASE("AMG 3D Poisson grid growth at requested 1e-8 is measured honestly") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    int counts[2] = {};
    for (int k = 0; k < 2; ++k) {
      Index size = k == 0 ? 32 : 64;
      auto a = poisson(size, size, size);
      auto matrix = upload(*dev, a);
      auto solver = make_solver(SolverKind::AmgPcg, *dev, a);
      // Representable manufactured solution exercises a stringent tolerance
      // without requiring arbitrary real-valued solutions to fit in float.
      std::vector<float> exact(a.pattern.n_rows);
      for (Index i = 0; i < a.pattern.n_rows; ++i)
        exact[i] = ((i % size + (i / size) % size + i / (size * size)) % 2) ? 1.0f : -1.0f;
      auto b = dev->upload<float>(rhs(a, exact));
      auto x = dev->alloc<float>(exact.size()); dev->fill(x, 0);
      auto s = solver->solve(*dev, matrix, x, b, {1e-8, 0, 400});
      auto xx = dev->download(x);
      auto rr = multiply(a, std::vector<double>(xx.begin(), xx.end()));
      auto bb = rhs(a, exact);
      double err = 0, real = 0;
      for (std::size_t i = 0; i < xx.size(); ++i) {
        err = std::max(err, std::abs(double(xx[i] - exact[i])));
        real += std::pow(double(bb[i]) - rr[i], 2);
      }
      CHECK(err < 2e-6);
      CHECK(std::sqrt(real) == doctest::Approx(s.final_residual).epsilon(1e-5));
      counts[k] = s.iterations;
      std::cout << "AMG_GRID backend=" << backend_name(dev->backend()) << " n=" << size
                << " iterations=" << s.iterations << " rel=" << s.final_residual / s.initial_residual
                << " converged=" << s.converged << " " << solver->describe() << '\n';
      CHECK(s.converged);
      CHECK(s.iterations <= 15);
    }
    CHECK(counts[1] < 1.5 * counts[0]);
  }
}

TEST_CASE("float representability limits never produce false convergence") {
  CsrMatrixHost a;
  a.pattern.n_rows = 1; a.pattern.row_ptr = {0, 1};
  a.pattern.col = {0}; a.pattern.diag = {0}; a.val = {7};
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto matrix = upload(*dev, a);
    auto b = dev->upload<float>(std::vector<float>{1});
    auto x = dev->alloc<float>(1); dev->fill(x, 0);
    auto solver = make_solver(SolverKind::PcgJacobi, *dev, a);
    auto stats = solver->solve(*dev, matrix, x, b, {1e-10, 0, 100});
    CHECK_FALSE(stats.converged);
    CHECK(stats.final_residual > 1e-10);
    const double xx = dev->download(x)[0];
    CHECK(stats.final_residual == doctest::Approx(std::abs(1 - 7 * xx)).epsilon(1e-6));
  }
}

TEST_CASE("AMG is grid robust for random forcing with low modes") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    int counts[2] = {};
    for (int k = 0; k < 2; ++k) {
      Index size = k == 0 ? 32 : 64;
      auto a = poisson(size, size, size);
      std::vector<float> bb(a.pattern.n_rows);
      for (Index i = 0; i < a.pattern.n_rows; ++i)
        bb[i] = float(int(hash(i + 53) % 2001) - 1000) / 1024.0f;
      auto b = dev->upload<float>(bb), x = dev->alloc<float>(bb.size());
      dev->fill(x, 0);
      auto matrix = upload(*dev, a);
      auto solver = make_solver(SolverKind::AmgPcg, *dev, a);
      auto stats = solver->solve(*dev, matrix, x, b, {1e-6, 0, 300});
      counts[k] = stats.iterations;
      std::cout << "AMG_RANDOM backend=" << backend_name(dev->backend()) << " n=" << size
                << " iterations=" << stats.iterations << " relative_residual="
                << stats.final_residual / stats.initial_residual << '\n';
      CHECK(stats.converged);
    }
    CHECK(counts[1] < 1.5 * counts[0]);
  }
}

TEST_CASE("a small reference term is preserved instead of projecting an invertible system") {
  CsrMatrixHost a;
  a.pattern.n_rows = 2; a.pattern.row_ptr = {0, 2, 4};
  a.pattern.col = {0, 1, 0, 1}; a.pattern.diag = {0, 3};
  a.val = {1 + std::ldexp(1.0, -20), -1, -1, 1 + std::ldexp(1.0, -20)};
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto matrix = upload(*dev, a);
    for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg}) {
      auto solver = make_solver(kind, *dev, a);
      CHECK(solver->describe().find("nullspace") == std::string::npos);
      auto x = dev->alloc<float>(2); dev->fill(x, 0);
      auto b = dev->upload<float>(std::vector<float>{1, 1});
      auto stats = solver->solve(*dev, matrix, x, b, {1e-6, 0, 100});
      CHECK(stats.converged);
      CHECK(dev->download(x)[0] == doctest::Approx(std::ldexp(1.0, 20)).epsilon(1e-6));
    }
  }
}

TEST_CASE("empty CSR and solvers are valid no-ops") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    CsrMatrixHost a;
    auto matrix = upload(*dev, a);
    auto x = dev->alloc<float>(0), b = dev->alloc<float>(0);
    for (auto kind : {SolverKind::Jacobi, SolverKind::PcgJacobi, SolverKind::BicgstabJacobi,
                      SolverKind::AmgPcg, SolverKind::AmgBicgstab}) {
      auto solver = make_solver(kind, *dev, a);
      auto stats = solver->solve(*dev, matrix, x, b, {});
      CHECK(stats.converged); CHECK(stats.iterations == 0); CHECK(stats.final_residual == 0);
    }
  }
}

TEST_CASE("strong upwind couplings coarsen in either direction and under permutation") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    for (double wind : {5.0, 20.0, 100.0}) {
      auto a = poisson(24, 16, 8, 1, false, false, wind);
      std::vector<Index> order(a.pattern.n_rows);
      std::iota(order.begin(), order.end(), 0);
      for (int ordering = 0; ordering < 3; ++ordering) {
        if (ordering == 1) std::reverse(order.begin(), order.end());
        if (ordering == 2) std::sort(order.begin(), order.end(), [](Index i, Index j) {
          return std::pair{hash(i), i} < std::pair{hash(j), j};
        });
        auto matrix = permute(a, order);
        auto solver = make_solver(SolverKind::AmgBicgstab, *dev, matrix);
        INFO("wind=" << wind << " order=" << ordering << " " << solver->describe());
        CHECK(solver->describe().find("1 levels") == std::string::npos);
      }
    }
  }
}

TEST_CASE("isolated diagonal rows do not prevent connected cells from coarsening") {
  auto a = poisson(8, 8);
  // Most cells are independent (e.g. masked solid cells), but the connected
  // fluid block must still get a coarse correction.
  for (int i = 0; i < 1024; ++i) {
    const Index row = a.pattern.n_rows++;
    a.pattern.diag.push_back(a.pattern.nnz()); a.pattern.col.push_back(row);
    a.val.push_back(2); a.pattern.row_ptr.push_back(a.pattern.nnz());
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto solver = make_solver(SolverKind::AmgPcg, *dev, a);
    CHECK(solver->describe().find("1 levels") == std::string::npos);
    check_solution(*dev, a, SolverKind::AmgPcg);
  }
}

TEST_CASE("default AMG avoids recursive host reductions and Krylov uses bounded storage") {
  InstrumentedDevice dev(Backend::Cpu);
  auto a = poisson(32, 32, 32);
  auto matrix = upload(dev, a);
  const auto before = dev.live_bytes;
  dev.peak_bytes = before;
  auto solver = make_solver(SolverKind::AmgPcg, dev, a);
  INFO("solver bytes/row=" << double(dev.live_bytes - before) / a.pattern.n_rows
       << " peak=" << double(dev.peak_bytes - before) / a.pattern.n_rows);
  CHECK(dev.live_bytes - before < 85 * std::size_t(a.pattern.n_rows));
  CHECK(dev.peak_bytes - before < 90 * std::size_t(a.pattern.n_rows));
  auto b = dev.upload<float>(random_solution(a.pattern.n_rows));
  auto x = dev.alloc<float>(b.size()); dev.fill(x, 0);
  dev.reductions = 0;
  auto s = solver->solve(dev, matrix, x, b, {0, 0, 3});
  CHECK(s.iterations == 3);
  CHECK(dev.reductions <= 3 * std::size_t(s.iterations) + 6);
}

TEST_CASE("Jacobi batches residual reductions and BiCGStab checks s only once") {
  for (auto kind : {SolverKind::Jacobi, SolverKind::BicgstabJacobi}) {
    InstrumentedDevice dev(Backend::Cpu);
    auto a = poisson(16, 12, 8, 1, false, false, 5);
    auto matrix = upload(dev, a);
    auto solver = make_solver(kind, dev, a);
    auto b = dev.upload<float>(random_solution(a.pattern.n_rows));
    auto x = dev.alloc<float>(b.size()); dev.fill(x, 0);
    dev.reductions = 0;
    auto s = solver->solve(dev, matrix, x, b, {0, 0, kind == SolverKind::Jacobi ? 17 : 4});
    INFO("reductions=" << dev.reductions << " iterations=" << s.iterations);
    if (kind == SolverKind::Jacobi) CHECK(dev.reductions <= 7);
    else CHECK(dev.reductions <= 6 * std::size_t(s.iterations) + 4);
  }
}

TEST_CASE("reused transport solves allocate nothing and avoid frozen reduction batches") {
  for(auto backend:Device::available()) {
    InstrumentedDevice dev(backend);
    auto a=poisson(8,6,4,1,false,false,5);
    for(auto d:a.pattern.diag) a.val[d]+=100;
    auto matrix=upload(dev,a);
    auto solver=make_solver(SolverKind::BicgstabJacobi,dev,a);
    auto reference=make_solver(SolverKind::BicgstabJacobi,dev,a);
    auto b=dev.upload<float>(random_solution(a.pattern.n_rows));
    auto x=dev.alloc<float>(b.size()),y=dev.alloc<float>(b.size());
    std::size_t first_reductions=0;
    for(int repeat=0;repeat<3;++repeat) {
      dev.fill(x,0);dev.fill(y,0);
      const auto allocations=dev.allocations;
      solver->refresh(dev,matrix);
      CHECK(dev.allocations==allocations);
      dev.device_reductions=0;
      const auto stats=solver->solve(dev,matrix,x,b,{1e-4,1e-12,100,8});
      REQUIRE(stats.converged);
      if(repeat==0) first_reductions=dev.device_reductions;
      else CHECK(dev.device_reductions<first_reductions);
      const auto oracle=reference->solve(dev,matrix,y,b,{1e-4,1e-12,100,1});
      CHECK(stats.iterations==oracle.iterations);
      CHECK(stats.final_residual==oracle.final_residual);
      CHECK(dev.download(x)==dev.download(y));
    }
    // Reusing validation storage must still reject bad refreshed coefficients.
    auto bad=std::vector<float>(a.val.begin(),a.val.end());
    bad.front()=std::numeric_limits<float>::quiet_NaN();
    dev.upload<float>(matrix.val,bad);
    CHECK_THROWS_AS(solver->refresh(dev,matrix),Error);
  }
}

TEST_CASE("Krylov exits do not reread an already verified current residual") {
  for(auto backend:Device::available()) for(auto kind:{SolverKind::PcgJacobi,SolverKind::BicgstabJacobi}) {
    InstrumentedDevice dev(backend);
    CsrMatrixHost a{{1,{0,1},{0},{0}},{2}};
    auto matrix=upload(dev,a);auto solver=make_solver(kind,dev,a);
    auto x=dev.alloc<float>(1),b=dev.alloc<float>(1);dev.fill(x,0);dev.fill(b,1);
    dev.reductions=0;
    const auto stats=solver->solve(dev,matrix,x,b,{1e-6,1e-12,20});
    REQUIRE(stats.converged);CHECK(stats.iterations==1);
    CHECK(stats.final_residual==0);CHECK(dev.download(x)==std::vector<float>{.5f});
    // Input b/x, initial defect, and one true b-A*x verification. A second
    // final read adds a GPU wait without observing a new solution.
    CHECK(dev.reductions==4);
  }
}

TEST_CASE("assembled random Neumann pressure uses the caller nullspace contract") {
  for (bool float_diagonal : {false, true}) {
    auto a = poisson(12, 10, 4, 1, true, true);
    // FV assembly sums float off-diagonals, introducing harmless row-sum noise.
    if (float_diagonal) for (Index i = 0; i < a.pattern.n_rows; ++i) {
      float d = 0;
      for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
        if (a.pattern.col[j] != i) d -= static_cast<float>(a.val[j]);
      a.val[a.pattern.diag[i]] = d;
    }
    KUTTA_FOR_EACH_BACKEND(dev) {
      auto matrix = upload(*dev, a);
      auto b = dev->alloc<float>(a.pattern.n_rows), x = dev->alloc<float>(a.pattern.n_rows);
      dev->fill(b, 0); dev->fill(x, 4);
      for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg}) {
        auto solver = make_solver(kind, *dev, a, {.nullspace = NullSpace::Constant});
        CHECK(solver->describe().find("nullspace") != std::string::npos);
        auto stats = solver->solve(*dev, matrix, x, b, {1e-5, 1e-7, 300});
        CHECK(stats.converged);
        CHECK(std::abs(dev->sum(x, x.size())) < 1e-5);
        std::vector<float> forcing(a.pattern.n_rows);
        forcing.front() = 1; forcing.back() = -1;
        dev->upload<float>(b, forcing);
        stats = solver->solve(*dev, matrix, x, b, {1e-5, 1e-7, 300});
        CHECK(stats.converged);
        const auto reference = reference_cg(a, forcing, true);
        const auto got = dev->download(x);
        double error = 0, scale = 0;
        for (std::size_t i = 0; i < got.size(); ++i) {
          error += std::pow(got[i] - reference[i], 2); scale += reference[i] * reference[i];
        }
        CHECK(std::sqrt(error / scale) < 2e-4);
        CHECK(std::abs(dev->sum(x, x.size())) < 1e-5);
        dev->fill(b, 1);
        CHECK_THROWS_AS(solver->solve(*dev, matrix, x, b, {}), Error);
        dev->fill(b, 0); dev->fill(x, 4);
      }
    }
  }
}

TEST_CASE("AMG rejects negative diagonals and Krylov restores its best verified iterate on divergence") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = poisson(16, 12, 4, 1, false, false, 20);
    for (Index i = 0; i < a.pattern.n_rows; i += 7)
      for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j) a.val[j] = -a.val[j];
    CHECK_THROWS_WITH_AS(make_solver(SolverKind::AmgBicgstab, *dev, a),
                         doctest::Contains("AMG: requires positive diagonals"), Error);
    auto positive = poisson(16, 12, 4, 1, false, false, 20);
    auto amg = make_solver(SolverKind::AmgBicgstab, *dev, positive);
    auto changed = upload(*dev, a);
    CHECK_THROWS_WITH_AS(amg->refresh(*dev, changed),
                         doctest::Contains("AMG: requires positive diagonals"), Error);
    auto matrix = upload(*dev, a);
    auto solver = make_solver(SolverKind::BicgstabJacobi, *dev, a);
    auto b = dev->upload<float>(rhs(a, random_solution(a.pattern.n_rows)));
    auto x = dev->alloc<float>(b.size()); dev->fill(x, 0);
    CHECK(solver->solve(*dev, matrix, x, b, {1e-5, 0, 300}).converged);
    // Deliberately violates PCG's SPD contract; both Krylov methods make a
    // huge first step. Returning it would amplify the initial residual 1e5x.
    CsrMatrixHost bad;
    bad.pattern.n_rows = 2; bad.pattern.row_ptr = {0, 2, 4};
    bad.pattern.col = {0, 1, 0, 1}; bad.pattern.diag = {0, 3};
    bad.val = {1, 999999, -999998, 1};
    matrix = upload(*dev, bad);
    for (auto kind : {SolverKind::PcgJacobi, SolverKind::BicgstabJacobi}) {
      solver = make_solver(kind, *dev, bad);
      b = dev->upload<float>(std::vector<float>{1, 1});
      x = dev->alloc<float>(2); dev->fill(x, 0);
      auto stats = solver->solve(*dev, matrix, x, b, {1e-8, 0, 8});
      CHECK_FALSE(stats.converged);
      CHECK(stats.final_residual <= stats.initial_residual);
      CHECK(dev->download(x) == std::vector<float>{0, 0});
      CHECK(stats.iterations == 1);
    }
  }
}

TEST_CASE("actual hierarchy Galerkin values match independent double PtAP before and after refresh") {
  for (auto backend : Device::available()) {
    InstrumentedDevice dev(backend);
    auto a = poisson(12, 10, 6, 1, true);
    auto matrix = upload(dev, a);
    auto solver = make_solver(SolverKind::AmgPcg, dev, a);
    auto b = dev.upload<float>(random_solution(a.pattern.n_rows));
    auto x = dev.alloc<float>(b.size());
    for (int refresh = 0; refresh < 2; ++refresh) {
      if (refresh) {
        for (Index i = 0; i < a.pattern.n_rows; ++i)
          for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
            a.val[j] = static_cast<float>(a.val[j] * (1.13 + 0.01 * (hash(i + a.pattern.col[j]) % 7)));
        // Keep SPD by adding a positive reaction larger than the perturbation.
        for (auto diagonal : a.pattern.diag) a.val[diagonal] += 10;
        auto changed = upload(dev, a); dev.copy(matrix.val, changed.val);
        const auto downloads = dev.downloads;
        solver->refresh(dev, matrix);
        CHECK(dev.downloads == downloads);
      }
      std::map<std::size_t, CsrMatrixHost> operators;
      std::map<std::size_t, std::vector<Index>> aggregates;
      dev.after_launch = [&](KernelId id, std::size_t n, std::span<const KernelArg> args) {
        if (id == KernelId::amg_coarse_solve) n = args[7].get<Index>();
        if ((id == KernelId::csr_residual || id == KernelId::csr_residual_guarded || id == KernelId::amg_smooth ||
             id == KernelId::amg_coarse_solve) && !operators.contains(n)) {
          auto& op = operators[n]; op.pattern.n_rows = static_cast<Index>(n);
          op.pattern.row_ptr = dev.read<Index>(*args[0].buffer);
          op.pattern.col = dev.read<Index>(*args[1].buffer);
          const auto values = dev.read<float>(*args[2].buffer);
          op.val.assign(values.begin(), values.end());
        }
        if ((id == KernelId::amg_prolong || id == KernelId::amg_prolong_guarded) && !aggregates.contains(n))
          aggregates[n] = dev.read<Index>(*args[0].buffer);
      };
      dev.fill(x, 0);
      solver->solve(dev, matrix, x, b, {0, 0, 1});
      dev.after_launch = {};
      REQUIRE(operators.size() >= 3);
      for (auto fine = operators.rbegin(); fine != operators.rend(); ++fine) {
        const auto coarse = std::next(fine);
        if (coarse == operators.rend()) break;
        const auto& op = fine->second;
        const auto& cp = coarse->second;
        REQUIRE(aggregates.contains(fine->first));
        const auto& map = aggregates.at(fine->first);
        // Independent row/column assembly in double: do not use the device
        // entry map or the implementation's accumulation order as the oracle.
        std::map<std::pair<Index, Index>, double> expected, magnitude;
        for (Index i = 0; i < op.pattern.n_rows; ++i) {
          if (map[i] < 0) continue;
          for (Index j = op.pattern.row_ptr[i]; j < op.pattern.row_ptr[i + 1]; ++j) {
            const Index col = map[op.pattern.col[j]];
            if (col < 0) continue;
            const auto key = std::pair{map[i], col};
            expected[key] += op.val[j]; magnitude[key] += std::abs(op.val[j]);
          }
        }
        CHECK(expected.size() == cp.val.size());
        for (Index i = 0; i < cp.pattern.n_rows; ++i)
          for (Index j = cp.pattern.row_ptr[i]; j < cp.pattern.row_ptr[i + 1]; ++j) {
            const auto key = std::pair{i, cp.pattern.col[j]};
            REQUIRE(expected.contains(key));
            CHECK(std::abs(cp.val[j] - expected.at(key)) <= 2e-7 * magnitude.at(key));
          }
      }
    }
  }
}

TEST_CASE("CPU and Metal solves have identical iteration counts residuals and solutions") {
  auto cpu = Device::create(Backend::Cpu);
  for (const auto& a : {poisson(4, 3), poisson(5, 4, 3),
                        poisson(5, 4, 1, 1, false, false, 20)}) {
    for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg,
                      SolverKind::BicgstabJacobi, SolverKind::AmgBicgstab}) {
      if (a.val[a.pattern.diag[0]] > 20 &&
          (kind == SolverKind::PcgJacobi || kind == SolverKind::AmgPcg)) continue;
      auto run = [&](Device& dev) {
        auto matrix = upload(dev, a); auto solver = make_solver(kind, dev, a);
        auto b = dev.upload<float>(rhs(a, random_solution(a.pattern.n_rows)));
        auto x = dev.alloc<float>(b.size()); dev.fill(x, 0);
        auto stats = solver->solve(dev, matrix, x, b, {1e-5, 0, 200});
        CHECK(stats.converged);
        return std::pair{stats, dev.download(x)};
      };
      const auto reference = run(*cpu);
      KUTTA_FOR_EACH_BACKEND(dev) {
        const auto result = run(*dev);
        CHECK(result.first.iterations == reference.first.iterations);
        CHECK(result.first.initial_residual == reference.first.initial_residual);
        CHECK(result.first.final_residual == reference.first.final_residual);
        CHECK(result.second == reference.second);
      }
    }
  }
}

TEST_CASE("Neumann PCG avoids duplicate fine projection and unused final preconditioning") {
  for (auto backend : Device::available()) {
    InstrumentedDevice dev(backend);
    const auto a = poisson(8, 8, 1, 1, false, true);
    auto matrix = upload(dev, a);
    auto solver = make_solver(SolverKind::AmgPcg, dev, a, {.nullspace = NullSpace::Constant});
    std::vector<float> forcing(a.pattern.n_rows);
    forcing.front() = 1; forcing.back() = -1;
    auto b = dev.upload<float>(forcing), x = dev.alloc<float>(b.size());
    dev.fill(x, 0); dev.reductions = 0;
    const auto stats = solver->solve(dev, matrix, x, b, {0, 0, 1});
    CHECK(stats.iterations == 1);
    CHECK_FALSE(stats.converged);
    INFO("Neumann reductions=" << dev.reductions);
    CHECK(dev.reductions <= 15);
  }
}

TEST_CASE("batched BiCGStab freezes on zero omega and zero shadow product") {
  CsrMatrixHost a;
  a.pattern.n_rows = 2; a.pattern.row_ptr = {0, 2, 4};
  a.pattern.col = {0, 1, 0, 1}; a.pattern.diag = {0, 3};
  // A*[1,1]=[4,0], alpha=1/2, s=[-1,1], A*s=[2,2].
  // Thus t.s=0 while t.t=8: an exact omega breakdown after one attempt.
  for (bool zero_omega : {true, false}) {
    a.val = zero_omega ? std::vector<double>{1, 3, -1, 1} :
                         std::vector<double>{1, 2, -4, 1}; // rhat.v=0
    for (auto backend : Device::available()) {
      InstrumentedDevice dev(backend);
      auto matrix = upload(dev, a);
      auto solver = make_solver(SolverKind::BicgstabJacobi, dev, a);
      auto b = dev.upload<float>(std::vector<float>{1, 1}), x = dev.alloc<float>(2);
      for (int interval : {1, 8, 31}) {
        dev.fill(x, 0); dev.downloads = dev.reductions = 0;
        auto stats = solver->solve(dev, matrix, x, b, {1e-8, 0, 31, interval});
        CHECK_FALSE(stats.converged);
        CHECK(stats.iterations == (zero_omega ? 1 : 0));
        CHECK(stats.final_residual == stats.initial_residual);
        CHECK(dev.downloads == 1); // one batch poll, even after frozen dispatches
        CHECK(dev.reductions == 4); // input norms, initial and actual final norms
        CHECK(dev.download(x) == std::vector<float>{0, 0});
      }
    }
  }
}

TEST_CASE("Krylov host checks are batched including K-cycle and Neumann projection") {
  for (auto backend : Device::available()) {
    InstrumentedDevice dev(backend);
    for (auto cycle : {AmgCycle::V, AmgCycle::K}) {
      for (bool neumann : {false, true}) {
        auto a = poisson(12, 10, 4, 1, false, neumann);
        auto matrix = upload(dev, a);
        AmgOptions options; options.cycle = cycle;
        options.nullspace = neumann ? NullSpace::Constant : NullSpace::None;
        for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg,
                          SolverKind::BicgstabJacobi, SolverKind::AmgBicgstab}) {
          if (neumann && (kind == SolverKind::BicgstabJacobi || kind == SolverKind::AmgBicgstab)) continue;
          auto solver = kind == SolverKind::PcgJacobi || kind == SolverKind::BicgstabJacobi ?
              make_solver(kind, dev, a, {.nullspace = options.nullspace}) :
              make_amg_solver(kind, dev, a, options);
          auto forcing = random_solution(a.pattern.n_rows);
          if (neumann) { std::fill(forcing.begin(), forcing.end(), 0); forcing.front() = 1; forcing.back() = -1; }
          auto b = dev.upload<float>(forcing), x = dev.alloc<float>(forcing.size());
          dev.fill(x, 0); dev.downloads = dev.reductions = dev.device_reductions = 0;
          auto stats = solver->solve(dev, matrix, x, b, {0, 0, 3, 8});
          CHECK(stats.iterations == 3);
          CHECK(dev.downloads == 1);
          CHECK(dev.reductions == (neumann ? 5 : 4));
          CHECK(dev.device_reductions > 0);
        }
      }
    }
  }
}
