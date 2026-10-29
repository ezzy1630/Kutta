#include <iostream>
#include <limits>

#include "../backend_helpers.hpp"
#include "kutta/linalg/amg.hpp"
#include "matrices.hpp"

using namespace kutta;
using namespace kutta::test_linalg;

namespace {
void check_solution(Device& dev, const CsrMatrixHost& a, SolverKind kind,
                    double tolerance = 2e-6, bool neumann = false, AmgCycle cycle = AmgCycle::K) {
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
  auto solver = kind == SolverKind::AmgPcg || kind == SolverKind::AmgBicgstab ?
      make_amg_solver(kind, dev, a, options) : make_solver(kind, dev, a);
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
                   SolverKind::AmgBicgstab, 2e-6, false, AmgCycle::V);
  }
}

TEST_CASE("constant Neumann nullspace is projected and incompatible RHS is rejected") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto a = poisson(8, 7, 1, 1, false, true);
    for (auto kind : {SolverKind::PcgJacobi, SolverKind::AmgPcg}) {
      check_solution(*dev, a, kind, 2e-6, true);
      auto solver = make_solver(kind, *dev, a);
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
      auto scaled = a; for (auto& v : scaled.val) v *= 4;
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
