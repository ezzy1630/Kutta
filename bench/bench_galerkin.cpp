#include <chrono>
#include <iomanip>
#include <iostream>

#include "kernels.hpp"
#include "../tests/linalg/matrices.hpp"

// Compare deterministic Galerkin refresh schemes on eight-cell aggregates of
// the 128^3 FV stencil. Both kernels use the same fine-index summation order.
int main() {
  using namespace kutta;
  using namespace kutta::test_linalg;
  constexpr Index side = 128, half = side / 2, nc = half * half * half;
  const auto a = poisson(side, side, side);
  std::vector<Index> map(a.pattern.n_rows), offsets(nc + 1), members(a.pattern.n_rows);
  for (Index i = 0; i < a.pattern.n_rows; ++i) {
    map[i] = (i % side) / 2 + ((i / side) % side) / 2 * half +
             (i / (side * side)) / 2 * half * half;
    ++offsets[map[i] + 1];
  }
  std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
  auto cursor = offsets;
  for (Index i = 0; i < a.pattern.n_rows; ++i) members[cursor[map[i]]++] = i;
  CsrPattern coarse; coarse.n_rows = nc;
  std::vector<Index> entry_offsets{0}, entries;
  std::vector<std::pair<Index, Index>> row;
  for (Index g = 0; g < nc; ++g) {
    row.clear();
    for (Index m = offsets[g]; m < offsets[g + 1]; ++m) {
      const auto i = members[m];
      for (Index j = a.pattern.row_ptr[i]; j < a.pattern.row_ptr[i + 1]; ++j)
        row.emplace_back(map[a.pattern.col[j]], j);
    }
    std::sort(row.begin(), row.end());
    Index last = -1;
    for (auto [col, source] : row) {
      if (col != last) {
        if (last != -1) entry_offsets.push_back(static_cast<Index>(entries.size()));
        coarse.col.push_back(col); last = col;
      }
      entries.push_back(source);
    }
    entry_offsets.push_back(static_cast<Index>(entries.size()));
    coarse.row_ptr.push_back(coarse.nnz());
  }
  for (auto backend : Device::available()) {
    auto dev = Device::create(backend);
    auto fine = upload(*dev, a);
    auto dm = dev->upload<Index>(map), dmo = dev->upload<Index>(offsets), dms = dev->upload<Index>(members);
    auto dcr = dev->upload<Index>(coarse.row_ptr), dcc = dev->upload<Index>(coarse.col);
    auto deo = dev->upload<Index>(entry_offsets), de = dev->upload<Index>(entries);
    auto mapped = dev->alloc<float>(coarse.nnz()), gathered = dev->alloc<float>(coarse.nnz());
    auto run = [&](bool gather) {
      if (gather) kern::amg_galerkin_rows(*dev, nc, dcr, dcc, dmo, dms,
                                         fine.row_ptr, fine.col, dm, fine.val, gathered);
      else kern::amg_galerkin(*dev, coarse.nnz(), deo, de, fine.val, mapped);
    };
    run(false); run(true); dev->synchronize();
    KUTTA_REQUIRE(dev->download(mapped) == dev->download(gathered), "Galerkin alternatives disagree");
    for (bool gather : {false, true}) {
      std::vector<double> samples;
      for (int sample = 0; sample < 5; ++sample) {
        dev->synchronize(); const auto start = std::chrono::steady_clock::now();
        for (int repeat = 0; repeat < 20; ++repeat) run(gather);
        dev->synchronize();
        samples.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count() / 20);
      }
      std::sort(samples.begin(), samples.end());
      std::cout << std::setprecision(9) << "{\"benchmark\":\"galerkin\",\"backend\":\""
                << backend_name(backend) << "\",\"scheme\":\"" << (gather ? "row_gather" : "entry_map")
                << "\",\"rows\":" << a.pattern.n_rows << ",\"median_s\":" << samples[2]
                << ",\"saved_bytes_per_row\":" <<
                   (gather ? double((entries.size() + entry_offsets.size()) * sizeof(Index)) / a.pattern.n_rows : 0)
                << "}" << std::endl;
    }
  }
}
