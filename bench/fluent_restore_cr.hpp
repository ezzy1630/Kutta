#pragma once
// Forensic restoration of a legacy Fluent .dat whose 0x0D bytes preceding
// 0x0A were deleted by a text-mode transfer. The adjoint cylinder tutorial's
// data file in ansys/example-data (commit f09067df48, ".gitattributes:
// * text=auto") lost them: its binary bodies are 15 bytes short and the file
// holds no CR LF pair where random doubles would hold about nine.
//
// kutta::fluent::read_data never repairs input; it rejects such a file and
// names the likely cause. This restoration is explicit and checked: each
// short 2300/3300 field body gets exactly as many CRs back as it lacks, each
// at the unique LF whose restoration realigns the following values, and the
// result must parse with every body length matching its header.
#include <cmath>
#include <cstdint>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace kutta::forensics {
struct CrInsertion {
  std::string section;      // the section header, e.g. "3300 (1 11 1 0 1 1 25740)"
  std::size_t offset = 0;   // byte offset of the restored CR in the output
  std::size_t value = 0;    // index of the value it belongs to in that section
};
struct CrRestoration {
  std::string bytes;
  std::vector<CrInsertion> inserted;
};

namespace detail_cr {
// Misaligned doubles decode to huge, tiny or non-finite values. A field value
// is "plausible" when zero or within [1e-30, max] in magnitude.
inline bool plausible(std::string_view bytes, std::size_t j, std::size_t width, double max) {
  if ((j + 1) * width > bytes.size())
    return true;
  double v = 0;
  if (width == 8) std::memcpy(&v, bytes.data() + j * 8, 8);
  else {float f = 0; std::memcpy(&f, bytes.data() + j * 4, 4); v = f;}
  return v == 0 || (std::isfinite(v) && std::abs(v) >= 1e-30 && std::abs(v) <= max);
}
inline std::string with_cr(std::string_view body, const std::vector<std::size_t>& before) {
  std::string out;out.reserve(body.size() + before.size());
  std::size_t from = 0;
  for (auto c : before) {out.append(body.substr(from, c - from));out.push_back('\r');from = c;}
  out.append(body.substr(from));
  return out;
}
// Positions in `body` before which a CR is restored, so that the body holds
// n values of `width` bytes. Little-endian bodies only (machine flag 60).
inline std::vector<std::size_t> locate(std::string_view body, std::size_t n, std::size_t width, double max,
                                       const std::string& name) {
  const std::size_t lost = n * width - body.size();
  constexpr std::size_t behind = 8, ahead = 64;
  std::vector<std::size_t> before;
  std::string restored(body);
  std::size_t j = 0;
  while (before.size() < lost) {
    while (j < n && plausible(restored, j, width, max)) ++j;
    if (j >= n) throw std::runtime_error(name + ": a lost byte cannot be located from the values");
    // Damaged-body offset of value j's first byte; all prior CRs lie before it.
    const std::size_t at = j * width - before.size();
    const std::size_t lo = std::max(before.empty() ? 0 : before.back() + 1, at > behind * width ? at - behind * width : 0);
    const std::size_t hi = std::min(at + width, body.size());
    std::vector<std::size_t> fits;
    for (std::size_t c = lo; c < hi; ++c) {
      if (body[c] != '\n') continue;
      auto trial = before;trial.push_back(c);
      const auto t = with_cr(body, trial);
      const std::size_t first = (c + before.size()) / width;
      bool ok = true;
      for (std::size_t q = first > 2 ? first - 2 : 0; q < std::min(n, first + ahead) && ok; ++q) ok = plausible(t, q, width, max);
      if (ok) fits.push_back(c);
    }
    if (fits.size() != 1)
      throw std::runtime_error(name + ": " + std::to_string(fits.size()) + " candidate CR positions near value " + std::to_string(j));
    before.push_back(fits.front());
    restored = with_cr(body, before);
    const std::size_t first = (fits.front() + before.size() - 1) / width;
    j = first > 2 ? first - 2 : 0;
  }
  for (std::size_t q = 0; q < n; ++q)
    if (!plausible(restored, q, width, max)) throw std::runtime_error(name + ": restored body still has implausible values");
  return before;
}
}  // namespace detail_cr

// `max` bounds plausible field magnitudes; 1e7 suits the cylinder's SI fields.
inline CrRestoration restore_stripped_cr(std::string_view in, double max = 1e7) {
  if (in.find("\r\n") != std::string_view::npos)
    throw std::runtime_error("file holds CR LF pairs; it was not stripped by a text-mode transfer");
  if (in.find("(4 (60 ") == std::string_view::npos)
    throw std::runtime_error("restoration supports little-endian (machine flag 60) files only");
  CrRestoration out;
  std::size_t copied = 0, p = 0;
  const std::string_view marker = "\nEnd of Binary Section";
  while ((p = in.find("\n(", p)) != std::string_view::npos) {
    p += 2;
    const bool wide = in.substr(p, 6) == "3300 (";
    const bool narrow = in.substr(p, 6) == "2300 (";
    if (!wide && !narrow) continue;
    const auto close = in.find(')', p);
    if (close == std::string_view::npos) throw std::runtime_error("truncated field header");
    const std::string section(in.substr(p, close + 1 - p));
    std::vector<std::uint64_t> h;
    for (std::size_t q = p + 6; q < close;) {
      while (q < close && in[q] == ' ') ++q;
      std::uint64_t v = 0;std::size_t digits = 0;
      while (q < close && in[q] >= '0' && in[q] <= '9') {v = 10 * v + std::uint64_t(in[q++] - '0');++digits;}
      if (!digits) {if (q < close) throw std::runtime_error("invalid field header " + section);break;}
      h.push_back(v);
    }
    if (h.size() != 7 || h[6] < h[5] || !h[2]) throw std::runtime_error("invalid field header " + section);
    const std::size_t width = wide ? 8 : 4, n = std::size_t(h[6] - h[5] + 1) * std::size_t(h[2]);
    const auto start = in.find('(', close + 1) + 1;
    const auto end = in.find(marker, start);
    if (start == 0 || end == std::string_view::npos) throw std::runtime_error("unterminated field " + section);
    const auto body = in.substr(start, end - start);
    if (body.size() > n * width) throw std::runtime_error("field body longer than its header: " + section);
    out.bytes.append(in.substr(copied, start - copied));
    if (body.size() < n * width) {
      const auto before = detail_cr::locate(body, n, width, max, section);
      for (std::size_t k = 0; k < before.size(); ++k) {
        const auto offset = out.bytes.size() + before[k] + k;
        out.inserted.push_back({section, offset, (before[k] + k) / width});
      }
      out.bytes.append(detail_cr::with_cr(body, before));
    } else out.bytes.append(body);
    copied = end;p = end;
  }
  out.bytes.append(in.substr(copied));
  return out;
}
}  // namespace kutta::forensics
