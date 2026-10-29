#include <bit>
#include <cctype>
#include <charconv>
#include <fstream>
#include <string_view>

#include "kutta/io/fluent_detail.hpp"

namespace kutta {
namespace {
class Parser {
 public:
  explicit Parser(std::string bytes) : bytes_(std::move(bytes)) {}
  bool end() {
    space();
    return pos_ == bytes_.size();
  }
  void space() {
    while (pos_ < bytes_.size() &&
           std::isspace(static_cast<unsigned char>(bytes_[pos_])))
      ++pos_;
  }
  bool next(char c) {
    space();
    return pos_ < bytes_.size() && bytes_[pos_] == c;
  }
  void expect(char c) {
    space();
    KUTTA_REQUIRE(pos_ < bytes_.size() && bytes_[pos_] == c,
                  "Fluent: expected '" + std::string(1, c) + "' at byte " +
                      std::to_string(pos_));
    ++pos_;
  }
  std::string token() {
    space();
    KUTTA_REQUIRE(pos_ < bytes_.size(), "Fluent: unexpected end of file");
    if (bytes_[pos_] == '"') {
      ++pos_;
      std::string out;
      bool closed = false;
      while (pos_ < bytes_.size()) {
        char c = bytes_[pos_++];
        if (c == '"') {
          closed = true;
          break;
        }
        if (c == '\\') {
          KUTTA_REQUIRE(pos_ < bytes_.size(),
                        "Fluent: truncated quoted escape");
          c = bytes_[pos_++];
        }
        out += c;
      }
      KUTTA_REQUIRE(closed, "Fluent: unterminated quoted string");
      return out;
    }
    auto start = pos_;
    while (pos_ < bytes_.size() &&
           !std::isspace(static_cast<unsigned char>(bytes_[pos_])) &&
           bytes_[pos_] != '(' && bytes_[pos_] != ')')
      ++pos_;
    KUTTA_REQUIRE(start != pos_,
                  "Fluent: expected token at byte " + std::to_string(pos_));
    return bytes_.substr(start, pos_ - start);
  }
  static Index integer(const std::string& s, int base = 16) {
    std::uint64_t v = 0;
    auto [p, e] = std::from_chars(s.data(), s.data() + s.size(), v, base);
    KUTTA_REQUIRE(
        e == std::errc{} && p == s.data() + s.size() &&
            v <= static_cast<std::uint64_t>(std::numeric_limits<Index>::max()),
        "Fluent: invalid/out-of-range integer '" + s + "'");
    return static_cast<Index>(v);
  }
  Index integer(int base = 16) { return integer(token(), base); }
  double real() {
    auto s = token();
    for (auto& c : s)
      if (c == 'D' || c == 'd') c = 'E';
    // from_chars intentionally rejects a leading '+'. Fluent accepts it.
    auto begin = s.data() + (s.starts_with('+') ? 1 : 0);
    double value = 0;
    auto [p, e] = std::from_chars(begin, s.data() + s.size(), value);
    KUTTA_REQUIRE(
        e == std::errc{} && p == s.data() + s.size() && std::isfinite(value),
        "Fluent: invalid coordinate '" + s + "'");
    return value;
  }
  std::vector<std::string> header() {
    expect('(');
    std::vector<std::string> h;
    while (!next(')')) h.push_back(token());
    expect(')');
    return h;
  }
  void rest() {  // remainder of an ASCII section, outer '(' already consumed
    int depth = 1;
    while (depth) {
      space();
      KUTTA_REQUIRE(pos_ < bytes_.size(), "Fluent: truncated section");
      auto c = bytes_[pos_];
      if (c == '(') {
        ++pos_;
        ++depth;
      } else if (c == ')') {
        ++pos_;
        --depth;
      } else
        token();
    }
  }
  std::uint64_t binary(unsigned width) {
    KUTTA_REQUIRE(width <= 8 && width <= bytes_.size() - pos_,
                  "Fluent: truncated binary payload");
    std::uint64_t v = 0;
    for (unsigned i = 0; i < width; ++i) {
      auto b = static_cast<unsigned char>(bytes_[pos_++]);
      v |= static_cast<std::uint64_t>(b) << (8 * (little_ ? i : width - 1 - i));
    }
    return v;
  }
  void skip_binary(std::uint64_t count, unsigned width) {
    KUTTA_REQUIRE(count <= (bytes_.size() - pos_) / width,
                  "Fluent: truncated binary payload");
    pos_ += static_cast<std::size_t>(count) * width;
  }
  Index binary_index() {
    auto v = binary(4);
    KUTTA_REQUIRE(
        v <= static_cast<std::uint64_t>(std::numeric_limits<Index>::max()),
        "Fluent: binary connectivity exceeds 32-bit index range");
    return static_cast<Index>(v);
  }
  double binary_real(bool wide) {
    double v = wide ? std::bit_cast<double>(binary(8))
                    : static_cast<double>(std::bit_cast<float>(
                          static_cast<std::uint32_t>(binary(4))));
    KUTTA_REQUIRE(
        std::isfinite(v),
        "Fluent: non-finite binary coordinate (check machine byte order)");
    return v;
  }
  void finish(bool binary_section, int section, bool body = true) {
    if (body) expect(')');
    if (binary_section) {
      KUTTA_REQUIRE(token() == "End" && token() == "of" &&
                        token() == "Binary" && token() == "Section",
                    "Fluent: missing binary section terminator");
      KUTTA_REQUIRE(integer(10) == section,
                    "Fluent: mismatched binary section terminator");
    }
    expect(')');
  }
  void byte_order(int flag) { little_ = flag == 60; }

 private:
  std::string bytes_;
  std::size_t pos_ = 0;
  bool little_ = true;
};

detail::FluentRaw parse(std::string bytes) {
  Parser p(std::move(bytes));
  detail::FluentRaw raw;
  std::map<int, detail::FluentZone> named;
  std::map<int, int> cell_types;
  Index declared_nodes = 0, declared_cells = 0, declared_faces = 0;
  bool dimension_seen = false;
  while (!p.end()) {
    p.expect('(');
    int section = p.integer(10);
    bool binary = section >= 2000;
    int code = binary ? section % 1000 : section;
    if (code == 2) {
      auto dimension = p.integer(10);
      KUTTA_REQUIRE(!dimension_seen || raw.dimension == dimension,
                    "Fluent: inconsistent grid dimensions");
      raw.dimension = dimension;
      dimension_seen = true;
      KUTTA_REQUIRE(raw.dimension == 2 || raw.dimension == 3,
                    "Fluent: dimension must be 2 or 3");
      p.expect(')');
      continue;
    }
    if (code == 4) {
      auto h = p.header();
      KUTTA_REQUIRE(!h.empty(), "Fluent: empty machine configuration");
      p.byte_order(Parser::integer(h[0], 10));
      p.rest();
      continue;
    }
    if (code == 39 || code == 45) {
      auto h = p.header();
      KUTTA_REQUIRE(h.size() >= 3, "Fluent: malformed zone metadata");
      int z = Parser::integer(h[0], 10);
      KUTTA_REQUIRE(named.emplace(z, detail::FluentZone{h[2], h[1]}).second,
                    "Fluent: duplicate named zone");
      p.rest();
      continue;
    }
    KUTTA_REQUIRE(code != 58 && code != 59 && code != 61,
                  "Fluent: adapted/hanging-node/nonconformal meshes are not "
                  "supported; export a conformal active-cell mesh");
    if (code == 41 && binary) {
      // Node flags: two hexadecimal endpoint ids, one 32-bit flag per node
      // in either precision. Flags are not connectivity and may use all bits.
      auto h = p.header();
      KUTTA_REQUIRE(h.size() == 2, "Fluent: malformed node flags header");
      auto first = Parser::integer(h[0]), last = Parser::integer(h[1]);
      KUTTA_REQUIRE(first > 0 && last >= first,
                    "Fluent: invalid node flags range");
      p.expect('(');
      p.skip_binary(static_cast<std::uint64_t>(last) - first + 1, 4);
      p.finish(true, section);
      continue;
    }
    if (code != 10 && code != 12 && code != 13 && code != 18) {
      KUTTA_REQUIRE(!binary, "Fluent: unsupported binary section " +
                                 std::to_string(section));
      p.rest();
      continue;
    }
    auto h = p.header();
    KUTTA_REQUIRE(h.size() >= 4,
                  "Fluent: grid section header requires at least four values");
    std::vector<Index> v;
    for (const auto& t : h) v.push_back(Parser::integer(t));
    if (code == 10 && v.size() > 4) {
      auto nd = v[4];
      KUTTA_REQUIRE(nd == 2 || nd == 3,
                    "Fluent: node dimension must be 2 or 3");
      KUTTA_REQUIRE(!dimension_seen || raw.dimension == nd,
                    "Fluent: inconsistent node dimension");
      raw.dimension = nd;
      dimension_seen = true;
    }
    if (code != 18 && v[0] == 0) {
      if (code == 10) declared_nodes = v[2];
      if (code == 12) declared_cells = v[2];
      if (code == 13) declared_faces = v[2];
      p.rest();
      continue;
    }
    auto first = code == 18 ? v[0] : v[1], last = code == 18 ? v[1] : v[2];
    KUTTA_REQUIRE(first > 0 && last >= first,
                  "Fluent: invalid section id range");
    auto number = [&]() { return binary ? p.binary_index() : p.integer(); };
    if (code == 10) {
      int nd = h.size() > 4 ? v[4] : raw.dimension;
      KUTTA_REQUIRE(nd == raw.dimension, "Fluent: inconsistent node dimension");
      p.expect('(');
      for (std::int64_t id = first; id <= last; ++id) {
        auto real = [&]() {
          return binary ? p.binary_real(section >= 3000) : p.real();
        };
        Vec3d point;
        point.x = real();
        point.y = real();
        if (nd == 3) point.z = real();
        KUTTA_REQUIRE(raw.nodes.emplace(static_cast<Index>(id), point).second,
                      "Fluent: overlapping node sections");
      }
      p.finish(binary, section);
    } else if (code == 12) {
      KUTTA_REQUIRE(v.size() >= 5, "Fluent: cell zone requires element type");
      KUTTA_REQUIRE(v[3] != 32,
                    "Fluent: inactive parent cells require unsupported "
                    "adapted-mesh handling");
      bool body = p.next('(');
      if (body) p.expect('(');
      KUTTA_REQUIRE(v[4] != 0 || body,
                    "Fluent: mixed cells require a type list");
      for (std::int64_t id = first; id <= last; ++id) {
        auto type = v[4] == 0 ? number() : v[4];
        KUTTA_REQUIRE(type >= 1 && type <= 7,
                      "Fluent: invalid cell element type");
        if (v[3] != 0) {
          KUTTA_REQUIRE(raw.cells.emplace(static_cast<Index>(id), v[0]).second,
                        "Fluent: overlapping cell sections");
          cell_types[v[0]] = v[3];
        }
      }
      // Uniform cells have no payload, but binary writers still emit the
      // section trailer. Do not scan raw binary bytes for parentheses.
      p.finish(binary, section, body);
    } else if (code == 13) {
      KUTTA_REQUIRE(v.size() >= 5, "Fluent: face zone requires face type");
      p.expect('(');
      raw.zones.try_emplace(v[0],
                            detail::FluentZone{"zone" + std::to_string(v[0]),
                                               detail::fluent_type(v[3])});
      for (std::int64_t id = first; id <= last; ++id) {
        auto n = (v[4] == 0 || v[4] == 5) ? number() : v[4];
        // Referenced coordinates may be in later sections. Check node ids
        // against the complete node table at assembly, not this partial table.
        KUTTA_REQUIRE(n >= 2, "Fluent: invalid face node count");
        detail::FluentFace f;
        f.zone = v[0];
        for (Index j = 0; j < n; ++j) f.nodes.push_back(number());
        f.c0 = number();
        f.c1 = number();
        KUTTA_REQUIRE(
            raw.faces.emplace(static_cast<Index>(id), std::move(f)).second,
            "Fluent: overlapping face sections");
      }
      p.finish(binary, section);
    } else {
      p.expect('(');
      for (std::int64_t i = first; i <= last; ++i) {
        auto a = number(), b = number();
        raw.periodic.push_back({a, b, v[2], v[3]});
      }
      p.finish(binary, section);
    }
  }
  for (const auto& [z, meta] : named) raw.zones[z] = meta;
  for (auto [z, type] : cell_types)
    if (!raw.zones.count(z))
      raw.zones[z] = {"zone" + std::to_string(z),
                      "fluent-cell-zone-type-" + std::to_string(type)};
  // Declarations count ids, not active rows (dead zones may be omitted).
  if (declared_nodes)
    KUTTA_REQUIRE(
        raw.nodes.empty() || raw.nodes.rbegin()->first <= declared_nodes,
        "Fluent: node exceeds declared range");
  if (declared_cells)
    KUTTA_REQUIRE(
        raw.cells.empty() || raw.cells.rbegin()->first <= declared_cells,
        "Fluent: cell exceeds declared range");
  if (declared_faces)
    KUTTA_REQUIRE(
        raw.faces.empty() || raw.faces.rbegin()->first <= declared_faces,
        "Fluent: face exceeds declared range");
  return raw;
}
}  // namespace
FluentMeshData read_fluent_msh_with_metadata(const std::filesystem::path& path,
                                             double thickness) {
  std::ifstream in(path, std::ios::binary);
  KUTTA_REQUIRE(in, "cannot open Fluent mesh: " + path.string());
  std::string bytes{std::istreambuf_iterator<char>(in), {}};
  return detail::assemble_fluent(parse(std::move(bytes)), thickness);
}
PolyMesh read_fluent_msh(const std::filesystem::path& path, double thickness) {
  return read_fluent_msh_with_metadata(path, thickness).mesh;
}
}  // namespace kutta
