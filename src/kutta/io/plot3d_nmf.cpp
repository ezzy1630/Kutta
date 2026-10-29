#include "kutta/io/plot3d_detail.hpp"

#include <cctype>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <sstream>

namespace kutta::detail {
namespace {
std::string canonical(std::string s) {
  for (auto& c : s) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (c == '-') c = '_';
  }
  return s;
}
PatchKind kind(const std::string& type) {
  if (type == "wall" || type == "viscous_solid" || type == "tangency")
    return PatchKind::Wall;
  if (type.starts_with("symmetry")) return PatchKind::Symmetry;
  if (type.starts_with("farfield")) return PatchKind::Farfield;
  if (type == "subsonic_inflow_pt") return PatchKind::PressureInlet;
  if (type == "inflow") return PatchKind::VelocityInlet;
  if (type == "back_pressure" || type == "outflow") return PatchKind::PressureOutlet;
  return PatchKind::Other;
}
struct Block {
  Plot3dBlock grid;
  std::array<int, 3> dims, axes;
};
struct Surface {
  int block, face;
  std::array<int, 2> start, end;
  int nu, nv;
  std::vector<Index> nodes;
  std::string suffix() const {
    return "_block" + std::to_string(block + 1) + "_f" + std::to_string(face) +
        "_" + std::to_string(start[0]) + "-" + std::to_string(end[0]) +
        "_" + std::to_string(start[1]) + "-" + std::to_string(end[1]);
  }
};
// NMF faces: kmin/kmax (i,j), imin/imax (j,k), jmin/jmax (k,i).
constexpr int fixed_axis[]{2, 2, 0, 0, 1, 1};
Surface surface(std::istream& in, const std::vector<Block>& blocks) {
  Surface s{};
  in >> s.block >> s.face >> s.start[0] >> s.end[0] >> s.start[1] >> s.end[1];
  KUTTA_REQUIRE(in && s.block >= 1 && s.block <= index(blocks.size()) &&
                    s.face >= 1 && s.face <= 6, "NMF: invalid block/face range");
  --s.block;
  const auto& b = blocks[s.block];
  const auto axis = fixed_axis[s.face - 1], u = (axis + 1) % 3, v = (axis + 2) % 3;
  for (int a = 0; a < 2; ++a) {
    const auto dim = b.dims[a == 0 ? u : v];
    KUTTA_REQUIRE(s.start[a] >= 1 && s.start[a] <= dim &&
                      s.end[a] >= 1 && s.end[a] <= dim && s.start[a] != s.end[a],
                  "NMF: point range outside block or has zero extent");
  }
  s.nu = std::abs(s.end[0] - s.start[0]) + 1;
  s.nv = std::abs(s.end[1] - s.start[1]) + 1;
  for (int j = 0; j < s.nv; ++j)
    for (int i = 0; i < s.nu; ++i) {
      std::array<int, 3> p{}, q{};
      p[axis] = s.face % 2 ? 0 : b.dims[axis] - 1;
      p[u] = s.start[0] - 1 + i * (s.end[0] > s.start[0] ? 1 : -1);
      p[v] = s.start[1] - 1 + j * (s.end[1] > s.start[1] ? 1 : -1);
      for (int a = 0; a < 3; ++a) q[b.axes[a]] = p[a];
      const auto& d = b.grid.dims;
      s.nodes.push_back(b.grid.point_offset + (q[2] * d[1] + q[1]) * d[0] + q[0]);
    }
  return s;
}
std::vector<Index> key(std::vector<Index> ns) {
  std::sort(ns.begin(), ns.end());
  return ns;
}
}  // namespace
void apply_neutral_map(const std::filesystem::path& path,
                       const std::vector<Plot3dBlock>& grids,
                       std::vector<Vec3d>& points, std::vector<Face>& faces,
                       std::vector<Patch>& patches) {
  std::ifstream file(path);
  KUTTA_REQUIRE(file, "NMF: cannot open " + path.string());
  int line_number = 0;
  auto line = [&]() {
    std::string text;
    while (std::getline(file, text)) {
      ++line_number;
      text.erase(text.find('#') == std::string::npos ? text.size() : text.find('#'));
      auto end = text.find_last_not_of(" \t\r");
      if (end != std::string::npos && text[end] == '\\') text.erase(end);
      if (text.find_first_not_of(" \t\r") != std::string::npos) return text;
    }
    return std::string{};
  };
  auto end_of_row = [](std::istream& in) {
    std::string extra;
    KUTTA_REQUIRE(!(in >> extra), "NMF: unexpected trailing fields");
  };
  try {
    std::istringstream header(line());
    int count = 0;
    header >> count;
    KUTTA_REQUIRE(header && count == index(grids.size()),
                  "NMF: block count disagrees with grid");
    end_of_row(header);
    std::vector<Block> blocks(grids.size());
    std::vector<bool> seen(grids.size());
    for (int n = 0; n < count; ++n) {
      std::istringstream row(line());
      int id = 0;
      std::array<int, 3> dims{};
      row >> id >> dims[0] >> dims[1] >> dims[2];
      KUTTA_REQUIRE(row && id >= 1 && id <= count && !seen[id - 1],
                    "NMF: invalid/duplicate block dimensions");
      end_of_row(row);
      seen[id - 1] = true;
      auto& b = blocks[id - 1];
      b = {grids[id - 1], dims, {0, 1, 2}};
      if (dims != b.grid.dims) {
        const std::array span_first{2, b.grid.dims[0], b.grid.dims[1]};
        KUTTA_REQUIRE(b.grid.extruded && dims == span_first,
                      "NMF: dimensions disagree with grid");
        b.axes = {2, 0, 1};
      }
    }
    std::map<std::vector<Index>, Index> boundary;
    for (Index f = 0; f < index(faces.size()); ++f)
      if (faces[f].neighbour < 0) boundary.emplace(key(faces[f].nodes), f);
    std::vector<bool> assigned(faces.size()), removed(faces.size());
    std::vector<Index> parent(points.size());
    std::iota(parent.begin(), parent.end(), 0);
    auto root = [&](Index p) {
      while (parent[p] != p) {
        parent[p] = parent[parent[p]];
        p = parent[p];
      }
      return p;
    };
    auto claim = [&](std::vector<Index> ns) {
      auto it = boundary.find(key(std::move(ns)));
      KUTTA_REQUIRE(it != boundary.end(), "NMF: range does not identify a boundary face");
      const auto f = it->second;
      KUTTA_REQUIRE(!assigned[f], "NMF: overlapping boundary/interface ranges");
      assigned[f] = true;
      return f;
    };
    patches.clear();
    std::map<std::string, Index> caps;
    bool has_connection = false;
    for (auto text = line(); !text.empty(); text = line()) {
      std::istringstream row(text);
      row >> std::ws;
      std::string type;
      const char quote = row.peek() == '"' ? '"' : '\'';
      row >> std::quoted(type, quote);
      const auto normalized = canonical(type);
      KUTTA_REQUIRE(!normalized.empty() &&
                        std::all_of(normalized.begin(), normalized.end(),
                                    [](unsigned char c) {
                                      return std::isalnum(c) || c == '_';
                                    }), "NMF: invalid type name");
      auto a = surface(row, blocks);
      const bool connection = normalized == "one_to_one";
      KUTTA_REQUIRE(normalized != "patched" && normalized != "collapsed" &&
                        !normalized.starts_with("pole_") && normalized != "unprocessed",
                    "NMF: unsupported topology type " + type);
      Surface b{};
      bool swap = false;
      Index patch = -1;
      if (connection) {
        b = surface(row, blocks);
        std::string flag;
        row >> flag;
        flag = canonical(flag);
        KUTTA_REQUIRE(flag == "true" || flag == "false",
                      "NMF: swap must be TRUE or FALSE");
        swap = flag == "true";
        KUTTA_REQUIRE(a.nu == (swap ? b.nv : b.nu) && a.nv == (swap ? b.nu : b.nv),
                      "NMF: one-to-one range sizes disagree");
        has_connection = true;
      } else {
        auto name = normalized + a.suffix();
        const auto& block = blocks[a.block];
        if (block.grid.extruded && block.axes[fixed_axis[a.face - 1]] == 2) {
          KUTTA_REQUIRE(kind(normalized) == PatchKind::Symmetry,
                        "NMF: extrusion caps must be symmetry boundaries");
          name = a.face % 2 ? "front" : "back";
          auto [it, added] = caps.emplace(name, index(patches.size()));
          patch = it->second;
          if (!added)
            KUTTA_REQUIRE(patches[patch].source_type == type,
                          "NMF: conflicting cap boundary types");
        }
        if (patch < 0 || patch == index(patches.size())) {
          patch = index(patches.size());
          patches.push_back({name, kind(normalized), type});
        }
      }
      end_of_row(row);
      auto other = [&](int i, int j) {
        return b.nodes[swap ? i * b.nu + j : j * b.nu + i];
      };
      if (connection)
        for (int j = 0; j < a.nv; ++j)
          for (int i = 0; i < a.nu; ++i) {
            auto x = a.nodes[j * a.nu + i], y = other(i, j);
            KUTTA_REQUIRE(norm(points[x] - points[y]) <=
                              1e-12 * std::max({1., norm(points[x]), norm(points[y])}),
                          "NMF: one-to-one points are not coincident");
            x = root(x);
            y = root(y);
            parent[std::max(x, y)] = std::min(x, y);
          }
      for (int j = 0; j < a.nv - 1; ++j)
        for (int i = 0; i < a.nu - 1; ++i) {
          const auto n = j * a.nu + i;
          auto f = claim({a.nodes[n], a.nodes[n + 1],
                          a.nodes[n + a.nu + 1], a.nodes[n + a.nu]});
          if (connection) {
            auto g = claim({other(i, j), other(i + 1, j),
                            other(i + 1, j + 1), other(i, j + 1)});
            KUTTA_REQUIRE(faces[f].owner != faces[g].owner,
                          "NMF: interface joins a cell to itself");
            faces[f].neighbour = faces[g].owner;
            faces[f].patch = -1;
            removed[g] = true;
          } else faces[f].patch = patch;
        }
    }
    for (const auto& [ns, f] : boundary) {
      (void)ns;
      KUTTA_REQUIRE(assigned[f], "NMF: boundary face has no BC/interface range");
    }
    if (has_connection) {
      std::vector<Vec3d> compact;
      std::vector<Index> remap(points.size(), -1);
      for (Index p = 0; p < index(points.size()); ++p)
        if (root(p) == p) {
          remap[p] = index(compact.size());
          compact.push_back(points[p]);
        }
      for (auto& f : faces)
        for (auto& p : f.nodes) p = remap[root(p)];
      points = std::move(compact);
      std::vector<Face> stitched;
      for (Index f = 0; f < index(faces.size()); ++f)
        if (!removed[f]) stitched.push_back(std::move(faces[f]));
      faces = std::move(stitched);
    }
  } catch (const Error& e) {
    throw Error(path.string() + ":" + std::to_string(line_number) + ": " + e.what());
  }
}
}  // namespace kutta::detail
