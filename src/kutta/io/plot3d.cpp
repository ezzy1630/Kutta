#include <cctype>
#include <fstream>
#include <locale>
#include <sstream>

#include "kutta/io/mesh_io.hpp"
#include "kutta/io/plot3d_detail.hpp"

namespace kutta {
namespace {
std::vector<int> integers(const std::string& line) {
  std::istringstream s(line);
  std::vector<int> v;
  int x;
  while (s >> x) v.push_back(x);
  return v;
}
}  // namespace
PolyMesh read_plot3d(const std::filesystem::path& path, double thickness) {
  auto nmf = path;
  nmf.replace_extension(".nmf");
  return read_plot3d(
      path, std::filesystem::exists(nmf) ? nmf : std::filesystem::path{}, thickness);
}
PolyMesh read_plot3d(const std::filesystem::path& path,
                     const std::filesystem::path& nmf, double thickness) {
  auto extension = path.extension().string();
  std::transform(
      extension.begin(), extension.end(), extension.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  KUTTA_REQUIRE(
      extension != ".gz",
      "Plot3D gzip input is not supported: gunzip the file before reading");
  std::ifstream in(path);
  KUTTA_REQUIRE(in, "cannot open Plot3D file: " + path.string());
  std::string first;
  do {
    KUTTA_REQUIRE(static_cast<bool>(std::getline(in, first)),
                  "Plot3D: missing header");
  } while (integers(first).empty());
  auto header = integers(first);
  int blocks = 1, nd = 0;
  std::string remainder;
  if (header.size() == 1) {
    blocks = header[0];
    KUTTA_REQUIRE(blocks > 0, "Plot3D: block count must be positive");
    std::string dims;
    do {
      KUTTA_REQUIRE(static_cast<bool>(std::getline(in, dims)),
                    "Plot3D: missing dimensions");
    } while (integers(dims).empty());
    auto row = integers(dims);
    const auto n = static_cast<std::size_t>(blocks);
    if (row.size() == 2 * n) nd = 2;
    else if (row.size() == 3 * n) nd = 3;
    else nd = extension == ".p2dfmt" || row.size() == 2 ? 2 : 3;
    remainder = dims + "\n";
  } else {
    KUTTA_REQUIRE(header.size() == 2 || header.size() == 3,
                  "Plot3D: expected block count or 2/3 grid dimensions");
    nd = static_cast<int>(header.size());
    remainder = first + "\n";
  }
  remainder.append(std::istreambuf_iterator<char>(in), {});
  for (auto& c : remainder)
    if (c == 'D' || c == 'd') c = 'E';
  std::istringstream data(remainder);
  data.imbue(std::locale::classic());
  struct Size {
    int ni = 0, nj = 0, nk = 1;
  };
  std::vector<Size> sizes(blocks);
  for (auto& s : sizes) {
    data >> s.ni >> s.nj;
    if (nd == 3) data >> s.nk;
    KUTTA_REQUIRE(data && s.ni >= 2 && s.nj >= 2 && s.nk >= 1,
                  "Plot3D: invalid/truncated grid dimensions");
    const auto cap =
        static_cast<std::uint64_t>(std::numeric_limits<Index>::max()) / 12;
    auto n = static_cast<std::uint64_t>(s.ni);
    KUTTA_REQUIRE(n <= cap && static_cast<std::uint64_t>(s.nj) <= cap / n,
                  "Plot3D: grid exceeds 32-bit mesh capacity");
    n *= s.nj;
    KUTTA_REQUIRE(static_cast<std::uint64_t>(s.nk) <= cap / n,
                  "Plot3D: grid exceeds 32-bit mesh capacity");
  }
  std::vector<Vec3d> all_points;
  std::vector<detail::Face> all_faces;
  std::vector<Patch> all_patches;
  Index all_cells = 0;
  std::vector<detail::Plot3dBlock> layouts;
  for (int block = 0; block < blocks; ++block) {
    auto [ni, nj, nk] = sizes[block];
    auto np = static_cast<std::size_t>(ni) * nj * nk;
    std::vector<Vec3d> points(np);
    for (int d = 0; d < nd; ++d)
      for (auto& p : points) {
        double value = 0;
        data >> value;
        KUTTA_REQUIRE(data && std::isfinite(value),
                      "Plot3D: truncated/non-finite coordinates in block " +
                          std::to_string(block + 1));
        if (d == 0)
          p.x = value;
        else if (d == 1)
          p.y = value;
        else
          p.z = value;
      }
    bool planar = nk == 1;
    if (planar) {
      KUTTA_REQUIRE(std::isfinite(thickness) && thickness > 0,
                    "Plot3D: extrusion thickness must be positive and finite");
      // One-plane 3D inputs must be xy planes; arbitrary plane extrusion would
      // silently choose an unintended physical direction.
      for (auto& p : points)
        KUTTA_REQUIRE(p.z == points[0].z,
                      "Plot3D: single k-plane must be parallel to xy");
      for (std::size_t i = 0; i < np; ++i) {
        auto p = points[i];
        p.z += thickness;
        points.push_back(p);
      }
      nk = 2;
    }
    auto node = [=](int i, int j, int k) {
      return static_cast<Index>((static_cast<std::int64_t>(k) * nj + j) * ni +
                                i);
    };
    std::vector<detail::CellFaces> cells;
    for (int k = 0; k < nk - 1; ++k)
      for (int j = 0; j < nj - 1; ++j)
        for (int i = 0; i < ni - 1; ++i) {
          auto a = node(i, j, k), b = node(i + 1, j, k),
               c = node(i + 1, j + 1, k), d = node(i, j + 1, k);
          auto e = node(i, j, k + 1), f = node(i + 1, j, k + 1),
               g = node(i + 1, j + 1, k + 1), h = node(i, j + 1, k + 1);
          cells.push_back({{a, b, c, d},
                           {e, h, g, f},
                           {a, e, f, b},
                           {b, f, g, c},
                           {c, g, h, d},
                           {d, h, e, a}});
        }
    std::vector<Patch> patches;
    auto prefix = "block" + std::to_string(block + 1) + "_";
    for (auto name : {"imin", "imax", "jmin", "jmax", "kmin", "kmax"})
      patches.push_back({prefix + name, PatchKind::Other, "plot3d"});
    if (planar) {
      patches[4] = {"front", PatchKind::Symmetry, "symmetry"};
      patches[5] = {"back", PatchKind::Symmetry, "symmetry"};
    }
    auto classify = [=](const std::vector<Index>& ns) -> Index {
      for (Index s = 0; s < 6; ++s) {
        bool all = true;
        for (auto v : ns) {
          auto i = v % ni, j = (v / ni) % nj, k = v / (ni * nj);
          all = all && (s == 0   ? i == 0
                        : s == 1 ? i == ni - 1
                        : s == 2 ? j == 0
                        : s == 3 ? j == nj - 1
                        : s == 4 ? k == 0
                                 : k == nk - 1);
        }
        if (all) return s;
      }
      throw Error("Plot3D: cannot classify block boundary");
    };
    auto m = detail::from_cells(std::move(points), cells, std::move(patches),
                                classify);
    auto po = detail::index(all_points.size());
    detail::index(all_points.size() + m.points.size());
    layouts.push_back({{ni, nj, nk}, po, planar});
    std::vector<Index> patch_ids;
    for (const auto& patch : m.patches) {
      auto found = std::find_if(all_patches.begin(), all_patches.end(),
                               [&](const Patch& p) { return p.name == patch.name; });
      patch_ids.push_back(static_cast<Index>(found - all_patches.begin()));
      if (found == all_patches.end()) all_patches.push_back(patch);
    }
    all_points.insert(all_points.end(), m.points.begin(), m.points.end());
    for (Index f = 0; f < m.n_faces(); ++f) {
      auto ns = detail::nodes(m, f);
      detail::Face face{
          {ns.begin(), ns.end()},
          m.owner[f] + all_cells,
          f < m.n_internal_faces() ? m.neighbour[f] + all_cells : -1,
          -1};
      for (auto& v : face.nodes) v += po;
      if (f >= m.n_internal_faces())
        for (Index p = 0; p < detail::index(m.patches.size()); ++p)
          if (f >= m.patches[p].start &&
              f < m.patches[p].start + m.patches[p].size) {
            face.patch = patch_ids[p];
            break;
          }
      all_faces.push_back(std::move(face));
    }
    all_cells = detail::index(static_cast<std::size_t>(all_cells) + m.n_cells);
  }
  std::string extra;
  KUTTA_REQUIRE(!(data >> extra),
                "Plot3D: unexpected trailing data (iblank/function files are "
                "unsupported)");
  if (!nmf.empty())
    detail::apply_neutral_map(nmf, layouts, all_points, all_faces, all_patches);
  auto m = detail::pack(std::move(all_points), all_cells, std::move(all_faces),
                        std::move(all_patches));
  validate(m);
  return m;
}
}  // namespace kutta
