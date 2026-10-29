#include <fstream>
#include <iomanip>
#include <locale>
#include <set>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta {
namespace {
std::string xml(const std::string& s) {
  std::string out;
  for (auto c : s) switch (c) {
      case '&':
        out += "&amp;";
        break;
      case '<':
        out += "&lt;";
        break;
      case '>':
        out += "&gt;";
        break;
      case '"':
        out += "&quot;";
        break;
      case '\'':
        out += "&apos;";
        break;
      default:
        KUTTA_REQUIRE(static_cast<unsigned char>(c) >= 32 || c == '\n' ||
                          c == '\t' || c == '\r',
                      "VTU name contains invalid XML control character");
        out += c;
    }
  return out;
}
// Recognize topological tet/hex connectivity, rather than assuming that sorted
// point ids happen to have the local ordering VTK requires.
std::vector<Index> fast_cell(const std::vector<std::vector<Index>>& faces,
                             const std::set<Index>& vertices) {
  if (faces.size() == 4 && vertices.size() == 4 &&
      std::all_of(faces.begin(), faces.end(),
                  [](const auto& f) { return f.size() == 3; })) {
    auto base = faces[0];
    std::reverse(base.begin(), base.end());
    for (auto v : vertices)
      if (std::find(base.begin(), base.end(), v) == base.end()) {
        base.push_back(v);
        break;
      }
    return base;
  }
  if (faces.size() != 6 || vertices.size() != 8 ||
      !std::all_of(faces.begin(), faces.end(),
                   [](const auto& f) { return f.size() == 4; }))
    return {};
  auto base = faces[0];
  std::reverse(base.begin(), base.end());
  std::vector<Index> top;
  for (auto v : base) {
    std::set<Index> connected;
    for (const auto& face : faces)
      for (std::size_t i = 0; i < 4; ++i) {
        auto a = face[i], b = face[(i + 1) % 4];
        if (b == v) std::swap(a, b);
        if (a == v && std::find(base.begin(), base.end(), b) == base.end())
          connected.insert(b);
      }
    if (connected.size() != 1) return {};
    top.push_back(*connected.begin());
  }
  if (std::set<Index>(top.begin(), top.end()).size() != 4) return {};
  base.insert(base.end(), top.begin(), top.end());
  return base;
}
template <class T>
void array(std::ostream& out, const char* type, const char* name,
           const std::vector<T>& values) {
  out << "<DataArray type=\"" << type << "\" Name=\"" << name
      << "\" format=\"ascii\">\n";
  for (auto v : values) out << v << ' ';
  out << "\n</DataArray>\n";
}
}  // namespace
void write_vtu(const std::filesystem::path& path, const PolyMesh& m,
               std::span<const CellField> fields) {
  validate(m);
  for (const auto& field : fields) {
    KUTTA_REQUIRE(field.components == 1 || field.components == 3,
                  "VTU field '" + field.name + "' requires 1 or 3 components");
    KUTTA_REQUIRE(
        field.values.size() ==
            static_cast<std::size_t>(m.n_cells) * field.components,
        "VTU field '" + field.name + "' size does not match cell count");
    xml(field.name);
  }
  auto a = build_cell_faces(m);
  std::vector<Index> connectivity, offsets, types, faces, faceoffsets;
  bool poly = false;
  for (Index c = 0; c < m.n_cells; ++c) {
    std::vector<std::vector<Index>> cell_faces;
    std::set<Index> vertices;
    for (auto i = a.offsets[c]; i < a.offsets[c + 1]; ++i) {
      auto ns = detail::nodes(m, a.faces[i]);
      std::vector<Index> f(ns.begin(), ns.end());
      if (a.signs[i] < 0) std::reverse(f.begin(), f.end());
      vertices.insert(f.begin(), f.end());
      cell_faces.push_back(std::move(f));
    }
    auto fast = fast_cell(cell_faces, vertices);
    if (!fast.empty()) {
      types.push_back(fast.size() == 4 ? 10 : 12);
      connectivity.insert(connectivity.end(), fast.begin(), fast.end());
      faceoffsets.push_back(-1);
    } else {
      poly = true;
      types.push_back(42);
      connectivity.insert(connectivity.end(), vertices.begin(), vertices.end());
      faces.push_back(detail::index(cell_faces.size()));
      for (const auto& f : cell_faces) {
        faces.push_back(detail::index(f.size()));
        faces.insert(faces.end(), f.begin(), f.end());
      }
      faceoffsets.push_back(detail::index(faces.size()));
    }
    offsets.push_back(detail::index(connectivity.size()));
  }
  std::ofstream out(path);
  KUTTA_REQUIRE(out, "cannot open VTU output: " + path.string());
  out.imbue(std::locale::classic());
  out << std::setprecision(17);
  out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" "
         "version=\"1.0\" "
         "byte_order=\"LittleEndian\">\n<UnstructuredGrid><Piece "
         "NumberOfPoints=\""
      << m.n_points() << "\" NumberOfCells=\"" << m.n_cells
      << "\">\n<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" "
         "format=\"ascii\">\n";
  for (auto p : m.points) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
  out << "</DataArray></Points>\n<Cells>\n";
  array(out, "Int32", "connectivity", connectivity);
  array(out, "Int32", "offsets", offsets);
  array(out, "UInt8", "types", types);
  if (poly) {
    array(out, "Int32", "faces", faces);
    array(out, "Int32", "faceoffsets", faceoffsets);
  }
  out << "</Cells><CellData>\n";
  for (const auto& f : fields) {
    out << "<DataArray type=\"Float64\" Name=\"" << xml(f.name)
        << "\" NumberOfComponents=\"" << f.components
        << "\" format=\"ascii\">\n";
    for (auto v : f.values) out << v << ' ';
    out << "\n</DataArray>\n";
  }
  out << "</CellData></Piece></UnstructuredGrid></VTKFile>\n";
  out.close();
  KUTTA_REQUIRE(out, "failed writing VTU: " + path.string());
}
}  // namespace kutta
