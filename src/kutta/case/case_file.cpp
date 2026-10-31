#include "kutta/case/fluent_case.hpp"
#include "kutta/case/hdf5.hpp"
#include "kutta/case/legacy.hpp"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <limits>
#include <sstream>

namespace kutta::fluent {
namespace {
void validate_topology(const std::vector<ZoneRange> &zs, std::optional<std::uint64_t> count) {
  std::uint64_t n = 0;
  std::set<std::int64_t> ids;
  for (std::size_t i = 0; i < zs.size(); ++i) {
    const auto &z = zs[i];
    if (!z.first || z.last < z.first || !ids.insert(z.zone_id).second ||
        (i && zs[i - 1].last >= z.first))
      throw std::runtime_error("invalid, duplicate or overlapping zone intervals");
    if (z.size() > std::numeric_limits<std::uint64_t>::max() - n)
      throw std::runtime_error("zone count overflow");
    n += z.size();
  }
  if (count && n != *count)
    throw std::runtime_error("mesh count disagrees with zone topology");
}
void check_thread_topology(FluentCase &c) {
  auto check = [&](const auto &ranges, const auto &threads, bool cell) {
    for (const auto &range : ranges) {
      const auto prefix = "mesh/zone/" + std::to_string(range.zone_id);
      auto thread = std::find_if(threads.begin(), threads.end(),
          [&](const auto &z) { return z.id == range.zone_id; });
      if (thread == threads.end()) {
        c.issues.push_back({prefix, "required: topology zone has no Thread Variables record"});
        continue;
      }
      if (range.name != thread->name)
        c.issues.push_back({prefix + "/name", "unsupported: topology and thread names disagree"});
      std::optional<int> expected;
      if (cell && thread->type == ZoneType::fluid)
        expected = 1;
      if (!cell) {
        switch (thread->type) {
        case ZoneType::interior: expected = 2; break;
        case ZoneType::wall: expected = 3; break;
        case ZoneType::pressure_outlet: expected = 5; break;
        case ZoneType::symmetry: expected = 7; break;
        case ZoneType::velocity_inlet: expected = 10; break;
        default: break;
        }
      }
      if (!expected || !range.type_code || range.type_code != expected)
        c.issues.push_back({prefix + "/type", "unsupported: topology and thread types disagree or are unverified"});
    }
  };
  check(c.mesh.cell_zones, c.cell_zones, true);
  check(c.mesh.face_zones, c.boundaries, false);
}
using legacy::balanced_end;
// Read only counts and ID intervals from mesh section headers. Bodies stay opaque;
// geometry/connectivity and PolyMesh construction belong exclusively to lane B.
void legacy_zone_header(std::string_view s, std::size_t begin, int section, MeshMetadata &mesh) {
  if (section >= 2000)
    section %= 1000;
  if (section != 10 && section != 12 && section != 13)
    return;
  while (begin < s.size() && std::isspace(static_cast<unsigned char>(s[begin])))
    ++begin;
  if (begin == s.size() || s[begin] != '(')
    throw std::runtime_error("missing legacy zone header");
  const auto end = balanced_end(s, begin);
  auto header = s.substr(begin + 1, end - begin - 2);
  std::vector<std::uint64_t> values;
  std::size_t p = 0;
  while (p < header.size()) {
    while (p < header.size() && std::isspace(static_cast<unsigned char>(header[p])))
      ++p;
    if (p == header.size())
      break;
    std::uint64_t value = 0;
    auto result = std::from_chars(header.data() + p, header.data() + header.size(), value, 16);
    if (result.ec != std::errc{} || result.ptr == header.data() + p)
      throw std::runtime_error("invalid hexadecimal zone header");
    values.push_back(value);
    p = static_cast<std::size_t>(result.ptr - header.data());
    if (values.size() > 8)
      throw std::runtime_error("oversized zone header");
  }
  if (values.size() < 4 || !values[1] || values[2] < values[1] ||
      values[0] > std::numeric_limits<int>::max() || values[3] > std::numeric_limits<int>::max())
    throw std::runtime_error("invalid legacy zone interval");
  if (values[0] == 0) {
    auto count = values[2] - values[1] + 1;
    if (section == 10)
      mesh.nodes = count;
    if (section == 12)
      mesh.cells = count;
    if (section == 13)
      mesh.faces = count;
  } else if (section == 12 || section == 13) {
    auto &zones = section == 12 ? mesh.cell_zones : mesh.face_zones;
    zones.push_back({static_cast<std::int64_t>(values[0]),
                     {},
                     values[1],
                     values[2],
                     static_cast<int>(values[3])});
  }
}
FluentCase legacy_case(const std::filesystem::path &path) {
  auto bytes = legacy::read_file(path, 512u * 1024 * 1024, "case");
  std::string_view s = bytes;
  std::map<std::string, std::vector<SExpr>> sections;
  std::optional<std::string> version, origin;
  std::optional<int> dimension;
  MeshMetadata metadata;
  std::size_t p = 0;
  while (p < s.size()) {
    while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p])))
      ++p;
    if (p == s.size())
      break;
    if (s[p] != '(')
      throw std::runtime_error("expected legacy section at byte " + std::to_string(p));
    auto begin = p++;
    while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p])))
      ++p;
    auto start = p;
    while (p < s.size() && std::isdigit(static_cast<unsigned char>(s[p])))
      ++p;
    int section = 0;
    auto conv = std::from_chars(s.data() + start, s.data() + p, section);
    if (conv.ec != std::errc{} || p == start)
      throw std::runtime_error("invalid legacy section number");
    legacy_zone_header(s, p, section, metadata);
    if (section >= 2000 && section < 4000) {
      p = legacy::skip_binary(s, p, section);
      continue;
    }
    p = balanced_end(s, begin);
    if (section != 0 && section != 1 && section != 2 && section != 37 && section != 38 &&
        section != 39 && section != 64)
      continue;
    auto e = parse_sexpr(s.substr(begin, p - begin));
    if (section == 37)
      sections["Rampant Variables"].push_back(std::move(e));
    else if (section == 38)
      sections["Cortex Variables"].push_back(std::move(e));
    else if (section == 64)
      sections["Domain Variables"].push_back(std::move(e));
    else if (section == 39)
      sections["Thread Variables"].push_back(std::move(e));
    else if (section == 2 && e.elements.size() == 2 && e.elements[1].kind == SExpr::Kind::integer)
      dimension = static_cast<int>(e.elements[1].integer);
    else if ((section == 0 || section == 1) && e.elements.size() > 1) {
      auto t = e.elements[1].atom_text();
      if (t && t->starts_with("fluent")) {
        auto end = t->find_first_not_of("0123456789.", 6);
        auto v = t->substr(6, end == t->npos ? end : end - 6);
        // CFF records major.minor; the legacy writer includes a patch component.
        auto dot = v.find('.');
        if (dot != v.npos) {
          auto second = v.find('.', dot + 1);
          v = v.substr(0, second);
        }
        version = std::string(v);
        origin = std::string(*t);
      }
    }
  }
  auto result = extract_case(std::move(sections), version, origin);
  if (dimension && result.dimension && dimension != result.dimension)
    throw std::runtime_error("legacy dimension disagrees with case-config");
  metadata.dimension = dimension;
  for (auto &z : metadata.cell_zones)
    for (auto &cell : result.cell_zones)
      if (z.zone_id == cell.id)
        z.name = cell.name;
  for (auto &z : metadata.face_zones)
    for (auto &boundary : result.boundaries)
      if (z.zone_id == boundary.id)
        z.name = boundary.name;
  auto order = [](auto &a, auto &b) { return a.first < b.first; };
  std::sort(metadata.cell_zones.begin(), metadata.cell_zones.end(), order);
  std::sort(metadata.face_zones.begin(), metadata.face_zones.end(), order);
  validate_topology(metadata.cell_zones, metadata.cells);
  validate_topology(metadata.face_zones, metadata.faces);
  result.mesh = std::move(metadata);
  return result;
}
#if KUTTA_HAS_HDF5
std::vector<std::string> split_names(std::string s) {
  std::vector<std::string> result;
  std::size_t p = 0;
  while (p < s.size()) {
    auto end = s.find(';', p);
    result.push_back(s.substr(p, end == s.npos ? end : end - p));
    if (end == s.npos)
      break;
    p = end + 1;
  }
  return result;
}
std::vector<ZoneRange> topology(hid_t f, const std::string &path) {
  if (!h5::exists(f, path))
    return {};
  auto ids = h5::numbers<std::int64_t>(f, path + "/id");
  auto first = h5::numbers<std::uint64_t>(f, path + "/minId"),
       last = h5::numbers<std::uint64_t>(f, path + "/maxId");
  auto names = split_names(h5::string(f, path + "/name"));
  auto types = h5::numbers<std::int64_t>(f, path + "/zoneType");
  if (ids.size() != first.size() || ids.size() != last.size() || ids.size() != names.size() ||
      ids.size() != types.size())
    throw std::runtime_error("inconsistent zone topology lengths");
  std::vector<ZoneRange> zones;
  std::set<std::int64_t> seen;
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] <= 0 || !first[i] || last[i] < first[i] || !seen.insert(ids[i]).second ||
        types[i] < 0 || types[i] > std::numeric_limits<int>::max())
      throw std::runtime_error("invalid zone topology range/ID");
    zones.push_back({ids[i], names[i], first[i], last[i], static_cast<int>(types[i])});
  }
  std::sort(zones.begin(), zones.end(), [](auto &a, auto &b) { return a.first < b.first; });
  for (std::size_t i = 1; i < zones.size(); ++i)
    if (zones[i - 1].last >= zones[i].first)
      throw std::runtime_error("overlapping zone topology ranges");
  return zones;
}

FluentCase cff(const std::filesystem::path &path) {
  auto file = h5::open(path);
  std::map<std::string, std::vector<SExpr>> sections;
  std::optional<std::string> version, origin;
  for (auto &name : h5::children(file, "/settings")) {
    auto s = h5::string(file, "/settings/" + name);
    if (name == "Version")
      version = s;
    else if (name == "Origin")
      origin = s;
    else if (name == "Rampant Variables" || name == "Thread Variables" ||
             name == "Domain Variables" || name == "Cortex Variables" ||
             (!s.empty() && s.front() == '('))
      sections.emplace(name, parse_sexprs(s));
    else {
      SExpr raw;
      raw.kind = SExpr::Kind::string;
      raw.text = s;
      sections[name].push_back(std::move(raw));
    }
  }
  auto c = extract_case(std::move(sections), version, origin);
  const std::string mesh = "/meshes/1";
  if (h5::exists(file, mesh)) {
    auto dim = h5::attribute(file, mesh, "dimension");
    if (dim && (*dim == 2 || *dim == 3))
      c.mesh.dimension = static_cast<int>(*dim);
    else
      throw std::runtime_error("missing or invalid mesh dimension");
    if (c.dimension && c.dimension != c.mesh.dimension)
      throw std::runtime_error("mesh dimension disagrees with case-config");
    c.mesh.cells = h5::attribute(file, mesh, "cellCount");
    c.mesh.faces = h5::attribute(file, mesh, "faceCount");
    c.mesh.nodes = h5::attribute(file, mesh, "nodeCount");
    c.mesh.cell_zones = topology(file, mesh + "/cells/zoneTopology");
    c.mesh.face_zones = topology(file, mesh + "/faces/zoneTopology");
    validate_topology(c.mesh.cell_zones, c.mesh.cells);
    validate_topology(c.mesh.face_zones, c.mesh.faces);
  } else
    c.issues.push_back({"mesh", "absent: no topology available"});
  return c;
}
#endif
} // namespace
FluentCase read_case(const std::filesystem::path &path) {
  auto c = [&]() {
    if (path.extension() == ".h5") {
#if KUTTA_HAS_HDF5
      return cff(path);
#else
      throw std::runtime_error(
          "HDF5 support is not built; configure with an available HDF5 C library");
#endif
    }
    return legacy_case(path);
  }();
  c.source_path = std::filesystem::absolute(path).lexically_normal();
  check_thread_topology(c);
  return c;
}
} // namespace kutta::fluent
