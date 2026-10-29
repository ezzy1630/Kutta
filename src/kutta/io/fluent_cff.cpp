#include "kutta/io/fluent_detail.hpp"
#if KUTTA_HAS_HDF5
#include <hdf5.h>
#include <zlib.h>

#include <cstring>
#include <set>
#include <sstream>

namespace kutta {
namespace {
// Small owning handles keep every HDF5 resource exception-safe. No HDF5 types
// escape this translation unit or become part of the mesh contract.
class H5 {
 public:
  H5(hid_t id, herr_t (*close)(hid_t), std::string context)
      : id_(id), close_(close) {
    KUTTA_REQUIRE(id >= 0, "CFF: cannot open " + context);
  }
  ~H5() { close_(id_); }
  H5(const H5&) = delete;
  H5& operator=(const H5&) = delete;
  operator hid_t() const { return id_; }

 private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};
void check(herr_t status, const std::string& context) {
  KUTTA_REQUIRE(status >= 0, "CFF: HDF5 failure reading " + context);
}
bool exists(hid_t file, const std::string& path) {
  htri_t value = 0;
  H5E_BEGIN_TRY { value = H5Lexists(file, path.c_str(), H5P_DEFAULT); }
  H5E_END_TRY;
  return value > 0;
}
std::int64_t attribute(hid_t object, const std::string& name) {
  KUTTA_REQUIRE(H5Aexists(object, name.c_str()) > 0,
                "CFF: missing attribute " + name);
  H5 a(H5Aopen(object, name.c_str(), H5P_DEFAULT), H5Aclose, name);
  H5 s(H5Aget_space(a), H5Sclose, name + " dataspace");
  KUTTA_REQUIRE(H5Sget_simple_extent_npoints(s) == 1,
                "CFF: attribute " + name + " must be scalar");
  std::int64_t v = 0;
  check(H5Aread(a, H5T_NATIVE_INT64, &v), name);
  return v;
}
Index id(std::int64_t n, bool zero = false) {
  KUTTA_REQUIRE(n >= (zero ? 0 : 1) && n <= std::numeric_limits<Index>::max(),
                "CFF: invalid/out-of-range id " + std::to_string(n));
  return static_cast<Index>(n);
}
std::pair<Index, Index> range(hid_t object) {
  auto first = id(attribute(object, "minId")),
       last = id(attribute(object, "maxId"));
  KUTTA_REQUIRE(last >= first, "CFF: reversed id range");
  return {first, last};
}
std::size_t count(hid_t dataset) {
  H5 s(H5Dget_space(dataset), H5Sclose, "dataset space");
  auto n = H5Sget_simple_extent_npoints(s);
  KUTTA_REQUIRE(n >= 0 && static_cast<std::uint64_t>(n) <=
                              static_cast<std::uint64_t>(
                                  std::numeric_limits<Index>::max()),
                "CFF: dataset exceeds mesh capacity");
  return static_cast<std::size_t>(n);
}
herr_t partial_chunk_error(unsigned, const H5E_error2_t* error, void* found) {
  if (error->desc && std::strstr(error->desc,
                                "chunk size is incorrect after being unfiltered"))
    *static_cast<bool*>(found) = true;
  return 0;
}
// HDF5 2.2 rejects older Fluent writers' unpadded final chunks (#6691).
// Only that diagnostic permits raw reads; every stream must still contain
// exactly its full chunk or valid final rows. Never pad missing mesh data.
template <class T>
void read_legacy_chunks(hid_t d, hid_t memory_type, std::vector<T>& out) {
  unsigned major = 0, minor = 0, release = 0;
  check(H5get_libversion(&major, &minor, &release), "library version");
  KUTTA_REQUIRE(major == 2 && minor == 2,
                "CFF: legacy chunk fallback requires HDF5 2.2");
  H5 space(H5Dget_space(d), H5Sclose, "legacy chunk space");
  auto rank = H5Sget_simple_extent_ndims(space);
  KUTTA_REQUIRE(rank == 1 || rank == 2, "CFF: unsupported legacy chunk rank");
  hsize_t dims[2]{1, 1}, chunk[2]{1, 1};
  check(H5Sget_simple_extent_dims(space, dims, nullptr), "chunk dimensions");
  H5 props(H5Dget_create_plist(d), H5Pclose, "chunk properties");
  KUTTA_REQUIRE(H5Pget_layout(props) == H5D_CHUNKED &&
                    H5Pget_chunk(props, rank, chunk) == rank && chunk[0] > 0 &&
                    (rank == 1 || chunk[1] == dims[1]),
                "CFF: unsupported legacy chunk layout");
  KUTTA_REQUIRE(H5Pget_nfilters(props) == 1 &&
                    H5Pget_filter2(props, 0, nullptr, nullptr, nullptr, 0,
                                   nullptr, nullptr) == H5Z_FILTER_DEFLATE,
                "CFF: unsupported legacy chunk filters");
  H5 stored_type(H5Dget_type(d), H5Tclose, "chunk numeric type");
  auto width = H5Tget_size(stored_type);
  KUTTA_REQUIRE((H5Tget_class(stored_type) == H5T_INTEGER ||
                 H5Tget_class(stored_type) == H5T_FLOAT) &&
                    width > 0 && width <= 8 && dims[1] > 0 &&
                    chunk[0] <= static_cast<hsize_t>(
                                    std::numeric_limits<Index>::max()) / dims[1],
                "CFF: invalid legacy numeric chunk");
  KUTTA_REQUIRE(dims[0] > 0 && dims[0] % chunk[0] != 0 &&
                    dims[0] <= out.size() / dims[1] &&
                    dims[0] * dims[1] == out.size(),
                "CFF: legacy fallback requires partial final rows within output");
  auto elements = static_cast<std::size_t>(chunk[0] * dims[1]);
  KUTTA_REQUIRE(elements <= std::numeric_limits<std::size_t>::max() /
                               std::max(width, sizeof(T)) &&
                    elements <= std::numeric_limits<uLong>::max() / width,
                "CFF: legacy chunk exceeds buffer capacity");
  auto full_bytes = elements * width;
  bool partial_final = false;
  for (hsize_t row = 0; row < dims[0]; row += chunk[0]) {
    hsize_t offset[2]{row, 0}, stored_bytes = 0;
    check(H5Dget_chunk_storage_size(d, offset, &stored_bytes), "chunk size");
    KUTTA_REQUIRE(stored_bytes > 0 && stored_bytes <= compressBound(full_bytes),
                  "CFF: invalid legacy compressed chunk size");
    std::vector<unsigned char> compressed(
        static_cast<std::size_t>(stored_bytes));
    std::uint32_t mask = 0;
#if H5_VERSION_GE(2, 0, 0)
    auto buffer_size = compressed.size();
    check(H5Dread_chunk2(d, H5P_DEFAULT, offset, &mask, compressed.data(),
                         &buffer_size), "raw legacy chunk");
    KUTTA_REQUIRE(buffer_size == compressed.size(), "CFF: changed chunk size");
#else
    check(H5Dread_chunk(d, H5P_DEFAULT, offset, &mask, compressed.data()),
          "raw legacy chunk");
#endif
    KUTTA_REQUIRE(mask == 0, "CFF: legacy deflate filter was skipped");
    // H5Tconvert operates in place and needs room for either numeric type.
    std::vector<unsigned char> decoded(elements * std::max(width, sizeof(T)));
    auto decoded_bytes = static_cast<uLongf>(full_bytes);
    auto consumed = static_cast<uLong>(compressed.size());
    auto status = uncompress2(decoded.data(), &decoded_bytes, compressed.data(),
                              &consumed);
    auto valid = static_cast<std::size_t>(
        std::min(chunk[0], dims[0] - row) * dims[1]);
    KUTTA_REQUIRE(status == Z_OK && consumed == compressed.size() &&
                      (decoded_bytes == full_bytes ||
                       decoded_bytes == valid * width),
                  "CFF: invalid/truncated legacy deflate chunk");
    partial_final = partial_final || decoded_bytes < full_bytes;
    check(H5Tconvert(stored_type, memory_type, valid, decoded.data(), nullptr,
                     H5P_DEFAULT), "legacy numeric conversion");
    std::memcpy(out.data() + static_cast<std::size_t>(row * dims[1]),
                decoded.data(), valid * sizeof(T));
  }
  KUTTA_REQUIRE(partial_final, "CFF: no legacy partial deflate chunk found");
}
template <class T>
std::vector<T> read(hid_t file, const std::string& path, hid_t type) {
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose, path);
  std::vector<T> out(count(d));
  if (!out.empty()) {
    herr_t status = 0;
    bool partial_error = false;
    herr_t walk_status = 0;
    H5E_BEGIN_TRY {
      status = H5Dread(d, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data());
      if (status < 0)
        walk_status = H5Ewalk2(H5E_DEFAULT, H5E_WALK_DOWNWARD,
                              partial_chunk_error, &partial_error);
    }
    H5E_END_TRY;
    check(walk_status, "error stack");
    if (status < 0) {
      KUTTA_REQUIRE(partial_error, "CFF: HDF5 failure reading " + path);
      try {
        read_legacy_chunks(d, type, out);
      } catch (const Error& error) {
        throw Error("CFF: cannot read " + path + ": " + error.what());
      }
    }
  }
  return out;
}
std::vector<std::int64_t> ints(hid_t file, const std::string& path) {
  return read<std::int64_t>(file, path, H5T_NATIVE_INT64);
}
std::vector<std::string> names(hid_t file, const std::string& path,
                               std::size_t n) {
  if (!exists(file, path)) return std::vector<std::string>(n);
  H5 d(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose, path);
  H5 t(H5Dget_type(d), H5Tclose, path + " type");
  KUTTA_REQUIRE(H5Tget_class(t) == H5T_STRING,
                "CFF: zone names must be HDF5 strings");
  auto size = count(d);
  std::vector<std::string> out;
  if (H5Tis_variable_str(t) > 0) {
    std::vector<char*> values(size, nullptr);
    check(H5Dread(d, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()), path);
    for (auto p : values) {
      out.emplace_back(p ? p : "");
      if (p) H5free_memory(p);
    }
  } else {
    auto width = H5Tget_size(t);
    KUTTA_REQUIRE(
        width > 0 && width <= std::numeric_limits<std::size_t>::max() /
                                  std::max<std::size_t>(size, 1),
        "CFF: invalid string width");
    std::vector<char> bytes(size * width);
    check(H5Dread(d, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, bytes.data()), path);
    for (std::size_t i = 0; i < size; ++i) {
      auto begin = bytes.begin() + i * width, end = begin + width;
      out.emplace_back(begin, std::find(begin, end, '\0'));
    }
  }
  if (out.size() == 1 && out[0].find(';') != std::string::npos) {
    std::vector<std::string> split;
    std::istringstream in(out[0]);
    std::string s;
    while (std::getline(in, s, ';')) split.push_back(s);
    out = std::move(split);
  }
  KUTTA_REQUIRE(out.size() == n,
                "CFF: zone name count does not match topology");
  return out;
}
struct Zones {
  std::vector<std::int64_t> ids, first, last, types;
  std::vector<std::string> labels;
};
Zones zones(hid_t file, const std::string& base, bool faces) {
  H5 group(H5Gopen2(file, base.c_str(), H5P_DEFAULT), H5Gclose, base);
  auto n = attribute(group, "nZones");
  KUTTA_REQUIRE(n > 0 && n <= std::numeric_limits<Index>::max(),
                "CFF: invalid zone count");
  Zones z;
  z.ids = ints(file, base + "/id");
  z.first = ints(file, base + "/minId");
  z.last = ints(file, base + "/maxId");
  // zoneType is the boundary-condition code. faceType is element topology
  // (triangle/quad/mixed); confusing these silently turns hex walls into
  // inlets.
  if (exists(file, base + "/zoneType"))
    z.types = ints(file, base + "/zoneType");
  else if (faces)
    throw Error(
        "CFF: faces/zoneTopology/zoneType is required for boundary types");
  else
    z.types.assign(static_cast<std::size_t>(n), 1);
  for (const auto* v : {&z.ids, &z.first, &z.last, &z.types})
    KUTTA_REQUIRE(v->size() == static_cast<std::size_t>(n),
                  "CFF: inconsistent zoneTopology column length at " + base);
  z.labels = names(file, base + "/name", static_cast<std::size_t>(n));
  for (std::size_t i = 0; i < z.ids.size(); ++i) {
    id(z.ids[i]);
    id(z.first[i]);
    id(z.last[i]);
    KUTTA_REQUIRE(z.last[i] >= z.first[i], "CFF: reversed zone range");
  }
  std::set<std::int64_t> ids(z.ids.begin(), z.ids.end());
  KUTTA_REQUIRE(ids.size() == z.ids.size(), "CFF: duplicate zone id");
  std::vector<std::pair<Index, Index>> ranges;
  for (std::size_t i = 0; i < z.ids.size(); ++i)
    ranges.emplace_back(id(z.first[i]), id(z.last[i]));
  std::sort(ranges.begin(), ranges.end());
  for (std::size_t i = 1; i < ranges.size(); ++i)
    KUTTA_REQUIRE(ranges[i - 1].second < ranges[i].first,
                  "CFF: overlapping zone ranges at " + base);
  return z;
}
template <class T>
void coverage(const Zones& z, const std::map<Index, T>& values) {
  std::uint64_t expected = 0;
  for (std::size_t i = 0; i < z.ids.size(); ++i) {
    expected += z.last[i] - z.first[i] + 1;
    auto it = values.find(id(z.first[i]));
    for (auto n = z.first[i]; n <= z.last[i]; ++n, ++it)
      KUTTA_REQUIRE(it != values.end() && it->first == n,
                    "CFF: missing section data in zone topology");
  }
  KUTTA_REQUIRE(values.size() == expected,
                "CFF: section data outside zone topology");
}
template <class F>
void sections(hid_t file, const std::string& base, bool datasets, F callback,
              bool numbered = true) {
  H5 g(H5Gopen2(file, base.c_str(), H5P_DEFAULT), H5Gclose, base);
  H5G_info_t info{};
  check(H5Gget_info(g, &info), base);
  auto n = numbered ? attribute(g, "nSections")
                    : static_cast<std::int64_t>(info.nlinks);
  KUTTA_REQUIRE(n >= 0 && n <= std::numeric_limits<Index>::max(),
                "CFF: invalid nSections");
  for (std::int64_t i = 1; i <= n; ++i) {
    auto name = std::to_string(i);
    if (!numbered) {
      auto length = H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC,
                                      i - 1, nullptr, 0, H5P_DEFAULT);
      KUTTA_REQUIRE(length >= 0, "CFF: cannot enumerate " + base);
      std::vector<char> buffer(static_cast<std::size_t>(length) + 1);
      check(H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC, i - 1,
                              buffer.data(), buffer.size(), H5P_DEFAULT), base);
      name = buffer.data();
    }
    auto path = base + "/" + name;
    H5 object(datasets ? H5Dopen2(file, path.c_str(), H5P_DEFAULT)
                       : H5Gopen2(file, path.c_str(), H5P_DEFAULT),
              datasets ? H5Dclose : H5Gclose, path);
    auto [a, b] = range(object);
    callback(path, a, b);
  }
}
}  // namespace
PolyMesh read_fluent_cff(const std::filesystem::path& path, double thickness) {
  H5 file(H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose,
          path.string());
  const std::string root = "/meshes/1";
  KUTTA_REQUIRE(exists(file, root), "CFF: file contains no /meshes/1 mesh");
  H5 mesh(H5Gopen2(file, root.c_str(), H5P_DEFAULT), H5Gclose, root);
  for (const auto& unsupported :
       {root + "/cells/tree", root + "/faces/tree", root + "/faces/interface",
        std::string("/special/Overset_DCI")})
    KUTTA_REQUIRE(!exists(file, unsupported),
                  "CFF: adapted/nonconformal/overset meshes are unsupported: " +
                      unsupported);
  detail::FluentRaw raw;
  auto dimension = attribute(mesh, "dimension");
  KUTTA_REQUIRE(dimension == 2 || dimension == 3,
                "CFF: dimension must be 2 or 3");
  raw.dimension = static_cast<int>(dimension);
  auto nz = zones(file, root + "/nodes/zoneTopology", false);
  // Fluent versions name coordinates by either zone id or section ordinal.
  // Names carry no placement information; the declared id ranges do.
  sections(file, root + "/nodes/coords", true,
           [&](const std::string& path_coords, Index a, Index b) {
             auto values = read<double>(file, path_coords, H5T_NATIVE_DOUBLE);
             auto n = static_cast<std::size_t>(b) - a + 1;
             KUTTA_REQUIRE(
                 values.size() == n * raw.dimension,
                 "CFF: node coordinate dataset length does not match "
                 "dimension/range");
             for (std::size_t i = 0; i < n; ++i) {
               auto k = i * raw.dimension;
               KUTTA_REQUIRE(raw.nodes
                                 .emplace(a + static_cast<Index>(i),
                                          Vec3d{values[k], values[k + 1],
                                                raw.dimension == 3
                                                    ? values[k + 2] : 0})
                                 .second,
                             "CFF: overlapping node zones");
             }
           }, false);
  coverage(nz, raw.nodes);
  auto cz = zones(file, root + "/cells/zoneTopology", false);
  for (std::size_t z = 0; z < cz.ids.size(); ++z) {
    KUTTA_REQUIRE(cz.types[z] != 32,
                  "CFF: inactive parent cells are unsupported");
    auto zone = id(cz.ids[z]);
    raw.zones[zone] = {
        cz.labels[z].empty() ? "zone" + std::to_string(zone) : cz.labels[z],
        "fluent-cell-zone-type-" + std::to_string(cz.types[z])};
    for (auto c = cz.first[z]; c <= cz.last[z]; ++c)
      KUTTA_REQUIRE(raw.cells.emplace(id(c), zone).second,
                    "CFF: overlapping cell zones");
  }
  std::map<Index, std::int64_t> cell_types;
  if (exists(file, root + "/cells/zoneTopology/cellType")) {
    auto types = ints(file, root + "/cells/zoneTopology/cellType");
    KUTTA_REQUIRE(types.size() == cz.ids.size(),
                  "CFF: cellType count does not match topology");
    for (std::size_t z = 0; z < types.size(); ++z) {
      KUTTA_REQUIRE(types[z] >= 0 && types[z] <= 7,
                    "CFF: invalid zone cell element type");
      cell_types.emplace(id(cz.ids[z]), types[z]);
    }
  }
  std::set<Index> typed;
  sections(file, root + "/cells/ctype", false,
           [&](const std::string& base, Index a, Index b) {
             H5 group(H5Gopen2(file, base.c_str(), H5P_DEFAULT), H5Gclose, base);
             auto type = attribute(group, "elementType");
             KUTTA_REQUIRE(type >= 0 && type <= 7,
                           "CFF: invalid cell element type");
             std::vector<std::int64_t> types;
             if (type == 0) {
               types = ints(file, base + "/cell-types");
               KUTTA_REQUIRE(types.size() == static_cast<std::size_t>(b) - a + 1,
                             "CFF: cell-types length/range mismatch");
             }
             for (std::int64_t c = a; c <= b; ++c) {
               auto it = raw.cells.find(id(c));
               KUTTA_REQUIRE(it != raw.cells.end() && typed.insert(id(c)).second,
                             "CFF: ctype section outside topology or overlapping");
               auto element = type ? type : types[static_cast<std::size_t>(c - a)];
               KUTTA_REQUIRE(element >= 1 && element <= 7,
                             "CFF: invalid cell element type");
               auto zone_type = cell_types.find(it->second);
               KUTTA_REQUIRE(zone_type == cell_types.end() || zone_type->second == 0 ||
                                 zone_type->second == element,
                             "CFF: ctype disagrees with zone cellType");
             }
           });
  KUTTA_REQUIRE(typed.size() == raw.cells.size(), "CFF: missing ctype data");
  auto fz = zones(file, root + "/faces/zoneTopology", true);
  for (std::size_t z = 0; z < fz.ids.size(); ++z) {
    auto zone = id(fz.ids[z]);
    auto type = id(fz.types[z], true);
    KUTTA_REQUIRE(!raw.zones.count(zone),
                  "CFF: cell and face zones share an id");
    raw.zones[zone] = {
        fz.labels[z].empty() ? "zone" + std::to_string(zone) : fz.labels[z],
        detail::fluent_type(type)};
    for (auto f = fz.first[z]; f <= fz.last[z]; ++f)
      KUTTA_REQUIRE(
          raw.faces.emplace(id(f), detail::FluentFace{{}, 0, 0, zone}).second,
          "CFF: overlapping face zones");
  }
  sections(file, root + "/faces/nodes", false,
           [&](const std::string& base, Index a, Index b) {
             auto counts = ints(file, base + "/nnodes"),
                  nodes = ints(file, base + "/nodes");
             KUTTA_REQUIRE(counts.size() == static_cast<std::size_t>(b) - a + 1,
                           "CFF: nnodes length/range mismatch");
             std::size_t ptr = 0;
             for (std::int64_t f = a; f <= b; ++f) {
               auto it = raw.faces.find(id(f));
               KUTTA_REQUIRE(it != raw.faces.end(),
                             "CFF: face nodes outside zone topology");
               auto& face = it->second;
               KUTTA_REQUIRE(face.nodes.empty(),
                             "CFF: overlapping face node sections");
               auto n = counts[static_cast<std::size_t>(f - a)];
               KUTTA_REQUIRE(n >= 2 && static_cast<std::uint64_t>(n) <=
                                           nodes.size() - ptr,
                             "CFF: invalid/truncated face node count");
               for (std::int64_t i = 0; i < n; ++i)
                 face.nodes.push_back(id(nodes[ptr++]));
             }
             KUTTA_REQUIRE(ptr == nodes.size(), "CFF: trailing face nodes");
           });
  for (bool c0 : {true, false}) {
    std::set<Index> seen;
    sections(
        file, root + (c0 ? "/faces/c0" : "/faces/c1"), true,
        [&](const std::string& base, Index a, Index b) {
          auto values = ints(file, base);
          KUTTA_REQUIRE(values.size() == static_cast<std::size_t>(b) - a + 1,
                        "CFF: c0/c1 dataset length/range mismatch");
          for (std::int64_t f = a; f <= b; ++f) {
            auto it = raw.faces.find(id(f));
            KUTTA_REQUIRE(it != raw.faces.end() && seen.insert(id(f)).second,
                          "CFF: c0/c1 section outside topology or overlapping");
            auto value = id(values[static_cast<std::size_t>(f - a)], true);
            if (c0)
              it->second.c0 = value;
            else
              it->second.c1 = value;
          }
        });
    if (c0)
      KUTTA_REQUIRE(seen.size() == raw.faces.size(),
                    "CFF: missing c0 connectivity");
  }
  for (const auto& [f, face] : raw.faces) {
    (void)f;
    KUTTA_REQUIRE(!face.nodes.empty(), "CFF: missing face node connectivity");
  }
  for (const auto& [name, size] :
       {std::pair{"nodeCount", raw.nodes.size()},
        std::pair{"cellCount", raw.cells.size()},
        std::pair{"faceCount", raw.faces.size()}}) {
    if (H5Aexists(mesh, name) > 0)
      KUTTA_REQUIRE(attribute(mesh, name) == static_cast<std::int64_t>(size),
                    std::string("CFF: ") + name + " disagrees with topology");
  }
  // CFF stores the Fluent file face winding, not the solver's in-memory normal.
  return detail::assemble_fluent(raw, thickness).mesh;
}
}  // namespace kutta
#else
namespace kutta {
PolyMesh read_fluent_cff(const std::filesystem::path&, double) {
  throw Error(
      "Fluent CFF reading requires HDF5 (KUTTA_HAS_HDF5=0). Install HDF5 and "
      "configure with -DKUTTA_ENABLE_HDF5=ON.");
}
}  // namespace kutta
#endif
