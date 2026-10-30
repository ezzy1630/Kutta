#pragma once
// Private HDF5 boundary shared by the case and data readers.
#if KUTTA_HAS_HDF5
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <hdf5.h>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
namespace kutta::fluent::h5 {
struct Handle {
  hid_t id = -1;
  herr_t (*close)(hid_t) = nullptr;
  Handle(hid_t value, herr_t (*closer)(hid_t)) : id(value), close(closer) {
    if (id < 0)
      throw std::runtime_error("HDF5 object open failed");
  }
  ~Handle() {
    if (id >= 0)
      close(id);
  }
  Handle(const Handle &) = delete;
  Handle &operator=(const Handle &) = delete;
  operator hid_t() const { return id; }
};
inline bool exists(hid_t file, const std::string &path) {
  htri_t result = 0;
  H5E_BEGIN_TRY { result = H5Lexists(file, path.c_str(), H5P_DEFAULT); }
  H5E_END_TRY;
  return result > 0;
}
inline Handle open(const std::filesystem::path &path) {
  hid_t id = -1;
  H5E_BEGIN_TRY { id = H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); }
  H5E_END_TRY;
  if (id < 0)
    throw std::runtime_error("cannot open HDF5 file: " + path.string());
  return {id, H5Fclose};
}
inline void check(herr_t code, const std::string &context) {
  if (code < 0)
    throw std::runtime_error("HDF5 read failed: " + context);
}
inline void local_link(hid_t file, const std::string &path) {
  // Refuse external/soft links rather than following references in untrusted files.
  std::size_t end = path.starts_with('/') ? 1 : 0;
  while (end < path.size()) {
    end = path.find('/', end);
    auto prefix = path.substr(0, end);
#if H5_VERSION_GE(1, 12, 0)
    H5L_info2_t info{};
    check(H5Lget_info2(file, prefix.c_str(), &info, H5P_DEFAULT), prefix);
#else  // Ubuntu 24.04 ships HDF5 1.10
    H5L_info_t info{};
    check(H5Lget_info(file, prefix.c_str(), &info, H5P_DEFAULT), prefix);
#endif
    if (info.type != H5L_TYPE_HARD)
      throw std::runtime_error("non-local HDF5 link: " + prefix);
    if (end == std::string::npos)
      break;
    ++end;
  }
}
inline void local_storage(hid_t dataset) {
  Handle properties(H5Dget_create_plist(dataset), H5Pclose);
  if (H5Pget_external_count(properties) != 0 || H5Pget_layout(properties) == H5D_VIRTUAL)
    throw std::runtime_error("external/virtual HDF5 dataset storage unsupported");
  auto n = H5Pget_nfilters(properties);
  if (n < 0)
    throw std::runtime_error("invalid HDF5 filters");
  for (int i = 0; i < n; ++i) {
    unsigned flags = 0;
    std::size_t count = 0;
    auto filter = H5Pget_filter2(properties, static_cast<unsigned>(i), &flags, &count, nullptr, 0,
                                 nullptr, nullptr);
    if (filter != H5Z_FILTER_DEFLATE && filter != H5Z_FILTER_SHUFFLE &&
        filter != H5Z_FILTER_FLETCHER32 && filter != H5Z_FILTER_SZIP && filter != H5Z_FILTER_NBIT &&
        filter != H5Z_FILTER_SCALEOFFSET)
      throw std::runtime_error("untrusted HDF5 plugin filter unsupported");
  }
}
inline std::vector<hsize_t> shape(hid_t ds) {
  Handle sp(H5Dget_space(ds), H5Sclose);
  int rank = H5Sget_simple_extent_ndims(sp);
  if (rank < 0 || rank > 2)
    throw std::runtime_error("unsupported dataset rank");
  std::vector<hsize_t> dims(static_cast<std::size_t>(rank));
  check(H5Sget_simple_extent_dims(sp, dims.data(), nullptr), "shape");
  return dims;
}
inline std::size_t count(const std::vector<hsize_t> &dims, std::size_t limit = 100'000'000) {
  std::size_t n = 1;
  for (auto d : dims) {
    if (d > limit || (d && n > limit / d))
      throw std::runtime_error("HDF5 allocation limit exceeded");
    n *= d;
  }
  return n;
}
inline std::vector<std::string> children(hid_t file, const std::string &path) {
  local_link(file, path);
  Handle g(H5Gopen2(file, path.c_str(), H5P_DEFAULT), H5Gclose);
  H5G_info_t info{};
  check(H5Gget_info(g, &info), path);
  if (info.nlinks > 100000)
    throw std::runtime_error("HDF5 group size limit exceeded");
  std::vector<std::string> result;
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    auto n = H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
    if (n < 0 || n > 65536)
      throw std::runtime_error("invalid HDF5 link name");
    std::string name(static_cast<std::size_t>(n) + 1, '\0');
    check(static_cast<herr_t>(H5Lget_name_by_idx(g, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(),
                                                 name.size(), H5P_DEFAULT)),
          path);
    name.resize(static_cast<std::size_t>(n));
    local_link(file, path + "/" + name);
    result.push_back(std::move(name));
  }
  return result;
}
inline std::string string(hid_t file, const std::string &path) {
  local_link(file, path);
  Handle ds(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  local_storage(ds);
  Handle ty(H5Dget_type(ds), H5Tclose);
  if (H5Tget_class(ty) != H5T_STRING || count(shape(ds)) != 1)
    throw std::runtime_error("expected one HDF5 string: " + path);
  if (H5Tis_variable_str(ty) > 0) {
    // Fixed strings are Fluent's observed layout; refuse VLEN before allocation.
    throw std::runtime_error("unbounded variable-length HDF5 string unsupported: " + path);
  }
  auto n = H5Tget_size(ty);
  if (n > 128 * 1024 * 1024)
    throw std::runtime_error("HDF5 string size limit exceeded");
  std::string value(n, '\0');
  check(H5Dread(ds, ty, H5S_ALL, H5S_ALL, H5P_DEFAULT, value.data()), path);
  auto nul = value.find('\0');
  if (nul != std::string::npos)
    value.resize(nul);
  return value;
}
template <class T> inline hid_t native();
template <> inline hid_t native<std::uint64_t>() { return H5T_NATIVE_UINT64; }
template <> inline hid_t native<std::int64_t>() { return H5T_NATIVE_INT64; }
template <> inline hid_t native<double>() { return H5T_NATIVE_DOUBLE; }
template <class T>
std::vector<T> numbers(hid_t file, const std::string &path,
                       std::vector<hsize_t> *dimensions = nullptr) {
  local_link(file, path);
  Handle ds(H5Dopen2(file, path.c_str(), H5P_DEFAULT), H5Dclose);
  local_storage(ds);
  auto dims = shape(ds);
  Handle ty(H5Dget_type(ds), H5Tclose);
  auto cls = H5Tget_class(ty);
  if (cls != H5T_INTEGER && cls != H5T_FLOAT)
    throw std::runtime_error("expected numeric dataset: " + path);
  if constexpr (std::is_integral_v<T>)
    if (cls != H5T_INTEGER)
      throw std::runtime_error("expected integer IDs: " + path);
  std::vector<T> out(count(dims));
  if (!out.empty())
    check(H5Dread(ds, native<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()), path);
  if (dimensions)
    *dimensions = std::move(dims);
  return out;
}
inline std::optional<std::uint64_t> attribute(hid_t file, const std::string &path,
                                              const char *name) {
  local_link(file, path);
  Handle obj(H5Oopen(file, path.c_str(), H5P_DEFAULT), H5Oclose);
  if (H5Aexists(obj, name) <= 0)
    return {};
  Handle a(H5Aopen(obj, name, H5P_DEFAULT), H5Aclose);
  Handle sp(H5Aget_space(a), H5Sclose);
  Handle ty(H5Aget_type(a), H5Tclose);
  if (H5Sget_simple_extent_npoints(sp) != 1 || H5Tget_class(ty) != H5T_INTEGER)
    throw std::runtime_error("expected scalar integer attribute");
  if (H5Tget_sign(ty) == H5T_SGN_2) {
    std::int64_t v = 0;
    check(H5Aread(a, H5T_NATIVE_INT64, &v), name);
    if (v < 0)
      throw std::runtime_error("negative ID/count attribute");
    return static_cast<std::uint64_t>(v);
  }
  std::uint64_t v = 0;
  check(H5Aread(a, H5T_NATIVE_UINT64, &v), name);
  return v;
}
} // namespace kutta::fluent::h5
#endif
