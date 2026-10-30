#pragma once

#include "kutta/core/types.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <span>
#include <type_traits>
#include <vector>

namespace kutta::binary {
// Explicit endian conversion; write arrays in blocks, never one syscall per
// scalar. Readers bound every allocation by the remaining file length.
template<class T> T little(T n) {
  static_assert(std::is_arithmetic_v<T>);
  if constexpr(std::endian::native==std::endian::big) {
    auto bytes=std::bit_cast<std::array<std::byte,sizeof(T)>>(n);
    std::reverse(bytes.begin(),bytes.end());return std::bit_cast<T>(bytes);
  }
  return n;
}
class Writer {
  std::ofstream out_;
 public:
  explicit Writer(const std::filesystem::path& p):out_(p,std::ios::binary) {KUTTA_REQUIRE(out_,"cannot open binary output: "+p.string());out_.exceptions(std::ios::failbit|std::ios::badbit);}
  void bytes(const void* data,std::size_t n) {out_.write(static_cast<const char*>(data),std::streamsize(n));}
  template<class T> void scalar(T n) {n=little(n);bytes(&n,sizeof(n));}
  template<class T> void array(std::span<const T> v) {
    if constexpr(std::endian::native==std::endian::little)bytes(v.data(),v.size_bytes());
    else {std::vector<T> block;block.reserve(16384);for(std::size_t i=0;i<v.size();){block.clear();for(int j=0;j<16384 && i<v.size();++j,++i)block.push_back(little(v[i]));bytes(block.data(),block.size()*sizeof(T));}}
  }
  void string(std::string_view s) {KUTTA_REQUIRE(s.size()<=UINT32_MAX,"binary string too large");scalar(std::uint32_t(s.size()));bytes(s.data(),s.size());}
  void close() {out_.close();}
};
class Reader {
  std::ifstream in_;
  std::uint64_t remaining_;
 public:
  explicit Reader(const std::filesystem::path& p):in_(p,std::ios::binary),remaining_(std::filesystem::file_size(p)) {KUTTA_REQUIRE(in_,"cannot open snapshot: "+p.string());in_.exceptions(std::ios::failbit|std::ios::badbit);}
  void bytes(void* data,std::size_t n) {KUTTA_REQUIRE(n<=remaining_,"truncated snapshot");in_.read(static_cast<char*>(data),std::streamsize(n));remaining_-=n;}
  template<class T> T scalar() {T n;bytes(&n,sizeof(n));return little(n);}
  template<class T> std::vector<T> array(std::size_t count) {
    KUTTA_REQUIRE(count<=remaining_/sizeof(T),"snapshot array exceeds file length");std::vector<T> out(count);bytes(out.data(),out.size()*sizeof(T));
    if constexpr(std::endian::native==std::endian::big)for(auto& n:out)n=little(n);
    return out;
  }
  std::string string() {const auto n=scalar<std::uint32_t>();KUTTA_REQUIRE(n<=remaining_,"snapshot string exceeds file length");std::string out(n,'\0');bytes(out.data(),n);return out;}
  std::uint64_t remaining() const {return remaining_;}
};
} // namespace kutta::binary
