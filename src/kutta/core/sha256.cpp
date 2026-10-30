#include "kutta/core/sha256.hpp"
#include "kutta/core/types.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace kutta {
namespace {
constexpr std::array<std::uint32_t,64> constants{
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
void compress(std::array<std::uint32_t,8>& h,const unsigned char* block) {
  std::array<std::uint32_t,64> w{};
  for(int i=0;i<16;++i) for(int j=0;j<4;++j) w[i]=(w[i]<<8)|block[4*i+j];
  for(int i=16;i<64;++i) {
    const auto x=w[i-15],y=w[i-2];
    w[i]=w[i-16]+(std::rotr(x,7)^std::rotr(x,18)^(x>>3))+w[i-7]+(std::rotr(y,17)^std::rotr(y,19)^(y>>10));
  }
  auto [a,b,c,d,e,f,g,t]=h;
  for(int i=0;i<64;++i) {
    const auto v=t+(std::rotr(e,6)^std::rotr(e,11)^std::rotr(e,25))+((e&f)^(~e&g))+constants[i]+w[i];
    const auto u=(std::rotr(a,2)^std::rotr(a,13)^std::rotr(a,22))+((a&b)^(a&c)^(b&c));
    t=g;g=f;f=e;e=d+v;d=c;c=b;b=a;a=v+u;
  }
  const std::array<std::uint32_t,8> state{a,b,c,d,e,f,g,t};
  for(int i=0;i<8;++i) h[i]+=state[i];
}
}
static std::string sha256_stream(std::istream&,const std::string&);
std::string sha256(std::string_view text) {
  std::istringstream in{std::string(text)};
  return sha256_stream(in,"text");
}
std::string sha256_file(const std::filesystem::path& path) {
  std::ifstream in(path,std::ios::binary);
  KUTTA_REQUIRE(in,"cannot hash input: "+path.string());
  return sha256_stream(in,path.string());
}
static std::string sha256_stream(std::istream& in,const std::string& name) {
  std::array<std::uint32_t,8> h{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  std::array<unsigned char,64> block{};
  std::uint64_t bytes=0;
  while(in.read(reinterpret_cast<char*>(block.data()),64)) {compress(h,block.data());bytes+=64;}
  KUTTA_REQUIRE(in.eof(),"failed reading input hash: "+name);
  const auto tail=static_cast<std::size_t>(in.gcount());bytes+=tail;
  block[tail]=0x80;
  std::fill(block.begin()+tail+1,block.end(),0);
  if(tail>=56) {compress(h,block.data());block.fill(0);}
  for(int i=0;i<8;++i) block[63-i]=static_cast<unsigned char>((bytes*8)>>(8*i));
  compress(h,block.data());
  std::ostringstream out;out<<std::hex<<std::setfill('0');
  for(auto word:h) out<<std::setw(8)<<word;
  return out.str();
}
}
