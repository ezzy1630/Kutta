#include "restart.hpp"

#include "kutta/io/binary.hpp"

#include <cmath>

namespace kutta::app {
namespace {
constexpr std::string_view magic="KUTRST01";

std::string_view layout(const fv::SolverSettings& settings) {
  if(settings.sst) return settings.sst->transition?"sst-intermittency":"sst";
  if(settings.realizable_ke) return "k-epsilon";
  KUTTA_REQUIRE(!settings.sa,"native restart does not support the SA closure");
  return "laminar";
}
}

setup::CellFields read_restart(const std::filesystem::path& path,std::string_view mesh_hash,
                               Index cells,const fv::SolverSettings& settings) {
  binary::Reader reader(path);
  std::array<char,magic.size()> header;
  reader.bytes(header.data(),header.size());
  KUTTA_REQUIRE(std::string_view(header.data(),header.size())==magic,"unsupported native restart format");
  const auto matches=[&](std::string_view expected) {
    if(reader.scalar<std::uint32_t>()!=expected.size()) return false;
    std::array<char,64> value;
    reader.bytes(value.data(),expected.size());
    return std::string_view(value.data(),expected.size())==expected;
  };
  KUTTA_REQUIRE(matches(mesh_hash),"native restart requires the identical mesh file (SHA-256 mismatch)");
  KUTTA_REQUIRE(matches(layout(settings)),"native restart turbulence/transition model mismatch");
  KUTTA_REQUIRE(reader.scalar<std::uint32_t>()==std::uint32_t(cells),"native restart cell count mismatch");
  const auto components=4+2*int(bool(settings.sst || settings.realizable_ke))+int(settings.sst && settings.sst->transition);
  KUTTA_REQUIRE(reader.remaining()==std::uint64_t(cells)*components*sizeof(float),"native restart field payload has the wrong length");
  const auto field=[&](auto valid) {
    auto values=reader.array<float>(std::size_t(cells));
    KUTTA_REQUIRE(std::all_of(values.begin(),values.end(),[&](float x){return std::isfinite(x) && valid(x);}),
                  "invalid value in native restart field");
    return values;
  };
  const auto finite=[](float){return true;};
  const auto nonnegative=[](float x){return x>=0;};
  const auto positive=[](float x){return x>0;};
  setup::CellFields result;
  result.velocity.resize(std::size_t(cells));
  for(auto component:{&Vec3d::x,&Vec3d::y,&Vec3d::z}) {
    const auto values=field(finite);
    for(Index c=0;c<cells;++c) result.velocity[c].*component=values[c];
  }
  const auto load=[&](std::vector<double>& into,auto valid) {
    const auto values=field(valid);
    into.assign(values.begin(),values.end());
  };
  load(result.pressure,finite);
  if(settings.sst || settings.realizable_ke) load(result.k,nonnegative);
  if(settings.sst) load(result.omega,positive);
  if(settings.realizable_ke) load(result.epsilon,positive);
  if(settings.sst && settings.sst->transition)
    load(result.intermittency,[](float x){return x>=0 && x<=1;});
  return result;
}

void write_restart(const std::filesystem::path& path,std::string_view mesh_hash,
                   Device& device,const fv::SteadySolver& solver,const fv::SolverSettings& settings) {
  auto temporary=path;
  temporary+=".part";
  try {
    binary::Writer writer(temporary);
    writer.bytes(magic.data(),magic.size());
    writer.string(mesh_hash);
    writer.string(layout(settings));
    writer.scalar(std::uint32_t(solver.mesh().n_cells));
    const auto& map=solver.permutation().cells_new_to_old;
    const auto field=[&](const Buffer<float>& buffer) {
      const auto values=device.download(buffer);
      std::vector<float> original(values.size());
      for(std::size_t c=0;c<values.size();++c) {
        KUTTA_REQUIRE(std::isfinite(values[c]),"non-finite native restart field");
        original[map[c]]=values[c];
      }
      writer.array<float>(original);
    };
    field(solver.velocity().x);
    field(solver.velocity().y);
    field(solver.velocity().z);
    field(solver.pressure());
    if(settings.sst || settings.realizable_ke) field(solver.kinetic_energy());
    if(settings.sst) field(solver.specific_dissipation());
    if(settings.realizable_ke) field(solver.dissipation());
    if(settings.sst && settings.sst->transition) field(solver.intermittency());
    writer.close();
    std::filesystem::rename(temporary,path);
  } catch(...) {
    std::error_code ec;
    std::filesystem::remove(temporary,ec);
    throw;
  }
}
}
