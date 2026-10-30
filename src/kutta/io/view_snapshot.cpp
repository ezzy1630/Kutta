#include "view.hpp"
#include "binary.hpp"
#include "kutta/core/parse.hpp"
#include "kutta/mesh/detail.hpp"
#include <cmath>
#include <limits>
#include <set>

namespace kutta::view {
namespace {
constexpr std::uint32_t magic=0x5254554b; // KUTR
void validate_fields(std::span<const Field> fields,std::size_t rows) {
  std::set<std::string> names;
  for(const auto& f:fields) {
    KUTTA_REQUIRE(names.insert(f.name).second,"duplicate snapshot field");
    KUTTA_REQUIRE((f.components==1 || f.components==3) && f.values.size()==rows*f.components,"snapshot field shape mismatch: "+f.name);
    for(float v:f.values)KUTTA_REQUIRE(!std::isinf(v),"infinite snapshot field: "+f.name);
  }
}
void write_fields(binary::Writer& w,std::span<const Field> fields) {
  w.scalar(std::uint32_t(fields.size()));
  for(const auto& f:fields){w.string(f.name);w.string(f.unit);w.scalar(f.components);w.array<float>(f.values);}
}
std::vector<Field> read_fields(binary::Reader& r,std::size_t rows) {
  const auto count=r.scalar<std::uint32_t>();KUTTA_REQUIRE(count<=1024,"too many snapshot fields");
  std::vector<Field> fields;
  for(std::uint32_t i=0;i<count;++i){Field f;f.name=r.string();f.unit=r.string();f.components=r.scalar<std::uint32_t>();KUTTA_REQUIRE(f.components==1 || f.components==3,"invalid field components");f.values=r.array<float>(rows*f.components);fields.push_back(std::move(f));}
  validate_fields(fields,rows);return fields;
}
}
Plane parse_plane(std::string_view text) {
  double v[6];
  for(int i=0;i<6;++i){const auto end=text.find(',');const auto n=parse_double(text.substr(0,end));KUTTA_REQUIRE(n && std::isfinite(*n),"--slice requires six finite numbers: ox,oy,oz,nx,ny,nz");v[i]=*n;KUTTA_REQUIRE((i==5)==(end==std::string_view::npos),"--slice requires exactly six comma-separated numbers");if(i<5)text.remove_prefix(end+1);}
  Plane p{{v[0],v[1],v[2]},{v[3],v[4],v[5]}};
  const double scale=std::max({std::abs(v[3]),std::abs(v[4]),std::abs(v[5])});
  KUTTA_REQUIRE(scale>0,"--slice normal must be nonzero");
  p.normal=detail::operator/(p.normal,scale);p.normal=detail::operator/(p.normal,detail::norm(p.normal));return p;
}
void write_snapshot(const std::filesystem::path& path,const PolyMesh& m,const Fields& s) {
  validate_fields(s.cells,m.n_cells);validate_fields(s.boundary,m.n_faces()-m.n_internal_faces());
  binary::Writer w(path);
  w.scalar(magic);w.scalar(std::uint32_t(2));w.scalar(std::uint32_t(s.reference));w.string(s.provenance);
  for(auto n:{m.n_points(),m.n_faces(),m.n_internal_faces(),m.n_cells,detail::index(m.face_nodes.size()),detail::index(m.patches.size())})w.scalar(n);
  // Double coordinates preserve translated meshes when cutting planes. View
  // payloads are float32; this private snapshot is never loaded by the UI.
  for(auto p:m.points){w.scalar(p.x);w.scalar(p.y);w.scalar(p.z);}
  w.array<Index>(m.face_offsets);w.array<Index>(m.face_nodes);w.array<Index>(m.owner);w.array<Index>(m.neighbour);
  for(const auto& p:m.patches){w.string(p.name);w.string(p.source_type);w.scalar(std::uint32_t(p.kind));w.scalar(std::int32_t(p.source_zone_id));w.scalar(p.start);w.scalar(p.size);}
  write_fields(w,s.cells);write_fields(w,s.boundary);
  for(const auto& n:{s.properties.density,s.properties.viscosity,s.properties.reference_density,s.properties.velocity,s.properties.pressure,s.properties.area,s.properties.length}) {
    w.scalar(std::uint32_t(n.has_value()));if(n){KUTTA_REQUIRE(std::isfinite(*n),"nonfinite snapshot reference");w.scalar(*n);}
  }
  w.close();
}
Snapshot read_snapshot(const std::filesystem::path& path) {
  binary::Reader r(path);KUTTA_REQUIRE(r.scalar<std::uint32_t>()==magic,"unknown run snapshot format");
  const auto version=r.scalar<std::uint32_t>();KUTTA_REQUIRE(version==1 || version==2,"unknown run snapshot version");
  Snapshot s;const auto reference=r.scalar<std::uint32_t>();KUTTA_REQUIRE(reference<=1,"invalid snapshot reference flag");s.reference=reference;s.provenance=r.string();
  auto& m=s.mesh;
  const auto np=r.scalar<Index>(),nf=r.scalar<Index>(),ni=r.scalar<Index>(),nc=r.scalar<Index>(),nn=r.scalar<Index>(),nz=r.scalar<Index>();
  KUTTA_REQUIRE(np>0 && nf>0 && ni>=0 && ni<=nf && nc>0 && nn>=0 && nz>=0,"invalid snapshot mesh sizes");
  KUTTA_REQUIRE(std::uint64_t(np)*24+std::uint64_t(nf)*8+std::uint64_t(ni)*4+std::uint64_t(nn)*4<=r.remaining(),"snapshot mesh exceeds file size");
  m.n_cells=nc;m.points.resize(np);
  for(auto& p:m.points){p={r.scalar<double>(),r.scalar<double>(),r.scalar<double>()};KUTTA_REQUIRE(detail::finite(p),"nonfinite snapshot point");}
  m.face_offsets=r.array<Index>(std::size_t(nf)+1);m.face_nodes=r.array<Index>(nn);m.owner=r.array<Index>(nf);m.neighbour=r.array<Index>(ni);
  KUTTA_REQUIRE(m.face_offsets.front()==0 && m.face_offsets.back()==nn,"invalid snapshot CSR");
  for(Index f=0;f<nf;++f){KUTTA_REQUIRE(m.face_offsets[f]>=0 && m.face_offsets[f+1]>=m.face_offsets[f]+3 && m.face_offsets[f+1]<=nn,"invalid snapshot face offsets");KUTTA_REQUIRE(m.owner[f]>=0 && m.owner[f]<nc,"invalid snapshot owner");if(f<ni)KUTTA_REQUIRE(m.neighbour[f]>m.owner[f] && m.neighbour[f]<nc,"invalid snapshot neighbour");}
  for(auto n:m.face_nodes)KUTTA_REQUIRE(n>=0 && n<np,"invalid snapshot point index");
  KUTTA_REQUIRE(nz<=nf,"invalid snapshot patch count");Index next=ni;std::set<std::string> names;
  for(Index i=0;i<nz;++i){Patch p;p.name=r.string();p.source_type=r.string();const auto kind=r.scalar<std::uint32_t>();KUTTA_REQUIRE(kind<=std::uint32_t(PatchKind::Other),"invalid patch kind");p.kind=PatchKind(kind);p.source_zone_id=r.scalar<std::int32_t>();p.start=r.scalar<Index>();p.size=r.scalar<Index>();KUTTA_REQUIRE(p.start==next && p.size>=0 && p.size<=nf-next && names.insert(p.name).second,"invalid snapshot patch range");next+=p.size;m.patches.push_back(std::move(p));}
  KUTTA_REQUIRE(next==nf,"incomplete snapshot patches");
  s.cells=read_fields(r,nc);s.boundary=read_fields(r,nf-ni);
  if(version==2)for(auto* value:{&s.properties.density,&s.properties.viscosity,&s.properties.reference_density,&s.properties.velocity,&s.properties.pressure,&s.properties.area,&s.properties.length}) {
    const auto present=r.scalar<std::uint32_t>();KUTTA_REQUIRE(present<=1,"invalid reference presence flag");
    if(present){const auto n=r.scalar<double>();KUTTA_REQUIRE(std::isfinite(n),"nonfinite snapshot reference");*value=n;}
  }
  KUTTA_REQUIRE(r.remaining()==0,"trailing snapshot bytes");return s;
}
} // namespace kutta::view
