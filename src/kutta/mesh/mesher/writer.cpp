#include "external_aero.hpp"
#include "kutta/mesh/detail.hpp"
#include <fstream>
#include <iomanip>

namespace kutta::mesher {
void write_msh(const std::filesystem::path& path,const PolyMesh& mesh) {
  std::ofstream out(path,std::ios::binary); out.imbue(std::locale::classic());
  KUTTA_REQUIRE(out,"mesh: cannot write "+path.string());
  out << "(0 \"Kutta external aero mesher v1\")\n(2 3)\n" << std::hex;
  out << "(10 (0 1 " << mesh.n_points() << " 0 3))\n(10 (1 1 " << mesh.n_points() << " 1 3)(\n" << std::setprecision(17);
  for(const auto p:mesh.points) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
  out << "))\n(12 (0 1 " << mesh.n_cells << " 0))\n(12 (2 1 " << mesh.n_cells << " 1 7))\n";
  out << "(13 (0 1 " << mesh.n_faces() << " 0))\n";
  auto faces=[&](int zone,Index first,Index count,int bc) {
    if(!count) return;
    out << "(13 (" << zone << ' ' << first+1 << ' ' << first+count << ' ' << bc << " 0)(\n";
    for(Index f=first;f<first+count;++f) {
      const auto ns=detail::nodes(mesh,f); out << ns.size();
      // Fluent's stored face normal points from c0 into the cell (the reader
      // reverses a boundary c0 face); reverse our outward owner ordering.
      for(auto it=ns.rbegin();it!=ns.rend();++it) out << ' ' << *it+1;
      out << ' ' << mesh.owner[f]+1 << ' ' << (f<mesh.n_internal_faces()?mesh.neighbour[f]+1:0) << '\n';
    }
    out << "))\n";
  };
  faces(3,0,mesh.n_internal_faces(),2);
  int zone=4;
  for(const auto& p:mesh.patches) faces(zone++,p.start,p.size,p.kind==PatchKind::VelocityInlet?10:p.kind==PatchKind::PressureOutlet?5:p.kind==PatchKind::Symmetry?7:3);
  out << std::dec << "(45 (2 fluid fluid)())\n";
  if(mesh.n_internal_faces()) out << "(45 (3 interior interior)())\n";
  zone=4;
  for(const auto& p:mesh.patches) {
    // Scheme quoted strings preserve STL names containing spaces and quotes.
    std::string name; for(char c:p.name) { if(c=='"' || c=='\\') name+='\\'; KUTTA_REQUIRE(static_cast<unsigned char>(c)>=32,"mesh: control character in zone name"); name+=c; }
    out << "(45 (" << zone++ << ' ' << p.source_type << " \"" << name << "\")())\n";
  }
  out.close(); KUTTA_REQUIRE(out,"mesh: write failed "+path.string());
}
}
