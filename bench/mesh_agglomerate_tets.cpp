// Convex tet-pair agglomeration for native-mesh validation. Every boundary
// face and every non-tet cell is retained; this is not Fluent Meshing output.
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <tuple>

using namespace kutta;
using namespace kutta::detail;

namespace {
bool is_tet(const PolyMesh& mesh,const kutta::CellFaces& cells,Index cell) {
  if(cells.offsets[cell+1]-cells.offsets[cell]!=4) return false;
  for(Index j=cells.offsets[cell];j<cells.offsets[cell+1];++j) {
    const auto face=cells.faces[j];
    if(mesh.face_offsets[face+1]-mesh.face_offsets[face]!=3) return false;
  }
  return true;
}

bool convex_pair(const PolyMesh& mesh,const MeshGeometry& geometry,
                 const kutta::CellFaces& cells,Index shared) {
  const auto owner=mesh.owner[shared],neighbour=mesh.neighbour[shared];
  std::set<Index> vertices;
  for(auto cell:{owner,neighbour})
    for(Index j=cells.offsets[cell];j<cells.offsets[cell+1];++j)
      for(auto vertex:nodes(mesh,cells.faces[j])) vertices.insert(vertex);
  for(auto cell:{owner,neighbour}) {
    for(Index j=cells.offsets[cell];j<cells.offsets[cell+1];++j) {
      const auto face=cells.faces[j];
      if(face==shared) continue;
      const auto normal=geometry.face_area[face]*double(cells.signs[j]);
      for(auto vertex:vertices) {
        const auto delta=mesh.points[vertex]-geometry.face_centroid[face];
        // Only convex unions have a cell centroid inside every face plane.
        if(dot(normal,delta)>1e-12*norm(normal)*std::max(norm(delta),1e-12)) return false;
      }
    }
  }
  return true;
}

struct BoundaryType {int code;const char* name;};
BoundaryType boundary_type(PatchKind kind) {
  switch(kind) {
    case PatchKind::Wall:return {3,"wall"};
    case PatchKind::VelocityInlet:return {10,"velocity-inlet"};
    case PatchKind::PressureOutlet:return {5,"pressure-outlet"};
    case PatchKind::Symmetry:return {7,"symmetry"};
    case PatchKind::PressureInlet:
    case PatchKind::Farfield:
    case PatchKind::Other:break;
  }
  KUTTA_REQUIRE(false,"tet agglomeration: unsupported boundary kind");
  return {};
}

int cell_type(const PolyMesh& mesh,const kutta::CellFaces& cells,Index cell) {
  Index triangles=0,quads=0;
  const auto count=cells.offsets[cell+1]-cells.offsets[cell];
  for(Index j=cells.offsets[cell];j<cells.offsets[cell+1];++j) {
    const auto face=cells.faces[j],size=mesh.face_offsets[face+1]-mesh.face_offsets[face];
    triangles+=size==3;quads+=size==4;
  }
  if(triangles==4 && count==4) return 2;
  if(triangles==2 && quads==3 && count==5) return 6;
  if(quads==6 && count==6) return 4;
  return 7;
}

void write_mesh(const PolyMesh& mesh,const std::string& fluid,const std::filesystem::path& path) {
  const auto cells=build_cell_faces(mesh);
  std::ofstream file(path);
  KUTTA_REQUIRE(bool(file),"tet agglomeration: cannot open output");
  file<<std::setprecision(17)<<"(0 \"Convex tet-pair agglomeration; original surface and prisms retained\")\n(2 3)\n";
  file<<"(10 (0 1 "<<std::hex<<mesh.points.size()<<" 0 3))\n(10 (1 1 "<<mesh.points.size()<<" 1 3)(\n";
  for(auto point:mesh.points) file<<std::dec<<point.x<<' '<<point.y<<' '<<point.z<<'\n';
  file<<"))\n(12 (0 1 "<<std::hex<<mesh.n_cells<<" 0))\n(12 (2 1 "<<mesh.n_cells<<" 1 0)(\n";
  for(Index cell=0;cell<mesh.n_cells;++cell) file<<cell_type(mesh,cells,cell)<<'\n';
  file<<"))\n(13 (0 1 "<<mesh.n_faces()<<" 0))\n";
  const auto zone=[&](int id,Index start,Index count,int kind) {
    file<<"(13 ("<<std::hex<<id<<' '<<start+1<<' '<<start+count<<' '<<kind<<" 0)(\n";
    for(Index face=start;face<start+count;++face) {
      file<<mesh.face_offsets[face+1]-mesh.face_offsets[face]<<' ';
      // Legacy Fluent winding points into c0; PolyMesh winding points out.
      for(auto j=mesh.face_offsets[face+1];j-->mesh.face_offsets[face];) file<<mesh.face_nodes[j]+1<<' ';
      file<<mesh.owner[face]+1<<' '<<(face<mesh.n_internal_faces()?mesh.neighbour[face]+1:0)<<'\n';
    }
    file<<"))\n";
  };
  zone(3,0,mesh.n_internal_faces(),2);
  int id=4;
  for(const auto& patch:mesh.patches) zone(id++,patch.start,patch.size,boundary_type(patch.kind).code);
  file<<"(45 (2 fluid "<<std::quoted(fluid)<<")())\n";
  id=4;
  for(const auto& patch:mesh.patches)
    file<<"(45 ("<<std::dec<<id++<<' '<<boundary_type(patch.kind).name<<' '<<std::quoted(patch.name)<<")())\n";
  file.close();
  KUTTA_REQUIRE(bool(file),"tet agglomeration: mesh write failed");
}

// Fluent sorts internal faces by owner/neighbour on import. Compare topology
// without interpreting that harmless ordering change as a connectivity change.
auto connectivity(const PolyMesh& mesh) {
  std::vector<std::tuple<Index,Index,std::vector<Index>>> faces;
  for(Index face=0;face<mesh.n_faces();++face) {
    const auto vertices=nodes(mesh,face);
    faces.emplace_back(mesh.owner[face],face<mesh.n_internal_faces()?mesh.neighbour[face]:-1,
                       std::vector<Index>(vertices.begin(),vertices.end()));
  }
  std::sort(faces.begin(),faces.end());
  return faces;
}
} // namespace

int main(int argc,char** argv) {
  try {
    KUTTA_REQUIRE(argc==3,"usage: mesh_agglomerate_tets input.msh new-output.msh");
    const std::filesystem::path output=argv[2];
    KUTTA_REQUIRE(!std::filesystem::exists(output),"tet agglomeration: output already exists");
    const auto mesh=read_fluent_msh(argv[1]);
    KUTTA_REQUIRE(mesh.cell_zones.size()<=1 && mesh.internal_face_zones.empty(),
                  "tet agglomeration: requires one fluid region without internal BCs");
    KUTTA_REQUIRE(mesh.cell_zones.empty() || mesh.cell_zones.front().source_type=="fluid",
                  "tet agglomeration: requires a fluid cell zone");
    for(const auto& patch:mesh.patches) (void)boundary_type(patch.kind);
    const auto geometry=compute_geometry(mesh);
    const auto cells=build_cell_faces(mesh);
    std::vector<bool> tet(mesh.n_cells),used(mesh.n_cells);
    std::vector<Index> candidates,representative(mesh.n_cells);
    for(Index cell=0;cell<mesh.n_cells;++cell) tet[cell]=is_tet(mesh,cells,cell);
    std::iota(representative.begin(),representative.end(),0);
    for(Index face=0;face<mesh.n_internal_faces();++face)
      if(tet[mesh.owner[face]] && tet[mesh.neighbour[face]]) candidates.push_back(face);
    std::stable_sort(candidates.begin(),candidates.end(),[&](Index a,Index b) {
      return dot(geometry.face_area[a],geometry.face_area[a])>dot(geometry.face_area[b],geometry.face_area[b]);
    });
    Index pairs=0;
    for(auto face:candidates) {
      const auto owner=mesh.owner[face],neighbour=mesh.neighbour[face];
      if(used[owner] || used[neighbour] || !convex_pair(mesh,geometry,cells,face)) continue;
      used[owner]=used[neighbour]=true;
      representative[neighbour]=owner;
      ++pairs;
    }
    PolyMesh merged;
    merged.points=mesh.points;
    std::vector<Index> ids(mesh.n_cells,-1),mapping(mesh.n_cells);
    for(Index cell=0;cell<mesh.n_cells;++cell)
      if(representative[cell]==cell) ids[cell]=merged.n_cells++;
    for(Index cell=0;cell<mesh.n_cells;++cell) mapping[cell]=ids[representative[cell]];
    const auto add=[&](Index face) {
      Index owner=mapping[mesh.owner[face]],neighbour=face<mesh.n_internal_faces()?mapping[mesh.neighbour[face]]:-1;
      if(owner==neighbour) return;
      const auto input=nodes(mesh,face);
      std::vector<Index> vertices(input.begin(),input.end());
      if(neighbour>=0 && owner>neighbour) {std::swap(owner,neighbour);std::reverse(vertices.begin(),vertices.end());}
      merged.face_nodes.insert(merged.face_nodes.end(),vertices.begin(),vertices.end());
      merged.face_offsets.push_back(Index(merged.face_nodes.size()));
      merged.owner.push_back(owner);
      if(neighbour>=0) merged.neighbour.push_back(neighbour);
    };
    for(Index face=0;face<mesh.n_internal_faces();++face) add(face);
    for(auto patch:mesh.patches) {
      const auto start=merged.n_faces();
      for(Index face=patch.start;face<patch.start+patch.size;++face) add(face);
      patch.start=start;
      merged.patches.push_back(patch);
    }
    validate(merged);
    const auto merged_geometry=compute_geometry(merged);
    std::vector<double> expected(merged.n_cells);
    for(Index cell=0;cell<mesh.n_cells;++cell) expected[mapping[cell]]+=geometry.cell_volume[cell];
    double volume_error=0;
    for(Index cell=0;cell<merged.n_cells;++cell)
      volume_error=std::max(volume_error,std::abs(merged_geometry.cell_volume[cell]-expected[cell])/expected[cell]);
    KUTTA_REQUIRE(volume_error<1e-10,"tet agglomeration: merged volumes changed");
    write_mesh(merged,mesh.cell_zones.empty()?"fluid":mesh.cell_zones.front().name,output);
    const auto check=read_fluent_msh(output);
    const auto checked_geometry=compute_geometry(check);
    KUTTA_REQUIRE(check.n_cells==merged.n_cells && connectivity(check)==connectivity(merged),
                  "tet agglomeration: round-trip connectivity changed");
    KUTTA_REQUIRE(check.points.size()==merged.points.size() && check.patches.size()==merged.patches.size(),
                  "tet agglomeration: round-trip points or patches changed");
    for(Index point=0;point<check.n_points();++point) {
      const auto a=check.points[point],b=merged.points[point];
      KUTTA_REQUIRE(a.x==b.x && a.y==b.y && a.z==b.z,"tet agglomeration: round-trip coordinates changed");
    }
    for(std::size_t i=0;i<check.patches.size();++i)
      KUTTA_REQUIRE(check.patches[i].name==merged.patches[i].name && check.patches[i].kind==merged.patches[i].kind &&
                    check.patches[i].size==merged.patches[i].size,"tet agglomeration: round-trip boundary changed");
    for(Index cell=0;cell<merged.n_cells;++cell)
      KUTTA_REQUIRE(std::abs(checked_geometry.cell_volume[cell]-expected[cell])/expected[cell]<1e-10,
                    "tet agglomeration: round-trip volumes changed");
    std::cout<<"original_cells="<<mesh.n_cells<<" cells="<<merged.n_cells<<" convex_pairs="<<pairs
             <<" relative_volume_error="<<volume_error
             <<" boundary_faces="<<merged.n_faces()-merged.n_internal_faces()<<'\n';
  } catch(const std::exception& error) {
    std::cerr<<error.what()<<'\n';
    return 1;
  }
}
