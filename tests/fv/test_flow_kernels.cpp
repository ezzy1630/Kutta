#include "../backend_helpers.hpp"
#include "kernels.hpp"
#include "kutta/fv/operators.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
using namespace kutta::detail;

TEST_CASE("FV momentum face interpolation is affine exact with unequal diagonals") {
  auto mesh=make_box_tet(3,3,3,{0,0,0},{1,1,1});
  const auto g=compute_geometry(mesh);
  const auto value=[](Vec3d p) {return float(1+.3*p.x-.2*p.y+.1*p.z);};
  std::vector<float> cell,boundary,diagonal;
  for(Index c=0;c<mesh.n_cells;++c) {cell.push_back(value(g.cell_centroid[c]));diagonal.push_back(float(1+c%7));}
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) boundary.push_back(value(g.face_centroid[f]));
  std::vector<float> baseline;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary),ap=dev->upload<float>(diagonal);
    auto gradient=fv::vector_field(*dev,mesh.n_cells);
    dev->fill(gradient.x,.3f);dev->fill(gradient.y,-.2f);dev->fill(gradient.z,.1f);
    auto face=fv::vector_field(*dev,mesh.n_faces());
    kern::fv_momentum_faces(*dev,m.n_faces,m.owner,m.neighbour,
      m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
      m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
      q,q,q,gradient.x,gradient.y,gradient.z,gradient.x,gradient.y,gradient.z,
      gradient.x,gradient.y,gradient.z,ap,bc,bc,bc,face.x,face.y,face.z,m.n_internal_faces);
    auto got=dev->download(face.x);
    for(Index f=0;f<m.n_faces;++f) CHECK(got[f]==doctest::Approx(value(g.face_centroid[f])).epsilon(1e-6));
    CHECK(dev->download(face.y)==got);CHECK(dev->download(face.z)==got);
    if(baseline.empty()) baseline=got;else CHECK(got==baseline);
  }
}

TEST_CASE("FV pressure symmetry boundary extrapolates only within its plane") {
  auto mesh=make_box_tet(2,2,2,{0,0,0},{1,1,1});
  const auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);
    auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,3);
    auto p=dev->alloc<float>(m.n_cells);dev->fill(p,2);
    auto gradient=fv::vector_field(*dev,m.n_cells);
    dev->fill(gradient.x,1);dev->fill(gradient.y,2);dev->fill(gradient.z,3);
    auto prescribed=dev->alloc<float>(m.n_boundary_faces()),bc=dev->alloc<float>(m.n_boundary_faces());
    dev->fill(prescribed,0);
    kern::fv_pressure_boundary(*dev,m.n_boundary_faces(),m.owner,kind,m.area.x,m.area.y,m.area.z,
      m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,p,gradient.x,gradient.y,gradient.z,
      prescribed,bc,m.n_internal_faces);
    auto got=dev->download(bc);
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) {
      const auto n=g.face_area[f]/norm(g.face_area[f]);
      auto r=g.face_centroid[f]-g.cell_centroid[mesh.owner[f]];
      r=r-n*dot(r,n);
      CHECK(got[f-m.n_internal_faces]==doctest::Approx(2+dot(Vec3d{1,2,3},r)).epsilon(1e-6));
    }
  }
}
