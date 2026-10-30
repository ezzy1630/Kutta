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
    auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);
    auto gradient=fv::vector_field(*dev,mesh.n_cells);
    dev->fill(gradient.x,.3f);dev->fill(gradient.y,-.2f);dev->fill(gradient.z,.1f);
    auto face=fv::vector_field(*dev,mesh.n_faces());
    kern::fv_momentum_faces(*dev,m.n_faces,m.owner,m.neighbour,kind,
      m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
      m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
      q,q,q,gradient.x,gradient.y,gradient.z,gradient.x,gradient.y,gradient.z,
      gradient.x,gradient.y,gradient.z,ap,bc,bc,bc,face.x,face.y,face.z,m.n_internal_faces);
    auto got=dev->download(face.x);
    for(Index f=0;f<m.n_faces;++f) CHECK(got[f]==doctest::Approx(value(g.face_centroid[f])).epsilon(1e-6));
    CHECK(dev->download(face.y)==got);CHECK(dev->download(face.z)==got);
    if(baseline.empty()) baseline=got;else CHECK(got==baseline);
    // Non-affine side gradients distinguish independent reconstructions from
    // the old weighted-point/mean-gradient formula. Outlets use the cell value.
    std::vector<float> gx(mesh.n_cells),gy(mesh.n_cells),gz(mesh.n_cells);
    for(Index c=0;c<mesh.n_cells;++c) {gx[c]=float(c%5)*.2f;gy[c]=float(c%3)*.3f;gz[c]=float(c%7)*.1f;}
    dev->upload<float>(gradient.x,gx);dev->upload<float>(gradient.y,gy);dev->upload<float>(gradient.z,gz);
    dev->fill(kind,2);
    kern::fv_momentum_faces(*dev,m.n_faces,m.owner,m.neighbour,kind,
      m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
      m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
      q,q,q,gradient.x,gradient.y,gradient.z,gradient.x,gradient.y,gradient.z,
      gradient.x,gradient.y,gradient.z,ap,bc,bc,bc,face.x,face.y,face.z,m.n_internal_faces);
    got=dev->download(face.x);
    for(Index f=0;f<m.n_faces;++f) {
      const auto c=mesh.owner[f];double expected=cell[c];
      if(f<m.n_internal_faces) {
        const auto n=mesh.neighbour[f];const double weight=diagonal[c]/(diagonal[c]+diagonal[n]);
        const auto rc=g.face_centroid[f]-g.cell_centroid[c],rn=g.face_centroid[f]-g.cell_centroid[n];
        expected=weight*(cell[c]+dot(Vec3d{gx[c],gy[c],gz[c]},rc))
          +(1-weight)*(cell[n]+dot(Vec3d{gx[n],gy[n],gz[n]},rn));
      }
      CHECK(got[f]==doctest::Approx(expected).epsilon(1e-6));
    }
  }
}

TEST_CASE("FV pressure symmetry boundary uses the full centroid displacement") {
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
      auto r=g.face_centroid[f]-g.cell_centroid[mesh.owner[f]];
      CHECK(got[f-m.n_internal_faces]==doctest::Approx(2+dot(Vec3d{1,2,3},r)).epsilon(1e-6));
    }
  }
}
