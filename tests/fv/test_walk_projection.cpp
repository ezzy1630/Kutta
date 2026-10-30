#include <array>
#include <iostream>
#include <limits>
#include <map>
#include "../backend_helpers.hpp"
#include "physics_meshes.hpp"
#include "kernels.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/mesh/topology.hpp"

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
struct Maximum {
  double value=0;
  void add(double e) {value=std::isfinite(e)?std::max(value,std::abs(e)):std::numeric_limits<double>::infinity();}
};
std::vector<float> legacy_faces(Device& dev,const DeviceMesh& m,const Buffer<Index>& kind,
    const Buffer<float>& q,const Buffer<float>& bc,const DeviceVec3& g,const Buffer<float>& ap,Buffer<float>& out) {
  auto spare1=dev.alloc<float>(m.n_faces),spare2=dev.alloc<float>(m.n_faces);
  kern::fv_momentum_faces(dev,m.n_faces,m.owner,m.neighbour,kind,
    m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
    q,q,q,g.x,g.y,g.z,g.x,g.y,g.z,g.x,g.y,g.z,ap,bc,bc,bc,out,spare1,spare2,m.n_internal_faces);
  return dev.download(out);
}
using Row=std::array<double,10>;
Row basis(Vec3d r) {return {r.x,r.y,r.z,.5*r.x*r.x,r.x*r.y,r.x*r.z,.5*r.y*r.y,r.y*r.z,.5*r.z*r.z,1};}
// Independent oracle: enumerate the walk-weighted cloud of one face and solve
// its weighted normal equations by Gaussian elimination in double.
Row walk_fit(const PolyMesh& mesh,const MeshGeometry& g,const kutta::CellFaces& cf,Index face,
    const std::vector<double>& cell,const std::vector<double>& boundary,double anchor) {
  const Index internal=mesh.n_internal_faces();const auto origin=g.face_centroid[face];
  std::array<std::array<double,11>,10> a{};
  for(Index start:{mesh.owner[face],face<internal?mesh.neighbour[face]:Index(-1)}) {
    if(start<0) continue;
    std::map<Index,double> walks{{start,1.}};
    for(int hop=0;hop<WalkFaceReconstruction::hops;++hop) {
      auto next=walks;
      for(auto [c,w]:walks) for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
        const auto f=cf.faces[j];if(f<internal) next[mesh.owner[f]==c?mesh.neighbour[f]:mesh.owner[f]]+=w;
      }
      walks=std::move(next);
    }
    for(auto [c,w]:walks) {
      const auto add=[&](Vec3d x,double q) {
        const auto b=basis(x-origin);
        for(int i=0;i<10;++i) {for(int j=0;j<10;++j) a[i][j]+=w*b[i]*b[j];a[i][10]+=w*b[i]*(q-anchor);}
      };
      add(g.cell_centroid[c],cell[c]);
      for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) if(cf.faces[j]>=internal) add(g.face_centroid[cf.faces[j]],boundary[cf.faces[j]-internal]);
    }
  }
  for(int i=0;i<10;++i) {
    int pivot=i;for(int r=i+1;r<10;++r) if(std::abs(a[r][i])>std::abs(a[pivot][i])) pivot=r;
    std::swap(a[i],a[pivot]);
    for(int r=0;r<10;++r) if(r!=i) {const double s=a[r][i]/a[i][i];for(int k=i;k<11;++k) a[r][k]-=s*a[i][k];}
  }
  Row beta{};for(int i=0;i<10;++i) beta[i]=a[i][10]/a[i][i];
  return beta;
}
}

TEST_CASE("FV walk projection preserves the normal stencil for quadratic fields") {
  KUTTA_FOR_EACH_BACKEND(dev) for(int variant=0;variant<6;++variant) {
    const auto mesh=variant==0?make_box_hex(4,4,4,{0,0,0},{1,1,1}):variant==1?make_box_tet(3,3,3,{0,0,0},{1,1,1}):
      test::physics_annulus(8<<std::min(variant-2,2),8*(8<<std::min(variant-2,2)),true);
    const auto geo=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
    const auto value=[](Vec3d x) {return x.x*x.x-2*x.x*x.y+3*x.y*x.y+.5*x.z*x.z+.3*x.x*x.z-.7*x.y*x.z+.2*x.x-.4*x.y+.6*x.z;};
    const auto derivative=[](Vec3d x) {return Vec3d{2*x.x-2*x.y+.3*x.z+.2,-2*x.x+6*x.y-.7*x.z-.4,x.z+.3*x.x-.7*x.y+.6};};
    std::vector<float> cell,boundary,diagonal,dx,dy,dz;
    for(auto x:geo.cell_centroid) {cell.push_back(float(value(x)));const auto g=derivative(x);dx.push_back(float(g.x));dy.push_back(float(g.y));dz.push_back(float(g.z));diagonal.push_back(1+.1f*(cell.size()%17));}
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(geo.face_centroid[f])));
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary),ap=dev->upload<float>(diagonal);
    // Variant 5 omits every boundary value, as wall functions do for walls.
    const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),variant==5?BoundaryObservation::Omitted:BoundaryObservation::FaceValue);
    auto reconstruction=std::make_shared<WalkFaceReconstruction>(*dev,mesh,geo,m,st,rows);
    INFO("backend "<<backend_name(dev->backend())<<" variant "<<variant);
    CHECK(reconstruction->fallback_faces()==0);
    if(variant==0) CHECK(reconstruction->corrected_faces()==0);
    else CHECK(reconstruction->corrected_faces()>0);
    auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);
    auto zero=dev->alloc<float>(m.n_cells),mass=dev->alloc<float>(m.n_faces),mu=dev->alloc<float>(m.n_faces);
    dev->fill(zero,0);dev->fill(mass,0);dev->fill(mu,1);
    TransportEquation equation(*dev,m,st,mesh,geo,{},reconstruction);TransportSettings controls;controls.limit=false;
    dev->upload<float>(equation.grad.x,dx);dev->upload<float>(equation.grad.y,dy);dev->upload<float>(equation.grad.z,dz);
    equation.assemble_reconstructed(q,bc,kind,mass,mu,zero,controls);
    auto face_values=dev->alloc<float>(m.n_faces);
    const auto predictor=legacy_faces(*dev,m,kind,q,bc,equation.grad,ap,face_values);
    // A zero cell limiter restores the original predictor exactly.
    auto theta=dev->alloc<float>(m.n_cells);dev->fill(theta,0);
    reconstruction->momentum_values(q,bc,ap,face_values,&theta);
    CHECK(dev->download(face_values)==predictor);
    reconstruction->momentum_values(q,bc,ap,face_values);
    const auto face=dev->download(face_values),correction=dev->download(equation.face_correction),diffusion=dev->download(st.diffusion);
    const auto mask=dev->download(reconstruction->mask());
    Maximum flux_error,value_error;
    for(Index f=0;f<m.n_faces;++f) {
      if(f>=m.n_internal_faces && variant==5) continue; // unobserved boundary faces keep the original stencil
      const auto p=mesh.owner[f];const auto area=geo.face_area[f],x=geo.face_centroid[f],n=area/norm(area);
      const double dp=dot(x-geo.cell_centroid[p],n),dn=f<m.n_internal_faces?dot(geo.cell_centroid[mesh.neighbour[f]]-x,n):0;
      const auto pf=x-n*dp,nf=x+n*dn;
      const double jump=(f<m.n_internal_faces?cell[mesh.neighbour[f]]:boundary[f-m.n_internal_faces])-double(cell[p]);
      flux_error.add((double(diffusion[f])*jump+correction[f])/norm(area)-(value(nf)-value(pf))/(dp+dn));
      if(f<m.n_internal_faces) value_error.add(face[f]-(dn*value(pf)+dp*value(nf))/(dp+dn));
    }
    std::cout<<"WALK_QUADRATIC backend="<<backend_name(dev->backend())<<" variant="<<variant<<" corrected="<<reconstruction->corrected_faces()
      <<" normal_stencil_Linf="<<flux_error.value<<" normal_interpolation_Linf="<<value_error.value<<std::endl;
    CHECK(flux_error.value<2e-4);CHECK(value_error.value<2e-5);
  }
}

TEST_CASE("FV walk projection moments match an explicit walk-weighted fit") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool poly:{false,true}) {
    const auto mesh=poly?test::physics_annulus(8,64,true):make_box_tet(3,3,3,{0,0,0},{1,1,1});
    const auto geo=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
    const auto cf=build_cell_faces(mesh);
    const auto value=[](Vec3d x) {return std::sin(x.x+2*x.y)+.3*x.x*x.z*x.z+std::exp(.5*x.y);};
    std::vector<float> cell,boundary,diagonal;std::vector<double> cell_d,boundary_d;
    for(auto x:geo.cell_centroid) {cell.push_back(float(value(x)));cell_d.push_back(cell.back());diagonal.push_back(1+.1f*(cell.size()%13));}
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) {boundary.push_back(float(value(geo.face_centroid[f])));boundary_d.push_back(boundary.back());}
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary),ap=dev->upload<float>(diagonal),out=dev->alloc<float>(m.n_faces);
    const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    WalkFaceReconstruction fit(*dev,mesh,geo,m,st,rows);
    fit.diffusion_correction(q,bc,out);const auto correction=dev->download(out);
    auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);auto grad=vector_field(*dev,m.n_cells);gradient(*dev,m,st,q,bc,grad);
    legacy_faces(*dev,m,kind,q,bc,grad,ap,out);fit.momentum_values(q,bc,ap,out);const auto face=dev->download(out);
    const auto mask=dev->download(fit.mask());
    Maximum correction_error,value_error;double correction_scale=0;
    for(Index f=0;f<m.n_faces;++f) {
      if(!mask[f]) {CHECK(correction[f]==0);continue;}
      const auto p=mesh.owner[f];const bool interior=f<m.n_internal_faces;
      const auto area=geo.face_area[f],x=geo.face_centroid[f],n=area/norm(area);
      const auto other=interior?geo.cell_centroid[mesh.neighbour[f]]:x;
      const double dp=dot(x-geo.cell_centroid[p],n),dn=dot(other-x,n),alpha=dot(area,area)/dot(area,other-geo.cell_centroid[p]);
      const double anchor=interior?.5*(cell_d[p]+cell_d[mesh.neighbour[f]]):cell_d[p];
      const auto beta=walk_fit(mesh,geo,cf,f,cell_d,boundary_d,anchor);
      const auto at=[&](Vec3d point) {const auto b=basis(point-x);double s=0;for(int k=0;k<10;++k) s+=beta[k]*b[k];return s;};
      const double expected=alpha*(at(x+n*dn)-at(other)-at(x-n*dp)+at(geo.cell_centroid[p]));
      correction_error.add(correction[f]-expected);correction_scale=std::max(correction_scale,std::abs(expected));
      if(interior) {
        const double w=double(diagonal[p])/(double(diagonal[p])+diagonal[mesh.neighbour[f]]),qn=cell_d[mesh.neighbour[f]];
        const double qstar=(dn*at(x-n*dp)+dp*at(x+n*dn))/(dp+dn);
        const double predicted=w*cell_d[p]+(1-w)*qn+(qstar-at(geo.cell_centroid[p]))-(1-w)*(at(other)-at(geo.cell_centroid[p]));
        value_error.add(face[f]-predicted);
      }
    }
    INFO("backend "<<backend_name(dev->backend())<<" poly "<<poly);
    std::cout<<"WALK_ORACLE backend="<<backend_name(dev->backend())<<" poly="<<poly<<" corrected="<<fit.corrected_faces()
      <<" correction_Linf="<<correction_error.value<<" correction_scale="<<correction_scale<<" value_Linf="<<value_error.value<<std::endl;
    CHECK(fit.corrected_faces()>0);
    CHECK(correction_error.value<1e-4*correction_scale+1e-6);CHECK(value_error.value<2e-6);
  }
}

TEST_CASE("FV walk projection is invariant under owner neighbour reversal") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool poly:{false,true}) {
    const auto original=poly?test::physics_annulus(8,64,true):make_box_tet(3,3,3,{0,0,0},{1,1,1});
    auto reversed=original;
    for(auto& c:reversed.owner) c=reversed.n_cells-1-c;
    for(auto& c:reversed.neighbour) c=reversed.n_cells-1-c;
    for(auto& zone:reversed.cell_zones) for(auto& c:zone.cells) c=reversed.n_cells-1-c;
    for(Index f=0;f<reversed.n_internal_faces();++f) {
      std::swap(reversed.owner[f],reversed.neighbour[f]);
      std::reverse(reversed.face_nodes.begin()+reversed.face_offsets[f],reversed.face_nodes.begin()+reversed.face_offsets[f+1]);
    }
    std::vector<float> corrections[2],faces[2];
    for(int side=0;side<2;++side) {
      const auto& mesh=side?reversed:original;const auto geo=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,geo);const auto st=make_stencil(*dev,mesh,geo);
      auto value=[](Vec3d x){return std::sin(x.x+2*x.y)+.3*x.x*x.z*x.z;};
      std::vector<float> cells,boundary,diagonal;
      for(auto x:geo.cell_centroid) {cells.push_back(float(value(x)));diagonal.push_back(1+.1f*((side?mesh.n_cells-Index(cells.size())+1:Index(cells.size()))%17));}
      for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(geo.face_centroid[f])));
      auto q=dev->upload<float>(cells),bc=dev->upload<float>(boundary),ap=dev->upload<float>(diagonal),c=dev->alloc<float>(m.n_faces),out=dev->alloc<float>(m.n_faces);
      auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);auto grad=vector_field(*dev,m.n_cells);gradient(*dev,m,st,q,bc,grad);
      const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
      WalkFaceReconstruction fit(*dev,mesh,geo,m,st,rows);fit.diffusion_correction(q,bc,c);
      legacy_faces(*dev,m,kind,q,bc,grad,ap,out);fit.momentum_values(q,bc,ap,out);
      corrections[side]=dev->download(c);faces[side]=dev->download(out);
    }
    Maximum correction_error,value_error;
    for(Index f=0;f<original.n_faces();++f) {
      correction_error.add(corrections[0][f]+(f<original.n_internal_faces()?1:-1)*corrections[1][f]);
      value_error.add(faces[0][f]-faces[1][f]);
    }
    INFO("backend "<<backend_name(dev->backend())<<" poly "<<poly);
    CHECK(correction_error.value<2e-6);CHECK(value_error.value<2e-6);
  }
}
