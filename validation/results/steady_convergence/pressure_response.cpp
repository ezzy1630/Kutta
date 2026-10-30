// Offline linear response of the existing pressure interpolation at a captured
// momentum diagonal. Host loops belong to this diagnostic, never the solve path.
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/topology.hpp"
#include "kernels.hpp"

int main(int argc,char** argv) {
  using namespace kutta;using namespace fv;
  try {
    KUTTA_REQUIRE(argc==4,"usage: pressure_response grid cells.csv output-dir");
    auto host=read_plot3d(argv[1],.01);renumber_cells_rcm(host);auto geom=compute_geometry(host);
    const int n=host.n_cells;KUTTA_REQUIRE(n>0 && n<=1024,"diagnostic is bounded to coarse grids");
    std::vector<float> diagonal,inverse,captured_pressure;
    std::ifstream cells(argv[2]);KUTTA_REQUIRE(cells,"cannot open captured cells");std::string line;std::getline(cells,line);
    while(std::getline(cells,line)) {
      std::istringstream in(line);std::vector<std::string> row;std::string value;
      while(std::getline(in,value,',')) row.push_back(value);
      KUTTA_REQUIRE(row.size()==6,"unexpected captured-cell schema");
      const auto c=diagonal.size();KUTTA_REQUIRE(c<std::size_t(n) && std::stoi(row[0])==int(c),"captured cell order mismatch");
      KUTTA_REQUIRE(std::abs(std::stod(row[1])-geom.cell_centroid[c].x)<1e-10 && std::abs(std::stod(row[2])-geom.cell_centroid[c].y)<1e-10,"captured geometry mismatch");
      diagonal.push_back(float(std::stod(row[4])));inverse.push_back(float(std::stod(row[5])));
      captured_pressure.push_back(float(std::stod(row[3])));
    }
    KUTTA_REQUIRE(diagonal.size()==std::size_t(n),"captured cell count mismatch");
    auto dev=Device::create(Backend::Cpu);auto mesh=upload_mesh(*dev,host,geom);auto stencil=make_stencil(*dev,host,geom);
    std::vector<Index> kinds(mesh.n_boundary_faces());std::vector<BoundaryObservation> observations(kinds.size());
    for(const auto& p:host.patches) for(int f=p.start;f<p.start+p.size;++f) {
      const int b=f-mesh.n_internal_faces;
      kinds[b]=p.kind==PatchKind::Wall?0:p.source_type=="subsonic_inflow_pt"?1:(p.source_type=="back_pressure" || p.source_type=="farfield_riem")?2:3;
      observations[b]=kinds[b]==2?BoundaryObservation::FaceValue:BoundaryObservation::Omitted;
    }
    auto fit=boundary_least_squares(*dev,host,geom,observations);
    auto kind=dev->upload<Index>(kinds);
    auto ap=dev->upload<float>(diagonal),ci=dev->upload<float>(inverse);
    auto zero=dev->alloc<float>(mesh.n_faces),one=dev->alloc<float>(mesh.n_faces),mass=dev->alloc<float>(mesh.n_faces),coefficient=dev->alloc<float>(mesh.n_faces);
    dev->fill(zero,0);dev->fill(one,1);
    // A unit pressure jump isolates the exact stored physical/correction face
    // coefficients from the generated CPU kernel.
    kern::fv_rhie_chow(*dev,mesh.n_faces,mesh.owner,mesh.neighbour,kind,mesh.area.x,mesh.area.y,mesh.area.z,
      stencil.diffusion,mesh.volume,ap,ci,zero,zero,zero,one,mass,coefficient,mesh.n_internal_faces,1);
    auto cr=dev->download(mass);auto matrix=upload(*dev,stencil.pattern);auto rhs=dev->alloc<float>(n);
    kern::fv_pressure_assemble(*dev,n,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,stencil.adjacent,
      stencil.matrix_entry,matrix.row_ptr,matrix.diag,coefficient,zero,matrix.val,rhs,-1);
    auto values=dev->download(matrix.val);
    const auto offsets=dev->download(mesh.cell_face_offsets),faces=dev->download(mesh.cell_faces),adjacent=dev->download(stencil.adjacent);
    const auto signs=dev->download(mesh.cell_face_sign),lsq=dev->download(fit.inverse);
    const auto dx=dev->download(fit.delta.x),dy=dev->download(fit.delta.y),dz=dev->download(fit.delta.z);
    const auto px=dev->download(mesh.owner_to_face.x),py=dev->download(mesh.owner_to_face.y),pz=dev->download(mesh.owner_to_face.z);
    const auto nx=dev->download(mesh.neighbour_to_face.x),ny=dev->download(mesh.neighbour_to_face.y),nz=dev->download(mesh.neighbour_to_face.z);
    std::vector<double> J(std::size_t(n)*n),L(J.size()),p(n),gx(n),gy(n),gz(n),flux(mesh.n_faces);
    for(int c=0;c<n;++c) for(int j=stencil.pattern.row_ptr[c];j<stencil.pattern.row_ptr[c+1];++j) L[std::size_t(c)*n+stencil.pattern.col[j]]=values[j];
    for(int col=0;col<n;++col) {
      std::fill(p.begin(),p.end(),0);p[col]=1;
      for(int c=0;c<n;++c) {
        double x=0,y=0,z=0;
        for(int j=offsets[c];j<offsets[c+1];++j) {const int a=adjacent[j];const double delta=(a>=0?p[a]:0)-p[c];x+=dx[j]*delta;y+=dy[j]*delta;z+=dz[j]*delta;}
        gx[c]=lsq[c]*x+lsq[c+n]*y+lsq[c+2*n]*z;
        gy[c]=lsq[c+n]*x+lsq[c+3*n]*y+lsq[c+4*n]*z;
        gz[c]=lsq[c+2*n]*x+lsq[c+4*n]*y+lsq[c+5*n]*z;
      }
      for(int f=0;f<mesh.n_faces;++f) {
        if(cr[f]==0) {flux[f]=0;continue;}
        const int c=host.owner[f];double jump=p[c]+gx[c]*px[f]+gy[c]*py[f]+gz[c]*pz[f];
        if(f<mesh.n_internal_faces) {
          const int a=host.neighbour[f];const double side=jump-(p[a]+gx[a]*nx[f]+gy[a]*ny[f]+gz[a]*nz[f]);
          const float x=px[f]-nx[f],y=py[f]-ny[f],z=pz[f]-nz[f];
          const float sx=.5f*(px[f]+nx[f]),sy=.5f*(py[f]+ny[f]),sz=.5f*(pz[f]+nz[f]);
          const float d2=x*x+y*y+z*z,s2=sx*sx+sy*sy+sz*sz;
          const float trust=std::clamp((4*d2-s2)/(3*d2),0.f,1.f);
          const double shared=p[c]-p[a]+.5*((gx[c]+gx[a])*x+(gy[c]+gy[a])*y+(gz[c]+gz[a])*z);
          jump=trust==1?side:shared+trust*(side-shared);
        }
        flux[f]=cr[f]*jump;
      }
      for(int c=0;c<n;++c) {double sum=0;for(int j=offsets[c];j<offsets[c+1];++j) sum+=signs[j]*flux[faces[j]];J[std::size_t(c)*n+col]=sum;}
    }
    const std::filesystem::path out=argv[3];std::filesystem::create_directories(out);
    std::vector<float> probe(n),boundary(kinds.size(),0);for(int c=0;c<n;++c) probe[c]=float(std::sin(.37*c));
    auto field=dev->upload<float>(probe),bc=dev->upload<float>(boundary);auto grad=vector_field(*dev,n);
    auto face=dev->alloc<float>(mesh.n_faces),jump=dev->alloc<float>(mesh.n_faces);
    gradient(*dev,mesh,stencil,fit,field,bc,grad);
    kern::fv_pressure_faces(*dev,mesh.n_faces,mesh.owner,mesh.neighbour,mesh.owner_to_face.x,mesh.owner_to_face.y,mesh.owner_to_face.z,
      mesh.neighbour_to_face.x,mesh.neighbour_to_face.y,mesh.neighbour_to_face.z,field,grad.x,grad.y,grad.z,bc,face,jump,mesh.n_internal_faces);
    kern::fv_rhie_chow(*dev,mesh.n_faces,mesh.owner,mesh.neighbour,kind,mesh.area.x,mesh.area.y,mesh.area.z,
      stencil.diffusion,mesh.volume,ap,ci,zero,zero,zero,jump,mass,coefficient,mesh.n_internal_faces,1);
    kern::fv_pressure_assemble(*dev,n,mesh.cell_face_offsets,mesh.cell_faces,mesh.cell_face_sign,stencil.adjacent,
      stencil.matrix_entry,matrix.row_ptr,matrix.diag,coefficient,mass,matrix.val,rhs,-1);
    const auto actual=dev->download(rhs);double error=0,norm=0;
    for(int i=0;i<n;++i) {double expected=0;for(int j=0;j<n;++j) expected+=J[std::size_t(i)*n+j]*probe[j];error+=std::pow(expected+actual[i],2);norm+=std::pow(expected,2);}
    const double relative=std::sqrt(error/norm);KUTTA_REQUIRE(relative<3e-5,"offline pressure response mismatches generated kernels");
    std::ofstream proof(out/"action-verification.json");proof<<std::setprecision(17)<<"{\"rows\":"<<n<<",\"generated_cpu_relative_l2\":"<<relative<<"}\n";
    auto write=[&](const char* name,const std::vector<double>& data) {std::ofstream stream(out/name);stream<<std::setprecision(17);for(int i=0;i<n;++i) {for(int j=0;j<n;++j) stream<<(j?",":"")<<data[std::size_t(i)*n+j];stream<<'\n';}};
    write("pressure_response.csv",J);write("pressure_correction.csv",L);
    std::ofstream alpha(out/"pressure_alpha.csv");alpha<<std::setprecision(17);for(int c=0;c<n;++c) alpha<<std::min(1.f,.7f*diagonal[c]*inverse[c])<<'\n';
    std::ofstream state(out/"pressure_state.csv");state<<std::setprecision(17);for(float value:captured_pressure) state<<value<<'\n';
    std::ofstream physical(out/"rhie_chow_physical.csv");physical<<std::setprecision(17);for(float value:cr) physical<<value<<'\n';
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
