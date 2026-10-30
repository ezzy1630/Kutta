// Host-side view-export benchmark. Fields are manufactured, never a CFD solve.
#include "kutta/io/view.hpp"
#include "kutta/core/json.hpp"
#include "kutta/mesh/detail.hpp"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace {
// Direct structured topology keeps a 10M-cell scalability probe bounded in
// memory; the general-purpose test generator builds maps of every face.
kutta::PolyMesh structured_10m() {
  using namespace kutta;
  constexpr Index nx=250,ny=250,nz=160;
  PolyMesh m;m.n_cells=nx*ny*nz;
  const auto node=[](Index i,Index j,Index k){return i+(nx+1)*(j+(ny+1)*k);};
  const auto cell=[](Index i,Index j,Index k){return i+nx*(j+ny*k);};
  m.points.reserve(std::size_t(nx+1)*(ny+1)*(nz+1));
  for(Index k=0;k<=nz;++k)for(Index j=0;j<=ny;++j)for(Index i=0;i<=nx;++i)m.points.push_back({double(i)/nx,double(j)/ny,double(k)/nz});
  const auto nf=std::size_t(nx+1)*ny*nz+std::size_t(ny+1)*nx*nz+std::size_t(nz+1)*nx*ny;
  m.face_offsets.reserve(nf+1);m.face_nodes.reserve(nf*4);m.owner.reserve(nf);m.neighbour.reserve(nf);
  const auto face=[&](std::array<Index,4> nodes,Index owner,Index neighbour){m.face_nodes.insert(m.face_nodes.end(),nodes.begin(),nodes.end());m.face_offsets.push_back(Index(m.face_nodes.size()));m.owner.push_back(owner);if(neighbour>=0)m.neighbour.push_back(neighbour);};
  for(Index k=0;k<nz;++k)for(Index j=0;j<ny;++j)for(Index i=1;i<nx;++i)face({node(i,j,k),node(i,j+1,k),node(i,j+1,k+1),node(i,j,k+1)},cell(i-1,j,k),cell(i,j,k));
  for(Index k=0;k<nz;++k)for(Index j=1;j<ny;++j)for(Index i=0;i<nx;++i)face({node(i,j,k),node(i,j,k+1),node(i+1,j,k+1),node(i+1,j,k)},cell(i,j-1,k),cell(i,j,k));
  for(Index k=1;k<nz;++k)for(Index j=0;j<ny;++j)for(Index i=0;i<nx;++i)face({node(i,j,k),node(i+1,j,k),node(i+1,j+1,k),node(i,j+1,k)},cell(i,j,k-1),cell(i,j,k));
  for(int side=0;side<6;++side) {
    const auto start=m.n_faces();
    if(side<2){const auto i=side==0?0:nx;for(Index k=0;k<nz;++k)for(Index j=0;j<ny;++j){std::array ns{node(i,j,k),node(i,j+1,k),node(i,j+1,k+1),node(i,j,k+1)};if(side==0)std::reverse(ns.begin(),ns.end());face(ns,cell(side==0?0:nx-1,j,k),-1);}}
    else if(side<4){const auto j=side==2?0:ny;for(Index k=0;k<nz;++k)for(Index i=0;i<nx;++i){std::array ns{node(i,j,k),node(i,j,k+1),node(i+1,j,k+1),node(i+1,j,k)};if(side==2)std::reverse(ns.begin(),ns.end());face(ns,cell(i,side==2?0:ny-1,k),-1);}}
    else {const auto k=side==4?0:nz;for(Index j=0;j<ny;++j)for(Index i=0;i<nx;++i){std::array ns{node(i,j,k),node(i+1,j,k),node(i+1,j+1,k),node(i,j+1,k)};if(side==4)std::reverse(ns.begin(),ns.end());face(ns,cell(i,j,side==4?0:nz-1),-1);}}
    m.patches.push_back({"side-"+std::to_string(side),PatchKind::Wall,"wall",side,start,m.n_faces()-start});
  }
  return m;
}
}
int main(int argc,char** argv) {
  using Clock=std::chrono::steady_clock;
  if(argc!=3){std::cerr<<"Usage: bench_export MESH|--10m OUTPUT_DIR\n";return 2;}
  try {
    const auto start=Clock::now();kutta::view::Snapshot s;s.mesh=std::string_view(argv[1])=="--10m"?structured_10m():kutta::read_mesh(argv[1]);
    const auto loaded=Clock::now();
    for(const auto& [name,unit,components]:std::vector<std::tuple<std::string,std::string,std::uint32_t>>{
        {"pressure","Pa",1},{"velocity","m/s",3},{"velocity_magnitude","m/s",1},{"k","m^2/s^2",1},{"omega","1/s",1},{"viscosity_ratio","1",1}}) {
      s.cells.push_back({name,unit,components,std::vector<float>(std::size_t(s.mesh.n_cells)*components,1)});
      if(name=="pressure" || name.starts_with("velocity"))s.boundary.push_back({name,unit,components,std::vector<float>(std::size_t(s.mesh.n_faces()-s.mesh.n_internal_faces())*components,1)});
    }
    const auto nb=s.mesh.n_faces()-s.mesh.n_internal_faces();
    s.boundary.push_back({"wall_shear","Pa",3,std::vector<float>(std::size_t(nb)*3,1)});
    s.boundary.push_back({"wall_shear_magnitude","Pa",1,std::vector<float>(nb,1)});
    s.boundary.push_back({"y_plus","1",1,std::vector<float>(nb,1)});
    const std::filesystem::path out=argv[2];std::filesystem::create_directories(out);
    const auto prepared=Clock::now();kutta::view::write_snapshot(out/"view-input.bin",s);const auto saved=Clock::now();
    kutta::Vec3d lo=s.mesh.points.front(),hi=lo;
    for(auto p:s.mesh.points){lo={std::min(lo.x,p.x),std::min(lo.y,p.y),std::min(lo.z,p.z)};hi={std::max(hi.x,p.x),std::max(hi.y,p.y),std::max(hi.z,p.z)};}
    const kutta::view::Plane plane{{(lo.x+hi.x)*.5,(lo.y+hi.y)*.5,(lo.z+hi.z)*.5},{0,0,1}};
    kutta::view::export_view(out,s,std::span(&plane,1));const auto exported=Clock::now();
    const auto seconds=[](auto a,auto b){return std::chrono::duration<double>(b-a).count();};
    std::uintmax_t bytes=0;for(const auto& entry:std::filesystem::directory_iterator(out/"view"))bytes+=entry.file_size();
    std::ofstream marker(out/"summary.json");marker<<"{\"v\":1,\"manufactured_benchmark_fields\":true}\n";marker.close();
    std::cout<<std::setprecision(9)<<"{\"v\":1,\"mesh\":"<<kutta::json_string(argv[1])<<",\"cells\":"<<s.mesh.n_cells<<",\"faces\":"<<s.mesh.n_faces()<<",\"nodes\":"<<s.mesh.n_points()
      <<",\"load_seconds\":"<<seconds(start,loaded)<<",\"prepare_seconds\":"<<seconds(loaded,prepared)<<",\"snapshot_write_seconds\":"<<seconds(prepared,saved)
      <<",\"export_seconds\":"<<seconds(saved,exported)<<",\"view_bytes\":"<<bytes<<",\"manufactured_fields\":true}\n";return 0;
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
