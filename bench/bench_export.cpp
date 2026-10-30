// Host-side view-export benchmark. Fields are manufactured, never a CFD solve.
#include "kutta/io/view.hpp"
#include "kutta/core/json.hpp"
#include "kutta/mesh/detail.hpp"
#include <chrono>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc,char** argv) {
  using Clock=std::chrono::steady_clock;
  if(argc!=3){std::cerr<<"Usage: bench_export MESH OUTPUT_DIR\n";return 2;}
  try {
    const auto start=Clock::now();kutta::view::Snapshot s;s.mesh=kutta::read_mesh(argv[1]);
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
