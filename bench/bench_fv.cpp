#include <array>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <set>
#include <type_traits>

#include "kernels.hpp"

#include "../validation/cases/cavity.hpp"

#ifndef KUTTA_BENCH_COMMIT
#define KUTTA_BENCH_COMMIT "unknown"
#endif

namespace {
using Clock=std::chrono::steady_clock;
struct KernelMeasurement {std::size_t calls=0,bytes=0;double seconds=0;};

// This adapter owns native float allocations for every four-byte buffer type.
// Integer payloads are copied as bits, and launch arguments recover the native
// backend handles. Disabled timing preserves asynchronous Device ordering.
class MeasuredDevice final:public kutta::Device {
  struct Storage final:kutta::BufferHandle {
    kutta::Buffer<float> native;
    explicit Storage(kutta::Buffer<float> buffer):native(std::move(buffer)) {}
    std::size_t bytes() const override {return native.size()*sizeof(float);}
  };
  std::unique_ptr<kutta::Device> native_;
  static Storage& storage(kutta::BufferHandle& b) {return static_cast<Storage&>(b);}
  static const Storage& storage(const kutta::BufferHandle& b) {return static_cast<const Storage&>(b);}
  template<class Operation>
  decltype(auto) measured(const std::string& name,std::size_t bytes,Operation operation) {
    if(!timing) return operation();
    native_->synchronize();const auto start=Clock::now();
    const auto record=[&] {
      native_->synchronize();auto& r=measurements[name];++r.calls;r.bytes+=bytes;
      r.seconds+=std::chrono::duration<double>(Clock::now()-start).count();
    };
    if constexpr(std::is_void_v<std::invoke_result_t<Operation>>) {operation();record();}
    else {auto result=operation();record();return result;}
  }
 public:
  bool timing=false;
  std::map<std::string,KernelMeasurement> measurements;
  explicit MeasuredDevice(kutta::Backend backend):native_(kutta::Device::create(backend)) {}
  kutta::Backend backend() const override {return native_->backend();}
  std::string name() const override {return native_->name();}
  void synchronize() override {native_->synchronize();}
  double sum(const kutta::Buffer<float>& x,std::size_t n) override {
    return measured("reduce.sum",4*n,[&]{return native_->sum(storage(*x.handle()).native,n);});
  }
  double dot(const kutta::Buffer<float>& x,const kutta::Buffer<float>& y,std::size_t n) override {
    return measured("reduce.dot",8*n,[&]{return native_->dot(storage(*x.handle()).native,storage(*y.handle()).native,n);});
  }
  double max_abs(const kutta::Buffer<float>& x,std::size_t n) override {
    return measured("reduce.max",4*n,[&]{return native_->max_abs(storage(*x.handle()).native,n);});
  }
  void launch(kutta::KernelId id,std::size_t n,std::span<const kutta::KernelArg> args) override {
    KUTTA_REQUIRE(args.size()<=32,"FV profiler argument capacity exceeded");
    std::array<kutta::KernelArg,32> translated;
    std::size_t bytes=0;
    for(std::size_t i=0;i<args.size();++i) {
      translated[i]=args[i];
      if(args[i].kind==kutta::KernelArg::Kind::Buffer) {
        if(timing) bytes+=args[i].buffer->bytes();
        translated[i].buffer=storage(*args[i].buffer).native.handle();
      }
    }
    measured(kutta::kernel_name(id),bytes,[&]{native_->launch(id,n,{translated.data(),args.size()});});
  }
 protected:
  void reduce_into(const kutta::Buffer<float>& x,const kutta::Buffer<float>& y,
      const kutta::Buffer<float>& u,const kutta::Buffer<float>& v,std::size_t n,
      kutta::Buffer<float>& dst,std::size_t offset,Reduction mode) override {
    check_reduction(x,y,u,v,n,dst,offset,mode);
    const auto& a=storage(*x.handle()).native;const auto& b=storage(*y.handle()).native;
    auto& out=storage(*dst.handle()).native;
    if(mode==Reduction::Sum) measured("reduce_into.sum",4*n,[&]{native_->sum_into(a,n,out,offset);});
    else if(mode==Reduction::Dot) measured("reduce_into.dot",8*n,[&]{native_->dot_into(a,b,n,out,offset);});
    else measured("reduce_into.dot2",16*n,[&]{native_->dot2_into(a,b,storage(*u.handle()).native,
      storage(*v.handle()).native,n,out,offset);});
  }
  std::unique_ptr<kutta::BufferHandle> alloc_bytes(std::size_t bytes) override {
    KUTTA_REQUIRE(bytes%4==0,"FV profiler requires four-byte elements");
    return std::make_unique<Storage>(native_->alloc<float>(bytes/4));
  }
  void upload_bytes(kutta::BufferHandle& dst,const void* source,std::size_t bytes) override {
    std::vector<float> bits(bytes/4);if(bytes) std::memcpy(bits.data(),source,bytes);
    native_->upload<float>(storage(dst).native,bits);
  }
  void download_bytes(const kutta::BufferHandle& source,void* dst,std::size_t bytes) override {
    auto bits=native_->download(storage(source).native);if(bytes) std::memcpy(dst,bits.data(),bytes);
  }
  void copy_bytes(kutta::BufferHandle& dst,const kutta::BufferHandle& src,std::size_t bytes) override {
    measured("runtime.copy",2*bytes,[&]{native_->copy(storage(dst).native,storage(src).native);});
  }
  void fill_bytes(kutta::BufferHandle& dst,const void* value,std::size_t width,std::size_t count) override {
    KUTTA_REQUIRE(width==4,"FV profiler fill width mismatch");
    float bits;std::memcpy(&bits,value,4);
    measured("runtime.fill",4*count,[&]{native_->fill(storage(dst).native,bits);});
  }
};

void performance(kutta::Backend backend,int size,int repeats) {
  MeasuredDevice dev(backend);
  auto mesh=kutta::make_box_hex(size,size,size,{0,0,0},{1,1,1});
  std::vector<kutta::fv::PatchCondition> boundaries;
  for(const auto& p:mesh.patches) {
    if(p.name=="zmin" || p.name=="zmax") boundaries.push_back({p.name,kutta::fv::Symmetry{}});
    else boundaries.push_back({p.name,kutta::fv::Wall{{p.name=="ymax"?1.:0.,0,0}}});
  }
  kutta::fv::SolverSettings settings;settings.momentum_relaxation=1;
  kutta::fv::SteadySolver solver(dev,std::move(mesh),boundaries,settings);
  solver.iterate(5);dev.synchronize();
  KUTTA_REQUIRE(std::isfinite(solver.report().maximum_velocity) && solver.report().maximum_velocity>0,
    "FV benchmark warmup did not produce moving flow");
  dev.measurements.clear();dev.timing=true;solver.iterate(1);dev.synchronize();dev.timing=false;
  const auto cells=solver.mesh().n_cells;
  std::size_t fv_bytes=0,total_bytes=0;
  std::cout<<std::setprecision(9);
  for(const auto& [name,r]:dev.measurements) {
    total_bytes+=r.bytes;if(name.starts_with("fv_")) fv_bytes+=r.bytes;
    std::cout<<"FV_KERNEL backend="<<kutta::backend_name(backend)<<" device=\""<<dev.name()<<"\" cells="<<cells
      <<" kernel="<<name<<" calls_per_outer="<<r.calls<<" us_per_call="<<1e6*r.seconds/r.calls
      <<" logical_bytes_per_cell_per_outer="<<double(r.bytes)/cells
      <<" logical_GB_s="<<double(r.bytes)/r.seconds/1e9<<std::endl;
  }
  const auto start=Clock::now();solver.iterate(repeats);dev.synchronize();
  const double seconds=std::chrono::duration<double>(Clock::now()-start).count();
  std::cout<<"FV_PERFORMANCE backend="<<kutta::backend_name(backend)<<" device=\""<<dev.name()<<"\" source="<<KUTTA_BENCH_COMMIT
    <<" cells="<<cells<<" measured_outer_iterations="<<repeats<<" outer_iterations_per_s="<<repeats/seconds
    <<" fv_logical_bytes_per_cell="<<double(fv_bytes)/cells
    <<" all_logical_bytes_per_cell="<<double(total_bytes)/cells
    <<" fv_ideal_ms_at_240GB_s="<<1e3*fv_bytes/240e9
    <<" last_pressure_iterations="<<solver.report().pressure.iterations
    <<" last_pressure_converged="<<solver.report().pressure.converged<<std::endl;
  if(size<=8) {
    auto native=kutta::Device::create(backend);
    kutta::fv::SteadySolver reference(*native,kutta::make_box_hex(size,size,size,{0,0,0},{1,1,1}),boundaries,settings);
    reference.iterate(6+repeats);
    KUTTA_REQUIRE(native->download(reference.velocity().x)==dev.download(solver.velocity().x) &&
      native->download(reference.velocity().y)==dev.download(solver.velocity().y) &&
      native->download(reference.velocity().z)==dev.download(solver.velocity().z) &&
      native->download(reference.pressure())==dev.download(solver.pressure()) &&
      native->download(reference.mass_flux())==dev.download(solver.mass_flux()),
      "FV profiling adapter changed the solution");
    std::cout<<"FV_PROFILER_ORACLE backend="<<kutta::backend_name(backend)<<" bitwise_equal=1"<<std::endl;
  }
}

void save_profiles(const kutta::validation::CavityResult& result,const std::string& stem) {
  std::ofstream output(stem+".csv");KUTTA_REQUIRE(output.good(),"Cannot write cavity profile: "+stem);
  output<<"component,coordinate,velocity\n"<<std::setprecision(12);
  for(auto p:result.u) output<<"u,"<<p.first<<","<<p.second<<"\n";
  for(auto p:result.v) output<<"v,"<<p.first<<","<<p.second<<"\n";
}

double profile_difference(const kutta::validation::Profile& a,const kutta::validation::Profile& b) {
  double error=0;
  for(auto p:a) error=std::max(error,std::abs(p.second-kutta::validation::interpolate(b,p.first)));
  for(auto p:b) error=std::max(error,std::abs(p.second-kutta::validation::interpolate(a,p.first)));
  return error;
}
}

int main(int argc,char** argv) {
  try {
    int reynolds=100,nx=129,ny=129,iterations=6000,size=100,repeats=3;
    bool run_cavity=false;
    std::string backend_filter;std::filesystem::path fluent;
    for(int i=1;i<argc;++i) {
      std::string arg=argv[i];
      KUTTA_REQUIRE(i+1<argc,"FV benchmark option needs a value: "+arg);
      const std::string value=argv[++i];
      if(arg=="--cavity") {reynolds=std::stoi(value);run_cavity=true;}
      else if(arg=="--nx") nx=std::stoi(value);
      else if(arg=="--ny") ny=std::stoi(value);
      else if(arg=="--iterations") iterations=std::stoi(value);
      else if(arg=="--backend") backend_filter=value;
      else if(arg=="--fluent") {fluent=value;run_cavity=true;}
      else if(arg=="--size") size=std::stoi(value);
      else if(arg=="--repeats") repeats=std::stoi(value);
      else throw kutta::Error("Unknown FV benchmark option: "+arg);
    }
    KUTTA_REQUIRE(reynolds==100 || reynolds==400 || reynolds==1000,"Cavity Re must be 100, 400 or 1000");
    KUTTA_REQUIRE(nx>1 && ny>1 && iterations>=100,"Invalid cavity dimensions/iterations");
    KUTTA_REQUIRE(size>1 && repeats>0 && repeats<10,"Invalid FV performance dimensions/repeats");
    const auto data=std::filesystem::path(__FILE__).parent_path().parent_path()/"validation/data";
    std::cout<<"FV source="<<KUTTA_BENCH_COMMIT<<std::endl;
    bool success=true,ran=false;
    for(auto backend:kutta::Device::available()) {
      if(!backend_filter.empty() && backend_filter!=kutta::backend_name(backend)) continue;
      ran=true;
      if(!run_cavity) {performance(backend,size,repeats);continue;}
      auto dev=kutta::Device::create(backend);
      std::cout<<"FV device="<<dev->name()<<std::endl;
      auto mesh=fluent.empty()?kutta::make_box_hex(nx,ny,1,{0,0,0},{1,1,1}):kutta::read_fluent_msh(fluent);
      std::vector<double> grid_x,grid_y;
      if(!fluent.empty()) {
        std::set<long long> x,y;
        for(auto p:kutta::compute_geometry(mesh).cell_centroid) {x.insert(std::llround(p.x*1e10));y.insert(std::llround(p.y*1e10));}
        nx=int(x.size());ny=int(y.size());
        KUTTA_REQUIRE(nx*ny==mesh.n_cells,"Fluent comparison expects a rectangular extruded grid");
        std::map<long long,double> nodes_x,nodes_y;
        for(auto p:mesh.points) {nodes_x.emplace(std::llround(p.x*1e10),p.x);nodes_y.emplace(std::llround(p.y*1e10),p.y);}
        for(auto [key,value]:nodes_x) {static_cast<void>(key);grid_x.push_back(value);}
        for(auto [key,value]:nodes_y) {static_cast<void>(key);grid_y.push_back(value);}
        KUTTA_REQUIRE(grid_x.size()==std::size_t(nx+1) && grid_y.size()==std::size_t(ny+1),"Fluent grid lines are not rectangular");
      }
      auto result=kutta::validation::cavity(*dev,std::move(mesh),reynolds,iterations);
      kutta::validation::compare_cavity(result,data,reynolds,backend,fluent.empty()?std::to_string(nx)+"x"+std::to_string(ny):fluent.string());
      success=success && result.converged;
      const std::string stem="build/fv-cavity-"+std::string(kutta::backend_name(backend))+"-Re"+std::to_string(reynolds)+
        "-"+(fluent.empty()?std::to_string(nx)+"x"+std::to_string(ny):"fluent");
      save_profiles(result,stem);
      if(!fluent.empty()) {
        auto uniform=kutta::validation::cavity(*dev,kutta::make_box_hex(nx,ny,1,{0,0,0},{1,1,1}),reynolds,iterations);
        kutta::validation::compare_cavity(uniform,data,reynolds,backend,"uniform-equal-resolution");
        std::cout<<"FV_V4_MATCH backend="<<kutta::backend_name(backend)<<" Re="<<reynolds<<" resolution="<<nx<<"x"<<ny
          <<" u_max="<<profile_difference(result.u,uniform.u)<<" v_max="<<profile_difference(result.v,uniform.v)<<std::endl;
        success=success && uniform.converged;save_profiles(uniform,stem+"-uniform");
        auto matched_mesh=kutta::make_box_hex(nx,ny,1,{0,0,0},{1,1,1});
        for(auto& p:matched_mesh.points) {p.x=grid_x[std::llround(p.x*nx)];p.y=grid_y[std::llround(p.y*ny)];}
        auto matched=kutta::validation::cavity(*dev,std::move(matched_mesh),reynolds,iterations);
        std::cout<<"FV_V4_MATCHED_GRID backend="<<kutta::backend_name(backend)<<" Re="<<reynolds
          <<" u_max="<<profile_difference(result.u,matched.u)<<" v_max="<<profile_difference(result.v,matched.v)
          <<" converged="<<matched.converged<<std::endl;
        success=success && matched.converged;save_profiles(matched,stem+"-matched-grid");
      }
    }
    KUTTA_REQUIRE(ran,"Requested FV backend is unavailable");
    return success?0:1;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
