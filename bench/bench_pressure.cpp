// Capture/replay the real FV pressure assembly, without changing solver gates.
#include <fstream>
#include <iomanip>
#include <iostream>
#include "../tests/linalg/instrumented_device.hpp"
#include "../validation/cases/tmr.hpp"
#include "kernels.hpp"

using namespace kutta;
namespace {
struct Assembly {
  int n=0,reference=-1,iteration=0;
  std::vector<Index> offsets,faces,adjacent,entry,rows,diag;
  std::vector<float> sign,coefficient,flux;
  CsrPattern pattern() const {
    KUTTA_REQUIRE(n>0 && offsets.size()==std::size_t(n)+1 && rows.size()==offsets.size() &&
      diag.size()==std::size_t(n) && offsets.front()==0 && rows.front()==0 &&
      offsets.back()>=0 && rows.back()>=n && faces.size()==std::size_t(offsets.back()) &&
      adjacent.size()==faces.size() && entry.size()==faces.size() && sign.size()==faces.size() &&
      coefficient.size()==flux.size() && reference>=-1 && reference<n,"invalid pressure replay dimensions");
    CsrPattern p{n,rows,std::vector<Index>(rows.back(),-1),diag};
    for(Index i=0;i<n;++i) {
      KUTTA_REQUIRE(rows[i]<=diag[i] && diag[i]<rows[i+1] && rows[i+1]<=rows.back() &&
        offsets[i]>=0 && offsets[i]<=offsets[i+1] && offsets[i+1]<=offsets.back(),"invalid replay row");
      p.col[diag[i]]=i;
      for(Index j=offsets[i];j<offsets[i+1];++j) {
        KUTTA_REQUIRE(faces[j]>=0 && std::size_t(faces[j])<coefficient.size(),"invalid replay face");
        if(adjacent[j]>=0) {
          KUTTA_REQUIRE(adjacent[j]<n && entry[j]>=rows[i] && entry[j]<rows[i+1],"invalid replay entry");
          p.col[entry[j]]=adjacent[j];
        }
      }
    }
    return p; // upload() validates the completed CSR pattern.
  }
};
template<class T> void write(std::ostream& out,const std::vector<T>& values) {
  out<<values.size();for(auto value:values) out<<' '<<value;out<<'\n';
}
template<class T> void read(std::istream& in,std::vector<T>& values) {
  std::size_t size=0;in>>size;
  KUTTA_REQUIRE(in && size<=std::size_t(std::numeric_limits<Index>::max()),"invalid pressure replay count");
  values.resize(size);for(auto& value:values) in>>value;
  KUTTA_REQUIRE(in,"truncated pressure replay");
}
void save(const std::string& path,const Assembly& a) {
  std::ofstream out(path);KUTTA_REQUIRE(out,"cannot write pressure assembly");
  out<<std::setprecision(17)<<a.n<<' '<<a.reference<<' '<<a.iteration<<'\n';
  write(out,a.offsets);write(out,a.faces);write(out,a.sign);write(out,a.adjacent);
  write(out,a.entry);write(out,a.rows);write(out,a.diag);write(out,a.coefficient);write(out,a.flux);
}
Assembly load(const std::string& path) {
  std::ifstream in(path);KUTTA_REQUIRE(in,"cannot read pressure assembly");Assembly a;
  in>>a.n>>a.reference>>a.iteration;
  read(in,a.offsets);read(in,a.faces);read(in,a.sign);read(in,a.adjacent);
  read(in,a.entry);read(in,a.rows);read(in,a.diag);read(in,a.coefficient);read(in,a.flux);
  return a;
}
void save_matrix(const std::string& path,const Assembly& a,const std::vector<float>& values) {
  const auto p=a.pattern();std::ofstream out(path);KUTTA_REQUIRE(out,"cannot write pressure matrix");
  out<<a.n<<' '<<a.iteration<<'\n'<<std::setprecision(17);
  for(Index i=0;i<a.n;++i) for(Index j=p.row_ptr[i];j<p.row_ptr[i+1];++j)
    out<<i<<' '<<p.col[j]<<' '<<values[j]<<'\n';
}
}
int main(int argc,char** argv) {
  try {
    KUTTA_REQUIRE(argc>=2,"usage: bench_pressure capture <tmr-root> <cpu|metal> <completed-step> <output-prefix> OR replay <assembly-file> <cpu|metal> <output-prefix>");
    const std::string mode=argv[1];
    KUTTA_REQUIRE((mode=="capture" && argc==6) || (mode=="replay" && argc==5),"invalid pressure diagnostic arguments");
    const std::string backend=argv[3],output=argv[argc-1];
    KUTTA_REQUIRE(backend=="cpu" || backend=="metal","unknown pressure backend");
    const auto kind=backend=="cpu"?Backend::Cpu:Backend::Metal;
    if(mode=="capture") {
      const int step=std::stoi(argv[4]);KUTTA_REQUIRE(step>=0,"negative capture step");
      test_linalg::InstrumentedDevice dev(kind);bool captured=false;int assembly=0;
      dev.after_launch=[&](KernelId id,std::size_t n,std::span<const KernelArg> args) {
        if(id!=KernelId::fv_pressure_assemble || assembly++!=step) return;
        Assembly a;a.n=Index(n);a.reference=args[11].get<Index>();a.iteration=step;
        a.offsets=dev.read<Index>(*args[0].buffer);a.faces=dev.read<Index>(*args[1].buffer);
        a.sign=dev.read<float>(*args[2].buffer);a.adjacent=dev.read<Index>(*args[3].buffer);
        a.entry=dev.read<Index>(*args[4].buffer);a.rows=dev.read<Index>(*args[5].buffer);
        a.diag=dev.read<Index>(*args[6].buffer);a.coefficient=dev.read<float>(*args[7].buffer);
        a.flux=dev.read<float>(*args[8].buffer);
        save(output+".assembly",a);save_matrix(output+".coo",a,dev.read<float>(*args[9].buffer));captured=true;
      };
      validation::tmr::Options options;options.max_iterations=step+1;
      const auto result=validation::tmr::run(dev,argv[2],validation::tmr::Case::Bump,0,options);
      KUTTA_REQUIRE(captured,"requested assembly was not reached: "+result.failure);
      std::cout<<"Captured pressure after "<<step<<" completed steps; later failure="<<result.failure<<'\n';
    } else {
      KUTTA_REQUIRE(mode=="replay","unknown pressure mode");
      const auto a=load(argv[2]);auto dev=Device::create(kind);auto matrix=upload(*dev,a.pattern());
      auto offsets=dev->upload<Index>(a.offsets),faces=dev->upload<Index>(a.faces),adjacent=dev->upload<Index>(a.adjacent),entry=dev->upload<Index>(a.entry);
      auto sign=dev->upload<float>(a.sign),coefficient=dev->upload<float>(a.coefficient),flux=dev->upload<float>(a.flux);
      auto rhs=dev->alloc<float>(a.n);
      kern::fv_pressure_assemble(*dev,a.n,offsets,faces,sign,adjacent,entry,matrix.row_ptr,matrix.diag,
        coefficient,flux,matrix.val,rhs,a.reference);
      const auto values=dev->download(matrix.val);save_matrix(output+".coo",a,values);
      CsrMatrixHost host{a.pattern(),{values.begin(),values.end()}};
      auto solver=make_solver(SolverKind::AmgPcg,*dev,host,{.structure=OperatorStructure::SymmetricMMatrix});
      solver->refresh(*dev,matrix);auto x=dev->alloc<float>(a.n);dev->fill(x,0);
      const auto stats=solver->solve(*dev,matrix,x,rhs,fv::SolverSettings{}.pressure_linear);
      std::cout<<std::setprecision(12)<<"pressure replay backend="<<backend<<" iterations="<<stats.iterations
        <<" converged="<<stats.converged<<" initial="<<stats.initial_residual<<" final="<<stats.final_residual
        <<" ratio="<<stats.final_residual/stats.initial_residual<<" max_correction="<<dev->max_abs(x,x.size())<<'\n';
    }
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
