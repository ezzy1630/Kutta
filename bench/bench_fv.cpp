#include <iostream>

#include "../validation/cases/cavity.hpp"

#ifndef KUTTA_BENCH_COMMIT
#define KUTTA_BENCH_COMMIT "unknown"
#endif

int main(int argc,char** argv) {
  try {
    int reynolds=100,nx=129,ny=129,iterations=6000;
    std::string backend_filter;std::filesystem::path fluent;
    for(int i=1;i<argc;++i) {
      std::string arg=argv[i];
      KUTTA_REQUIRE(i+1<argc,"FV benchmark option needs a value: "+arg);
      const std::string value=argv[++i];
      if(arg=="--cavity") reynolds=std::stoi(value);
      else if(arg=="--nx") nx=std::stoi(value);
      else if(arg=="--ny") ny=std::stoi(value);
      else if(arg=="--iterations") iterations=std::stoi(value);
      else if(arg=="--backend") backend_filter=value;
      else if(arg=="--fluent") fluent=value;
      else throw kutta::Error("Unknown FV benchmark option: "+arg);
    }
    KUTTA_REQUIRE(reynolds==100 || reynolds==400 || reynolds==1000,"Cavity Re must be 100, 400 or 1000");
    KUTTA_REQUIRE(nx>1 && ny>1 && iterations>=100,"Invalid cavity dimensions/iterations");
    const auto data=std::filesystem::path(__FILE__).parent_path().parent_path()/"validation/data";
    std::cout<<"FV source="<<KUTTA_BENCH_COMMIT<<std::endl;
    bool success=true;
    for(auto backend:kutta::Device::available()) {
      if(!backend_filter.empty() && backend_filter!=kutta::backend_name(backend)) continue;
      auto dev=kutta::Device::create(backend);
      std::cout<<"FV device="<<dev->name()<<std::endl;
      auto mesh=fluent.empty()?kutta::make_box_hex(nx,ny,1,{0,0,0},{1,1,1}):kutta::read_fluent_msh(fluent);
      auto result=kutta::validation::cavity(*dev,std::move(mesh),reynolds,iterations);
      kutta::validation::compare_cavity(result,data,reynolds,backend,fluent.empty()?std::to_string(nx)+"x"+std::to_string(ny):fluent.string());
      success=success && result.converged;
    }
    return success?0:1;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
