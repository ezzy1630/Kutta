#include "export.hpp"
#include "protocol.hpp"
#include "kutta/post/aero.hpp"
#include "kutta/core/json.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <sstream>

namespace kutta::app {
namespace {
const Json& at(const Json& j,const std::string& key){return std::get<Json::Object>(j.value).at(key);}
std::string number(const Json& j,int precision=4) {
  const auto* n=std::get_if<double>(&j.value);if(!n)return "-";
  std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(precision)<<*n;return s.str();
}
void print(const std::filesystem::path& file,const Json& result) {
  std::cout<<std::left<<std::setw(24)<<"row"<<std::right<<std::setw(12)<<"drag N"<<std::setw(12)<<"downforce N"<<std::setw(10)<<"cd"<<std::setw(10)<<"cl"<<std::setw(10)<<"share"<<'\n';
  for(const auto& row:std::get<Json::Array>(at(result,"rows").value)) {
    const auto& coefficients=at(row,"coefficients");
    std::cout<<std::left<<std::setw(24)<<std::get<std::string>(at(row,"name").value)<<std::right<<std::setw(12)<<number(at(at(at(row,"force"),"total"),"drag"))
      <<std::setw(12)<<number(at(row,"downforce"))<<std::setw(10)<<number(at(coefficients,"cd"))<<std::setw(10)<<number(at(coefficients,"cl"))<<std::setw(10)<<number(at(row,"downforce_share"))<<'\n';
  }
  const auto& balance=at(result,"balance");
  std::cout<<"Downforce/drag: "<<number(at(result,"lift_to_drag"))<<"; balance: "<<std::get<std::string>(at(balance,"status").value);
  if(at(balance,"front_load")!=Json(nullptr))std::cout<<", front "<<number(at(balance,"front_load"))<<" N ("<<number(at(balance,"front_percent"))<<"%), rear "<<number(at(balance,"rear_load"))<<" N, CoP "<<number(at(balance,"cop_from_front_axle"))<<" m behind the front axle";
  std::cout<<'\n';
  if(const auto& p=std::get<std::string>(at(result,"provenance").value);!p.empty())std::cout<<p<<'\n';
  std::cout<<"Wrote "<<file.string()<<'\n';
}
}
int aero(int argc,char** argv) {
  bool json=false;std::filesystem::path directory,spec_file;bool usage=false;
  for(int i=2;i<argc;++i) {
    const std::string_view arg=argv[i];
    if(arg=="--json" && !json)json=true;
    else if(arg=="--aero" && spec_file.empty() && i+1<argc && *argv[i+1])spec_file=argv[++i];
    else if(!arg.starts_with("--") && directory.empty() && !arg.empty())directory=argv[i];
    else usage=true;
  }
  const auto fail=[&](const std::string& message,int code){if(json)error_event(message);std::cerr<<"kutta aero: "<<message<<'\n';return code;};
  if(usage || directory.empty())return fail("Usage: kutta aero DIR [--aero spec.json] [--json]",2);
  post::AeroSpec spec;
  try {if(!spec_file.empty())spec=post::read_aero_spec(spec_file);}
  catch(const std::exception& e){return fail(e.what(),2);}
  try {
    const auto loads_file=directory/"loads.json";
    if(!std::filesystem::is_regular_file(loads_file))throw Error("run has no loads.json; it was solved before aero loads were written, so solve it again");
    std::ifstream in(loads_file,std::ios::binary);if(!in)throw Error("cannot read "+loads_file.string());
    const std::string text{std::istreambuf_iterator<char>(in),{}};
    // A trusted run output; its size scales with the number of wall zones.
    auto result=post::aero(post::parse_loads(parse_json(text,std::size_t(256)*1024*1024)),spec);
    const auto target=directory/"aero.json";
    const auto stage=directory/(".aero-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count())+".json");
    try {
      std::ofstream out;out.exceptions(std::ios::failbit|std::ios::badbit);out.open(stage,std::ios::binary);
      out<<dump_json(result)<<'\n';out.close();
      std::filesystem::rename(stage,target);
    }catch(...){std::error_code ec;std::filesystem::remove(stage,ec);throw;}
    if(json) {
      auto& event=std::get<Json::Object>(result.value);event.emplace("event","aero");
      std::cout<<dump_json(result)<<std::endl;
    } else print(target,result);
    return 0;
  }
  catch(const post::InvalidAeroSpec& e){return fail(e.what(),2);}
  catch(const std::exception& e){return fail(e.what(),1);}
}
} // namespace kutta::app
