#include "protocol.hpp"

#include "kutta/core/json.hpp"
#include "kutta/io/view.hpp"
#include "kutta/core/parse.hpp"
#include <cmath>

#include <iostream>

namespace kutta::app {
int mesh_section(int argc,char** argv) {
  try {
    KUTTA_REQUIRE(argc==5 && std::string_view(argv[3])=="--sample",
                  "Usage: kutta mesh-section MESH --sample x0,y0,z0,x1,y1,z1");
    std::string_view sample=argv[4];
    double values[6];
    for(int i=0;i<6;++i) {
      const auto comma=sample.find(',');
      KUTTA_REQUIRE((i==5)==(comma==std::string_view::npos),"section requires six comma-separated coordinates");
      const auto value=parse_double(sample.substr(0,comma));
      KUTTA_REQUIRE(value && std::isfinite(*value),"section coordinates must be finite numbers");
      values[i]=*value;
      if(i<5)sample.remove_prefix(comma+1);
    }
    const auto mesh=read_mesh(argv[2]);
    std::cout<<dump_json(view::mesh_gap_section(mesh,{values[0],values[1],values[2]},
                                                  {values[3],values[4],values[5]}))<<'\n';
    return 0;
  } catch(const std::exception& e) {
    error_event(e.what());
    std::cerr<<"kutta mesh-section: "<<e.what()<<'\n';
    return 1;
  }
}
}
