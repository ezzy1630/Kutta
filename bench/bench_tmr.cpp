#include <charconv>
#include <iostream>
#include "../validation/cases/tmr.hpp"

int main(int argc,char** argv) {
  using namespace kutta;using namespace validation::tmr;
  try {
    KUTTA_REQUIRE(argc>=5,"usage: bench_tmr <data-root> <plate|bump> <cpu|metal|cuda> <level 0..2> [iterations] [output-dir] [default|2003|1994|vorticity|no-wall-production|no-cross-clip|symmetry-top]");
    auto integer=[](const char* s){int value=0;const auto [end,error]=std::from_chars(s,s+std::char_traits<char>::length(s),value);KUTTA_REQUIRE(error==std::errc{} && *end=='\0' && value>0,"invalid positive integer");return value;};
    const std::string_view name=argv[2],backend=argv[3];
    KUTTA_REQUIRE(name=="plate" || name=="bump","unknown TMR case");KUTTA_REQUIRE(backend=="cpu" || backend=="metal" || backend=="cuda","unknown backend");
    const std::string_view level_text=argv[4];KUTTA_REQUIRE(level_text=="0" || level_text=="1" || level_text=="2","level must be 0..2");const int level=argv[4][0]-'0';
    Options o;if(argc>5) o.max_iterations=integer(argv[5]);if(argc>6) o.output=argv[6];
    const std::string_view mode=argc>7?argv[7]:"default";
    if(mode=="2003") o.model.gamma_set=SstGammaSet::Menter2003;
    else if(mode=="1994") {o.model.viscosity_limiter=SstViscosityLimiter::Vorticity;o.model.production_clip=20;o.model.limit_omega_production=false;o.model.cdkw_floor=1e-20f;}
    else if(mode=="vorticity") o.model.vorticity_production=true;
    else if(mode=="no-wall-production") o.model.wall_production=false;
    else if(mode=="no-cross-clip") o.model.limit_cross_diffusion=false;
    else if(mode=="symmetry-top") o.symmetry_top=true;
    else KUTTA_REQUIRE(mode=="default","unknown model/BC diagnostic");
    auto dev=Device::create(backend=="cpu"?Backend::Cpu:(backend=="metal"?Backend::Metal:Backend::Cuda));
    const auto r=run(*dev,argv[1],name=="plate"?Case::FlatPlate:Case::Bump,level,o);
    std::cout<<std::setprecision(12)<<"RESULT source="<<KUTTA_BENCH_COMMIT<<" device="<<dev->name()<<" mode="<<mode<<" cells="<<r.cells<<" iterations="<<r.iterations<<" setup_s="<<r.setup_seconds<<" solve_s="<<r.solve_seconds<<" converged="<<r.converged<<" defect="<<defect(r.residual)<<" change="<<r.load_change<<" continuity_raw="<<r.residual.continuity.raw<<" continuity_scale="<<r.residual.continuity.fluent_normalization<<" corrected_mass="<<r.residual.corrected_mass_imbalance<<" ru="<<r.residual.u.scaled<<" rv="<<r.residual.v.scaled<<" rk="<<r.residual.k.scaled<<" rw="<<r.residual.omega.scaled<<" drag="<<r.values.drag<<" ymaxplus="<<r.values.max_yplus;
    for(auto q:r.values.cf) std::cout<<" cf="<<q;for(auto q:r.values.cp) std::cout<<" cp="<<q;std::cout<<'\n';
    if(!r.failure.empty()) std::cout<<"FAILURE completed_iterations="<<r.iterations<<" reason="<<r.failure<<'\n';
    return r.converged?0:2;
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
