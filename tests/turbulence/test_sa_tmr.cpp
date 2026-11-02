#include "../backend_helpers.hpp"
#include "support.hpp"
#include "kutta/fv/turbulence/sa.hpp"

using namespace kutta;
using namespace sst_test;

TEST_CASE("TMR SA production retains raw modified strain when its ratio is safeguarded") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=box(*dev,1);State state(*dev,mesh);
    dev->fill(mesh.wall_distance,.1f);dev->fill(state.rho,1);dev->fill(state.mu,.001f);dev->fill(state.k,.003f);
    for(auto* g:{&state.gk.x,&state.gk.y,&state.gk.z}) dev->fill(*g,0);
    auto strain=filled(*dev,1,1),rotation=filled(*dev,1,1);
    SaOptions options;options.formulation=SaFormulation::TmrStandard;options.tmr_clip_production=false;
    SaModel model(*dev,mesh,options);model.update({state.rho,state.mu,state.k,state.gk,strain,rotation});
    const double q=double(float(.003)),nu=double(float(.001)),d=double(float(.1)),kp=options.kappa;
    const double chi=q/nu,cv=options.cv1,f1=std::pow(chi,3)/(std::pow(chi,3)+std::pow(cv,3));
    const double f2=1-chi/(1+chi*f1),raw_rate=1+q*f2/(kp*kp*d*d);
    const double ratio=std::min(q/(std::max(raw_rate,.3)*kp*kp*d*d),10.);
    const double g=ratio+options.cw2*(std::pow(ratio,6)-ratio),cw3=options.cw3;
    const double fw=g*std::pow((1+std::pow(cw3,6))/(std::pow(g,6)+std::pow(cw3,6)),1./6);
    const double ft2=1.2*std::exp(-.5*chi*chi),cw1=options.cb1/(kp*kp)+(1+options.cb2)/options.sigma;
    // TMR's transport equation uses raw modified strain in production.
    // Note 1(b) protects only the denominator used to form r and hence fw.
    const double expected=options.cb1*(1-ft2)*raw_rate*q-(cw1*fw-options.cb1*ft2/(kp*kp))*q*q/(d*d);
    const auto& fields=model.fields();const double source=dev->download(fields.source)[0],sink=dev->download(fields.sink)[0];
    CHECK(std::abs(double(dev->download(fields.modified_rate)[0])-raw_rate)<1e-6);
    CHECK(std::abs(source-sink*q-expected)<1e-7);
    CHECK(source>=0);CHECK(sink>=0);
  }
}
