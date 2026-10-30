#include <doctest/doctest.h>
#include "kutta/core/sha256.hpp"
#include "../io/test_support.hpp"
#include <cstdlib>
#include <map>
#if KUTTA_HAS_HDF5
#include "../io/cff_fixture.hpp"
#include "kutta/case/hdf5.hpp"
namespace {
using namespace kutta;
using namespace kutta::fluent;
using namespace cff_fixture;
void text(hid_t f,const std::string& path,const std::string& value) {
  h5::Handle space(H5Screate(H5S_SCALAR),H5Sclose),type(H5Tcopy(H5T_C_S1),H5Tclose);
  h5::check(H5Tset_size(type,value.size()+1),path);
  h5::Handle d(H5Dcreate2(f,path.c_str(),type,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
  h5::check(H5Dwrite(d,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,value.c_str()),path);
}
std::string variables(bool unsupported) {
  std::string rp=R"((37 (
  (case-config ((rp-3d? . #t) (rp-double? . #t) (rp-seg? . #t) (rp-unsteady? . #f)
    (rp-lam? . #f) (rp-ke? . #f) (rp-kw? . #t) (rf-energy? . #f) (rp-des? . #f)))
  (kw-geko-on? #f) (kw-easm-on? #f) (kw-wj-bsl-earsm-on? #f) (kw-sst-on? #t) (kw-bsl-on? #f) (kw-std-on? #f) (gravity? #f)
  (kw-not-trans-sst-production-limiter? #t) (sst-pk-factor 10) (kw-low-re-mod? #f)
  (curvature-correction? #f) (curvature-correction-ccurv 1) (kw-sst-strain? #f)
  (kw-not-trans-sst-kato-launder? #f) (kw-vorticity-pk? #f) (kw-shear-mod? #t)
  (turb-compress-mod? #f) (omega/min-limit 1e-20) (k/min-limit 1e-14)
  (turb-viscosity/max-limit 1e5) (sst-limit-cross-diff-to-n-times-w-destruct 1e6)
  (sst-linearize-neg-cross-diff? #t) (turbulence/wall-omega-treatment "correlation")
  (pressure/scheme 12) (mom/scheme 1) (k/scheme 1) (omega/scheme 1) (flow/scheme 21)
  (recon/cell-lsq? #t) (recon/node-lsq? #f) (recon/lim/type 0) (recon/limit? #t)
  (recon/lim/dx-to-face? #t) (recon/lim/filter/on? #f) (recon/cell-lsf? #f)
  (recon/relax/relax? #f) (recon/relax/steady-urf 1)
  (pbs/flux-index 1) (pbs/flux-auto-select? #f)
  (wall-distance/method 2) (wall-distance/quad-correction? #t)
  (pseudo-transient-formulation? #t) (pseudo-auto-time-step? #f) (pseudo-time-step .01)
  (pseudo-time-method/segregated-pbns/dt-method 1) (pseudo-time-method/coupled-pbns/dt-method 1)
  (pseudo-transient-mode 0) (pseudo-auto-time-step-scale-factor 1)
  (pseudo/autotime-lscale-option 1) (pseudo-auto-time-step .01) (pressure-coupled/courant-number 5)
  (mom/pseudo-relax .7) (pressure/pseudo-relax 1) (density/pseudo-relax 1)
  (body-force/pseudo-relax 1) (turb-viscosity/pseudo-relax 1) (k/pseudo-relax .75) (omega/pseudo-relax .75)
  (materials ((water fluid (density (constant . 1)) (viscosity (constant . .01)))))
  (operating-pressure 101325) (reference-area 3) (reference-length 1) (reference-density 1)
  (reference-velocity 2) (reference-pressure 0) (reference-viscosity .01)
  (pressure-reference/x 0) (pressure-reference/y 0) (pressure-reference/z 0)
  (hyb-init? #f) (fmg-init/enable? #f) (x-velocity/default .1) (y-velocity/default 0)
  (z-velocity/default 0) (pressure/default 0) (k/default .01) (omega/default 10)
  (monitor/report-definitions (
    ((name "drag") (report-definition force-drag (zone-names "shell") (force-vector 1 0 0) (scaled? #t)))
    ((name "lift") (report-definition force-lift (zone-names "shell") (force-vector 0 1 0) (scaled? #f)))
    ((name "force") (report-definition force (zone-names "shell") (force-vector 0 0 1) (scaled? #f)))))
  )))";
  if(unsupported) {
    auto p=rp.find("(recon/cell-lsf? #f)");rp.replace(p,std::string("(recon/cell-lsf? #f)").size(),"(recon/cell-lsf? #t)");
    p=rp.find("(curvature-correction? #f)");rp.replace(p,std::string("(curvature-correction? #f)").size(),"(curvature-correction? #t)");
  }
  return rp;
}
void case_file(const std::filesystem::path& path,bool unsupported) {
  // Asymmetric dimensions and multiple cells exercise mesh/field ordering.
  synthetic(path,make_box_hex(2,2,1,{0,0,0},{2,1,1}),false);
  h5::Handle f(H5Fopen(path.string().c_str(),H5F_ACC_RDWR,H5P_DEFAULT),H5Fclose);
  group(f,"/settings");text(f,"/settings/Rampant Variables",variables(unsupported));
  text(f,"/settings/Thread Variables",R"(
    (39 (2 interior interior 1) ())
    (39 (20 fluid water 1) ((material . water) (mrf-motion? . #f) (porous? . #f) (laminar? . #f)
      (mrf-omega (constant . 0)) (mrf-relative-to-thread . -1) (motion-spec . 0)
      (mrf-grid-x-vel (constant . 0)) (mrf-grid-y-vel (constant . 0)) (mrf-grid-z-vel (constant . 0))
      (mrf-x-origin (constant . 0)) (mrf-y-origin (constant . 0)) (mrf-z-origin (constant . 0))
      (mrf-ai (constant . 0)) (mrf-aj (constant . 0)) (mrf-ak (constant . 1)) (porosity (constant . 1))))
    (39 (10 wall shell 1) ((motion-bc . 0) (relative? . #f) (specified-shear? . #f)
      (shear-bc . 0) (rough-bc . 0) (shear-x (constant . 0)) (shear-y (constant . 0)) (shear-z (constant . 0))
      (roughness-height (constant . 0)) (roughness-const (constant . .5)))))");
  text(f,"/settings/Version","24.1.0");
}
void data_file(const std::filesystem::path& path) {
  h5::Handle f(H5Fcreate(path.string().c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT),H5Fclose);
  group(f,"/settings");text(f,"/settings/Data Variables","(37 ((number-of-iterations 10)))");
  for(const auto& [name,value]:std::map<std::string,double>{{"SV_U",0},{"SV_V",0},{"SV_W",0},{"SV_P",0},{"SV_K",.2},{"SV_O",10}}) {
    const auto base="/results/1/phase-1/cells/"+name;group(f,base);
    const std::vector<double> values{value,value,value,value};const hsize_t count=4;
    h5::Handle space(H5Screate_simple(1,&count,nullptr),H5Sclose);
    const auto ds=base+"/1";
    h5::Handle d(H5Dcreate2(f,ds.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT),H5Dclose);
    h5::check(H5Dwrite(d,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,values.data()),ds);
    attr(f,ds,"minId",1);attr(f,ds,"maxId",4);
  }
}
std::string quote(const std::filesystem::path& path) {
#ifdef _WIN32
  return '"'+std::filesystem::path(path).make_preferred().string()+'"';
#else
  std::string out="'";for(char c:path.string()) out+=c=='\''?"'\\''":std::string(1,c);return out+"'";
#endif
}
}
TEST_CASE("Solve CLI runs generated CFF cases and data, rejects settings, and emits parseable outputs") {
  MeshTemp t;
  case_file(t.path("tiny.cas.h5"),false);case_file(t.path("unsupported.cas.h5"),true);
  data_file(t.path("tiny.dat.h5"));
  std::string command=quote(KUTTA_PYTHON_PATH)+" "+quote(std::filesystem::path(KUTTA_SOURCE_DIR)/"tests/app/check_solve.py")+" "+quote(KUTTA_CLI_PATH)+" "+quote(t.directory);
#ifdef _WIN32
  command='"'+command+'"';
#endif
  CHECK(std::system(command.c_str())==0);
}
#endif
TEST_CASE("SHA256 input fingerprints match standard vectors across padding boundaries") {
  MeshTemp t;
  for(const auto& [text,hash]:std::vector<std::pair<std::string,std::string>>{
    {"","e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"},
    {"abc","ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
    {"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq","248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"},
    {std::string(1000000,'a'),"cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"}}) {
    mesh_text(t.path("hash"),text);CHECK(kutta::sha256_file(t.path("hash"))==hash);
  }
}
