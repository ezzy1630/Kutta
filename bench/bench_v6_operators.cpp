// V6 operator probe. Seed the mixing elbow with Fluent's saved fields and,
// before any step, report the saved-state transport residuals and compare
// Kutta's mass flux with Fluent's saved SV_FLUX by face skewness (flux rows
// from tools/forensics/dump_interior_flux.py; pass - to skip). With
// "converge", solve with V6's stopping rule and attribute the field error
// versus Fluent to wall-adjacent cells and to interior cells by worst skew.
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/setup/fluent.hpp"

using namespace kutta;
using namespace kutta::detail;
using namespace kutta::fv;
namespace {
struct Bin {
  std::vector<double> relative;double square=0,reference=0;
  void add(double ours,double saved) {relative.push_back(std::abs(ours-saved)/std::max(std::abs(saved),1e-30));square+=(ours-saved)*(ours-saved);reference+=saved*saved;}
  double quantile(double q) {if(relative.empty()) return 0;std::sort(relative.begin(),relative.end());return relative[std::size_t(q*(relative.size()-1))];}
};
}
int main(int argc,char** argv) {
  KUTTA_REQUIRE(argc>=3,"Usage: bench_v6_operators fluent_dir interior_flux.txt|- [cpu|metal] [converge]");
  const std::filesystem::path root=argv[1];
  const auto c=fluent::read_case(root/"mixing_elbow.cas.h5");const auto mesh=read_fluent_cff(root/"mixing_elbow.cas.h5");
  const auto d=fluent::read_data(root/"mixing_elbow.dat.h5",c);
  const auto saved=setup::fluent_fields(c,d,mesh);
  const auto setup=setup::fluent_setup(c,mesh,setup::parse_approximations("wfgc"));
  std::map<std::pair<Index,Index>,double> fluent_flux;std::map<std::pair<Index,Index>,int> multiplicity;
  const bool compare_flux=std::string_view(argv[2])!="-";
  if(compare_flux) {
    std::ifstream in(argv[2]);std::size_t n=0;Index cells=0;in>>n>>cells;KUTTA_REQUIRE(in && cells==mesh.n_cells,"flux dump cell count");
    for(std::size_t j=0;j<n;++j) {Index a,b;double f;in>>a>>b>>f;const auto key=std::minmax(a,b);fluent_flux[key]=a<b?f:-f;++multiplicity[key];}
  }
  const auto g=compute_geometry(mesh);
  std::vector<int> faces_per_cell(mesh.n_cells,0);
  for(Index f=0;f<mesh.n_faces();++f) {++faces_per_cell[mesh.owner[f]];if(f<mesh.n_internal_faces()) ++faces_per_cell[mesh.neighbour[f]];}
  std::map<int,int> histogram;for(auto k:faces_per_cell) ++histogram[k];
  std::cout<<"V6_MESH cells="<<mesh.n_cells<<" faces="<<mesh.n_faces()<<" internal="<<mesh.n_internal_faces()<<" faces_per_cell";
  for(auto [k,count]:histogram) std::cout<<' '<<k<<':'<<count;
  std::cout<<'\n';
  const char* backend=argc>3?argv[3]:"metal";
  for(auto b:Device::available()) if(std::string_view(backend_name(b))==backend) {
    auto dev=Device::create(b);
    auto settings=setup.settings;settings.monitor_interval=50;
    SteadySolver solver(*dev,mesh,setup.boundaries,settings);
    setup::initialize(solver,setup::Initialization{saved});
    const auto r=solver.monitor();
    const auto& sm=solver.mesh();const auto& sg=solver.geometry();const auto& perm=solver.permutation().cells_new_to_old;
    const auto mass=dev->download(solver.mass_flux());
    // Skew is the tangential offset of both centroids from the face normal
    // line, relative to the normal centroid distance.
    Bin bins[4];const char* names[]{"orthogonal","skew<0.05","skew<0.2","skew>=0.2"};Bin all;double divergence=0;
    std::vector<double> div(sm.n_cells,0);
    for(Index f=0;f<sm.n_faces();++f) {div[sm.owner[f]]+=mass[f];if(f<sm.n_internal_faces()) div[sm.neighbour[f]]-=mass[f];}
    for(auto v:div) divergence+=std::abs(v);
    for(Index f=0;compare_flux && f<sm.n_internal_faces();++f) {
      const Index o=perm[sm.owner[f]],n=perm[sm.neighbour[f]];const auto key=std::minmax(o,n);
      if(multiplicity[key]!=1) continue;
      const double reference=o<n?fluent_flux.at(key):-fluent_flux.at(key);
      const auto s=sg.face_area[f];const auto unit=s/norm(s);
      const auto rp=sg.face_centroid[f]-sg.cell_centroid[sm.owner[f]],rn=sg.cell_centroid[sm.neighbour[f]]-sg.face_centroid[f];
      const double dp=dot(rp,unit),dn=dot(rn,unit);
      const double skew=(norm(rp-unit*dp)+norm(rn-unit*dn))/(dp+dn);
      const int bin=skew<1e-6?0:(skew<.05?1:(skew<.2?2:3));
      bins[bin].add(mass[f],reference);all.add(mass[f],reference);
    }
    std::cout<<std::setprecision(7)<<"V6_PROBE backend="<<backend<<" u="<<r.u.scaled<<" v="<<r.v.scaled<<" w="<<r.w.scaled<<" k="<<r.k.scaled<<" omega="<<r.omega.scaled
      <<" continuity_raw="<<r.continuity.raw<<" sum_abs_div="<<divergence<<'\n';
    for(int k=0;compare_flux && k<5;++k) {
      auto& bin=k<4?bins[k]:all;
      std::cout<<"V6_FLUX bin="<<(k<4?names[k]:"all")<<" faces="<<bin.relative.size()
        <<" relative_L2="<<(bin.reference>0?std::sqrt(bin.square/bin.reference):0)
        <<" median="<<bin.quantile(.5)<<" p90="<<bin.quantile(.9)<<'\n';
    }
    if(argc<=4) continue;
    // Converge with V6's own stopping rule, then attribute the squared field
    // error to cell classes: wall-adjacent or not, and worst face skew.
    const double continuity_target=d.residuals.at("continuity").rows.back()[0];
    for(int it=0;it<4000;it+=50) {
      solver.iterate(50);const auto q=solver.monitor();
      if(std::max({q.u.scaled,q.v.scaled,q.w.scaled,q.k.scaled,q.omega.scaled})<1e-6 && q.continuity.raw<continuity_target) break;
    }
    std::vector<double> worst(sm.n_cells,0);std::vector<int> wall(sm.n_cells,0);
    for(Index f=0;f<sm.n_faces();++f) {
      const auto s=sg.face_area[f];const auto unit=s/norm(s);
      const auto rp=sg.face_centroid[f]-sg.cell_centroid[sm.owner[f]];const double dp=dot(rp,unit);
      double skew=norm(rp-unit*dp)/dp;
      if(f<sm.n_internal_faces()) {
        const auto rn=sg.cell_centroid[sm.neighbour[f]]-sg.face_centroid[f];const double dn=dot(rn,unit);
        skew=(norm(rp-unit*dp)+norm(rn-unit*dn))/(dp+dn);worst[sm.neighbour[f]]=std::max(worst[sm.neighbour[f]],skew);
      } else {
        const auto& patch=*std::find_if(sm.patches.begin(),sm.patches.end(),[&](const auto& p){return f>=p.start && f<p.start+p.size;});
        if(patch.kind==PatchKind::Wall) wall[sm.owner[f]]=1;
      }
      worst[sm.owner[f]]=std::max(worst[sm.owner[f]],skew);
    }
    const auto u=dev->download(solver.velocity().x),v=dev->download(solver.velocity().y),w=dev->download(solver.velocity().z);
    const auto p=dev->download(solver.pressure()),k=dev->download(solver.kinetic_energy()),o=dev->download(solver.specific_dissipation());
    struct Sum {double error=0,reference=0;int cells=0;};
    const char* classes[]{"wall","interior,skew<0.1","interior,skew<0.3","interior,skew>=0.3"};
    for(const char* field:{"velocity","pressure","k","omega"}) {
      Sum sums[4];double total=0,reference=0;
      for(Index cell=0;cell<sm.n_cells;++cell) {
        const auto old=perm[cell];double e=0,r=0;
        if(field[0]=='v') {const auto s=saved.velocity[old];e=std::pow(u[cell]-s.x,2)+std::pow(v[cell]-s.y,2)+std::pow(w[cell]-s.z,2);r=dot(s,s);}
        if(field[0]=='p') {e=std::pow(p[cell]-saved.pressure[old],2);r=saved.pressure[old]*saved.pressure[old];}
        if(field[0]=='k') {e=std::pow(k[cell]-saved.k[old],2);r=saved.k[old]*saved.k[old];}
        if(field[0]=='o') {e=std::pow(o[cell]-saved.omega[old],2);r=saved.omega[old]*saved.omega[old];}
        const int c=wall[cell]?0:(worst[cell]<.1?1:(worst[cell]<.3?2:3));
        sums[c].error+=e;sums[c].reference+=r;++sums[c].cells;total+=e;reference+=r;
      }
      std::cout<<"V6_LOCAL field="<<field<<" relative_L2="<<std::sqrt(total/reference);
      for(int c=0;c<4;++c) std::cout<<" ["<<classes[c]<<" cells="<<sums[c].cells<<" share="<<sums[c].error/total
        <<" class_L2="<<std::sqrt(sums[c].error/std::max(sums[c].reference,1e-300))<<']';
      std::cout<<'\n';
    }
  }
}
