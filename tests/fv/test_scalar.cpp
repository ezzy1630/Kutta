#include <cmath>
#include <iostream>
#include <numbers>

#include "../backend_helpers.hpp"
#include "kutta/fv/transport.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
using namespace kutta::detail;

namespace {
double manufactured(Device& dev,int n,bool tet,fv::ConvectionScheme scheme,double speed) {
  auto mesh=tet?make_box_tet(n,n,n,{0,0,0},{1,1,1}):make_box_hex(n,n,n,{0,0,0},{1,1,1});
  auto g=compute_geometry(mesh); auto m=upload_mesh(dev,mesh,g); auto s=fv::make_stencil(dev,mesh,g);
  fv::TransportEquation equation(dev,m,s,mesh,g);
  const double pi=std::numbers::pi;
  const auto exact=[pi](Vec3d p){return std::sin(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z);};
  std::vector<float> source,boundary,flux;
  for(auto p:g.cell_centroid) source.push_back(float(3*pi*pi*exact(p)+speed*pi*std::cos(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z)));
  for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(exact(g.face_centroid[f])));
  for(auto a:g.face_area) flux.push_back(float(speed*a.x));
  auto q=dev.alloc<float>(m.n_cells);dev.fill(q,0);
  auto mu=dev.alloc<float>(m.n_faces);dev.fill(mu,1);
  auto kinds=dev.alloc<Index>(m.n_boundary_faces());dev.fill(kinds,0);
  auto bc=dev.upload<float>(boundary),src=dev.upload<float>(source),mass=dev.upload<float>(flux);
  fv::TransportSettings settings;settings.convection=scheme;settings.limit=true;
  auto old=dev.alloc<float>(m.n_cells);
  int iterations=0;
  for(;iterations<200;++iterations) {
    dev.copy(old,q); equation.assemble(q,bc,kinds,mass,mu,src,settings);
    auto stats=equation.solve(q,{1e-5,1e-10,2000});
    if(!stats.converged && stats.final_residual>2e-6) {
      INFO(stats.final_residual);CHECK(stats.converged);
    }
    // Only test monitors download; the operator/outer step never does.
    auto got=dev.download(q),before=dev.download(old);
    double update=0;for(Index c=0;c<m.n_cells;++c) update=std::max(update,std::abs(double(got[c]-before[c])));
    if(update<2e-7) break;
  }
  CHECK(iterations<200);
  auto got=dev.download(q); double error=0,vol=0;
  for(Index c=0;c<m.n_cells;++c){error+=g.cell_volume[c]*std::pow(got[c]-exact(g.cell_centroid[c]),2);vol+=g.cell_volume[c];}
  error=std::sqrt(error/vol);
  std::cout<<"FV_V1 backend="<<backend_name(dev.backend())<<" tet="<<tet<<" n="<<n
    <<" speed="<<speed<<" scheme="<<(scheme==fv::ConvectionScheme::FirstOrderUpwind?1:2)
    <<" iterations="<<iterations+1<<" L2="<<error<<'\n';
  return error;
}
}

TEST_CASE("FV V1 manufactured scalar diffusion and upwind refinement" * doctest::test_suite("slow")) {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    for(int mode=0;mode<3;++mode) {
      const auto scheme=mode==1?fv::ConvectionScheme::FirstOrderUpwind:fv::ConvectionScheme::SecondOrderUpwind;
      const double speed=mode==0?0:4;
      double previous=0;
      for(int n:{8,16,32}) {
        double error=manufactured(*dev,n,tet,scheme,speed);
        if(previous>0) {
          double slope=std::log2(previous/error);
          std::cout<<"FV_V1_ORDER backend="<<backend_name(dev->backend())<<" tet="<<tet<<" mode="<<mode<<" n="<<n<<" slope="<<slope<<'\n';
          if(n==32) CHECK(slope>=(mode==0?1.9:(mode==2?1.8:.8)));
        }
        previous=error;
      }
    }
  }
}

TEST_CASE("FV transport kernels preserve an affine solution on skew tetrahedra") {
  auto mesh=make_box_tet(3,3,3,{0,0,0},{1,1,1});
  auto g=compute_geometry(mesh);
  const auto exact=[](Vec3d p) {return float(2+.3*p.x-.2*p.y+.1*p.z);};
  std::vector<float> field,boundary,flux,baseline;
  for(auto c:g.cell_centroid) field.push_back(exact(c));
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) boundary.push_back(exact(g.face_centroid[f]));
  for(auto a:g.face_area) flux.push_back(float(.8*a.x+.4*a.y-.2*a.z));
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);
    fv::TransportEquation equation(*dev,m,s,mesh,g);
    auto q=dev->upload<float>(field),bc=dev->upload<float>(boundary),mass=dev->upload<float>(flux);
    auto source=dev->alloc<float>(m.n_cells),mu=dev->alloc<float>(m.n_faces),r=dev->alloc<float>(m.n_cells);
    auto kind=dev->alloc<Index>(m.n_boundary_faces());dev->fill(kind,0);dev->fill(mu,.7f);dev->fill(source,.14f);
    fv::TransportSettings settings;settings.limit=false;
    equation.assemble(q,bc,kind,mass,mu,source,settings);
    residual(*dev,equation.matrix,q,equation.rhs,r);
    CHECK(dev->max_abs(r,r.size())<2e-6);
    auto got=dev->download(equation.matrix.val),rhs=dev->download(equation.rhs);
    got.insert(got.end(),rhs.begin(),rhs.end());
    if(baseline.empty()) baseline=got;
    else for(std::size_t i=0;i<got.size();++i) CHECK(got[i]==doctest::Approx(baseline[i]).epsilon(1e-6));
  }
}

TEST_CASE("FV equations share one uploaded CSR pattern and value workspace") {
  auto mesh=make_box_hex(2,2,2,{0,0,0},{1,1,1});auto g=compute_geometry(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto s=fv::make_stencil(*dev,mesh,g);
    auto shared=std::make_shared<CsrMatrix>(upload(*dev,s.pattern));
    fv::TransportEquation a(*dev,m,s,mesh,g,shared),b(*dev,m,s,mesh,g,shared);
    CHECK(a.matrix.row_ptr.handle()==b.matrix.row_ptr.handle());
    CHECK(a.matrix.col.handle()==b.matrix.col.handle());
    CHECK(a.matrix.diag.handle()==b.matrix.diag.handle());
    CHECK(a.matrix.val.handle()==b.matrix.val.handle());
  }
}

TEST_CASE("FV scalar preconditioning resolves an anisotropic operator and refreshes its values") {
  auto mesh=make_box_hex(64,16,1,{0,0,0},{1,.02,1});
  const auto geometry=compute_geometry(mesh);
  std::vector<Index> kinds;std::vector<float> boundary;
  for(Index f=mesh.n_internal_faces();f<mesh.n_faces();++f) {
    kinds.push_back(std::abs(geometry.face_area[f].x)>0?0:1);
    boundary.push_back(float(geometry.face_centroid[f].x));
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    const auto m=upload_mesh(*dev,mesh,geometry);
    const auto stencil=fv::make_stencil(*dev,mesh,geometry);
    auto workspace=std::make_shared<KrylovWorkspace>(m.n_cells);
    fv::TransportEquation equation(*dev,m,stencil,mesh,geometry,{}, {},workspace);
    equation.set_preconditioner(fv::TransportPreconditioner::Amg);
    equation.enable_compensated_solve();
    auto field=dev->alloc<float>(m.n_cells),mass=dev->alloc<float>(m.n_faces);
    auto mu=dev->alloc<float>(m.n_faces),source=dev->alloc<float>(m.n_cells);
    auto bc=dev->upload<float>(boundary);auto kind=dev->upload<Index>(kinds);
    dev->fill(mass,0);dev->fill(source,0);
    fv::TransportSettings settings;settings.limit=false;
    // The exact steady solution is x, with zero transverse derivative. A new
    // diffusivity exercises value refresh on the retained hierarchy.
    for(float diffusivity:{1.f,2.f}) {
      dev->fill(field,0);dev->fill(mu,diffusivity);
      for(int correction=0;correction<3;++correction) {
        equation.assemble(field,bc,kind,mass,mu,source,settings);
        // Acceptance is the physical solution after bounded defect corrections;
        // a float inner solve can return its best iterate at a residual floor.
        equation.solve(field,{1e-3,1e-20,100});
      }
      const auto values=dev->download(field);double error=0;
      for(Index c=0;c<m.n_cells;++c)
        error=std::max(error,std::abs(double(values[c])-geometry.cell_centroid[c].x));
      CHECK(error<2e-4);
    }
  }
}
