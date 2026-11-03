#include <cmath>

#include "../backend_helpers.hpp"
#include "kutta/fv/face_projection.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;

TEST_CASE("FV integrated boundary diffusion resolves a cubic wall flux on hexes and tetrahedra") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    const auto mesh=tet?make_box_tet(4,4,4,{0,0,0},{1,1,1}):make_box_hex(4,4,4,{0,0,0},{1,1,1});
    const auto g=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,g);const auto st=make_stencil(*dev,mesh,g);
    const auto value=[](Vec3d p){return p.x*p.y*p.y;};
    std::vector<float> cell,boundary;
    for(auto p:g.cell_centroid) cell.push_back(float(value(p)));
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(g.face_centroid[f])));
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary),correction=dev->alloc<float>(m.n_faces);
    const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    ProjectedFaceReconstruction fit(*dev,mesh,g,rows,{},ProjectedFaceReconstruction::Mode::Integrated);
    REQUIRE(fit.fallback_faces()==0);
    fit.diffusion_correction(q,bc,correction);
    const auto secondary=dev->download(correction),primary=dev->download(st.diffusion);
    double flux=0;
    for(const auto& patch:mesh.patches) if(patch.name=="xmin")
      for(Index f=patch.start;f<patch.start+patch.size;++f)
        flux+=primary[f]*(boundary[f-m.n_internal_faces]-cell[mesh.owner[f]])+secondary[f];
    // Integral of -d(x*y*y)/dx over x=0 is -integral_0^1 y*y dy.
    CHECK(std::abs(flux+1./3)<1e-6);
  }
}

TEST_CASE("FV automatic unlimited laminar reconstruction resolves quadratic wall shear") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    auto mesh=tet?make_box_tet(3,3,3,{0,0,0},{1,1,1}):make_box_hex(3,3,3,{0,0,0},{1,1,1});
    const auto velocity=[](Vec3d p){return Vec3d{0,p.x+p.x*p.x,0};};
    std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches)
      if(p.name=="xmin") bc.push_back({p.name,Wall{}});
      else bc.push_back({p.name,VelocityInlet{velocity}});
    SolverSettings settings;settings.viscosity=1;settings.limit=false;settings.laminar_reconstruction=LaminarReconstruction::Automatic;
    SteadySolver solver(*dev,std::move(mesh),bc,settings);
    solver.initialize(velocity,[](Vec3d){return 0.;});
    const std::vector<std::string> walls{"xmin"};
    const auto f=solver.forces(walls);
    CHECK(norm(f.viscous-Vec3d{0,1,0})<1e-5);
    CHECK(norm(f.pressure)<1e-7);
  }
}

TEST_CASE("FV automatic laminar reconstruction preserves rotated single-layer meshes") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=make_box_hex(4,4,1,{0,0,0},{1,1,.25});
    const double c=std::cos(.4),s=std::sin(.4);
    for(auto& p:mesh.points) p={c*p.x-s*p.z,p.y,s*p.x+c*p.z};
    const auto velocity=[=](Vec3d p){const double x=c*p.x+s*p.z;return Vec3d{0,x+x*x,0};};
    std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches)
      if(p.name=="xmin") bc.push_back({p.name,Wall{}});
      else bc.push_back({p.name,VelocityInlet{velocity}});
    ForceReport previous;
    for(auto mode:{LaminarReconstruction::Projected,LaminarReconstruction::Automatic}) {
      SolverSettings settings;settings.viscosity=1;settings.limit=false;settings.laminar_reconstruction=mode;
      SteadySolver solver(*dev,mesh,bc,settings);
      solver.initialize(velocity,[](Vec3d){return 0.;});
      const std::vector<std::string> walls{"xmin"};const auto f=solver.forces(walls);
      if(mode==LaminarReconstruction::Automatic) {
        CHECK(f.viscous.x==previous.viscous.x);
        CHECK(f.viscous.y==previous.viscous.y);
        CHECK(f.viscous.z==previous.viscous.z);
      }
      previous=f;
    }
  }
}

TEST_CASE("FV integrated reconstruction rejects unsupported gradient wall-function and limiter contracts") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=make_box_hex(2,2,2,{0,0,0},{1,1,1});std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches) bc.push_back({p.name,Wall{}});
    SolverSettings settings;settings.limit=false;settings.laminar_reconstruction=LaminarReconstruction::Integrated;
    settings.gradient=GradientScheme::GreenGaussCell;
    CHECK_THROWS_WITH(SteadySolver(*dev,mesh,bc,settings),doctest::Contains("integrated laminar reconstruction requires"));
    settings.gradient=GradientScheme::LeastSquares;settings.sst=SstSettings{};
    CHECK_THROWS_WITH(SteadySolver(*dev,mesh,bc,settings),doctest::Contains("integrated laminar reconstruction requires"));
    settings.sst.reset();settings.limit=true;
    CHECK_THROWS_WITH(SteadySolver(*dev,mesh,bc,settings),doctest::Contains("requires unlimited convection"));
  }
}

TEST_CASE("FV integrated unlimited tet channel converges and limited automatic keeps the original operator") {
  const auto mesh=make_box_tet(8,4,1,{0,0,0},{2,1,.2});
  std::vector<PatchCondition> bc;
  const auto inlet=[](Vec3d p){return Vec3d{6*p.y*(1-p.y),0,0};};
  for(const auto& p:mesh.patches) {
    if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{inlet}});
    else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0}});
    else if(p.name=="zmin" || p.name=="zmax") bc.push_back({p.name,Symmetry{}});
    else bc.push_back({p.name,Wall{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SolverSettings settings;settings.viscosity=.05;settings.limit=false;settings.monitor_interval=50;
    SteadySolver integrated(*dev,mesh,bc,settings);
    integrated.initialize([](Vec3d){return Vec3d{};},[](Vec3d){return 0.;});
    integrated.iterate(400);const auto r=integrated.monitor();
    CHECK(std::max(r.u.scaled,r.v.scaled)<2e-7);
    CHECK(r.continuity.raw<1e-6);
    std::vector<float> u,v,p;
    settings.limit=true;
    for(auto mode:{LaminarReconstruction::Projected,LaminarReconstruction::Automatic}) {
      settings.laminar_reconstruction=mode;
      SteadySolver limited(*dev,mesh,bc,settings);
      limited.initialize([](Vec3d){return Vec3d{};},[](Vec3d){return 0.;});limited.iterate(20);
      const auto current_u=dev->download(limited.velocity().x),current_v=dev->download(limited.velocity().y),current_p=dev->download(limited.pressure());
      if(mode==LaminarReconstruction::Automatic) {CHECK(current_u==u);CHECK(current_v==v);CHECK(current_p==p);}
      u=current_u;v=current_v;p=current_p;
    }
  }
}
