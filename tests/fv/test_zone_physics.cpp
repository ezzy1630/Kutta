#include <iostream>
#include <numbers>
#include "../backend_helpers.hpp"
#include "physics_meshes.hpp"
#include "kutta/fv/boundary_physics.hpp"
#include "kutta/setup/fluent_physics.hpp"

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;
namespace {
SolverSettings controls(double mu=1) {
  SolverSettings s;s.viscosity=mu;s.limit=false;s.pseudo_time_step=.1;
  s.momentum_linear={1e-5,1e-11,200};s.pressure_linear={1e-5,1e-12,300};return s;
}
double velocity_error(Device& dev,const SteadySolver& s,const std::function<Vec3d(Vec3d)>& exact) {
  const auto u=dev.download(s.velocity().x),v=dev.download(s.velocity().y),w=dev.download(s.velocity().z);
  double square=0,volume=0;
  for(Index c=0;c<s.mesh().n_cells;++c) {auto d=Vec3d{u[c],v[c],w[c]}-exact(s.geometry().cell_centroid[c]);
    square+=s.geometry().cell_volume[c]*dot(d,d);volume+=s.geometry().cell_volume[c];}
  return std::sqrt(square/volume);
}
double rotational_pressure(double r,bool solid) {
  if(solid) return .5*r*r;
  constexpr double a=-1./3,b=4./3;
  return .5*a*a*r*r+2*a*b*std::log(r)-.5*b*b/(r*r);
}
double rotational_pressure_error(Device& dev,const SteadySolver& s,bool solid) {
  const auto p=dev.download(s.pressure());double mean=0,volume=0,square=0;
  for(Index c=0;c<s.mesh().n_cells;++c) {
    const auto x=s.geometry().cell_centroid[c];const auto v=s.geometry().cell_volume[c];
    mean+=v*(p[c]-rotational_pressure(std::hypot(x.x,x.y),solid));volume+=v;
  }
  mean/=volume;
  for(Index c=0;c<s.mesh().n_cells;++c) {
    const auto x=s.geometry().cell_centroid[c];const auto d=p[c]-rotational_pressure(std::hypot(x.x,x.y),solid)-mean;
    square+=s.geometry().cell_volume[c]*d*d;
  }
  return std::sqrt(square/volume)/std::abs(rotational_pressure(2,solid)-rotational_pressure(1,solid));
}
}

TEST_CASE("FSAE tensor source and MRF kernels match the CPU oracle") {
  auto mesh=make_box_hex(3,1,1,{0,0,0},{1,1,.1});test::physics_zone(mesh);
  const auto g=compute_geometry(mesh);
  CellZonePhysics z{7,RigidFrame{{.2,-.1,0},{{.1,.2,0},{0,0,2}}},PorousResistance{{2,7,4},{3,1,5},{2,1,0},{1,3,1},.4,true}};
  std::vector<float> oracle;
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto dm=upload_mesh(*dev,mesh,g);ZonePhysics physics(*dev,mesh,g,{&z,1},1.3,.7);
    auto u=vector_field(*dev,mesh.n_cells);dev->fill(u.x,1.2f);dev->fill(u.y,-.3f);dev->fill(u.z,.4f);
    auto diagonal=dev->alloc<float>(mesh.n_cells),rhs=dev->alloc<float>(mesh.n_cells);
    const auto basis=porous_basis(*z.porous);
    for(int k=0;k<3;++k) {
      dev->fill(diagonal,0);dev->fill(rhs,0);physics.momentum(k,u,dm.volume,diagonal,rhs);
      const auto a=dev->download(diagonal),b=dev->download(rhs);
      for(Index c=0;c<mesh.n_cells;++c) {
        const Vec3d absolute{1.2,-.3,.4};const auto relative=absolute-z.frame->velocity(g.cell_centroid[c]);
        Vec3d force=cross(z.frame->rotation.angular_velocity,absolute-z.frame->translation)*(-1.3);
        const double d[]{2,7,4},inertia[]{3,1,5};
        for(int j=0;j<3;++j) force=force-basis[j]*((.7*d[j]+.65*norm(relative)*inertia[j])*dot(basis[j],relative));
        const double q[]{1.2,-.3,.4},f[]{force.x,force.y,force.z};
        CHECK((b[c]-a[c]*q[k])/g.cell_volume[c]==doctest::Approx(f[k]).epsilon(2e-6));
        if(dev->backend()==Backend::Cpu) oracle.push_back(b[c]);
        else CHECK(b[c]==doctest::Approx(oracle[k*mesh.n_cells+c]).epsilon(2e-6));
      }
    }
  }
}
TEST_CASE("FSAE invalid porous bases and MRF interfaces fail closed") {
  CHECK_THROWS(porous_basis(PorousResistance{{},{},{1,0,0},{2,0,0}}));
  auto mesh=make_box_hex(2,1,1,{0,0,0},{1,1,1});mesh.cell_zones={{"left","fluid",7,{0}},{"right","fluid",8,{1}}};
  CellZonePhysics z{7,RigidFrame{{1,0,0},{}},{}};
  auto dev=Device::create(Backend::Cpu);
  CHECK_THROWS_WITH(ZonePhysics(*dev,mesh,compute_geometry(mesh),{&z,1},1,1),doctest::Contains("zero normal frame velocity"));
}
TEST_CASE("FSAE clamped fan curve evaluates bounded speed without limiting flow") {
  auto mesh=make_box_hex(2,1,1,{0,0,0},{1,1,.1});const auto g=compute_geometry(mesh);
  mesh.internal_face_zones={{"fan","fan",71,{{0,1}}}};
  const InteriorCondition fan{71,FanCurve{{5,-1,-1},false,std::array<double,2>{.2,.5}}};
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto m=upload_mesh(*dev,mesh,g);auto stencil=make_stencil(*dev,mesh,g);
    InteriorJumps jump(*dev,mesh,g,{&fan,1},2,.1);
    std::vector<float> initial(mesh.n_faces(),0);initial[0]=.2f;
    auto mass=dev->upload<float>(initial),coefficient=dev->alloc<float>(mesh.n_faces()),diagonal=dev->alloc<float>(mesh.n_cells);
    dev->fill(coefficient,1);dev->fill(diagonal,1);
    jump.update(mass);jump.flux(m,stencil,diagonal,mass,coefficient);
    // |v|=1 is evaluated at the curve bound .5, giving 4.25 Pa. The
    // constitutive flux receives that head; its velocity itself is not clipped.
    CHECK(dev->download(mass)[0]==doctest::Approx(.285).epsilon(2e-6));
    CHECK(dev->download(coefficient)[0]==1);
  }
}

TEST_CASE("FSAE V7.1 translating-wall Couette hex and poly" * doctest::test_suite("slow")) {
  for(bool poly:{false,true}) {
    auto mesh=poly?test::physics_poly_box(4,12):make_box_hex(4,12,1,{0,0,0},{1,1,.1});
    const auto exact=[](Vec3d p){return Vec3d{p.y,0,0};};std::vector<PatchCondition> bc;
    for(auto p:mesh.patches) {
      if(p.name=="ymax") bc.push_back({p.name,Wall{{1,4,0}}});
      else if(p.name=="ymin") bc.push_back({p.name,Wall{}});
      else if(p.name=="xmin" || p.name=="xmax") bc.push_back({p.name,VelocityInlet{std::function<Vec3d(Vec3d)>(exact)}});
      else bc.push_back({p.name,Symmetry{}});
    }
    std::vector<float> oracle;
    KUTTA_FOR_EACH_BACKEND(dev) {
      SteadySolver s(*dev,mesh,bc,controls());s.iterate(400);
      const auto e=velocity_error(*dev,s,exact);const auto r=s.monitor();
      std::cout<<"FSAE V7.1 backend="<<backend_name(dev->backend())<<" poly="<<poly<<" velocity_L2="<<e<<" momentum="<<r.u.scaled<<std::endl;
      CHECK(e<2e-4);CHECK(r.u.scaled<1e-5);
      const auto result=dev->download(s.velocity().x);
      if(oracle.empty()) oracle=result;else for(std::size_t i=0;i<result.size();++i) CHECK(std::abs(result[i]-oracle[i])<2e-5);
    }
  }
}
TEST_CASE("FSAE V7.2 rotating-wall Taylor Couette hex and poly" * doctest::test_suite("slow")) {
  for(bool poly:{false,true}) {
    auto mesh=test::physics_annulus(poly?32:8,poly?256:64,poly);std::vector<PatchCondition> bc;
    for(auto p:mesh.patches) {
      if(p.name=="inner") bc.push_back({p.name,Wall{{},Rotation{{},{0,0,1}}}});
      else if(p.name=="outer") bc.push_back({p.name,Wall{}});
      else bc.push_back({p.name,Symmetry{}});
    }
    const auto exact=[](Vec3d p) {const double r2=p.x*p.x+p.y*p.y,angular=-1./3+4./(3*r2);return Vec3d{-angular*p.y,angular*p.x,0};};
    std::vector<float> oracle;
    KUTTA_FOR_EACH_BACKEND(dev) {
      SteadySolver s(*dev,mesh,bc,controls());
      // Refinement reduces the spatial error, while the unchanged residual gate
      // checks that the analytic seed has relaxed to this mesh's fixed point.
      s.initialize(exact,[](Vec3d p){return rotational_pressure(std::hypot(p.x,p.y),false);});s.iterate(1200);
      const auto e=velocity_error(*dev,s,exact);const auto r=s.monitor();
      const auto ep=rotational_pressure_error(*dev,s,false);
      std::cout<<"FSAE V7.2 backend="<<backend_name(dev->backend())<<" poly="<<poly<<" velocity_L2="<<e<<" pressure_relative_L2="<<ep<<" momentum="<<r.u.scaled<<std::endl;
      CHECK(e<.01);CHECK(ep<.01);CHECK(r.u.scaled<1e-5);
      const auto result=dev->download(s.velocity().x);
      if(oracle.empty()) oracle=result;else for(std::size_t i=0;i<result.size();++i) CHECK(std::abs(result[i]-oracle[i])<2e-5);
    }
  }
}
TEST_CASE("FSAE V8 Darcy Forchheimer pressure drop and oblique directions" * doctest::test_suite("slow")) {
  for(double angle:{0.,.6}) for(double porosity:{.35,1.}) {
    auto mesh=make_box_hex(24,3,1,{0,0,0},{1,1,.1});test::physics_zone(mesh);
    const Vec3d d{std::cos(angle),std::sin(angle),0},t{-d.y,d.x,0};
    for(auto& p:mesh.points) p={d.x*p.x+t.x*p.y,d.y*p.x+t.y*p.y,p.z};
    std::vector<PatchCondition> bc;
    for(auto p:mesh.patches) {
      if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{d}});
      else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{}});
      else bc.push_back({p.name,Symmetry{}});
    }
    auto opt=controls(.1);opt.density=2;opt.cell_zone_physics={{7,{},PorousResistance{{30,200,400},{4,20,40},d,d+t,porosity,true}}};
    KUTTA_FOR_EACH_BACKEND(dev) {
      SteadySolver s(*dev,mesh,bc,opt);s.iterate(400);
      const auto p=dev->download(s.pressure());double error=0;
      for(Index c=0;c<mesh.n_cells;++c) error=std::max(error,std::abs(p[c]-7*(1-dot(d,s.geometry().cell_centroid[c]))));
      const auto e=velocity_error(*dev,s,[d](Vec3d){return d;});
      std::cout<<"FSAE V8 backend="<<backend_name(dev->backend())<<" angle="<<angle<<" porosity="<<porosity<<" pressure_Linf="<<error<<" velocity_L2="<<e<<std::endl;
      CHECK(error<.007);CHECK(e<1e-4);CHECK(s.monitor().u.scaled<1e-5);
    }
  }
}

TEST_CASE("FSAE V9 absolute-velocity MRF analytic solutions and frame invariance" * doctest::test_suite("slow")) {
  for(bool solid:{true,false}) {
    auto mesh=test::physics_annulus(8,64,false);test::physics_zone(mesh);
    const auto exact=[solid](Vec3d p) {const auto a=solid?1.:(-1./3+4./(3*(p.x*p.x+p.y*p.y)));return Vec3d{-a*p.y,a*p.x,0};};
    std::vector<PatchCondition> bc;
    for(auto p:mesh.patches) {
      if(p.name=="inner" || (solid && p.name=="outer")) bc.push_back({p.name,Wall{{},Rotation{{},{0,0,1}}}});
      else if(p.name=="outer") bc.push_back({p.name,Wall{}});
      else bc.push_back({p.name,Symmetry{}});
    }
    KUTTA_FOR_EACH_BACKEND(dev) {
      std::vector<float> baseline;
      for(double omega:{0.,.4,1.}) {
        auto opt=controls();if(omega) opt.cell_zone_physics={{7,RigidFrame{{},{{},{0,0,omega}}},{}}};
        SteadySolver s(*dev,mesh,bc,opt);s.iterate(700);
        const auto error=velocity_error(*dev,s,exact);const auto u=dev->download(s.velocity().x);
        const auto ep=rotational_pressure_error(*dev,s,solid);
        double difference=0;
        if(baseline.empty()) baseline=u;
        else for(std::size_t i=0;i<u.size();++i) difference=std::max(difference,std::abs(double(u[i]-baseline[i])));
        std::cout<<"FSAE V9 backend="<<backend_name(dev->backend())<<" solid="<<solid<<" frame_omega="<<omega<<" velocity_L2="<<error<<" pressure_relative_L2="<<ep<<" frame_Linf="<<difference<<std::endl;
        CHECK(error<.01);CHECK(ep<.01);CHECK(difference<.003);CHECK(s.monitor().u.scaled<1e-5);
      }
    }
  }
}

TEST_CASE("FSAE V10 incompressible pressure and mass-flow inlets" * doctest::test_suite("slow")) {
  for(bool pressure:{false,true}) {
    auto mesh=make_box_hex(16,3,1,{0,0,0},{1,1,.1});std::vector<PatchCondition> bc;
    for(auto p:mesh.patches) {
      if(p.name=="xmin") {
        if(pressure) bc.push_back({p.name,PressureInlet{1}});
        else bc.push_back({p.name,MassFlowInlet{.2}});
      } else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{}});
      else bc.push_back({p.name,Symmetry{}});
    }
    auto opt=controls(.1);opt.density=2;
    KUTTA_FOR_EACH_BACKEND(dev) {
      SteadySolver s(*dev,mesh,bc,opt);s.iterate(500);
      const auto error=velocity_error(*dev,s,[](Vec3d){return Vec3d{1,0,0};});
      std::cout<<"FSAE V10 backend="<<backend_name(dev->backend())<<" pressure="<<pressure<<" velocity_L2="<<error<<std::endl;
      CHECK(error<1e-4);CHECK(dev->max_abs(s.pressure(),mesh.n_cells)<1e-4);
      CHECK(s.monitor().corrected_mass_imbalance<1e-6);
    }
  }
}

TEST_CASE("FSAE V11 internal fan polynomial and porous jump pressure discontinuity" * doctest::test_suite("slow")) {
  for(int mode:{0,1,2,3}) {
    auto mesh=make_box_hex(24,3,1,{0,0,0},{1,1,.1});const auto g=compute_geometry(mesh);
    InternalFaceZone zone{"sheet","fan",71,{}};
    for(Index f=0;f<mesh.n_internal_faces();++f) if(std::abs(g.face_centroid[f].x-.5)<1e-9)
      zone.faces.push_back({f,g.face_area[f].x>0?1:-1});
    mesh.internal_face_zones.push_back(zone);
    std::vector<PatchCondition> bc;
    const double speed=mode==3?-1:1;
    for(auto p:mesh.patches) {
      if(p.name=="xmin") bc.push_back({p.name,VelocityInlet{Vec3d{speed,0,0}}});
      else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,Vec3d{speed,0,0}}});
      else bc.push_back({p.name,Symmetry{}});
    }
    auto opt=controls(.1);opt.density=2;
    if(mode==0 || mode==1) opt.interior_conditions={{71,FanCurve{{5,-1,-1},mode==1}}};
    else opt.interior_conditions={{71,PorousJump{.01,.2,4}}};
    const double jump=mode==0?3:mode==1?-3:mode==2?-2.8:2.8;
    std::vector<float> oracle;
    KUTTA_FOR_EACH_BACKEND(dev) {
      SteadySolver s(*dev,mesh,bc,opt);s.iterate(600);
      const auto p=dev->download(s.pressure());double error=0;
      for(Index c=0;c<mesh.n_cells;++c) error=std::max(error,std::abs(p[c]-(s.geometry().cell_centroid[c].x<.5?-jump:0)));
      const auto velocity=velocity_error(*dev,s,[speed](Vec3d){return Vec3d{speed,0,0};});
      const auto r=s.monitor();
      std::cout<<"FSAE V11 backend="<<backend_name(dev->backend())<<" mode="<<mode<<" expected_jump="<<jump<<" pressure_Linf="<<error<<" velocity_L2="<<velocity<<" momentum="<<r.u.scaled<<std::endl;
      CHECK(error<.003);CHECK(velocity<1e-4);CHECK(r.u.scaled<1e-5);
      if(oracle.empty()) oracle=p;else for(std::size_t i=0;i<p.size();++i) CHECK(std::abs(p[i]-oracle[i])<1e-4);
    }
  }
}

TEST_CASE("FSAE internal face-zone source normals survive RCM") {
  auto mesh=make_box_hex(4,3,1,{0,0,0},{1,1,.1});const auto before=compute_geometry(mesh);
  InternalFaceZone zone{"sheet","fan",71,{}};Vec3d source_normal{};
  for(Index f=0;f<mesh.n_internal_faces();++f) if(std::abs(before.face_centroid[f].x-.5)<1e-9) {
    const auto direction=before.face_area[f].x>0?1:-1;zone.faces.push_back({f,direction});source_normal=source_normal+before.face_area[f]*double(direction);
  }
  mesh.internal_face_zones={zone};renumber_cells_rcm(mesh);const auto after=compute_geometry(mesh);Vec3d normal{};
  REQUIRE(mesh.internal_face_zones.size()==1);
  for(auto [f,direction]:mesh.internal_face_zones[0].faces) normal=normal+after.face_area[f]*double(direction);
  CHECK(norm(normal-source_normal)<1e-12);
}

TEST_CASE("FSAE V7.3 SST moving-wall and translating-frame invariance" * doctest::test_suite("slow")) {
  auto mesh=make_box_hex(4,12,1,{0,0,0},{1,1,.1});test::physics_zone(mesh);
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::vector<float> base_u,base_k,base_omega;
    for(double shift:{0.,2.}) {
      std::vector<PatchCondition> bc;
      for(auto p:mesh.patches) {
        if(p.name=="ymax") bc.push_back({p.name,Wall{{1+shift,0,0}}});
        else if(p.name=="ymin") bc.push_back({p.name,Wall{{shift,0,0}}});
        else if(p.name=="xmin" || p.name=="xmax") bc.push_back({p.name,VelocityInlet{
          std::function<Vec3d(Vec3d)>([shift](Vec3d x){return Vec3d{x.y+shift,0,0};}),TurbulenceValues{1e-4,2}}});
        else bc.push_back({p.name,Symmetry{}});
      }
      auto opt=controls();opt.sst=SstSettings{};
      if(shift) opt.cell_zone_physics={{7,RigidFrame{{shift,0,0},{}},{}}};
      SteadySolver s(*dev,mesh,bc,opt);s.iterate(400);
      const auto u=dev->download(s.velocity().x),k=dev->download(s.kinetic_energy()),omega=dev->download(s.specific_dissipation());
      double du=0,dk=0,dw=0;
      if(!shift) {base_u=u;base_k=k;base_omega=omega;}
      for(std::size_t i=0;i<u.size();++i) {
        du=std::max(du,std::abs(u[i]-shift-base_u[i]));dk=std::max(dk,std::abs(double(k[i]-base_k[i])));
        dw=std::max(dw,std::abs(double(omega[i]-base_omega[i]))/base_omega[i]);
      }
      std::cout<<"FSAE V7.3 backend="<<backend_name(dev->backend())<<" shift="<<shift<<" velocity_Linf="<<du<<" k_Linf="<<dk<<" omega_relative_Linf="<<dw<<std::endl;
      CHECK(du<2e-4);CHECK(dk<2e-7);CHECK(dw<2e-4);
      CHECK(std::isfinite(s.monitor().u.scaled));
    }
  }
}

TEST_CASE("FSAE V9.2 conformal MRF interface preserves absolute Taylor Couette" * doctest::test_suite("slow")) {
  auto mesh=test::physics_annulus(8,64,true);const auto g=compute_geometry(mesh);
  CellZone inner{"inner","fluid",7,{}},outer{"outer","fluid",8,{}};
  for(Index c=0;c<mesh.n_cells;++c) (std::hypot(g.cell_centroid[c].x,g.cell_centroid[c].y)<1.5?inner:outer).cells.push_back(c);
  mesh.cell_zones={inner,outer};auto opt=controls();opt.cell_zone_physics={{7,RigidFrame{{},{{},{0,0,1}}},{}}};
  std::vector<PatchCondition> bc;
  for(auto p:mesh.patches) {
    if(p.name=="inner") bc.push_back({p.name,Wall{{},{},true}}); // stationary RELATIVE to rotor
    else if(p.name=="outer") bc.push_back({p.name,Wall{}});
    else bc.push_back({p.name,Symmetry{}});
  }
  KUTTA_FOR_EACH_BACKEND(dev) {
    SteadySolver s(*dev,mesh,bc,opt);s.iterate(700);
    const auto e=velocity_error(*dev,s,[](Vec3d p){const auto a=-1./3+4./(3*(p.x*p.x+p.y*p.y));return Vec3d{-a*p.y,a*p.x,0};});
    std::cout<<"FSAE V9.2 backend="<<backend_name(dev->backend())<<" velocity_L2="<<e<<" momentum="<<s.monitor().u.scaled<<std::endl;
    CHECK(e<.01);CHECK(s.report().u.scaled<1e-5);
  }
}

TEST_CASE("FSAE V10.2 SST pressure and mass inlets match velocity-inlet transport" * doctest::test_suite("slow")) {
  const auto mesh=make_box_hex(12,2,1,{0,0,0},{1,1,.1});
  KUTTA_FOR_EACH_BACKEND(dev) {
    std::vector<float> base_k,base_omega;
    for(int mode:{0,1,2}) {
      std::vector<PatchCondition> bc;const TurbulenceCondition turbulence=IntensityLength{.05,.1};
      for(auto p:mesh.patches) {
        if(p.name=="xmin") {
          if(mode==0) bc.push_back({p.name,VelocityInlet{Vec3d{1,0,0},turbulence}});
          if(mode==1) bc.push_back({p.name,MassFlowInlet{.2,{},turbulence}});
          if(mode==2) bc.push_back({p.name,PressureInlet{1,{},turbulence}});
        } else if(p.name=="xmax") bc.push_back({p.name,PressureOutlet{0,NormalBackflow{},BackflowPressure::Static,turbulence}});
        else bc.push_back({p.name,Symmetry{}});
      }
      auto opt=controls(.01);opt.density=2;opt.sst=SstSettings{};
      SteadySolver s(*dev,mesh,bc,opt);
      const auto values=setup::turbulence_values(turbulence,1,opt.density,opt.viscosity);
      s.initialize_fields(std::vector<Vec3d>(mesh.n_cells,{1,0,0}),std::vector<double>(mesh.n_cells,0),
        std::vector<double>(mesh.n_cells,values.k),std::vector<double>(mesh.n_cells,values.omega));s.iterate(300);
      const auto k=dev->download(s.kinetic_energy()),omega=dev->download(s.specific_dissipation());
      if(mode==0) {base_k=k;base_omega=omega;}
      double dk=0,dw=0;
      for(std::size_t i=0;i<k.size();++i) {
        dk=std::max(dk,std::abs(double(k[i]-base_k[i]))/base_k[i]);dw=std::max(dw,std::abs(double(omega[i]-base_omega[i]))/base_omega[i]);
      }
      std::cout<<"FSAE V10.2 backend="<<backend_name(dev->backend())<<" mode="<<mode<<" k_relative_Linf="<<dk<<" omega_relative_Linf="<<dw<<std::endl;
      CHECK(dk<1e-4);CHECK(dw<1e-4);
    }
  }
}

TEST_CASE("FSAE V9.3 rotation and translation preserve Taylor Couette pressure" * doctest::test_suite("slow")) {
  auto mesh=test::physics_annulus(8,64,false);test::physics_zone(mesh);const Vec3d shift{.7,-.2,0};
  std::vector<PatchCondition> bc;
  for(auto p:mesh.patches) {
    if(p.name=="inner") bc.push_back({p.name,Wall{shift,Rotation{{},{0,0,1}}}});
    else if(p.name=="outer") bc.push_back({p.name,Wall{shift}});
    else bc.push_back({p.name,Symmetry{}});
  }
  auto opt=controls();opt.cell_zone_physics={{7,RigidFrame{shift,{{},{0,0,.4}}},{}}};
  KUTTA_FOR_EACH_BACKEND(dev) {
    SteadySolver s(*dev,mesh,bc,opt);s.iterate(700);
    const auto eu=velocity_error(*dev,s,[shift](Vec3d p){const auto a=-1./3+4./(3*(p.x*p.x+p.y*p.y));return shift+Vec3d{-a*p.y,a*p.x,0};});
    const auto ep=rotational_pressure_error(*dev,s,false);
    std::cout<<"FSAE V9.3 backend="<<backend_name(dev->backend())<<" velocity_L2="<<eu<<" pressure_relative_L2="<<ep<<std::endl;
    CHECK(eu<.01);CHECK(ep<.01);CHECK(s.monitor().u.scaled<1e-5);
  }
}
