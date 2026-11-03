#include <cmath>

#include "../backend_helpers.hpp"
#include "kutta/fv/face_projection.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "kutta/mesh/detail.hpp"

using namespace kutta;
using namespace kutta::fv;
using namespace kutta::detail;

TEST_CASE("FV cubic face diffusion integrates a cubic field on hexes and tetrahedra") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    const auto mesh=tet?make_box_tet(4,4,4,{0,0,0},{1,1,1}):make_box_hex(4,4,4,{0,0,0},{1,1,1});
    const auto g=compute_geometry(mesh);const auto m=upload_mesh(*dev,mesh,g);const auto st=make_stencil(*dev,mesh,g);
    const auto value=[](Vec3d p){return p.x*p.x*p.x;};
    std::vector<float> cell,boundary;
    for(auto p:g.cell_centroid) cell.push_back(float(value(p)));
    for(Index f=m.n_internal_faces;f<m.n_faces;++f) boundary.push_back(float(value(g.face_centroid[f])));
    auto q=dev->upload<float>(cell),bc=dev->upload<float>(boundary),correction=dev->alloc<float>(m.n_faces);
    const std::vector<BoundaryObservation> rows(m.n_boundary_faces(),BoundaryObservation::FaceValue);
    ProjectedFaceReconstruction fit(*dev,mesh,g,rows,{},ProjectedFaceReconstruction::Mode::Cubic);
    REQUIRE(fit.fallback_faces()==0);
    fit.diffusion_correction(q,bc,correction);
    const auto secondary=dev->download(correction),primary=dev->download(st.diffusion);
    std::vector<double> divergence(m.n_cells);
    for(Index f=0;f<m.n_faces;++f) {
      const double other=f<m.n_internal_faces?cell[mesh.neighbour[f]]:boundary[f-m.n_internal_faces];
      const double flux=primary[f]*(other-cell[mesh.owner[f]])+secondary[f];
      divergence[mesh.owner[f]]+=flux;
      if(f<m.n_internal_faces) divergence[mesh.neighbour[f]]-=flux;
    }
    double maximum=0;
    for(Index c=0;c<m.n_cells;++c)
      maximum=std::max(maximum,std::abs(divergence[c]/g.cell_volume[c]-6*g.cell_centroid[c].x));
    CHECK(maximum<1e-4);
  }
}

TEST_CASE("FV automatic volumetric laminar reconstruction resolves quadratic wall shear") {
  KUTTA_FOR_EACH_BACKEND(dev) for(bool tet:{false,true}) {
    auto mesh=tet?make_box_tet(3,3,3,{0,0,0},{1,1,1}):make_box_hex(3,3,3,{0,0,0},{1,1,1});
    const auto velocity=[](Vec3d p){return Vec3d{0,p.x+p.x*p.x,0};};
    std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches)
      if(p.name=="xmin") bc.push_back({p.name,Wall{}});
      else bc.push_back({p.name,VelocityInlet{velocity}});
    SolverSettings settings;settings.viscosity=1;settings.laminar_reconstruction=LaminarReconstruction::Automatic;
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
      SolverSettings settings;settings.viscosity=1;settings.laminar_reconstruction=mode;
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

TEST_CASE("FV explicit cubic reconstruction rejects unsupported gradient and wall-function contracts") {
  KUTTA_FOR_EACH_BACKEND(dev) {
    auto mesh=make_box_hex(2,2,2,{0,0,0},{1,1,1});std::vector<PatchCondition> bc;
    for(const auto& p:mesh.patches) bc.push_back({p.name,Wall{}});
    SolverSettings settings;settings.laminar_reconstruction=LaminarReconstruction::Cubic;
    settings.gradient=GradientScheme::GreenGaussCell;
    CHECK_THROWS_WITH(SteadySolver(*dev,mesh,bc,settings),doctest::Contains("cubic laminar reconstruction requires"));
    settings.gradient=GradientScheme::LeastSquares;settings.sst=SstSettings{};
    CHECK_THROWS_WITH(SteadySolver(*dev,mesh,bc,settings),doctest::Contains("cubic laminar reconstruction requires"));
  }
}
