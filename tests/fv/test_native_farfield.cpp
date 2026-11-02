#include "../backend_helpers.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <utility>

using namespace kutta;
using namespace kutta::fv;

TEST_CASE("Normal-foot pressure boundary gradients recover uniform flow on native tetrahedra") {
  // Eleven unchanged far-field tets from the generated native Ahmed mesh.
  // Slip faces or walls translating at U, plus cut inflow/outflow faces,
  // admit exact uniform flow. A persistent cycle cannot be wake shedding.
  const auto path=std::filesystem::path(KUTTA_SOURCE_DIR)/"tests/data/meshes/native_farfield_tets.msh";
  const auto mesh=read_fluent_msh(path);
  REQUIRE(mesh.n_cells==11);
  const auto geometry=compute_geometry(mesh);
  const TurbulenceCondition inflow=TurbulenceValues{1e-10,1};
  const std::vector<PatchCondition> base_bc{
    {"inlet",VelocityInlet{Vec3d{1,0,0},inflow}},
    {"outlet",PressureOutlet{0,NormalBackflow{},BackflowPressure::Static,inflow}},
    {"slip",Symmetry{}}};
  SolverSettings settings;
  settings.pressure_boundary_gradient=PressureBoundaryGradient::NormalFoot;
  settings.viscosity=1.044/2.78e6;
  settings.sst=SstSettings{};
  std::vector<Vec3d> velocity(mesh.n_cells,Vec3d{1,0,0});
  std::vector<double> pressure,kinetic(mesh.n_cells,1e-10),omega(mesh.n_cells,1);
  for(auto p:geometry.cell_centroid) pressure.push_back(.01*std::sin(37*p.x+13*p.y+23*p.z));
  KUTTA_FOR_EACH_BACKEND(dev) {
    for(bool moving_wall:{false,true}) {
      INFO(std::string(moving_wall?"moving wall":"symmetry"));
      auto bc=base_bc;
      if(moving_wall) bc.back().condition=Wall{Vec3d{1,0,0}};
      SteadySolver solver(*dev,mesh,bc,settings);
      solver.initialize_fields(velocity,pressure,kinetic,omega);
      solver.iterate(300);
      const auto report=solver.monitor();
      double error=0;
      const auto& u=solver.velocity();
      for(const auto& [component,expected]:{std::pair{&u.x,1.f},std::pair{&u.y,0.f},std::pair{&u.z,0.f}})
        for(auto value:dev->download(*component)) error=std::max(error,std::abs(double(value-expected)));
      CHECK(report.continuity.raw<1e-6);
      CHECK(error<1e-5);
      CHECK(dev->max_abs(solver.pressure(),mesh.n_cells)<1e-5);
    }

    auto unsupported=settings;
    unsupported.sst.reset();
    unsupported.gradient=GradientScheme::GreenGaussCell;
    CHECK_THROWS_WITH(SteadySolver(*dev,mesh,base_bc,unsupported),
      doctest::Contains("normal-foot pressure boundary gradients require least-squares gradients"));
  }
}
