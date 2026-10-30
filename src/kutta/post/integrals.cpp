#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <cmath>

namespace kutta::post {
using namespace detail;
namespace {
Json vector(Vec3d p) { return Json::Array{p.x, p.y, p.z}; }
Vec3d vec(const Sample &s) {
  KUTTA_REQUIRE(s.components == 3, "expected vector field");
  return {s.value[0], s.value[1], s.value[2]};
}
Json values(const std::array<double, 3> &v, std::uint32_t components,
            double denominator) {
  if (denominator <= 0)
    return nullptr;
  if (components == 1)
    return v[0] / denominator;
  return Json::Array{v[0] / denominator, v[1] / denominator,
                     v[2] / denominator};
}
struct Weighted {
  FieldInfo field;
  std::array<double, 3> area_sum{}, mass_sum{};
  double area = 0, mass = 0;
};
struct Result {
  double area = 0, mass = 0, absolute_mass = 0;
  Vec3d pressure, viscous, pressure_moment, viscous_moment;
  std::vector<Weighted> weighted;
  bool forces_valid = true;
};
Json pressure_center(const Result &r, const Integral &spec, double scale) {
  const double force2 = dot(r.pressure, r.pressure);
  if (force2 <= 0)
    return Json::Object{{"point", nullptr},
                        {"status", "zero_pressure_force"},
                        {"residual_moment", vector(r.pressure_moment)}};
  Vec3d location = spec.origin + cross(r.pressure, r.pressure_moment) / force2;
  if (spec.pressure_plane) {
    const auto &plane = *spec.pressure_plane;
    const double denominator = dot(plane.normal, r.pressure),
                 distance = dot(plane.normal, plane.origin - location);
    if (std::abs(denominator) <= 1e-12 * std::sqrt(force2)) {
      if (std::abs(distance) > scale * 1e-10)
        return Json::Object{{"point", nullptr},
                            {"status", "resultant_parallel_to_reference_plane"},
                            {"residual_moment", vector(r.pressure_moment)}};
    } else
      location = location + r.pressure * (distance / denominator);
  }
  const auto residual =
      r.pressure_moment - cross(location - spec.origin, r.pressure);
  return Json::Object{
      {"point", vector(location)},
      {"status", "defined"},
      {"constraint", spec.pressure_plane ? "reference_plane"
                                         : "closest_point_to_moment_origin"},
      {"residual_moment", vector(residual)},
      {"zero_moment",
       norm(residual) <= 1e-10 * std::max(1.0, norm(r.pressure_moment))}};
}
Json report(const Result &r, const Integral &spec,
            const view::Reference &reference, std::string name, int zone_id,
            double scale, bool solver_flux) {
  Json::Object averages;
  for (const auto &w : r.weighted)
    averages.emplace(
        w.field.name,
        Json::Object{
            {"unit", w.field.unit},
            {"area_weighted_average",
             values(w.area_sum, w.field.components, w.area)},
            {"mass_flow_weighted_average",
             values(w.mass_sum, w.field.components, w.mass)},
            {"area_integral", values(w.area_sum, w.field.components, 1)},
            {"area_coverage", r.area > 0 ? w.area / r.area : 0},
            {"mass_coverage",
             r.absolute_mass > 0 ? w.mass / r.absolute_mass : 0}});
  const auto total = r.pressure + r.viscous,
             moment = r.pressure_moment + r.viscous_moment;
  Json coefficients = nullptr, moment_coefficients = nullptr;
  if (reference.reference_density && *reference.reference_density > 0 &&
      reference.velocity && *reference.velocity > 0 && reference.area &&
      *reference.area > 0 && r.forces_valid) {
    const double scale_force = .5 * *reference.reference_density *
                               *reference.velocity * *reference.velocity *
                               *reference.area;
    const auto projections = [&](Vec3d f) {
      return Json::Object{{"drag", dot(f, spec.drag) / scale_force},
                          {"lift", dot(f, spec.lift) / scale_force},
                          {"side", dot(f, spec.side) / scale_force}};
    };
    coefficients = Json::Object{{"pressure", projections(r.pressure)},
                                {"viscous", projections(r.viscous)},
                                {"total", projections(total)}};
    if (reference.length && *reference.length > 0)
      moment_coefficients = Json::Object{
          {"pressure",
           vector(r.pressure_moment / (scale_force * *reference.length))},
          {"viscous",
           vector(r.viscous_moment / (scale_force * *reference.length))},
          {"total", vector(moment / (scale_force * *reference.length))}};
  }
  return Json::Object{
      {"zone", std::move(name)},
      {"zone_id", zone_id},
      {"area", r.area},
      {"mass_flow", r.mass},
      {"absolute_mass_flow", r.absolute_mass},
      {"mass_flow_source",
       solver_flux ? "solver_face_flux"
                   : "density_times_reconstructed_boundary_velocity"},
      {"averages", std::move(averages)},
      {"force", r.forces_valid
                    ? Json(Json::Object{{"pressure", vector(r.pressure)},
                                        {"viscous", vector(r.viscous)},
                                        {"total", vector(total)}})
                    : Json(nullptr)},
      {"moment", r.forces_valid ? Json(Json::Object{
                                      {"pressure", vector(r.pressure_moment)},
                                      {"viscous", vector(r.viscous_moment)},
                                      {"total", vector(moment)}})
                                : Json(nullptr)},
      {"coefficients", std::move(coefficients)},
      {"moment_coefficients", std::move(moment_coefficients)},
      {"centre_of_pressure",
       r.forces_valid ? pressure_center(r, spec, scale) : Json(nullptr)}};
}
} // namespace
Json integrate(const Volume &volume, const Integral &spec) {
  const auto &snapshot = volume.snapshot();
  const auto &mesh = snapshot.mesh;
  const auto &reference = snapshot.properties;
  const auto flux =
      std::find_if(snapshot.boundary.begin(), snapshot.boundary.end(),
                   [](const auto &f) { return f.name == "mass_flux"; });
  const bool solver_flux = flux != snapshot.boundary.end();
  KUTTA_REQUIRE(solver_flux || (reference.density && *reference.density > 0),
                "mass-flow reports need saved density or solver face flux");
  KUTTA_REQUIRE(reference.pressure, "force reports need saved reference "
                                    "pressure (or post reference override)");
  std::vector<const Patch *> zones;
  if (spec.zones.empty())
    for (const auto &p : mesh.patches)
      zones.push_back(&p);
  else
    for (const auto &name : spec.zones) {
      const auto p =
          std::find_if(mesh.patches.begin(), mesh.patches.end(),
                       [&](const auto &p) { return p.name == name; });
      KUTTA_REQUIRE(p != mesh.patches.end(), "unknown report zone: " + name);
      zones.push_back(&*p);
    }
  Result total;
  for (const auto &name : spec.fields)
    total.weighted.push_back({volume.info(name), {}, {}, 0, 0});
  Json::Array reports;
  for (const auto *zone : zones) {
    Result r;
    for (const auto &w : total.weighted)
      r.weighted.push_back({w.field, {}, {}, 0, 0});
    for (Index face = zone->start; face < zone->start + zone->size; ++face) {
      const auto area_vector = volume.geometry().face_area[face],
                 center = volume.geometry().face_centroid[face];
      const double area = norm(area_vector);
      r.area += area;
      const auto velocity = vec(volume.boundary("velocity", face, center));
      const double mass = solver_flux
                              ? flux->values[face - mesh.n_internal_faces()]
                              : *reference.density * dot(velocity, area_vector);
      KUTTA_REQUIRE(std::isfinite(mass),
                    "nonfinite mass flux in report zone: " + zone->name);
      r.mass += mass;
      r.absolute_mass += std::abs(mass);
      for (auto &w : r.weighted) {
        const auto s = volume.boundary(w.field.name, face, center);
        bool valid = true;
        for (std::uint32_t k = 0; k < s.components; ++k)
          valid &= std::isfinite(s.value[k]);
        if (valid) {
          w.area += area;
          w.mass += std::abs(mass);
          for (std::uint32_t k = 0; k < s.components; ++k) {
            w.area_sum[k] += s.value[k] * area;
            w.mass_sum[k] += s.value[k] * std::abs(mass);
          }
        }
      }
      const auto pressure = volume.boundary("pressure", face, center).value[0];
      Vec3d viscous;
      if (zone->kind == PatchKind::Wall) {
        const auto shear = volume.boundary("wall_shear", face, center);
        viscous = vec(shear) *
                  (-area); // stored stress is on fluid; report load on body
      } else {
        KUTTA_REQUIRE(reference.viscosity && *reference.viscosity > 0,
                      "non-wall viscous traction requires molecular viscosity");
        const auto gradient = volume.velocity_gradient(mesh.owner[face]);
        const auto normal = area_vector / area;
        double mu = *reference.viscosity;
        const auto ratio = std::find_if(
            snapshot.cells.begin(), snapshot.cells.end(),
            [](const auto &f) { return f.name == "viscosity_ratio"; });
        if (ratio != snapshot.cells.end())
          mu *= 1 + ratio->values[mesh.owner[face]];
        const double divergence = gradient[0].x + gradient[1].y + gradient[2].z;
        const Vec3d traction{2 * gradient[0].x * normal.x +
                                 (gradient[0].y + gradient[1].x) * normal.y +
                                 (gradient[0].z + gradient[2].x) * normal.z,
                             (gradient[1].x + gradient[0].y) * normal.x +
                                 2 * gradient[1].y * normal.y +
                                 (gradient[1].z + gradient[2].y) * normal.z,
                             (gradient[2].x + gradient[0].z) * normal.x +
                                 (gradient[2].y + gradient[1].z) * normal.y +
                                 2 * gradient[2].z * normal.z};
        viscous = (traction - normal * (2 * divergence / 3)) * (-mu * area);
      }
      const auto force = area_vector * (pressure - *reference.pressure);
      r.forces_valid &= std::isfinite(pressure) && finite(viscous);
      r.pressure = r.pressure + force;
      r.viscous = r.viscous + viscous;
      r.pressure_moment =
          r.pressure_moment + cross(center - spec.origin, force);
      r.viscous_moment =
          r.viscous_moment + cross(center - spec.origin, viscous);
    }
    reports.push_back(report(r, spec, reference, zone->name,
                             zone->source_zone_id, volume.scale(),
                             solver_flux));
    total.area += r.area;
    total.mass += r.mass;
    total.absolute_mass += r.absolute_mass;
    total.pressure = total.pressure + r.pressure;
    total.viscous = total.viscous + r.viscous;
    total.pressure_moment = total.pressure_moment + r.pressure_moment;
    total.viscous_moment = total.viscous_moment + r.viscous_moment;
    total.forces_valid &= r.forces_valid;
    for (std::size_t i = 0; i < r.weighted.size(); ++i) {
      auto &to = total.weighted[i];
      const auto &from = r.weighted[i];
      to.area += from.area;
      to.mass += from.mass;
      for (std::uint32_t k = 0; k < to.field.components; ++k) {
        to.area_sum[k] += from.area_sum[k];
        to.mass_sum[k] += from.mass_sum[k];
      }
    }
  }
  return Json::Object{{"name", spec.name},
                      {"origin", vector(spec.origin)},
                      {"directions", Json::Object{{"drag", vector(spec.drag)},
                                                  {"lift", vector(spec.lift)},
                                                  {"side", vector(spec.side)}}},
                      {"units", Json::Object{{"area", "m^2"},
                                             {"mass_flow", "kg/s"},
                                             {"force", "N"},
                                             {"moment", "N m"},
                                             {"centre_of_pressure", "m"}}},
                      {"zones", std::move(reports)},
                      {"total", report(total, spec, reference, "selected-zones",
                                       -1, volume.scale(), solver_flux)}};
}
} // namespace kutta::post
