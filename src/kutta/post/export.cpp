#include "kutta/core/json.hpp"
#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <set>

namespace kutta::post {
using namespace detail;
namespace {
std::string csv(std::string_view s) {
  std::string out = "\"";
  for (char c : s)
    out += c == '"' ? "\"\"" : std::string(1, c);
  return out + '"';
}
void probe_csv(const std::filesystem::path &path, const view::Geometry &g) {
  std::ofstream out(path);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out.imbue(std::locale::classic());
  out << std::setprecision(9) << "sample,x,y,z";
  for (const auto &f : g.fields)
    for (std::uint32_t c = 0; c < f.components; ++c)
      out << ','
          << csv(f.name +
                 (f.components == 3 ? std::string(".") + "xyz"[c] : ""));
  out << '\n';
  for (std::size_t i = 0; i < g.positions.size() / 3; ++i) {
    out << i << ',' << g.positions[i * 3] << ',' << g.positions[i * 3 + 1]
        << ',' << g.positions[i * 3 + 2];
    for (const auto &f : g.fields)
      for (std::uint32_t c = 0; c < f.components; ++c) {
        const float v = f.values[i * f.components + c];
        out << ',';
        if (std::isfinite(v))
          out << v;
      }
    out << '\n';
  }
  out.close();
}
} // namespace
void enrich(view::Snapshot &snapshot, const Request &request) {
  for (const auto &[key, value] : request.reference) {
    auto &r = snapshot.properties;
    if (key == "density") {
      KUTTA_REQUIRE(!r.density || *r.density == value,
                    "density is a solved material property; change it with "
                    "solve --setup");
      r.density = value;
    } else if (key == "viscosity") {
      KUTTA_REQUIRE(!r.viscosity || *r.viscosity == value,
                    "viscosity is a solved material property; change it with "
                    "solve --setup");
      r.viscosity = value;
    } else if (key == "reference_density")
      r.reference_density = value;
    else if (key == "velocity")
      r.velocity = value;
    else if (key == "pressure")
      r.pressure = value;
    else if (key == "area")
      r.area = value;
    else if (key == "length")
      r.length = value;
  }
  if (!request.reference.empty()) {
    if (!snapshot.provenance.empty())
      snapshot.provenance += "; ";
    snapshot.provenance += "POST REFERENCE OVERRIDES: ";
    for (const auto &[key, value] : request.reference)
      snapshot.provenance += key + "=" + std::to_string(value) + " ";
  }
  const Volume volume(snapshot);
  std::vector<std::string> names = request.fields;
  if (names.empty())
    for (auto suffix : {"", "@fluent"})
      for (auto base :
           {"pressure_coefficient", "total_pressure",
            "total_pressure_coefficient", "velocity_x", "velocity_y",
            "velocity_z", "vorticity", "vorticity_magnitude", "q_criterion",
            "helicity", "skin_friction", "skin_friction_vector"})
        try {
          const auto name = std::string(base) + suffix;
          (void)volume.info(name);
          names.push_back(name);
        } catch (const Error &) {
        }
  std::vector<view::Field> cells, boundary;
  for (const auto &name : names) {
    const auto info = volume.info(name);
    if (!info.boundary_only) {
      view::Field field{name, info.unit, info.components, {}};
      field.values.reserve(std::size_t(snapshot.mesh.n_cells) *
                           info.components);
      for (Index c = 0; c < snapshot.mesh.n_cells; ++c) {
        const auto s =
            volume.sample(name, c, volume.geometry().cell_centroid[c]);
        for (std::uint32_t k = 0; k < s.components; ++k)
          field.values.push_back(float(s.value[k]));
      }
      cells.push_back(std::move(field));
    }
    const auto nb = snapshot.mesh.n_faces() - snapshot.mesh.n_internal_faces();
    view::Field field{
        name, info.unit, info.components,
        std::vector<float>(std::size_t(nb) * info.components, NAN)};
    for (const auto &patch : snapshot.mesh.patches) {
      if (info.boundary_only && patch.kind != PatchKind::Wall)
        continue;
      for (Index f = patch.start; f < patch.start + patch.size; ++f) {
        const auto s =
            volume.boundary(name, f, volume.geometry().face_centroid[f]);
        for (std::uint32_t k = 0; k < s.components; ++k)
          field.values[std::size_t(f - snapshot.mesh.n_internal_faces()) *
                           info.components +
                       k] = float(s.value[k]);
      }
    }
    boundary.push_back(std::move(field));
  }
  const auto replace = [](std::vector<view::Field> &target,
                          std::vector<view::Field> source) {
    for (auto &field : source) {
      auto old = std::find_if(target.begin(), target.end(), [&](const auto &f) {
        return f.name == field.name;
      });
      if (old == target.end())
        target.push_back(std::move(field));
      else
        *old = std::move(field);
    }
  };
  replace(snapshot.cells, std::move(cells));
  replace(snapshot.boundary, std::move(boundary));
}
Outputs write_outputs(const std::filesystem::path &stage,
                      const view::Snapshot &snapshot, const Request &request) {
  const Volume volume(snapshot);
  Outputs out;
  Json::Array traces, isos, probes, reports, catalog;
  const auto write = [&](const std::string &file,
                         const view::Geometry &geometry) {
    auto artifact = view::write_geometry(stage / file, geometry);
    out.fields.insert(out.fields.end(), artifact.fields.begin(),
                      artifact.fields.end());
    return std::move(artifact.descriptor);
  };
  for (std::size_t i = 0; i < request.traces.size(); ++i) {
    const auto &spec = request.traces[i];
    const auto result = trace(volume, spec);
    auto descriptor =
        write("streamlines-" + std::to_string(i) + ".bin", result.geometry);
    descriptor.emplace("name", spec.name);
    descriptor.emplace("mode", spec.pathline ? "pathline" : "streamline");
    descriptor.emplace("velocity_time_model",
                       spec.frames.empty() ? "steady"
                                           : "linear_between_saved_frames");
    descriptor.emplace("lines", result.lines);
    traces.emplace_back(std::move(descriptor));
  }
  for (std::size_t i = 0; i < request.isosurfaces.size(); ++i) {
    const auto &spec = request.isosurfaces[i];
    auto descriptor =
        write("iso-" + std::to_string(i) + ".bin", isosurface(volume, spec));
    descriptor.emplace("name", spec.name);
    descriptor.emplace("field", spec.field);
    descriptor.emplace("value", spec.value);
    descriptor.emplace("component", spec.component);
    isos.emplace_back(std::move(descriptor));
  }
  for (std::size_t i = 0; i < request.probes.size(); ++i) {
    const auto result = probe(volume, request.probes[i]);
    auto descriptor =
        write("probe-" + std::to_string(i) + ".bin", result.geometry);
    descriptor.emplace("name", std::visit([](const auto &s) { return s.name; },
                                          request.probes[i]));
    descriptor.emplace("type",
                       std::holds_alternative<LineProbe>(request.probes[i])
                           ? "line"
                           : "section");
    descriptor.emplace("lines", result.lines);
    const auto csv_file = "probe-" + std::to_string(i) + ".csv";
    probe_csv(stage / csv_file, result.geometry);
    descriptor.emplace("csv", csv_file);
    probes.emplace_back(std::move(descriptor));
  }
  for (const auto &spec : request.integrals)
    reports.push_back(integrate(volume, spec));
  if (!reports.empty()) {
    std::ofstream file(stage / "reports.json");
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << dump_json(Json::Object{{"v", 1}, {"reports", reports}}) << '\n';
    file.close();
    out.manifest.emplace("reports_file", "reports.json");
  }
  for (const auto &field : volume.available())
    catalog.emplace_back(Json::Object{
        {"name", field.name},
        {"unit", field.unit},
        {"components", int(field.components)},
        {"location", field.boundary_only ? "boundary" : "volume"}});
  out.manifest.emplace("streamlines", std::move(traces));
  out.manifest.emplace("isosurfaces", std::move(isos));
  out.manifest.emplace("probes", std::move(probes));
  out.manifest.emplace("reports", std::move(reports));
  out.manifest.emplace("available_fields", std::move(catalog));
  const auto &r = snapshot.properties;
  Json::Object reference;
  for (const auto &[key, value] :
       std::vector<std::pair<std::string, std::optional<double>>>{
           {"density", r.density},
           {"viscosity", r.viscosity},
           {"reference_density", r.reference_density},
           {"velocity", r.velocity},
           {"pressure", r.pressure},
           {"area", r.area},
           {"length", r.length}})
    reference.emplace(key, value ? Json(*value) : Json(nullptr));
  out.manifest.emplace("reference_values", std::move(reference));
  return out;
}
} // namespace kutta::post
