#include "kutta/core/json.hpp"
#include "kutta/core/path.hpp"
#include "kutta/core/sha256.hpp"
#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <chrono>
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
std::string text(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  KUTTA_REQUIRE(in, "cannot read " + path.string());
  return {std::istreambuf_iterator<char>(in), {}};
}
void write_text(const std::filesystem::path &path, const std::string &s) {
  std::ofstream out(path, std::ios::binary);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out << s;
  out.close();
}
Json::Object skip(const std::string &name, const std::string *operation) {
  return {{"name", name},
          {"operation", operation ? Json(*operation) : Json(nullptr)}};
}
// Strict keeps every name; failures surface when the operation samples it.
// The available policy drops names this operation cannot sample.
std::vector<std::string> usable(const Volume &volume,
                                const std::vector<std::string> &names,
                                bool lenient, bool volume_only,
                                const std::string &operation,
                                Json::Array &skipped) {
  if (!lenient)
    return names;
  std::vector<std::string> out;
  for (const auto &name : names) {
    bool ok;
    try {
      ok = !volume.info(name).boundary_only || !volume_only;
    } catch (const Error &) {
      ok = false;
    }
    if (ok)
      out.push_back(name);
    else
      skipped.emplace_back(skip(name, &operation));
  }
  return out;
}

// One cached unit of work: a directory of output files under canonical names
// (data.bin, data.csv, zone files) and unit.json with descriptors/metadata.
// Entries are written to a temporary directory, then renamed into place.
class Cache {
  std::filesystem::path root_;
  std::set<std::string> used_;

public:
  explicit Cache(const std::filesystem::path &run) : root_(run / "post-cache") {
    std::filesystem::create_directories(root_);
  }
  std::filesystem::path directory(const std::string &key) const {
    return root_ / key;
  }
  std::optional<Json::Object> load(const std::string &key) {
    used_.insert(key);
    const auto file = root_ / key / "unit.json";
    if (!std::filesystem::exists(file))
      return {};
    try {
      // Trusted, self-written data: traces can carry megabytes of line
      // metadata, beyond the request-file limit.
      auto unit = parse_json(text(file), std::size_t(1) << 34);
      if (auto *o = std::get_if<Json::Object>(&unit.value))
        return std::move(*o);
    } catch (const std::exception &) {
    }
    // A damaged entry is recomputed, never trusted.
    std::error_code ec;
    std::filesystem::remove_all(root_ / key, ec);
    return {};
  }
  std::filesystem::path begin(const std::string &key) const {
    for (unsigned n = 0;; ++n) {
      const auto tmp = root_ / (".tmp-" + key + "-" + std::to_string(n));
      if (std::filesystem::create_directory(tmp))
        return tmp;
    }
  }
  void commit(const std::string &key, const std::filesystem::path &tmp,
              const Json::Object &unit) {
    write_text(tmp / "unit.json", dump_json(unit));
    std::error_code ec;
    std::filesystem::rename(tmp, root_ / key, ec);
    if (ec)
      std::filesystem::remove_all(tmp, ec);
  }
  // Bounded disk: keep exactly the entries the published view used.
  void prune() const {
    std::vector<std::filesystem::path> stale;
    for (const auto &entry : std::filesystem::directory_iterator(root_))
      if (!used_.contains(entry.path().filename().string()))
        stale.push_back(entry.path());
    for (const auto &path : stale) {
      std::error_code ec;
      std::filesystem::remove_all(path, ec);
    }
  }
};
// Hard links make reuse free for multi-GB outputs; both copies are immutable.
void place(const std::filesystem::path &from, const std::filesystem::path &to) {
  std::error_code ec;
  std::filesystem::create_hard_link(from, to, ec);
  if (ec)
    std::filesystem::copy_file(from, to);
}
// A run's solution identity without hashing a multi-GB snapshot: its size and
// mtime, plus the run's completion summary.
Json snapshot_identity(const std::filesystem::path &run) {
  const auto input = run / "view-input.bin";
  return Json::Object{
      {"path", utf8_path(std::filesystem::absolute(run))},
      {"size", std::to_string(std::filesystem::file_size(input))},
      {"mtime", std::to_string((long long)std::filesystem::last_write_time(input)
                                   .time_since_epoch()
                                   .count())},
      {"summary", sha256(text(run / "summary.json"))}};
}

Json::Object trace_unit(const Volume &volume, Trace spec, bool lenient,
                        const std::filesystem::path &dir) {
  Json::Array skipped;
  spec.fields = usable(volume, spec.fields, lenient, true, spec.name, skipped);
  const auto result = trace(volume, spec);
  auto d = view::write_geometry(dir / "data.bin", result.geometry).descriptor;
  d.emplace("name", spec.name);
  d.emplace("mode", spec.pathline ? "pathline" : "streamline");
  d.emplace("velocity_time_model",
            spec.frames.empty() ? "steady" : "linear_between_saved_frames");
  d.emplace("lines", result.lines);
  return {{"descriptor", std::move(d)}, {"unavailable", std::move(skipped)}};
}
Json::Object iso_unit(const Volume &volume, Iso spec, bool lenient,
                      const std::filesystem::path &dir) {
  Json::Array skipped;
  spec.fields = usable(volume, spec.fields, lenient, true, spec.name, skipped);
  auto d =
      view::write_geometry(dir / "data.bin", isosurface(volume, spec)).descriptor;
  d.emplace("name", spec.name);
  d.emplace("field", spec.field);
  d.emplace("value", spec.value);
  d.emplace("component", spec.component);
  return {{"descriptor", std::move(d)}, {"unavailable", std::move(skipped)}};
}
Json::Object probe_unit(const Volume &volume, Probe spec, bool lenient,
                        const std::filesystem::path &dir) {
  Json::Array skipped;
  const bool line = std::holds_alternative<LineProbe>(spec);
  std::visit(
      [&](auto &p) {
        p.fields = usable(volume, p.fields, lenient, line, p.name, skipped);
      },
      spec);
  const auto result = probe(volume, spec);
  auto d = view::write_geometry(dir / "data.bin", result.geometry).descriptor;
  d.emplace("name", std::visit([](const auto &s) { return s.name; }, spec));
  d.emplace("type", line ? "line" : "section");
  d.emplace("lines", result.lines);
  probe_csv(dir / "data.csv", result.geometry);
  d.emplace("csv", "data.csv");
  return {{"descriptor", std::move(d)}, {"unavailable", std::move(skipped)}};
}
Json::Object integral_unit(const Volume &volume, Integral spec, bool lenient) {
  Json::Array skipped;
  spec.fields = usable(volume, spec.fields, lenient, false, spec.name, skipped);
  return {{"report", integrate(volume, spec)},
          {"unavailable", std::move(skipped)}};
}
} // namespace
void apply_reference(view::Snapshot &snapshot, const Request &request) {
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
}
std::vector<std::string> surface_fields(const Volume &volume,
                                        const Request &request,
                                        Json::Array *skipped) {
  std::vector<std::string> names;
  if (request.fields.empty()) {
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
    return names;
  }
  for (const auto &name : request.fields) {
    if (!request.available_fields) {
      names.push_back(name);
      continue;
    }
    try {
      (void)volume.info(name);
      names.push_back(name);
    } catch (const Error &) {
      if (skipped)
        skipped->emplace_back(skip(name, nullptr));
    }
  }
  return names;
}
Derived derive(const Volume &volume, std::span<const std::string> names,
               bool cells, bool boundary) {
  const auto &mesh = volume.snapshot().mesh;
  Derived out;
  for (const auto &name : names) {
    const auto info = volume.info(name);
    if (cells && !info.boundary_only) {
      view::Field field{name, info.unit, info.components, {}};
      field.values.reserve(std::size_t(mesh.n_cells) * info.components);
      for (Index c = 0; c < mesh.n_cells; ++c) {
        const auto s =
            volume.sample(name, c, volume.geometry().cell_centroid[c]);
        for (std::uint32_t k = 0; k < s.components; ++k)
          field.values.push_back(float(s.value[k]));
      }
      out.cells.push_back(std::move(field));
    }
    if (!boundary)
      continue;
    const auto nb = mesh.n_faces() - mesh.n_internal_faces();
    view::Field field{
        name, info.unit, info.components,
        std::vector<float>(std::size_t(nb) * info.components, NAN)};
    for (const auto &patch : mesh.patches) {
      if (info.boundary_only && patch.kind != PatchKind::Wall)
        continue;
      for (Index f = patch.start; f < patch.start + patch.size; ++f) {
        const auto s =
            volume.boundary(name, f, volume.geometry().face_centroid[f]);
        for (std::uint32_t k = 0; k < s.components; ++k)
          field.values[std::size_t(f - mesh.n_internal_faces()) *
                           info.components +
                       k] = float(s.value[k]);
      }
    }
    out.boundary.push_back(std::move(field));
  }
  return out;
}
void attach(view::Snapshot &snapshot, Derived derived) {
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
  replace(snapshot.cells, std::move(derived.cells));
  replace(snapshot.boundary, std::move(derived.boundary));
}
void enrich(view::Snapshot &snapshot, const Request &request) {
  apply_reference(snapshot, request);
  Derived derived;
  {
    const Volume volume(snapshot);
    derived = derive(volume, surface_fields(volume, request), true, true);
  }
  attach(snapshot, std::move(derived));
}
CacheStats export_view(const std::filesystem::path &run,
                       std::span<const view::Plane> planes,
                       const Request &request, std::string_view version) {
  const auto start = std::chrono::steady_clock::now();
  CacheStats stats;
  // Every key shares the solution, references, policy and build identity.
  Json::Object reference;
  for (const auto &[k, v] : request.reference)
    reference.emplace(k, v);
  const Json common = Json::Object{
      {"format", "kutta-post-cache-1"},
      {"version", std::string(version)},
      {"snapshot", snapshot_identity(run)},
      {"reference", std::move(reference)},
      {"fields_policy", request.available_fields ? "available" : "strict"}};
  const auto key = [&](std::string_view kind, Json spec) {
    return sha256(dump_json(Json::Object{
        {"common", common}, {"kind", std::string(kind)}, {"spec", std::move(spec)}}));
  };
  Json::Array global;
  for (const auto &f : request.fields)
    global.emplace_back(f);
  struct Op {
    std::string kind, key, name;
    std::optional<Json::Object> unit;
  };
  std::vector<Op> ops;
  const auto add = [&](std::string kind, const std::string &identity,
                       const std::string &name, Json::Array frames = {}) {
    KUTTA_REQUIRE(!identity.empty(),
                  "post operation has no canonical identity: " + name);
    ops.push_back({kind,
                   key(kind, Json::Object{{"operation", identity},
                                          {"frames", std::move(frames)}}),
                   name,
                   {}});
  };
  for (const auto &t : request.traces) {
    // A pathline also depends on each saved frame's solution.
    Json::Array frames;
    for (const auto &f : t.frames)
      frames.emplace_back(snapshot_identity(f.run));
    add("streamline", t.identity, t.name, std::move(frames));
  }
  for (const auto &i : request.isosurfaces)
    add("iso", i.identity, i.name);
  for (const auto &p : request.probes)
    std::visit([&](const auto &s) { add("probe", s.identity, s.name); }, p);
  for (const auto &i : request.integrals)
    add("integral", i.identity, i.name);

  Cache cache(run);
  const auto base_key = key("base", global);
  auto base = cache.load(base_key);
  std::vector<std::string> slice_keys;
  std::vector<std::optional<Json::Object>> slices;
  for (const auto &p : planes) {
    slice_keys.push_back(
        key("slice", Json::Object{{"fields", global},
                                  {"origin", Json::Array{p.origin.x, p.origin.y, p.origin.z}},
                                  {"normal", Json::Array{p.normal.x, p.normal.y, p.normal.z}}}));
    slices.push_back(cache.load(slice_keys.back()));
  }
  for (auto &op : ops)
    op.unit = cache.load(op.key);
  const auto hit = [&](const auto &unit) { stats.reused += unit.has_value(); };
  hit(base);
  for (const auto &s : slices)
    hit(s);
  for (const auto &op : ops)
    hit(op.unit);
  const bool slices_missing = std::any_of(
      slices.begin(), slices.end(), [](const auto &s) { return !s; });
  const bool ops_missing = std::any_of(ops.begin(), ops.end(),
                                       [](const auto &o) { return !o.unit; });

  // Misses read the snapshot and construct one Volume; full hits do neither.
  if (!base || slices_missing || ops_missing) {
    auto snapshot = view::read_snapshot(run / "view-input.bin");
    apply_reference(snapshot, request);
    Derived derived;
    Json::Array global_skipped, available;
    {
      const Volume volume(snapshot);
      const auto compute = [&](Op &op, auto &&work) {
        const auto tmp = cache.begin(op.key);
        try {
          op.unit = work(tmp);
          cache.commit(op.key, tmp, *op.unit);
        } catch (const std::exception &e) {
          std::error_code ec;
          std::filesystem::remove_all(tmp, ec);
          throw Error("post operation \"" + op.name + "\": " + e.what());
        }
        ++stats.computed;
      };
      std::size_t n = 0;
      const bool lenient = request.available_fields;
      for (const auto &t : request.traces)
        if (auto &op = ops[n++]; !op.unit)
          compute(op, [&](const auto &dir) {
            return trace_unit(volume, t, lenient, dir);
          });
      for (const auto &i : request.isosurfaces)
        if (auto &op = ops[n++]; !op.unit)
          compute(op, [&](const auto &dir) {
            return iso_unit(volume, i, lenient, dir);
          });
      for (const auto &p : request.probes)
        if (auto &op = ops[n++]; !op.unit)
          compute(op, [&](const auto &dir) {
            return probe_unit(volume, p, lenient, dir);
          });
      for (const auto &i : request.integrals)
        if (auto &op = ops[n++]; !op.unit)
          compute(op, [&](const auto &) {
            return integral_unit(volume, i, lenient);
          });
      if (!base || slices_missing)
        derived = derive(volume, surface_fields(volume, request, &global_skipped),
                         slices_missing, !base);
      if (!base)
        for (const auto &field : volume.available())
          available.emplace_back(Json::Object{
              {"name", field.name},
              {"unit", field.unit},
              {"components", int(field.components)},
              {"location", field.boundary_only ? "boundary" : "volume"}});
    } // Volume caches gradients by name: attach only after it is gone.
    attach(snapshot, std::move(derived));
    if (!base) {
      const auto tmp = cache.begin(base_key);
      try {
        std::vector<Index> local(snapshot.mesh.points.size(), -1);
        Json::Array zones;
        for (std::size_t i = 0; i < snapshot.mesh.patches.size(); ++i)
          zones.emplace_back(
              view::write_zone(tmp, i, snapshot.mesh, snapshot, local, true));
        const auto &r = snapshot.properties;
        Json::Object values;
        for (const auto &[k, v] :
             std::vector<std::pair<std::string, std::optional<double>>>{
                 {"density", r.density},
                 {"viscosity", r.viscosity},
                 {"reference_density", r.reference_density},
                 {"velocity", r.velocity},
                 {"pressure", r.pressure},
                 {"area", r.area},
                 {"length", r.length}})
          values.emplace(k, v ? Json(*v) : Json(nullptr));
        base = Json::Object{
            {"header", view::manifest_header(snapshot.reference,
                                             snapshot.provenance,
                                             snapshot.mesh.points)},
            {"zones", std::move(zones)},
            {"available_fields", std::move(available)},
            {"reference_values", std::move(values)},
            {"unavailable", std::move(global_skipped)}};
        cache.commit(base_key, tmp, *base);
      } catch (...) {
        std::error_code ec;
        std::filesystem::remove_all(tmp, ec);
        throw;
      }
      ++stats.computed;
    }
    if (slices_missing) {
      const auto cf = build_cell_faces(snapshot.mesh);
      const auto eps = view::slice_tolerance(snapshot.mesh);
      for (std::size_t i = 0; i < planes.size(); ++i)
        if (!slices[i]) {
          const auto tmp = cache.begin(slice_keys[i]);
          try {
            slices[i] = Json::Object{
                {"descriptor", view::write_slice(tmp / "data.bin", snapshot,
                                                 planes[i], cf, eps)}};
            cache.commit(slice_keys[i], tmp, *slices[i]);
          } catch (...) {
            std::error_code ec;
            std::filesystem::remove_all(tmp, ec);
            throw;
          }
          ++stats.computed;
        }
    }
  }

  // Assemble the publication. Cached and computed parts take the same path,
  // so a warm export equals a cold one apart from post_cache.
  const auto stage = view::staging(run);
  try {
    const auto get = [](const Json::Object &o, const char *k) -> const Json & {
      return o.at(k);
    };
    auto manifest = std::get<Json::Object>(get(*base, "header").value);
    const auto base_dir = cache.directory(base_key);
    Json::Array zones = std::get<Json::Array>(get(*base, "zones").value);
    for (const auto &z : zones) {
      const auto &zone = std::get<Json::Object>(z.value);
      for (const auto *d : {&zone, zone.contains("edges")
                                       ? &std::get<Json::Object>(zone.at("edges").value)
                                       : nullptr})
        if (d) {
          const auto file = std::get<std::string>(d->at("file").value);
          place(base_dir / file, stage / file);
        }
    }
    Json::Array slice_entries, traces, isos, probes, reports;
    Json::Array unavailable = std::get<Json::Array>(get(*base, "unavailable").value);
    for (std::size_t i = 0; i < planes.size(); ++i) {
      auto d = std::get<Json::Object>(get(*slices[i], "descriptor").value);
      const auto file = "slice-" + std::to_string(i) + ".bin";
      place(cache.directory(slice_keys[i]) / "data.bin", stage / file);
      d.insert_or_assign("file", file);
      slice_entries.emplace_back(std::move(d));
    }
    std::map<std::string, std::size_t> ordinals;
    for (const auto &op : ops) {
      const auto &unit = *op.unit;
      for (const auto &s : std::get<Json::Array>(unit.at("unavailable").value))
        unavailable.push_back(s);
      if (op.kind == "integral") {
        reports.push_back(unit.at("report"));
        continue;
      }
      const auto i = std::to_string(ordinals[op.kind]++);
      const auto dir = cache.directory(op.key);
      auto d = std::get<Json::Object>(unit.at("descriptor").value);
      const auto prefix = op.kind == "streamline" ? "streamlines-"
                          : op.kind == "iso"      ? "iso-"
                                                  : "probe-";
      place(dir / "data.bin", stage / (prefix + i + ".bin"));
      d.insert_or_assign("file", prefix + i + ".bin");
      if (op.kind == "probe") {
        place(dir / "data.csv", stage / ("probe-" + i + ".csv"));
        d.insert_or_assign("csv", "probe-" + i + ".csv");
      }
      (op.kind == "streamline" ? traces : op.kind == "iso" ? isos : probes)
          .emplace_back(std::move(d));
    }
    if (!reports.empty()) {
      write_text(stage / "reports.json",
                 dump_json(Json::Object{{"v", 1}, {"reports", reports}}) + "\n");
      manifest.emplace("reports_file", "reports.json");
    }
    manifest.emplace("zones", std::move(zones));
    manifest.emplace("slices", std::move(slice_entries));
    manifest.emplace("streamlines", std::move(traces));
    manifest.emplace("isosurfaces", std::move(isos));
    manifest.emplace("probes", std::move(probes));
    manifest.emplace("reports", std::move(reports));
    manifest.emplace("available_fields", get(*base, "available_fields"));
    manifest.emplace("reference_values", get(*base, "reference_values"));
    manifest.emplace("unavailable_fields", std::move(unavailable));
    stats.seconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - start)
                        .count();
    manifest.emplace("post_cache",
                     Json::Object{{"reused", double(stats.reused)},
                                  {"computed", double(stats.computed)},
                                  {"seconds", stats.seconds}});
    view::write_manifest(stage / "manifest.json", std::move(manifest));
    view::publish(stage, run);
  } catch (...) {
    std::error_code ec;
    std::filesystem::remove_all(stage, ec);
    throw;
  }
  cache.prune();
  return stats;
}
} // namespace kutta::post
