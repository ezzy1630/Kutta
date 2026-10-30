#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <cmath>
#include <random>

namespace kutta::post {
using namespace detail;
namespace {
Vec3d vec(const Sample &s) {
  KUTTA_REQUIRE(s.components == 3, "velocity must have three components");
  return {s.value[0], s.value[1], s.value[2]};
}
const Patch &patch(const PolyMesh &m, const std::string &name) {
  auto it = std::find_if(m.patches.begin(), m.patches.end(),
                         [&](const auto &p) { return p.name == name; });
  KUTTA_REQUIRE(it != m.patches.end(), "unknown surface zone: " + name);
  return *it;
}
double radical(std::uint64_t n, unsigned base) {
  double out = 0, f = 1;
  while (n) {
    f /= base;
    out += f * (n % base);
    n /= base;
  }
  return out;
}
void same_mesh(const PolyMesh &a, const PolyMesh &b) {
  KUTTA_REQUIRE(a.n_cells == b.n_cells && a.points.size() == b.points.size() &&
                    a.face_offsets == b.face_offsets &&
                    a.face_nodes == b.face_nodes && a.owner == b.owner &&
                    a.neighbour == b.neighbour,
                "pathline frame topology differs from run mesh");
  for (std::size_t i = 0; i < a.points.size(); ++i)
    KUTTA_REQUIRE(a.points[i].x == b.points[i].x &&
                      a.points[i].y == b.points[i].y &&
                      a.points[i].z == b.points[i].z,
                  "pathline frame coordinates differ from run mesh");
}
// Only two frames are resident. Geometry/adjacency is shared after exact mesh
// identity checks; velocity gradients belong to each frame's saved fields.
class Series {
  const Volume &base_;
  const Trace &spec_;
  std::span<const Volume *const> supplied_;
  struct Cached {
    std::size_t id = SIZE_MAX;
    std::unique_ptr<view::Snapshot> snapshot;
    std::unique_ptr<Volume> volume;
  };
  mutable std::array<Cached, 2> cache_;
  const Volume &frame(std::size_t id) const {
    if (!supplied_.empty())
      return *supplied_[id];
    for (const auto &c : cache_)
      if (c.id == id)
        return *c.volume;
    std::size_t slot = 0;
    if (cache_[0].id != SIZE_MAX && cache_[0].id + 1 == id)
      slot = 1;
    else if (cache_[1].id != SIZE_MAX &&
             (cache_[1].id + 1 == id || cache_[1].id == id + 1))
      slot = 0;
    auto &c = cache_[slot];
    c.volume.reset();
    KUTTA_REQUIRE(
        std::filesystem::exists(spec_.frames[id].run / "summary.json"),
        "pathline frame is not a completed run");
    c.snapshot = std::make_unique<view::Snapshot>(
        view::read_snapshot(spec_.frames[id].run / "view-input.bin"));
    same_mesh(base_.snapshot().mesh, c.snapshot->mesh);
    c.volume = std::make_unique<Volume>(*c.snapshot, &base_);
    c.id = id;
    return *c.volume;
  }

public:
  Series(const Volume &v, const Trace &t, std::span<const Volume *const> frames)
      : base_(v), spec_(t), supplied_(frames) {
    KUTTA_REQUIRE(frames.empty() || frames.size() == t.frames.size(),
                  "pathline frame count mismatch");
    if (!frames.empty())
      for (const auto *f : frames)
        same_mesh(v.snapshot().mesh, f->snapshot().mesh);
  }
  Sample value(std::string_view field, Index cell, Vec3d p, double time) const {
    if (spec_.frames.empty())
      return base_.sample(field, cell, p);
    const auto end =
        std::upper_bound(spec_.frames.begin(), spec_.frames.end(), time,
                         [](double t, const Frame &f) { return t < f.time; });
    const std::size_t hi = end == spec_.frames.end()
                               ? spec_.frames.size() - 1
                               : std::size_t(end - spec_.frames.begin()),
                      lo = hi == 0 ? 0 : hi - 1;
    const auto a = frame(lo).sample(field, cell, p);
    if (lo == hi)
      return a;
    const auto b = frame(hi).sample(field, cell, p);
    const double f =
        std::clamp((time - spec_.frames[lo].time) /
                       (spec_.frames[hi].time - spec_.frames[lo].time),
                   0.0, 1.0);
    auto out = a;
    for (std::uint32_t i = 0; i < out.components; ++i)
      out.value[i] = (1 - f) * a.value[i] + f * b.value[i];
    return out;
  }
};
struct Trial {
  Vec3d point{};
  Index cell = -1;
  std::optional<Walk> boundary;
};
Trial rk4(const Volume &volume, const Series &series, Vec3d p, Index cell,
          double time, double dt) {
  std::optional<Walk> boundary;
  const auto velocity = [&](Vec3d x, double t) {
    const auto walked = volume.walk(p, cell, x);
    if (walked.boundary) {
      boundary = walked;
      return Vec3d{};
    }
    if (walked.cell < 0)
      return Vec3d{NAN, NAN, NAN};
    return vec(series.value("velocity", walked.cell, x, t));
  };
  const auto a = velocity(p, time),
             b = velocity(p + a * (dt * .5), time + dt * .5),
             c = velocity(p + b * (dt * .5), time + dt * .5),
             d = velocity(p + c * dt, time + dt);
  if (boundary)
    return {p, cell, boundary};
  const auto result = p + (a + b * 2 + c * 2 + d) * (dt / 6);
  if (!finite(result))
    return {{NAN, NAN, NAN}, -1, {}};
  const auto walked = volume.walk(p, cell, result);
  return {walked.point, walked.cell,
          walked.boundary ? std::optional<Walk>(walked) : std::nullopt};
}
} // namespace
std::vector<Vec3d> seeds(const Volume &volume, const Seeds &specification) {
  std::vector<Vec3d> out;
  std::visit(
      [&](const auto &s) {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T, Rake>) {
          for (int i = 0; i < s.count; ++i)
            out.push_back(s.start +
                          (s.end - s.start) *
                              (s.count == 1 ? .5 : double(i) / (s.count - 1)));
        } else if constexpr (std::is_same_v<T, PlaneGrid>) {
          for (int j = 0; j < s.nv; ++j)
            for (int i = 0; i < s.nu; ++i)
              out.push_back(s.origin +
                            s.u * (s.nu == 1 ? .5 : double(i) / (s.nu - 1)) +
                            s.v * (s.nv == 1 ? .5 : double(j) / (s.nv - 1)));
        } else {
          struct Triangle {
            Vec3d a, b, c, normal;
            double cumulative = 0;
          };
          std::vector<Triangle> triangles;
          double area = 0;
          const auto &m = volume.snapshot().mesh;
          for (const auto &name : s.zones) {
            const auto &zone = patch(m, name);
            for (Index f = zone.start; f < zone.start + zone.size; ++f) {
              std::vector<Vec3d> p;
              for (auto n : nodes(m, f))
                p.push_back(m.points[n]);
              const auto indices = volume.face_triangles(f);
              for (std::size_t j = 0; j < indices.size(); j += 3) {
                const auto a = p[indices[j]], b = p[indices[j + 1]],
                           c = p[indices[j + 2]], normal = cross(b - a, c - a);
                const auto size = norm(normal) * .5;
                if (size > 0) {
                  area += size;
                  triangles.push_back({a, b, c, normal / norm(normal), area});
                }
              }
            }
          }
          KUTTA_REQUIRE(area > 0, "surface seed zones have no area");
          std::mt19937_64 rng(s.random_seed);
          const auto random = [&]() {
            return std::generate_canonical<double, 53>(rng);
          };
          for (int i = 0; i < s.count; ++i) {
            const double target =
                (s.random ? random() : (i + .5) / s.count) * area;
            const auto &t = *std::lower_bound(
                triangles.begin(), triangles.end(), target,
                [](const auto &tri, double x) { return tri.cumulative < x; });
            const double r = std::sqrt(s.random ? random() : radical(i + 1, 2)),
                         q = s.random ? random() : radical(i + 1, 3);
            out.push_back(t.a * (1 - r) + t.b * (r * (1 - q)) + t.c * (r * q) -
                          t.normal * s.offset);
          }
        }
      },
      specification);
  return out;
}
Curves trace(const Volume &volume, const Trace &spec,
             std::span<const Volume *const> frames) {
  Curves out;
  auto &g = out.geometry;
  g.line_offsets.push_back(0);
  for (const auto &field : spec.fields) {
    const auto info = volume.info(field);
    KUTTA_REQUIRE(!info.boundary_only,
                  "trace field is boundary-only: " + field);
    g.fields.push_back({field, info.unit, info.components, {}});
  }
  g.fields.push_back({"arc_length", "m", 1, {}});
  g.fields.push_back({"time", "s", 1, {}});
  const Series series(volume, spec, frames);
  const auto points = seeds(volume, spec.seeds);
  const double tolerance =
      spec.absolute_tolerance > 0
          ? spec.absolute_tolerance
          : std::max(volume.scale() * 1e-8, volume.tolerance());
  for (std::size_t seed = 0; seed < points.size(); ++seed)
    for (int direction : {-1, 1})
      if (spec.direction == 0 || spec.direction == direction) {
        Vec3d p = points[seed];
        Index cell = volume.locate(p);
        double length = 0, time = spec.start_time;
        int steps = 0, accepted_steps = 0;
        std::string reason = "max_steps";
        const auto append = [&](Vec3d x, Index c, double at, double distance) {
          const auto id = std::uint32_t(g.positions.size() / 3);
          KUTTA_REQUIRE(g.positions.size() / 3 < UINT32_MAX,
                        "too many streamline vertices");
          g.positions.insert(g.positions.end(),
                             {float(x.x), float(x.y), float(x.z)});
          g.indices.push_back(id);
          for (std::size_t i = 0; i < spec.fields.size(); ++i) {
            Sample s;
            s.components = g.fields[i].components;
            s.value = {NAN, NAN, NAN};
            if (c >= 0)
              s = series.value(spec.fields[i], c, x, at);
            for (std::uint32_t j = 0; j < s.components; ++j)
              g.fields[i].values.push_back(float(s.value[j]));
          }
          g.fields[g.fields.size() - 2].values.push_back(float(distance));
          g.fields.back().values.push_back(float(at));
        };
        append(p, cell, time, length);
        if (cell < 0)
          reason = "outside";
        else {
          double dt = 0;
          for (; steps < spec.max_steps; ++steps) {
            if (spec.max_length - length <= tolerance) {
              reason = "max_length";
              break;
            }
            const auto u = vec(series.value("velocity", cell, p, time));
            const double speed = norm(u);
            if (!std::isfinite(speed)) {
              reason = "missing_velocity";
              break;
            }
            if (speed <= 1e-14 && spec.frames.empty()) {
              reason = "stagnation";
              break;
            }
            const double distance =
                spec.step_length > 0
                    ? spec.step_length
                    : std::cbrt(volume.geometry().cell_volume[cell]) * .4;
            const double limit = std::min(
                {speed > 0 ? distance / speed : INFINITY,
                 speed > 0 ? (spec.max_length - length) / speed : INFINITY,
                 spec.max_time - std::abs(time - spec.start_time)});
            if (limit <= 0) {
              reason = length >= spec.max_length ? "max_length" : "max_time";
              break;
            }
            double available = INFINITY;
            if (!spec.frames.empty()) {
              available = direction > 0 ? spec.frames.back().time - time
                                        : time - spec.frames.front().time;
              if (available <= 1e-14) {
                reason = "time_range";
                break;
              }
              const double epsilon = 1e-13 * std::max(1.0, std::abs(time));
              for (const auto &frame : spec.frames) {
                const double distance_time = direction * (frame.time - time);
                if (distance_time > epsilon)
                  available = std::min(available, distance_time);
              }
            }
            dt = std::min({dt > 0 ? dt : limit, limit, available});
            bool accepted = false, terminal = false;
            Trial result;
            double error = 0;
            for (int retry = 0; retry < 48; ++retry) {
              const auto full =
                             rk4(volume, series, p, cell, time, direction * dt),
                         half = rk4(volume, series, p, cell, time,
                                    direction * dt * .5);
              if (full.boundary || half.boundary) {
                const auto hit =
                    half.boundary ? *half.boundary : *full.boundary;
                if (norm(hit.point - p) <= tolerance * 2 ||
                    dt <= 1e-14 * std::max(1.0, std::abs(time))) {
                  result = {hit.point, hit.cell, {}};
                  accepted = true;
                  terminal = true;
                  break;
                }
                dt *= .5;
                continue;
              }
              if (full.cell < 0 || half.cell < 0) {
                dt *= .5;
                continue;
              }
              const auto end =
                  rk4(volume, series, half.point, half.cell,
                      time + direction * dt * .5, direction * dt * .5);
              if (end.boundary) {
                if (norm(end.boundary->point - p) <= tolerance * 2) {
                  result = {end.boundary->point, end.boundary->cell, {}};
                  accepted = true;
                  terminal = true;
                  break;
                }
                dt *= .5;
                continue;
              }
              if (end.cell < 0) {
                dt *= .5;
                continue;
              }
              // Speed can change within a step, especially when an unsteady
              // particle starts at rest. Retry the integration instead of
              // clipping its position and assigning an incorrect flight time.
              const double travel = norm(end.point - p);
              const double max_segment =
                  std::min(distance, spec.max_length - length);
              if (travel > max_segment + tolerance) {
                dt *= std::clamp(max_segment / travel, .1, .9);
                continue;
              }
              error = norm(end.point - full.point) / 15;
              const double allowed =
                  tolerance + spec.relative_tolerance * norm(end.point - p);
              if (error <= allowed) {
                result = end;
                accepted = true;
                break;
              }
              dt *= std::clamp(.9 * std::pow(allowed / error, .2), .1, .5);
            }
            if (!accepted) {
              reason = "step_underflow";
              break;
            }
            double segment = norm(result.point - p);
            if (segment > spec.max_length - length) {
              const double f = (spec.max_length - length) / segment;
              result.point = p + (result.point - p) * f;
              dt *= f;
              segment = spec.max_length - length;
              terminal = true;
              reason = "max_length";
            }
            if (segment <= volume.tolerance() * .01 && !terminal &&
                spec.frames.empty()) {
              reason = "step_underflow";
              break;
            }
            length += segment;
            time += direction * dt;
            p = result.point;
            cell = result.cell;
            append(p, cell, time, length);
            ++accepted_steps;
            if (terminal) {
              if (reason != "max_length")
                reason = "boundary";
              break;
            }
            if (length >= spec.max_length * (1 - 1e-12)) {
              reason = "max_length";
              break;
            }
            if (std::abs(time - spec.start_time) >=
                spec.max_time * (1 - 1e-12)) {
              reason = "max_time";
              break;
            }
            const double allowed =
                tolerance + spec.relative_tolerance * segment;
            dt *= error > 0
                      ? std::clamp(.9 * std::pow(allowed / error, .2), .5, 2.0)
                      : 2.0;
          }
        }
        g.line_offsets.push_back(std::uint32_t(g.indices.size()));
        out.lines.emplace_back(Json::Object{{"seed", double(seed)},
                                            {"direction", direction},
                                            {"stop_reason", reason},
                                            {"steps", accepted_steps},
                                            {"length", length},
                                            {"time", time}});
      }
  return out;
}
} // namespace kutta::post
