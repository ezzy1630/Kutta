#include "kutta/mesh/detail.hpp"
#include "post.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

namespace kutta::post {
using namespace detail;
namespace {
struct Point {
  Vec3d position;
  std::vector<Sample> fields;
};
void fields(view::Geometry &g, const Volume &volume,
            const std::vector<std::string> &names) {
  for (const auto &name : names) {
    const auto f = volume.info(name);
    g.fields.push_back({name, f.unit, f.components, {}});
  }
  g.fields.push_back({"arc_length", "m", 1, {}});
  g.fields.push_back({"valid", "1", 1, {}});
}
void append(view::Geometry &g, const Point &p, double arc, bool valid) {
  g.positions.insert(
      g.positions.end(),
      {float(p.position.x), float(p.position.y), float(p.position.z)});
  for (std::size_t i = 0; i < p.fields.size(); ++i)
    for (std::uint32_t j = 0; j < p.fields[i].components; ++j)
      g.fields[i].values.push_back(float(p.fields[i].value[j]));
  g.fields[g.fields.size() - 2].values.push_back(float(arc));
  g.fields.back().values.push_back(valid ? 1 : 0);
}
Curves line_probe(const Volume &volume, const LineProbe &spec) {
  Curves out;
  auto &g = out.geometry;
  fields(g, volume, spec.fields);
  for (const auto &name : spec.fields)
    KUTTA_REQUIRE(!volume.info(name).boundary_only,
                  "line probe field is boundary-only: " + name);
  g.line_offsets.push_back(0);
  bool active = false;
  Index hint = -1;
  Vec3d previous;
  for (int i = 0; i < spec.line.count; ++i) {
    const double t =
        spec.line.count == 1 ? .5 : double(i) / (spec.line.count - 1);
    const auto p = spec.line.start + (spec.line.end - spec.line.start) * t;
    Index cell = -1;
    if (hint >= 0) {
      const auto walked = volume.walk(previous, hint, p);
      if (!walked.boundary)
        cell = walked.cell;
    }
    if (cell < 0)
      cell = volume.locate(p);
    Point point{p, {}};
    for (const auto &name : spec.fields) {
      Sample sample;
      sample.components = volume.info(name).components;
      sample.value = {NAN, NAN, NAN};
      if (cell >= 0)
        sample = volume.sample(name, cell, p);
      point.fields.push_back(sample);
    }
    append(g, point, norm(p - spec.line.start), cell >= 0);
    if (cell >= 0) {
      g.indices.push_back(std::uint32_t(i));
      active = true;
    } else if (active) {
      g.line_offsets.push_back(std::uint32_t(g.indices.size()));
      active = false;
    }
    hint = cell;
    previous = p;
  }
  if (active)
    g.line_offsets.push_back(std::uint32_t(g.indices.size()));
  for (std::size_t i = 1; i < g.line_offsets.size(); ++i)
    out.lines.emplace_back(Json::Object{
        {"closed", false},
        {"samples", int(g.line_offsets[i] - g.line_offsets[i - 1])}});
  return out;
}
Curves section_probe(const Volume &volume, const SectionProbe &spec) {
  Curves out;
  auto &g = out.geometry;
  fields(g, volume, spec.fields);
  g.line_offsets.push_back(0);
  const auto &m = volume.snapshot().mesh;
  for (const auto &zone_name : spec.zones) {
    const auto zone =
        std::find_if(m.patches.begin(), m.patches.end(),
                     [&](const auto &p) { return p.name == zone_name; });
    KUTTA_REQUIRE(zone != m.patches.end(),
                  "unknown section zone: " + zone_name);
    std::map<Index, std::size_t> node_ids;
    std::vector<Point> original;
    std::vector<double> weights;
    std::vector<std::vector<bool>> missing;
    for (Index f = zone->start; f < zone->start + zone->size; ++f) {
      const auto center = volume.geometry().face_centroid[f];
      const double area = norm(volume.geometry().face_area[f]);
      std::vector<Sample> samples;
      for (const auto &name : spec.fields)
        samples.push_back(volume.boundary(name, f, center));
      for (auto n : nodes(m, f)) {
        auto [it, added] = node_ids.emplace(n, original.size());
        if (added) {
          Point p{m.points[n], {}};
          for (const auto &name : spec.fields) {
            Sample s;
            s.components = volume.info(name).components;
            p.fields.push_back(s);
          }
          original.push_back(std::move(p));
          weights.push_back(0);
          missing.emplace_back(spec.fields.size(), false);
        }
        const auto id = it->second;
        weights[id] += area;
        for (std::size_t k = 0; k < samples.size(); ++k)
          for (std::uint32_t c = 0; c < samples[k].components; ++c) {
            original[id].fields[k].value[c] += area * samples[k].value[c];
            missing[id][k] =
                missing[id][k] || !std::isfinite(samples[k].value[c]);
          }
      }
    }
    for (std::size_t i = 0; i < original.size(); ++i)
      for (std::size_t k = 0; k < spec.fields.size(); ++k)
        for (std::uint32_t c = 0; c < original[i].fields[k].components; ++c)
          original[i].fields[k].value[c] =
              missing[i][k] ? NAN : original[i].fields[k].value[c] / weights[i];
    std::vector<Point> cuts;
    std::map<std::array<std::int64_t, 3>, std::vector<std::uint32_t>> buckets;
    const auto origin = m.points.front();
    const double epsilon = volume.tolerance();
    const auto weld = [&](Point p) {
      const auto delta = (p.position - origin) / epsilon;
      const std::array<std::int64_t, 3> key{
          std::llround(delta.x), std::llround(delta.y), std::llround(delta.z)};
      for (int x = -1; x <= 1; ++x)
        for (int y = -1; y <= 1; ++y)
          for (int z = -1; z <= 1; ++z) {
            const auto it = buckets.find({key[0] + x, key[1] + y, key[2] + z});
            if (it != buckets.end())
              for (auto id : it->second)
                if (norm(cuts[id].position - p.position) <= epsilon)
                  return id;
          }
      const auto id = std::uint32_t(cuts.size());
      cuts.push_back(std::move(p));
      buckets[key].push_back(id);
      return id;
    };
    struct Edge {
      int coplanar = 0;
      bool ordinary = false;
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>, Edge> edges;
    const auto edge = [&](std::uint32_t a, std::uint32_t b, bool coplanar) {
      if (a == b)
        return;
      const auto key = std::minmax(a, b);
      auto &e = edges[{key.first, key.second}];
      if (coplanar)
        ++e.coplanar;
      else
        e.ordinary = true;
    };
    const auto interpolate = [&](const Point &a, const Point &b, double t) {
      Point p{a.position + (b.position - a.position) * t, a.fields};
      for (std::size_t i = 0; i < p.fields.size(); ++i)
        for (std::uint32_t j = 0; j < p.fields[i].components; ++j)
          p.fields[i].value[j] =
              (1 - t) * a.fields[i].value[j] + t * b.fields[i].value[j];
      return p;
    };
    for (Index f = zone->start; f < zone->start + zone->size; ++f) {
      const auto ns = nodes(m, f);
      std::vector<Vec3d> p;
      for (auto n : ns)
        p.push_back(m.points[n]);
      const auto triangles = view::triangulate(p);
      for (std::size_t t = 0; t < triangles.size(); t += 3) {
        const Point *v[3] = {&original[node_ids.at(ns[triangles[t]])],
                             &original[node_ids.at(ns[triangles[t + 1]])],
                             &original[node_ids.at(ns[triangles[t + 2]])]};
        double d[3];
        for (int i = 0; i < 3; ++i)
          d[i] = dot(v[i]->position - spec.plane.origin, spec.plane.normal);
        if (std::abs(d[0]) <= epsilon && std::abs(d[1]) <= epsilon &&
            std::abs(d[2]) <= epsilon) {
          const auto a = weld(*v[0]), b = weld(*v[1]), c = weld(*v[2]);
          edge(a, b, true);
          edge(b, c, true);
          edge(c, a, true);
          continue;
        }
        std::vector<std::uint32_t> hits;
        const auto add = [&](Point point) {
          const auto id = weld(std::move(point));
          if (std::find(hits.begin(), hits.end(), id) == hits.end())
            hits.push_back(id);
        };
        for (int i = 0; i < 3; ++i) {
          const auto j = (i + 1) % 3;
          if (std::abs(d[i]) <= epsilon)
            add(*v[i]);
          if ((d[i] > epsilon && d[j] < -epsilon) ||
              (d[i] < -epsilon && d[j] > epsilon))
            add(interpolate(*v[i], *v[j], d[i] / (d[i] - d[j])));
        }
        if (hits.size() == 2)
          edge(hits[0], hits[1], false);
      }
    }
    std::vector<std::array<std::uint32_t, 2>> segments;
    for (const auto &[key, e] : edges)
      if (e.ordinary || e.coplanar % 2 == 1)
        segments.push_back({key.first, key.second});
    std::vector<std::vector<std::size_t>> adjacency(cuts.size());
    for (std::size_t i = 0; i < segments.size(); ++i) {
      adjacency[segments[i][0]].push_back(i);
      adjacency[segments[i][1]].push_back(i);
    }
    std::vector<bool> used(segments.size());
    const auto curve = [&](std::uint32_t start, std::size_t first) {
      std::vector<std::uint32_t> ids{start};
      auto current = start;
      auto next = first;
      for (;;) {
        used[next] = true;
        current = segments[next][0] == current ? segments[next][1]
                                               : segments[next][0];
        ids.push_back(current);
        if (current == start || adjacency[current].size() != 2)
          break;
        const auto &a = adjacency[current];
        next = used[a[0]] ? a[1] : a[0];
        if (used[next])
          break;
      }
      const bool closed = ids.back() == start;
      if (closed) {
        ids.pop_back();
        Vec3d area;
        for (std::size_t i = 0; i < ids.size(); ++i)
          area = area + cross(cuts[ids[i]].position - spec.plane.origin,
                              cuts[ids[(i + 1) % ids.size()]].position -
                                  spec.plane.origin);
        if (dot(area, spec.plane.normal) < 0)
          std::reverse(ids.begin(), ids.end());
        const auto min =
            std::min_element(ids.begin(), ids.end(), [&](auto a, auto b) {
              const auto p = cuts[a].position, q = cuts[b].position;
              return std::tie(p.x, p.y, p.z) < std::tie(q.x, q.y, q.z);
            });
        std::rotate(ids.begin(), min, ids.end());
        ids.push_back(ids.front());
      }
      double arc = 0;
      for (std::size_t i = 0; i < ids.size(); ++i) {
        if (i)
          arc += norm(cuts[ids[i]].position - cuts[ids[i - 1]].position);
        const auto vertex = std::uint32_t(g.positions.size() / 3);
        append(g, cuts[ids[i]], arc, true);
        g.indices.push_back(vertex);
      }
      g.line_offsets.push_back(std::uint32_t(g.indices.size()));
      out.lines.emplace_back(Json::Object{{"zone", zone_name},
                                          {"closed", closed},
                                          {"length", arc},
                                          {"samples", double(ids.size())}});
    };
    for (std::size_t p = 0; p < adjacency.size(); ++p)
      if (adjacency[p].size() != 2)
        for (auto edge : adjacency[p])
          if (!used[edge])
            curve(std::uint32_t(p), edge);
    for (std::size_t e = 0; e < segments.size(); ++e)
      if (!used[e])
        curve(segments[e][0], e);
  }
  return out;
}
} // namespace
Curves probe(const Volume &volume, const Probe &spec) {
  return std::visit(
      [&](const auto &p) -> Curves {
        if constexpr (std::is_same_v<std::decay_t<decltype(p)>, LineProbe>)
          return line_probe(volume, p);
        else
          return section_probe(volume, p);
      },
      spec);
}
} // namespace kutta::post
