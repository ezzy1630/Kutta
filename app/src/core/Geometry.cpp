#include "core/Geometry.h"

#include <cmath>
#include <limits>

namespace kapp {

double triangleArea(const float* a, const float* b, const float* c)
{
    const QVector3D ab(b[0] - a[0], b[1] - a[1], b[2] - a[2]), ac(c[0] - a[0], c[1] - a[1], c[2] - a[2]);
    return 0.5 * double(QVector3D::crossProduct(ab, ac).length());
}

namespace {

double sampleValue(const MeshData::Field& f, std::uint32_t v, int component)
{
    if (f.desc.components == 1) return f.values[v];
    if (component >= 0) return f.values[3 * std::size_t(v) + std::size_t(component)];
    const double x = f.values[3 * std::size_t(v)], y = f.values[3 * std::size_t(v) + 1], z = f.values[3 * std::size_t(v) + 2];
    return std::sqrt(x * x + y * y + z * z);
}

} // namespace

std::vector<SectionPoint> planeSection(const MeshData& mesh, const QString& fieldName, int component,
                                       const QVector3D& origin, const QVector3D& normalIn)
{
    std::vector<SectionPoint> out;
    const MeshData::Field* f = mesh.field(fieldName);
    if (!f || mesh.polylines) return out;
    const QVector3D n = normalIn.normalized();
    const auto& P = mesh.positions;
    auto pos = [&](std::uint32_t v) { return QVector3D(P[3 * v], P[3 * v + 1], P[3 * v + 2]); };
    for (std::size_t t = 0; t + 2 < mesh.indices.size(); t += 3) {
        const std::uint32_t idx[3] = {mesh.indices[t], mesh.indices[t + 1], mesh.indices[t + 2]};
        double d[3];
        for (int k = 0; k < 3; ++k) d[k] = double(QVector3D::dotProduct(pos(idx[k]) - origin, n));
        SectionPoint hits[2];
        int nh = 0;
        for (int k = 0; k < 3 && nh < 2; ++k) {
            const int j = (k + 1) % 3;
            // Half-open test so a vertex exactly on the plane is counted once.
            if ((d[k] < 0) == (d[j] < 0)) continue;
            const double s = d[k] / (d[k] - d[j]);
            const QVector3D p = pos(idx[k]) + (pos(idx[j]) - pos(idx[k])) * float(s);
            const double va = sampleValue(*f, idx[k], component), vb = sampleValue(*f, idx[j], component);
            hits[nh++] = {p, va + (vb - va) * s};
        }
        if (nh == 2) {
            out.push_back(hits[0]);
            out.push_back(hits[1]);
        }
    }
    return out;
}

std::optional<SurfaceIntegral> integrate(const MeshData& mesh, const QString& fieldName, IntegralKind kind)
{
    const MeshData::Field* f = fieldName.isEmpty() ? nullptr : mesh.field(fieldName);
    if (kind != IntegralKind::Area && !f) return std::nullopt;
    SurfaceIntegral r;
    double sum = 0, lo = std::numeric_limits<double>::infinity(), hi = -lo;
    const auto& P = mesh.positions;
    for (std::size_t t = 0; t + 2 < mesh.indices.size(); t += 3) {
        const std::uint32_t a = mesh.indices[t], b = mesh.indices[t + 1], c = mesh.indices[t + 2];
        const double area = triangleArea(&P[3 * a], &P[3 * b], &P[3 * c]);
        if (kind == IntegralKind::Area) { r.area += area; continue; }
        const double va = sampleValue(*f, a, -1), vb = sampleValue(*f, b, -1), vc = sampleValue(*f, c, -1);
        if (!std::isfinite(va) || !std::isfinite(vb) || !std::isfinite(vc)) continue;
        const double mean = (va + vb + vc) / 3.0;
        r.area += area;
        sum += mean * area;
        lo = std::min({lo, va, vb, vc});
        hi = std::max({hi, va, vb, vc});
        ++r.samples;
    }
    switch (kind) {
    case IntegralKind::Area: r.value = r.area; break;
    case IntegralKind::Integral: r.value = sum; break;
    case IntegralKind::AreaWeightedAverage:
        if (r.area <= 0) return std::nullopt;
        r.value = sum / r.area;
        break;
    case IntegralKind::Minimum: if (!r.samples) return std::nullopt; r.value = lo; break;
    case IntegralKind::Maximum: if (!r.samples) return std::nullopt; r.value = hi; break;
    }
    return r;
}

} // namespace kapp
