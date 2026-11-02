#include "core/Geometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace kapp {

std::optional<WheelGeometry> estimateWheel(const MeshData& mesh, int dimension)
{
    if (mesh.polylines || mesh.indices.empty()) return std::nullopt;
    const auto point = [&mesh](std::uint32_t i) {
        return QVector3D(mesh.positions[3 * i], mesh.positions[3 * i + 1], mesh.positions[3 * i + 2]);
    };
    QVector3D axes[3]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    if (dimension == 3) {
        // A zone can contain both tyres on an axle. Its overall narrowest
        // extent is then radial, so determine the axle on a connected part.
        std::vector<std::uint32_t> parent(mesh.vertexCount());
        std::iota(parent.begin(), parent.end(), 0);
        auto root = [&parent](std::uint32_t i) {
            while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; }
            return i;
        };
        for (std::size_t t = 0; t < mesh.indices.size(); t += 3) {
            const auto a = root(mesh.indices[t]);
            parent[root(mesh.indices[t + 1])] = a;
            parent[root(mesh.indices[t + 2])] = a;
        }
        std::vector<double> areas(mesh.vertexCount());
        for (std::size_t t = 0; t < mesh.indices.size(); t += 3)
            areas[root(mesh.indices[t])] += triangleArea(&mesh.positions[3 * mesh.indices[t]],
                &mesh.positions[3 * mesh.indices[t + 1]], &mesh.positions[3 * mesh.indices[t + 2]]);
        const auto largest = std::uint32_t(std::max_element(areas.begin(), areas.end()) - areas.begin());
        if (!(areas[largest] > 0)) return std::nullopt;
        // Area weighting avoids letting a finely meshed patch dictate the
        // principal axes. Subtract a local origin to retain small dimensions.
        const auto origin = point(largest);
        double mean[3]{}, moment[3][3]{};
        for (std::size_t t = 0; t < mesh.indices.size(); t += 3) {
            if (root(mesh.indices[t]) != largest) continue;
            QVector3D p[3]{point(mesh.indices[t]) - origin, point(mesh.indices[t + 1]) - origin, point(mesh.indices[t + 2]) - origin};
            const double weight = triangleArea(&mesh.positions[3 * mesh.indices[t]],
                &mesh.positions[3 * mesh.indices[t + 1]], &mesh.positions[3 * mesh.indices[t + 2]]) / areas[largest];
            for (int i = 0; i < 3; ++i) {
                for (int k = 0; k < 3; ++k) mean[i] += weight * p[k][i] / 3;
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        for (int l = 0; l < 3; ++l)
                            moment[i][j] += weight * double(p[k][i]) * p[l][j] / (k == l ? 6 : 12);
            }
        }
        double cov[3][3], vectors[3][3]{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) cov[i][j] = moment[i][j] - mean[i] * mean[j];
        // Jacobi diagonalization of the symmetric 3x3 surface covariance.
        for (int step = 0; step < 32; ++step) {
            int p = 0, q = 1;
            for (int i = 0; i < 3; ++i)
                for (int j = i + 1; j < 3; ++j)
                    if (std::abs(cov[i][j]) > std::abs(cov[p][q])) { p = i; q = j; }
            if (std::abs(cov[p][q]) <= 1e-12 * std::max({std::abs(cov[0][0]), std::abs(cov[1][1]), std::abs(cov[2][2])})) break;
            const double angle = 0.5 * std::atan2(2 * cov[p][q], cov[q][q] - cov[p][p]);
            const double c = std::cos(angle), s = std::sin(angle);
            for (int i = 0; i < 3; ++i) {
                const double a = cov[i][p], b = cov[i][q];
                cov[i][p] = c * a - s * b; cov[i][q] = s * a + c * b;
                const double x = vectors[i][p], y = vectors[i][q];
                vectors[i][p] = c * x - s * y; vectors[i][q] = s * x + c * y;
            }
            for (int j = 0; j < 3; ++j) {
                const double a = cov[p][j], b = cov[q][j];
                cov[p][j] = c * a - s * b; cov[q][j] = s * a + c * b;
            }
        }
        int thin = 0;
        for (int i = 1; i < 3; ++i) if (cov[i][i] < cov[thin][thin]) thin = i;
        for (int i = 0; i < 3; ++i) {
            const int column = (thin + 1 + i) % 3;  // axle last
            axes[i] = QVector3D(float(vectors[0][column]), float(vectors[1][column]), float(vectors[2][column])).normalized();
        }
    }
    const double infinity = std::numeric_limits<double>::infinity();
    double lo[3]{infinity, infinity, infinity}, hi[3]{-infinity, -infinity, -infinity};
    for (auto index : mesh.indices) {
        const auto p = point(index);
        for (int i = 0; i < 3; ++i) {
            const double v = QVector3D::dotProduct(p, axes[i]);
            lo[i] = std::min(lo[i], v); hi[i] = std::max(hi[i], v);
        }
    }
    WheelGeometry wheel;
    wheel.radius = std::max(hi[0] - lo[0], hi[1] - lo[1]) * 0.5;
    if (!(wheel.radius > 0) || !std::isfinite(wheel.radius)) return std::nullopt;
    // Planar meshes are exported with an artificial unit-depth extrusion;
    // its midpoint and thickness are not the original wheel geometry.
    QVector3D centre;
    for (int i = 0; i < (dimension == 2 ? 2 : 3); ++i) centre += axes[i] * float((lo[i] + hi[i]) * 0.5);
    for (int i = 0; i < 3; ++i) { wheel.centre[i] = centre[i]; wheel.axis[i] = axes[2][i]; }
    wheel.axis = wheel.rollingAxis(dimension);
    // The export uses float32 positions; showing more than six significant
    // digits suggests precision the surface cannot supply. Remove float noise
    // around zero too, so an aligned axle reads (0, -1, 0).
    auto rounded = [](double x) {
        if (x == 0) return 0.0;
        const double scale = std::pow(10.0, 5 - std::floor(std::log10(std::abs(x))));
        return std::round(x * scale) / scale;
    };
    const double noise = wheel.radius * 1e-7;
    wheel.radius = rounded(wheel.radius);
    for (int i = 0; i < 3; ++i) {
        wheel.centre[i] = std::abs(wheel.centre[i]) < noise ? 0 : rounded(wheel.centre[i]);
        wheel.axis[i] = std::abs(wheel.axis[i]) < 1e-7 ? 0 : rounded(wheel.axis[i]);
    }
    return wheel;
}

} // namespace kapp
