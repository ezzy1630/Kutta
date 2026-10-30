#pragma once

// Geometry on exported view meshes: plane sections (XY plots along a wing
// section) and surface integrals. Used by the fixture CLI and tests; the real
// CLI owns these computations for real runs (docs/app-protocol-requests.md).

#include "core/ViewFormat.h"

#include <QVector3D>

#include <optional>
#include <vector>

namespace kapp {

struct SectionPoint { QVector3D position; double value; };

// Intersect a triangle mesh with a plane; returns segment end points with the
// field interpolated along edges (pairs: [2k], [2k+1] form one segment).
std::vector<SectionPoint> planeSection(const MeshData& mesh, const QString& field, int component,
                                       const QVector3D& origin, const QVector3D& normal);

enum class IntegralKind { Area, Integral, AreaWeightedAverage, Minimum, Maximum };
struct SurfaceIntegral { double area = 0; double value = 0; std::size_t samples = 0; };
// Area-weighted triangle integration of a vertex field (vector fields use magnitude).
std::optional<SurfaceIntegral> integrate(const MeshData& mesh, const QString& field, IntegralKind kind);

double triangleArea(const float* a, const float* b, const float* c);

} // namespace kapp
