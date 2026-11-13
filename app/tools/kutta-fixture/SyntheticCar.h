#pragma once

// A SYNTHETIC FSAE car for demonstrating and testing the app. Geometry is
// procedural; fields come from a potential-flow + vortex + wake model with
// invented constants. Nothing here is a CFD result or evidence about any car:
// every file the fixture CLI writes from it carries "synthetic": true.

#include "core/ViewFormat.h"

#include <QJsonObject>
#include <QMap>
#include <QString>
#include <QVector3D>

#include <vector>

namespace fixture {

// A part of the car moved and scaled about its own centre, as a swapped STL
// solid of a different size would be (fixture meshes from `mesh`).
struct SolidTransform {
    QVector3D pivot, offset;
    float scale = 1;
    QVector3D apply(const QVector3D& p) const { return pivot + (p - pivot) * scale + offset; }
};

struct CarParams {
    double speed = 16.0;          // inlet velocity magnitude, m/s
    double density = 1.225;       // kg/m^3
    double viscosity = 1.7894e-5; // Pa s
    double refArea = 1.0, refDensity = 1.225, refVelocity = 16.0, refPressure = 0.0, refLength = 1.53;
    double convergence = 1.0;     // 0..1: live snapshots blend from the initial field
    int detail = 1;               // tessellation multiplier (~d^2 triangles)
    QMap<QString, SolidTransform> transforms;  // by wall zone
};

// A zone's surface alone (no fields) at a tessellation resolution (1 = the
// views'), e.g. coarse for sample STL files.
kapp::MeshData carZoneSurface(const QString& zone, double resolution);

// The car's wall zones (the solids an STL of it holds) and their bounds.
QStringList carWallZones();
kapp::Bounds carZoneBounds(const QString& zone);
// The transform that maps a zone's own bounds onto `target`: uniform scale
// from the diagonals, then the centre shift.
SolidTransform fitTransform(const kapp::Bounds& zone, const kapp::Bounds& target);

struct CarZone { QString name, type; int id; bool hasGeometry; };
const std::vector<CarZone>& carZones();

// One boundary zone: geometry plus fields (Kutta, and Fluent/delta when `reference`).
kapp::MeshData carZoneMesh(const QString& zone, const CarParams& p, bool reference);
kapp::MeshData carSlice(const QVector3D& origin, const QVector3D& normal, const CarParams& p, bool reference);
kapp::MeshData carStreamlines(const QVector3D& from, const QVector3D& to, int count, double maxLength,
                              const CarParams& p);
kapp::MeshData carIsoSurface(const QString& field, double value, const CarParams& p);

struct ProbeSeries { QString name, unit; std::vector<double> y; };
struct ProbeResult { std::vector<double> distance; std::vector<QVector3D> points; std::vector<ProbeSeries> series; };
ProbeResult carLineProbe(const QVector3D& from, const QVector3D& to, int samples, const QStringList& fields,
                         const CarParams& p);

// Force and moment about the origin (SYNTHETIC: invented fields, not CFD).
struct ZoneForce { QString zone; QVector3D pressure, viscous, pressureMoment, viscousMoment; double area = 0; };
std::vector<ZoneForce> carForces(const CarParams& p);

// Synthetic case summary in the inspect --json shape.
QJsonObject carInspect(const CarParams& p, bool withData, const QJsonArray& rawSettings);

} // namespace fixture
