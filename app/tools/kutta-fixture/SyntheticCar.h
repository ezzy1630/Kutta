#pragma once

// A SYNTHETIC FSAE car for demonstrating and testing the app. Geometry is
// procedural; fields come from a potential-flow + vortex + wake model with
// invented constants. Nothing here is a CFD result or evidence about any car:
// every file the fixture CLI writes from it carries "synthetic": true.

#include "core/ViewFormat.h"

#include <QJsonObject>
#include <QString>
#include <QVector3D>

#include <vector>

namespace fixture {

struct CarParams {
    double speed = 16.0;          // inlet velocity magnitude, m/s
    double density = 1.225;       // kg/m^3
    double viscosity = 1.7894e-5; // Pa s
    double refArea = 1.0, refDensity = 1.225, refVelocity = 16.0, refPressure = 0.0;
    double convergence = 1.0;     // 0..1: live snapshots blend from the initial field
    int detail = 1;               // tessellation multiplier (~d^2 triangles)
};

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

struct ZoneForce { QString zone; QVector3D pressure, viscous; };
std::vector<ZoneForce> carForces(const CarParams& p);

// Synthetic case summary in the inspect --json shape.
QJsonObject carInspect(const CarParams& p, bool withData, const QJsonArray& rawSettings);

} // namespace fixture
