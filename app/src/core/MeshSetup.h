#pragma once

// What `kutta mesh` is asked to build (mesh.json, v1), the wind tunnel it
// implies around the geometry, and a planning estimate of its size. Lengths
// are metres here; mesh.json is written in the project's unit.

#include "core/Stl.h"
#include "core/ViewFormat.h"

#include <QJsonObject>
#include <QString>
#include <QStringList>

#include <optional>
#include <vector>

namespace kapp {

struct RefinementBox {
    QString name;
    QVector3D min, max;  // m
    double size = 0.02;  // target cell size inside, m
};

struct TunnelSetup {
    // Distances in body lengths (the geometry's extent along the flow).
    double upstream = 3, downstream = 6, side = 2, top = 2;
    bool ground = true;     // the floor is a moving ground at the geometry's lowest point
    QString symmetry;       // "" or "y": half model, symmetry plane at y = 0
};

struct MeshSetup {
    LengthUnit units = LengthUnit::Metre;  // the STL files' unit, and mesh.json's
    QString flow = "+x", up = "+z";
    TunnelSetup tunnel;
    double baseSize = 0.5;                 // largest (far-field) cell, m
    double minSize = 0.002, maxSize = 0.04;  // surface cells, m
    double curvatureDegrees = 18;          // max normal angle per surface cell
    int proximityCells = 2;                // cells across a gap
    std::vector<RefinementBox> boxes;
    int layers = 5;
    double firstHeight = 0.0005;           // m
    double growth = 1.2;
    QStringList layerZones{"*"};
    qint64 maxCells = 20'000'000;

    QJsonObject toJson() const;  // mesh.json v1, lengths in `units`
    static MeshSetup fromJson(const QJsonObject& o);
    // Problems that would make the mesher refuse the setup; empty when valid.
    QStringList problems() const;
    bool operator==(const MeshSetup& o) const { return toJson() == o.toJson(); }
};

// Settings that suit the geometry: sizes from its length, a wake box behind
// it and a box around it.
MeshSetup defaultMeshSetup(const Bounds& geometry, LengthUnit units);

// Axis index (0 x, 1 y, 2 z) and sign of "+x", "-z", ...
std::pair<int, int> parseAxis(const QString& axis);

// The wind tunnel around the geometry, in metres.
Bounds tunnelBounds(const MeshSetup& s, const Bounds& geometry);
// Names of the tunnel's faces as `kutta mesh` names them (inlet, outlet,
// ground, sides, top, symmetry), with their roles, for display.
struct TunnelFace { QString name, role; int axis = 0; bool atMax = false; };
std::vector<TunnelFace> tunnelFaces(const MeshSetup& s);

// Rough size of the mesh before it is built. A planning figure: the surface
// band of an octree graded 2:1 from the surface size, prism layers on the
// surface, refinement boxes and the far field.
struct CellEstimate {
    double surfaceFaces = 0, octreeCells = 0, layerCells = 0, boxCells = 0, farFieldCells = 0;
    double total() const { return octreeCells + layerCells + boxCells + farFieldCells; }
    double solverBytes() const;  // device memory to solve it
    bool overLimit = false;      // above max_cells
};
CellEstimate estimateCells(const MeshSetup& s, const Bounds& geometry, double surfaceArea);
// Device bytes per cell the planning estimate assumes for a steady SST solve.
constexpr double kPlanningBytesPerCell = 1500;

// First cell height for a target y+ at a speed, flat-plate correlation over
// the body length (Schlichting skin friction), in metres.
double firstHeightForYPlus(double yPlus, double speed, double length, double density = 1.225, double viscosity = 1.7894e-5);
double yPlusForFirstHeight(double height, double speed, double length, double density = 1.225, double viscosity = 1.7894e-5);

} // namespace kapp
