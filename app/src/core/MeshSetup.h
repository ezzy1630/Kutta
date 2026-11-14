#pragma once

// What `kutta mesh` is asked to build (mesh.json, v1), the wind tunnel it
// implies around the geometry, and a planning estimate of its size. Lengths
// are metres here and in mesh.json; `units` there scales STL coordinates only.

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
    bool ground = true;     // the floor is the ground plane z = 0 (geometry must not go below it)
    QString symmetry;       // "" or "y": half model, symmetry plane at y = 0
};

// Defaults are kutta mesh's own (a setup that omits a key gets these).
struct MeshSetup {
    LengthUnit units = LengthUnit::Metre;  // the STL files' unit (mesh.json lengths are always m)
    QString flow = "+x", up = "+z";        // kutta mesh v1 supports only +x and +z
    TunnelSetup tunnel;
    double baseSize = 0.25;                // largest (far-field) cell, m
    double minSize = 0.006, maxSize = 0.03;  // surface cells, m
    double curvatureDegrees = 20;          // max normal angle per surface cell
    int proximityCells = 3;                // cells across a gap (1 to 20)
    double proximityMinSize = 0;           // m; gaps may refine below minSize to this; 0: the mesher's default (minSize / 4)
    std::vector<RefinementBox> boxes;
    int layers = 4;
    double firstHeight = 0.002;            // m
    double growth = 1.3;
    QStringList layerZones{"*"};
    qint64 maxCells = 10'000'000;

    QJsonObject toJson() const;  // mesh.json v1, lengths in metres
    // `lengthsInUnits`: geometry.json from before the real mesher, whose
    // lengths were written in `units`.
    static MeshSetup fromJson(const QJsonObject& o, bool lengthsInUnits = false);
    // Problems that would make the mesher refuse the setup; empty when valid.
    QStringList problems() const;
    bool operator==(const MeshSetup& o) const { return toJson() == o.toJson(); }
};

// Settings that suit the geometry: sizes from its length, a wake box behind
// it and a box around it.
MeshSetup defaultMeshSetup(const Bounds& geometry, LengthUnit units);

// Axis index (0 x, 1 y, 2 z) and sign of "+x", "-z", ...
std::pair<int, int> parseAxis(const QString& axis);

// The wind tunnel around the geometry, in metres, as kutta mesh builds it:
// distances in body x lengths from the body's bounds, the ground at z = 0
// (else `side` lengths below the body), upper limits rounded out to whole
// base cells.
Bounds tunnelBounds(const MeshSetup& s, const Bounds& geometry);
// How the ground (the floor at z = 0) moves in the case: with the air, as a
// car on the road sees it (new-case's default, moving-ground), or not at all,
// the fixed floor of a wind tunnel such as the Ahmed body's (stationary-ground).
// Neither is in the force reports. Not part of mesh.json: the mesh is the same.
enum class GroundMotion { Moving, Stationary };
QString groundMotionKey(GroundMotion g);         // "moving", "stationary"
GroundMotion groundMotionFrom(const QString& key);  // anything else: Moving
// Names of the tunnel's faces as `kutta mesh` names them (inlet, outlet,
// ground or bottom, top, side-low or symmetry, side-high), with their
// new-case roles; the ground's follows `ground`.
struct TunnelFace { QString name, role; int axis = 0; bool atMax = false; };
std::vector<TunnelFace> tunnelFaces(const MeshSetup& s, GroundMotion ground = GroundMotion::Moving);
// Whether a new-case role is one of the ground roles.
bool isGroundRole(const QString& role);
// Names kutta mesh reserves for the tunnel; a body zone cannot use them.
const QStringList& reservedZoneNames();

// Rough size of the mesh before it is built. A planning figure: the surface
// band of an octree graded 2:1 from the surface size, prism layers on the
// surface, refinement boxes and the far field.
struct CellEstimate {
    double surfaceFaces = 0, octreeCells = 0, layerCells = 0, boxCells = 0, farFieldCells = 0;
    double total() const { return octreeCells + layerCells + boxCells + farFieldCells; }
    double solverBytes() const;  // device memory to solve it (kSolverBytesPerCell)
    bool overLimit = false;      // above max_cells
};
CellEstimate estimateCells(const MeshSetup& s, const Bounds& geometry, double surfaceArea);

// First cell height for a target y+ at a speed, flat-plate correlation over
// the body length (Schlichting skin friction), in metres.
double firstHeightForYPlus(double yPlus, double speed, double length, double density = 1.225, double viscosity = 1.7894e-5);
double yPlusForFirstHeight(double height, double speed, double length, double density = 1.225, double viscosity = 1.7894e-5);

} // namespace kapp
