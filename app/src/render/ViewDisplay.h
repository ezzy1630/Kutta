#pragma once

#include "core/Colormap.h"
#include "core/Fields.h"
#include "render/Camera.h"

#include <QJsonObject>
#include <QStringList>

namespace kapp {

// Everything one viewport shows of its dataset. Plain data: the inspector edits
// it, the viewport renders it, saved views and run history store it.
struct ViewDisplay {
    FieldRef field{"pressure", FieldRole::Kutta, -1};
    bool contours = true;          // false: solid grey shading
    ColormapId colormap = ColormapId::Viridis;
    bool colormapFollowsRole = true;  // deltas switch to the diverging map
    int bands = 0;                 // 0 smooth; otherwise banded with contour lines
    RangeMode rangeMode = RangeMode::Auto;
    double fixedMin = 0, fixedMax = 1;
    bool edges = false;
    bool zoneColors = false;       // surfaces without a field take their zone's colour
    bool lighting = true;
    QStringList visible;           // item ids; empty means the dataset's defaults
    bool visibleInitialized = false;

    // Vectors (velocity on slices and flow boundaries; wall shear on walls).
    bool vectors = false;
    QString vectorItem;            // item id carrying the vector field
    float vectorLength = 0;        // world units; 0 = automatic
    int vectorStride = 4;

    // Pathlines and streamlines: width in logical pixels, arrowheads along
    // them, and a solid colour instead of the viewport's field.
    float ribbonWidth = 3.0f;
    bool lineArrows = false;
    bool lineSolid = false;

    Camera camera;
    bool cameraInitialized = false;

    ColormapId effectiveColormap() const;
    bool symmetricRange() const { return field.role == FieldRole::Delta || field.role == FieldRole::Run; }

    QJsonObject toJson() const;
    static ViewDisplay fromJson(const QJsonObject& o);
};

} // namespace kapp
