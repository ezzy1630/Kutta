#pragma once

#include "core/ViewFormat.h"

#include <QJsonObject>
#include <QMarginsF>
#include <QMatrix4x4>
#include <QPointF>
#include <QQuaternion>
#include <QSizeF>
#include <QString>
#include <QVector3D>

namespace kapp {

// Turntable camera about a target. Z is up (Fluent cases for cars use Z up;
// the up axis is configurable for cases that do not).
struct Camera {
    QVector3D target{0, 0, 0};
    float distance = 5;
    float yaw = -35;     // degrees about the up axis
    float pitch = 25;    // degrees above the horizon
    float fovY = 30;     // degrees; 0 means orthographic
    QVector3D up{0, 0, 1};

    QVector3D forward() const;  // from eye toward target
    QVector3D eye() const { return target - forward() * distance; }
    QMatrix4x4 view() const;
    QMatrix4x4 projection(float aspect, float nearZ, float farZ) const;
    // Frame bounds in a viewport of the given aspect.
    void fit(const Bounds& b, float aspect, QMarginsF reserve = {});
    // Ray through a point in normalized device coordinates (-1..1, y up).
    void ray(QPointF ndc, float aspect, QVector3D* origin, QVector3D* dir) const;

    void orbit(float dYaw, float dPitch);
    void pan(QPointF ndcDelta, float aspect);
    // Zoom by factor (<1 closer) toward the point under the cursor on the focal plane.
    void zoomAt(QPointF ndc, float aspect, float factor);

    enum class Preset { Front, Back, Left, Right, Top, Bottom, Isometric };
    void preset(Preset p);
    static QString presetName(Preset p);

    QJsonObject toJson() const;
    static Camera fromJson(const QJsonObject& o);
    bool operator==(const Camera&) const = default;
};

} // namespace kapp
