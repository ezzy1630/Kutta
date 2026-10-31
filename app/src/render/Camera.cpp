#include "render/Camera.h"

#include <QJsonArray>
#include <QtMath>

#include <algorithm>

namespace kapp {

namespace {
// Basis with `up` as world vertical.
void basis(const QVector3D& up, QVector3D* ex, QVector3D* ey, QVector3D* ez)
{
    *ez = up.normalized();
    const QVector3D helper = std::abs(ez->x()) < 0.9f ? QVector3D(1, 0, 0) : QVector3D(0, 1, 0);
    *ey = QVector3D::crossProduct(*ez, helper).normalized();
    *ex = QVector3D::crossProduct(*ey, *ez);
}
} // namespace

QVector3D Camera::forward() const
{
    QVector3D ex, ey, ez;
    basis(up, &ex, &ey, &ez);
    const float cy = std::cos(qDegreesToRadians(yaw)), sy = std::sin(qDegreesToRadians(yaw));
    const float cp = std::cos(qDegreesToRadians(pitch)), sp = std::sin(qDegreesToRadians(pitch));
    // Eye sits at -forward*distance: yaw 0, pitch 0 looks along +ex.
    return (ex * (cy * cp) + ey * (sy * cp) - ez * sp).normalized();
}

QMatrix4x4 Camera::view() const
{
    QMatrix4x4 m;
    QVector3D f = forward();
    QVector3D u = up.normalized();
    if (std::abs(QVector3D::dotProduct(f, u)) > 0.999f) {
        QVector3D ex, ey, ez;
        basis(up, &ex, &ey, &ez);
        u = QVector3D::crossProduct(ey, f).normalized();
        if (u.isNull()) u = ex;
    }
    m.lookAt(eye(), target, u);
    return m;
}

QMatrix4x4 Camera::projection(float aspect, float nearZ, float farZ) const
{
    QMatrix4x4 p;
    if (fovY <= 0) {
        const float h = distance * std::tan(qDegreesToRadians(15.0f));
        p.ortho(-h * aspect, h * aspect, -h, h, nearZ, farZ);
    } else {
        p.perspective(fovY, aspect, nearZ, farZ);
    }
    return p;
}

void Camera::fit(const Bounds& b, float aspect, QMarginsF reserve)
{
    if (!b.valid) return;
    target = b.center();
    // Fit the box's projected corners, not its bounding sphere: long, flat
    // shapes (a car, a slice) then fill the view instead of floating in it.
    const QVector3D f = forward();
    QVector3D right = QVector3D::crossProduct(f, up);
    if (right.lengthSquared() < 1e-8f) right = QVector3D::crossProduct(f, QVector3D(0, 1, 0));
    right.normalize();
    const QVector3D upv = QVector3D::crossProduct(right, f).normalized();
    const float margin = 0.9f;
    // reserve: fractions of the viewport kept clear for overlays (title, legend).
    const float fx = float(std::clamp(1.0 - reserve.left() - reserve.right(), 0.3, 1.0));
    const float fy = float(std::clamp(1.0 - reserve.top() - reserve.bottom(), 0.3, 1.0));
    const float tanFull = std::tan(qDegreesToRadians((fovY > 0 ? fovY : 30.0f) * 0.5f));
    const float tanY = tanFull * margin * fy;
    const float tanX = tanFull * margin * fx * std::max(aspect, 0.2f);
    float need = 1e-6f;
    for (int i = 0; i < 8; ++i) {
        const QVector3D c((i & 1) ? b.max.x() : b.min.x(), (i & 2) ? b.max.y() : b.min.y(), (i & 4) ? b.max.z() : b.min.z());
        const QVector3D d = c - target;
        const float x = std::abs(QVector3D::dotProduct(d, right));
        const float y = std::abs(QVector3D::dotProduct(d, upv));
        const float z = QVector3D::dotProduct(d, f);  // > 0: beyond the target
        // Orthographic views scale with distance alone; perspective adds depth.
        const float depth = fovY > 0 ? z : 0.0f;
        need = std::max({need, x / tanX - depth, y / tanY - depth});
    }
    distance = std::max(need, b.radius() * 0.05f);
    // Centre the box in the free area: shift the target opposite to the offset.
    const float cx = float(reserve.left() - reserve.right()), cy = float(reserve.bottom() - reserve.top());
    target -= right * (cx * distance * tanFull * std::max(aspect, 0.2f)) + upv * (cy * distance * tanFull);
}

void Camera::fit(const std::vector<QVector3D>& points, float aspect, QMarginsF reserve, float fill)
{
    if (points.empty()) return;
    const QVector3D f = forward();
    QVector3D right = QVector3D::crossProduct(f, up);
    if (right.lengthSquared() < 1e-8f) right = QVector3D::crossProduct(f, QVector3D(0, 1, 0));
    right.normalize();
    const QVector3D upv = QVector3D::crossProduct(right, f).normalized();
    // Extent in the view plane; the target moves to its middle.
    QVector3D c;
    for (const QVector3D& p : points) c += p;
    c /= float(points.size());
    float x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY, z0 = INFINITY, z1 = -INFINITY;
    for (const QVector3D& p : points) {
        const QVector3D d = p - c;
        const float x = QVector3D::dotProduct(d, right), y = QVector3D::dotProduct(d, upv), z = QVector3D::dotProduct(d, f);
        x0 = std::min(x0, x), x1 = std::max(x1, x), y0 = std::min(y0, y), y1 = std::max(y1, y), z0 = std::min(z0, z), z1 = std::max(z1, z);
    }
    const float mx = (x0 + x1) * 0.5f, my = (y0 + y1) * 0.5f;
    target = c + right * mx + upv * my;
    const float fx = float(std::clamp(1.0 - reserve.left() - reserve.right(), 0.3, 1.0));
    const float fy = float(std::clamp(1.0 - reserve.top() - reserve.bottom(), 0.3, 1.0));
    const float tanFull = std::tan(qDegreesToRadians((fovY > 0 ? fovY : 30.0f) * 0.5f));
    const float tanY = tanFull * fill * fy;
    const float tanX = tanFull * fill * fx * std::max(aspect, 0.05f);
    // Perspective: a point nearer the eye than the target needs more room.
    // Orthographic views scale with distance alone.
    float need = 1e-6f;
    for (const QVector3D& p : points) {
        const QVector3D d = p - target;
        const float x = std::abs(QVector3D::dotProduct(d, right)), y = std::abs(QVector3D::dotProduct(d, upv));
        const float depth = fovY > 0 ? QVector3D::dotProduct(d, f) : 0.0f;
        need = std::max({need, x / tanX - depth, y / tanY - depth});
    }
    const float size = std::max({x1 - x0, y1 - y0, z1 - z0});
    distance = std::max(need, size * 0.02f);
    const float cx = float(reserve.left() - reserve.right()), cy = float(reserve.bottom() - reserve.top());
    target -= right * (cx * distance * tanFull * std::max(aspect, 0.05f)) + upv * (cy * distance * tanFull);
}

void Camera::ray(QPointF ndc, float aspect, QVector3D* origin, QVector3D* dir) const
{
    const QMatrix4x4 inv = (projection(aspect, distance * 0.01f, distance * 100.f) * view()).inverted();
    const QVector3D a = inv.map(QVector3D(float(ndc.x()), float(ndc.y()), -1));
    const QVector3D b = inv.map(QVector3D(float(ndc.x()), float(ndc.y()), 1));
    *origin = a;
    *dir = (b - a).normalized();
}

void Camera::orbit(float dYaw, float dPitch)
{
    yaw = std::fmod(yaw + dYaw, 360.0f);
    pitch = std::clamp(pitch + dPitch, -89.5f, 89.5f);
}

void Camera::pan(QPointF d, float aspect)
{
    const QMatrix4x4 v = view();
    const QVector3D right(v(0, 0), v(0, 1), v(0, 2)), upv(v(1, 0), v(1, 1), v(1, 2));
    const float h = distance * std::tan(qDegreesToRadians((fovY > 0 ? fovY : 30.0f) * 0.5f));
    target -= right * float(d.x()) * h * aspect + upv * float(d.y()) * h;
}

void Camera::zoomAt(QPointF ndc, float aspect, float factor)
{
    QVector3D o, dir;
    ray(ndc, aspect, &o, &dir);
    const QVector3D n = forward();
    const float denom = QVector3D::dotProduct(dir, n);
    QVector3D p = target;
    if (std::abs(denom) > 1e-6f) p = o + dir * (QVector3D::dotProduct(target - o, n) / denom);
    target += (p - target) * (1 - factor);
    distance *= factor;
}

void Camera::planView()
{
    // With y as the vertical, yaw 90 at the horizon looks along -z.
    up = QVector3D(0, 1, 0);
    yaw = 90;
    pitch = 0;
}

void Camera::preset(Preset p)
{
    switch (p) {
    case Preset::Front: yaw = 0; pitch = 0; break;       // looking downstream (+x)
    case Preset::Back: yaw = 180; pitch = 0; break;
    case Preset::Left: yaw = 90; pitch = 0; break;
    case Preset::Right: yaw = -90; pitch = 0; break;
    // Plan views keep +x to the right: cars and most ducts are long in x.
    case Preset::Top: yaw = 90; pitch = 89.5f; break;
    case Preset::Bottom: yaw = 90; pitch = -89.5f; break;
    case Preset::Isometric: yaw = -35; pitch = 25; break;
    }
}

QString Camera::presetName(Preset p)
{
    switch (p) {
    case Preset::Front: return QStringLiteral("Front (+X)");
    case Preset::Back: return QStringLiteral("Back (−X)");
    case Preset::Left: return QStringLiteral("Left (+Y)");
    case Preset::Right: return QStringLiteral("Right (−Y)");
    case Preset::Top: return QStringLiteral("Top (−Z)");
    case Preset::Bottom: return QStringLiteral("Bottom (+Z)");
    case Preset::Isometric: return QStringLiteral("Isometric");
    }
    return {};
}

QJsonObject Camera::toJson() const
{
    return {{"target", QJsonArray{target.x(), target.y(), target.z()}},
            {"distance", distance},
            {"yaw", yaw},
            {"pitch", pitch},
            {"fov", fovY},
            {"up", QJsonArray{up.x(), up.y(), up.z()}}};
}

Camera Camera::fromJson(const QJsonObject& o)
{
    Camera c;
    auto v3 = [](const QJsonValue& v, QVector3D d) {
        const QJsonArray a = v.toArray();
        return a.size() == 3 ? QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble())) : d;
    };
    c.target = v3(o.value("target"), c.target);
    c.distance = float(o.value("distance").toDouble(c.distance));
    c.yaw = float(o.value("yaw").toDouble(c.yaw));
    c.pitch = float(o.value("pitch").toDouble(c.pitch));
    c.fovY = float(o.value("fov").toDouble(c.fovY));
    c.up = v3(o.value("up"), c.up);
    return c;
}

} // namespace kapp
