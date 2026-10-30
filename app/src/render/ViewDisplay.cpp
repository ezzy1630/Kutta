#include "render/ViewDisplay.h"

#include <QJsonArray>

namespace kapp {

ColormapId ViewDisplay::effectiveColormap() const
{
    if (colormapFollowsRole && symmetricRange() && !kapp::colormap(colormap).diverging) return defaultColormap(true);
    return colormap;
}

QJsonObject ViewDisplay::toJson() const
{
    QJsonArray vis;
    for (const QString& v : visible) vis.append(v);
    return {{"field", field.key()},
            {"contours", contours},
            {"colormap", kapp::colormap(colormap).key},
            {"colormap_follows_role", colormapFollowsRole},
            {"bands", bands},
            {"range", rangeMode == RangeMode::Auto ? "auto" : "fixed"},
            {"min", fixedMin},
            {"max", fixedMax},
            {"edges", edges},
            {"lighting", lighting},
            {"visible", vis},
            {"vectors", vectors},
            {"vector_item", vectorItem},
            {"vector_length", vectorLength},
            {"vector_stride", vectorStride},
            {"ribbon_width", ribbonWidth},
            {"camera", camera.toJson()}};
}

ViewDisplay ViewDisplay::fromJson(const QJsonObject& o)
{
    ViewDisplay d;
    if (o.contains("field")) d.field = FieldRef::fromKey(o.value("field").toString());
    d.contours = o.value("contours").toBool(true);
    d.colormap = colormapFromKey(o.value("colormap").toString()).value_or(ColormapId::Viridis);
    d.colormapFollowsRole = o.value("colormap_follows_role").toBool(true);
    d.bands = o.value("bands").toInt(0);
    d.rangeMode = o.value("range").toString() == "fixed" ? RangeMode::Fixed : RangeMode::Auto;
    d.fixedMin = o.value("min").toDouble(0);
    d.fixedMax = o.value("max").toDouble(1);
    d.edges = o.value("edges").toBool(false);
    d.lighting = o.value("lighting").toBool(true);
    for (const QJsonValue& v : o.value("visible").toArray()) d.visible << v.toString();
    d.visibleInitialized = o.contains("visible");
    d.vectors = o.value("vectors").toBool(false);
    d.vectorItem = o.value("vector_item").toString();
    d.vectorLength = float(o.value("vector_length").toDouble(0));
    d.vectorStride = o.value("vector_stride").toInt(4);
    d.ribbonWidth = float(o.value("ribbon_width").toDouble(2.5));
    if (o.contains("camera")) {
        d.camera = Camera::fromJson(o.value("camera").toObject());
        d.cameraInitialized = true;
    }
    return d;
}

} // namespace kapp
