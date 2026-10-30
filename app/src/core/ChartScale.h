#pragma once

#include <QString>

#include <vector>

namespace kapp {

struct Tick { double value; bool major; QString label; };

// One chart axis: maps data to pixels and produces readable ticks.
class Axis {
public:
    enum class Scale { Linear, Log };

    void setScale(Scale s) { m_scale = s; }
    Scale scale() const { return m_scale; }
    // Data range; log axes clamp to positive values.
    void setRange(double lo, double hi);
    double lo() const { return m_lo; }
    double hi() const { return m_hi; }
    void setPixels(double p0, double p1) { m_p0 = p0; m_p1 = p1; }

    double toPixel(double v) const;
    double fromPixel(double p) const;

    // Ticks for roughly `target` labels across the axis.
    std::vector<Tick> ticks(int target = 6) const;

    // Range covering data with nice end points (linear: 1/2/5 steps; log: decades).
    static void niceRange(Scale s, double& lo, double& hi);
    // Zoom about a data value by factor (<1 zooms in), in the axis's own scale.
    void zoom(double about, double factor);
    void pan(double pixels);

private:
    double t(double v) const;  // value -> transformed (log10 or identity)
    double inv(double u) const;
    Scale m_scale = Scale::Linear;
    double m_lo = 0, m_hi = 1, m_p0 = 0, m_p1 = 1;
};

double niceStep(double rough);
// Compact engineering formatting: "0.0012", "1.2e-05", "1500".
QString formatNumber(double v, int significant = 4);
QString formatTick(double v, double step);

} // namespace kapp
