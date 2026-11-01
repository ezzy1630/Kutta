#include "core/ChartScale.h"

#include <QLocale>

#include <algorithm>
#include <cmath>

namespace kapp {

namespace {
constexpr double kLogFloor = 1e-300;
}

void Axis::setRange(double lo, double hi)
{
    if (!std::isfinite(lo) || !std::isfinite(hi)) return;
    if (hi < lo) std::swap(lo, hi);
    if (m_scale == Scale::Log) {
        lo = std::max(lo, kLogFloor);
        hi = std::max(hi, lo * 10);
    } else if (hi == lo) {
        const double pad = lo == 0 ? 1 : std::abs(lo) * 0.05;
        lo -= pad;
        hi += pad;
    }
    m_lo = lo;
    m_hi = hi;
}

double Axis::t(double v) const { return m_scale == Scale::Log ? std::log10(std::max(v, kLogFloor)) : v; }
double Axis::inv(double u) const { return m_scale == Scale::Log ? std::pow(10.0, u) : u; }

double Axis::toPixel(double v) const
{
    const double a = t(m_lo), b = t(m_hi);
    return m_p0 + (t(v) - a) / (b - a) * (m_p1 - m_p0);
}

double Axis::fromPixel(double p) const
{
    const double a = t(m_lo), b = t(m_hi);
    return inv(a + (p - m_p0) / (m_p1 - m_p0) * (b - a));
}

double niceStep(double rough)
{
    if (!(rough > 0) || !std::isfinite(rough)) return 1;
    const double e = std::pow(10.0, std::floor(std::log10(rough)));
    const double f = rough / e;
    const double n = f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10;
    return n * e;
}

std::vector<Tick> Axis::ticks(int target) const
{
    std::vector<Tick> out;
    if (m_scale == Scale::Log) {
        const int d0 = int(std::floor(std::log10(m_lo))), d1 = int(std::ceil(std::log10(m_hi)));
        const int decades = std::max(1, d1 - d0);
        // Label every few decades, counted down from the top. Prefer a step
        // that divides the span, so a residual history from 1e0 to 1e-10 is
        // labelled at both ends (1e0, 1e-5, 1e-10), not left unlabelled at
        // the floor where it ends.
        const int fewest = std::max(1, int(std::ceil(double(decades) / std::max(target, 1))));
        int every = fewest;
        for (int e = fewest; e <= std::min(decades, 2 * fewest); ++e)
            if (decades % e == 0) { every = e; break; }
        for (int d = d0; d <= d1; ++d) {
            const double v = std::pow(10.0, d);
            const bool major = (d1 - d) % every == 0;
            if (v >= m_lo * (1 - 1e-9) && v <= m_hi * (1 + 1e-9))
                out.push_back({v, major, major ? QStringLiteral("1e%1").arg(d) : QString()});
            if (decades <= 6)
                for (int m = 2; m <= 9; ++m) {
                    const double mv = m * v;
                    if (mv > m_lo && mv < m_hi) out.push_back({mv, false, {}});
                }
        }
        return out;
    }
    const double step = niceStep((m_hi - m_lo) / std::max(target, 1));
    const double start = std::ceil(m_lo / step - 1e-9) * step;
    for (double v = start; v <= m_hi + step * 1e-9; v += step) {
        const double clean = std::abs(v) < step * 1e-9 ? 0.0 : v;
        out.push_back({clean, true, formatTick(clean, step)});
        if (out.size() > 200) break;
    }
    return out;
}

void Axis::niceRange(Scale s, double& lo, double& hi)
{
    if (!std::isfinite(lo) || !std::isfinite(hi)) { lo = s == Scale::Log ? 1e-6 : 0; hi = 1; return; }
    if (s == Scale::Log) {
        lo = std::pow(10.0, std::floor(std::log10(std::max(lo, kLogFloor))));
        hi = std::pow(10.0, std::ceil(std::log10(std::max(hi, kLogFloor))));
        if (hi <= lo) hi = lo * 10;
        return;
    }
    if (hi == lo) { const double pad = lo == 0 ? 1 : std::abs(lo) * 0.05; lo -= pad; hi += pad; }
    const double step = niceStep((hi - lo) / 5);
    lo = std::floor(lo / step) * step;
    hi = std::ceil(hi / step) * step;
}

void Axis::zoom(double about, double factor)
{
    const double c = t(about), a = t(m_lo), b = t(m_hi);
    const double na = c + (a - c) * factor, nb = c + (b - c) * factor;
    if (nb - na < 1e-12 * std::max(1.0, std::abs(c))) return;
    m_lo = inv(na);
    m_hi = inv(nb);
}

void Axis::pan(double pixels)
{
    const double a = t(m_lo), b = t(m_hi);
    const double du = -pixels / (m_p1 - m_p0) * (b - a);
    m_lo = inv(a + du);
    m_hi = inv(b + du);
}

QString formatNumber(double v, int sig)
{
    if (!std::isfinite(v)) return QStringLiteral("—");
    if (v == 0) return QStringLiteral("0");
    const double a = std::abs(v);
    const QLocale c = QLocale::c();
    if (a >= 1e-3 && a < 1e6) {
        const int decimals = std::clamp(sig - 1 - int(std::floor(std::log10(a))), 0, 8);
        QString s = c.toString(v, 'f', decimals);
        if (s.contains('.')) {
            while (s.endsWith('0')) s.chop(1);
            if (s.endsWith('.')) s.chop(1);
        }
        return s;
    }
    QString s = c.toString(v, 'e', std::max(0, sig - 1));
    // "1.200e-05" -> "1.2e-05"
    const qsizetype e = s.indexOf('e');
    QString mant = s.left(e), exp = s.mid(e);
    if (mant.contains('.')) {
        while (mant.endsWith('0')) mant.chop(1);
        if (mant.endsWith('.')) mant.chop(1);
    }
    return mant + exp;
}

QString formatScientific(double v, int sig)
{
    if (!std::isfinite(v)) return QStringLiteral("—");
    if (v == 0) return QStringLiteral("0");
    const QString s = QLocale::c().toString(v, 'e', std::max(0, sig - 1));
    const qsizetype e = s.indexOf('e');
    QString mant = s.left(e);
    if (mant.contains('.')) {
        while (mant.endsWith('0')) mant.chop(1);
        if (mant.endsWith('.')) mant.chop(1);
    }
    return mant + 'e' + QString::number(s.mid(e + 1).toInt());
}

QString formatTick(double v, double step)
{
    if (v == 0) return QStringLiteral("0");
    const double a = std::max(std::abs(v), step);
    if (a >= 1e5 || step < 1e-4) return formatNumber(v, 3);
    const int decimals = std::max(0, -int(std::floor(std::log10(step) + 1e-9)));
    return QLocale::c().toString(v, 'f', decimals);
}

} // namespace kapp
