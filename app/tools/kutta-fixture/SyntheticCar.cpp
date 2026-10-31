#include "SyntheticCar.h"

#include "core/Geometry.h"

#include <QJsonArray>
#include <QtMath>

#include <algorithm>
#include <cmath>
#include <functional>
#include <unordered_map>

namespace fixture {

using kapp::MeshData;

namespace {

constexpr double kPi = 3.14159265358979323846;

double smoothstep(double a, double b, double x)
{
    const double t = std::clamp((x - a) / (b - a), 0.0, 1.0);
    return t * t * (3 - 2 * t);
}

double spow(double v, double e) { return std::copysign(std::pow(std::abs(v), e), v); }

// ---- body shapes -------------------------------------------------------------------------

struct Section { double w, h, c, hump; };

constexpr double kNoseX = -1.25, kTailX = 1.30;

Section bodySection(double x)
{
    const double s = std::clamp((x - kNoseX) / (kTailX - kNoseX), 0.0, 1.0);
    const double cap = std::sqrt(std::max(0.0, 1 - std::pow(std::abs(2 * s - 1), 7)));
    Section q;
    q.w = (0.05 + 0.24 * smoothstep(0, 0.42, s) - 0.11 * smoothstep(0.72, 1, s)) * cap;
    q.h = (0.05 + 0.17 * smoothstep(0, 0.40, s) + 0.03 * smoothstep(0.4, 0.6, s) - 0.10 * smoothstep(0.78, 1, s)) * cap;
    q.c = 0.20 + 0.11 * smoothstep(0, 0.45, s);
    q.hump = 0.13 * std::exp(-std::pow((s - 0.66) / 0.10, 2)) * cap;  // engine cover / roll hoop fairing
    return q;
}

constexpr double kPodX0 = -0.05, kPodX1 = 0.95, kPodY = 0.43;

Section podSection(double x)
{
    const double s = std::clamp((x - kPodX0) / (kPodX1 - kPodX0), 0.0, 1.0);
    const double cap = std::sqrt(std::max(0.0, 1 - std::pow(std::abs(2 * s - 1), 6)));
    return {0.16 * cap * (1 - 0.35 * smoothstep(0.5, 1, s)), 0.13 * cap * (1 - 0.3 * smoothstep(0.55, 1, s)),
            0.24 - 0.04 * smoothstep(0.5, 1, s), 0};
}

struct Wing {
    double xLE, zLE, chord, aoaDeg, camber, thick, y0, y1;
    double cl;   // section lift coefficient (negative: downforce)
    bool front;
};

const std::vector<Wing>& wings()
{
    static const std::vector<Wing> w{
        {-1.66, 0.065, 0.32, 3, 0.07, 0.11, -0.72, 0.72, -1.7, true},
        {-1.39, 0.105, 0.18, 18, 0.06, 0.10, -0.72, -0.22, -1.3, true},
        {-1.39, 0.105, 0.18, 18, 0.06, 0.10, 0.22, 0.72, -1.3, true},
        {-1.25, 0.170, 0.13, 32, 0.05, 0.10, -0.72, -0.24, -1.1, true},
        {-1.25, 0.170, 0.13, 32, 0.05, 0.10, 0.24, 0.72, -1.1, true},
        {1.42, 0.92, 0.42, 4, 0.08, 0.12, -0.50, 0.50, -1.6, false},
        {1.80, 1.00, 0.25, 22, 0.07, 0.11, -0.50, 0.50, -1.3, false},
        {2.03, 1.12, 0.18, 38, 0.06, 0.10, -0.50, 0.50, -1.1, false},
    };
    return w;
}

struct Box { QVector3D lo, hi; };
const std::vector<Box>& endplates()
{
    static const std::vector<Box> b{
        {{-1.70f, -0.735f, 0.025f}, {-1.08f, -0.725f, 0.30f}}, {{-1.70f, 0.725f, 0.025f}, {-1.08f, 0.735f, 0.30f}},
        {{1.38f, -0.515f, 0.80f}, {2.26f, -0.505f, 1.32f}}, {{1.38f, 0.505f, 0.80f}, {2.26f, 0.515f, 1.32f}}};
    return b;
}

struct Tyre { QVector3D centre; };
const std::vector<Tyre>& tyres()
{
    static const std::vector<Tyre> t{{{-0.78f, -0.62f, 0.23f}}, {{-0.78f, 0.62f, 0.23f}},
                                     {{0.80f, -0.62f, 0.23f}}, {{0.80f, 0.62f, 0.23f}}};
    return t;
}
constexpr double kTyreR = 0.23, kTyreHalfW = 0.10;

// Airfoil (inverted cambered NACA-style) in chord coordinates.
void airfoil(double beta, double camber, double thick, double* yc, double* yt)
{
    *yt = 5 * thick * (0.2969 * std::sqrt(beta) - 0.126 * beta - 0.3516 * beta * beta + 0.2843 * beta * beta * beta -
                       0.1036 * beta * beta * beta * beta);
    *yc = -camber * 4 * beta * (1 - beta);  // camber toward the ground: downforce
}

// Local chord coordinates of a point relative to a wing element.
void toChord(const Wing& w, const QVector3D& p, double* s, double* n)
{
    const double a = qDegreesToRadians(w.aoaDeg);
    const double dx = p.x() - w.xLE, dz = p.z() - w.zLE;
    // Trailing edge up for a downforce element: rotate by +aoa about y.
    *s = (dx * std::cos(a) + dz * std::sin(a)) / w.chord;
    *n = (-dx * std::sin(a) + dz * std::cos(a)) / w.chord;
}

QVector3D fromChord(const Wing& w, double s, double n, double y)
{
    const double a = qDegreesToRadians(w.aoaDeg);
    return QVector3D(float(w.xLE + w.chord * (s * std::cos(a) - n * std::sin(a))), float(y),
                     float(w.zLE + w.chord * (s * std::sin(a) + n * std::cos(a))));
}

bool insideSolid(const QVector3D& p)
{
    const double x = p.x(), y = p.y(), z = p.z();
    if (z < 0) return true;
    if (x > kNoseX && x < kTailX) {
        const Section q = bodySection(x);
        const double zt = z - q.c;
        const double hh = zt > 0 ? q.h + q.hump : q.h;
        if (q.w > 1e-4 && hh > 1e-4 && std::pow(std::abs(y) / q.w, 3) + std::pow(std::abs(zt) / hh, 3) < 1) return true;
    }
    if (x > kPodX0 && x < kPodX1) {
        const Section q = podSection(x);
        if (q.w > 1e-4 && std::pow(std::abs(std::abs(y) - kPodY) / q.w, 3) + std::pow(std::abs(z - q.c) / q.h, 3) < 1)
            return true;
    }
    for (const Tyre& t : tyres()) {
        const double r = std::hypot(x - t.centre.x(), z - t.centre.z());
        if (r < kTyreR && std::abs(y - t.centre.y()) < kTyreHalfW) return true;
    }
    for (const Wing& w : wings()) {
        if (y < w.y0 || y > w.y1) continue;
        double s, n;
        toChord(w, p, &s, &n);
        if (s < 0 || s > 1) continue;
        double yc, yt;
        airfoil(s, w.camber, w.thick, &yc, &yt);
        if (std::abs(n - yc) < yt) return true;
    }
    for (const Box& b : endplates())
        if (x > b.lo.x() && x < b.hi.x() && y > b.lo.y() && y < b.hi.y() && z > b.lo.z() && z < b.hi.z()) return true;
    return false;
}

// ---- flow model ----------------------------------------------------------------------------

struct Source { QVector3D pos; double m; };       // point source strength (m^3/s)
struct Doublet { QVector3D pos; double mu; };      // x-oriented doublet
struct Segment { QVector3D a, b; double gamma; };  // straight vortex segment

struct FlowModel {
    double U = 16;
    std::vector<Source> sources;
    std::vector<Doublet> doublets;
    std::vector<Segment> segments;
    struct TipLine { QVector3D start; double strength; };
    std::vector<TipLine> tips;  // for total-pressure losses in vortex cores

    explicit FlowModel(double speed) : U(speed)
    {
        auto addLine = [&](double x0, double x1, int n, auto section, double y) {
            const double dx = (x1 - x0) / n;
            for (int i = 0; i < n; ++i) {
                const double xa = x0 + i * dx, xb = xa + dx, xm = 0.5 * (xa + xb);
                const Section sa = section(xa), sb = section(xb), sm = section(xm);
                const double Aa = kPi * sa.w * (sa.h + 0.5 * sa.hump), Ab = kPi * sb.w * (sb.h + 0.5 * sb.hump);
                const double m = U * (Ab - Aa);
                sources.push_back({QVector3D(float(xm), float(y), float(sm.c + 0.3 * sm.hump)), m});
            }
        };
        addLine(kNoseX, kTailX, 48, bodySection, 0.0);
        addLine(kPodX0, kPodX1, 24, podSection, kPodY);
        addLine(kPodX0, kPodX1, 24, podSection, -kPodY);
        for (const Tyre& t : tyres()) doublets.push_back({t.centre, U * std::pow(0.19, 3) / 2});
        // Wings: horseshoe vortices at chordwise stations, thin-airfoil loading.
        constexpr int stations = 5;
        for (const Wing& w : wings()) {
            const double gammaTotal = 0.5 * U * w.chord * w.cl * (w.front ? 1.25 : 1.0);
            double weights[stations], wsum = 0;
            for (int k = 0; k < stations; ++k) {
                const double s = (k + 0.5) / stations;
                weights[k] = std::sqrt((1 - s) / s);
                wsum += weights[k];
            }
            for (int k = 0; k < stations; ++k) {
                const double s = (k + 0.25) / stations;
                const double g = gammaTotal * weights[k] / wsum;
                const QVector3D a = fromChord(w, s, 0, w.y0), b = fromChord(w, s, 0, w.y1);
                const QVector3D farA(a.x() + 30, a.y(), a.z()), farB(b.x() + 30, b.y(), b.z());
                // Horseshoe: farA -> a -> b -> farB (bound vortex along +y).
                segments.push_back({farA, a, g});
                segments.push_back({a, b, g});
                segments.push_back({b, farB, g});
            }
            tips.push_back({fromChord(w, 1.0, 0, w.y0), std::abs(gammaTotal)});
            tips.push_back({fromChord(w, 1.0, 0, w.y1), std::abs(gammaTotal)});
        }
        // Ground images (z -> -z). Image vortices take the opposite circulation.
        const std::size_t ns = sources.size(), nd = doublets.size(), nv = segments.size();
        for (std::size_t i = 0; i < ns; ++i) {
            Source s = sources[i];
            s.pos.setZ(-s.pos.z());
            sources.push_back(s);
        }
        for (std::size_t i = 0; i < nd; ++i) {
            Doublet d = doublets[i];
            d.pos.setZ(-d.pos.z());
            doublets.push_back(d);
        }
        for (std::size_t i = 0; i < nv; ++i) {
            Segment s = segments[i];
            s.a.setZ(-s.a.z());
            s.b.setZ(-s.b.z());
            s.gamma = -s.gamma;
            segments.push_back(s);
        }
    }

    QVector3D potential(const QVector3D& p) const
    {
        double ux = U, uy = 0, uz = 0;
        for (const Source& s : sources) {
            const double rx = p.x() - s.pos.x(), ry = p.y() - s.pos.y(), rz = p.z() - s.pos.z();
            const double r2 = rx * rx + ry * ry + rz * rz + 4e-4;
            const double k = s.m / (4 * kPi * r2 * std::sqrt(r2));
            ux += k * rx; uy += k * ry; uz += k * rz;
        }
        for (const Doublet& d : doublets) {
            const double rx = p.x() - d.pos.x(), ry = p.y() - d.pos.y(), rz = p.z() - d.pos.z();
            const double r2 = rx * rx + ry * ry + rz * rz + 4e-4, r = std::sqrt(r2);
            const double xr = rx / r;
            const double k = d.mu / (r2 * r);
            ux += k * (3 * xr * rx / r - 1); uy += k * 3 * xr * ry / r; uz += k * 3 * xr * rz / r;
        }
        for (const Segment& s : segments) {
            // Biot-Savart for a finite segment with a smoothed core.
            const QVector3D r1 = p - s.a, r2v = p - s.b, r0 = s.b - s.a;
            const QVector3D c = QVector3D::crossProduct(r1, r2v);
            const double c2 = c.lengthSquared();
            const double l1 = r1.length(), l2 = r2v.length();
            if (l1 < 1e-6 || l2 < 1e-6) continue;
            const double core = 0.03 + 0.02 * std::max(0.0, double(p.x() - s.a.x())) * 0.1;
            const double d2 = c2 / std::max(double(r0.lengthSquared()), 1e-12);
            const double k = s.gamma / (4 * kPi) * QVector3D::dotProduct(r0, r1 / float(l1) - r2v / float(l2)) /
                             (c2 + core * core * r0.lengthSquared());
            (void)d2;
            ux += k * c.x(); uy += k * c.y(); uz += k * c.z();
        }
        return QVector3D(float(ux), float(uy), float(uz));
    }

    // Non-potential wake deficits (fraction of U) and total-pressure losses.
    double wakeDeficit(const QVector3D& p) const
    {
        const double x = p.x(), y = p.y(), z = p.z();
        double d = 0;
        if (x > 1.1) {
            const double dx = x - 1.1, sy = 0.24 + 0.07 * dx, sz = 0.20 + 0.05 * dx;
            d += 0.38 * std::exp(-dx / 5.0) * std::exp(-y * y / (sy * sy) - std::pow(z - 0.36, 2) / (sz * sz)) * smoothstep(1.1, 1.5, x);
        }
        for (const Tyre& t : tyres()) {
            const double dx = x - t.centre.x() - 0.15;
            if (dx < -0.15) continue;
            const double s = 0.14 + 0.05 * std::max(0.0, dx);
            d += 0.62 * std::exp(-std::max(0.0, dx) / 2.5) * smoothstep(-0.15, 0.1, dx) *
                 std::exp(-std::pow(y - t.centre.y(), 2) / (s * s) - std::pow(z - 0.20, 2) / (s * s));
        }
        if (x > 2.1) {
            const double dx = x - 2.1, sz = 0.05 + 0.04 * dx, zc = 1.16 + 0.10 * dx;
            d += 0.30 * std::exp(-dx / 3.0) * std::exp(-std::pow(z - zc, 2) / (sz * sz)) * smoothstep(0.55, 0.40, std::abs(y));
        }
        return std::min(d, 0.95);
    }

    double coreLoss(const QVector3D& p) const
    {
        double loss = 0;
        for (const TipLine& t : tips) {
            const double dx = p.x() - t.start.x();
            if (dx < 0) continue;
            // Tip vortices drift inboard/up slightly as they convect.
            const double yc = t.start.y() * (1 - 0.04 * std::min(dx, 6.0)), zc = t.start.z() + 0.02 * dx;
            const double rc = 0.035 + 0.012 * dx;
            const double r2 = std::pow(p.y() - yc, 2) + std::pow(p.z() - zc, 2);
            loss += 1.6 * std::min(1.0, t.strength / 4.0) * std::exp(-r2 / (rc * rc)) * std::exp(-dx / 7.0);
        }
        return loss;
    }

    QVector3D velocity(const QVector3D& p) const
    {
        QVector3D u = potential(p);
        const double d = wakeDeficit(p);
        u.setX(u.x() * float(1 - d));
        return u;
    }

    double totalPressureCoefficient(const QVector3D& p) const
    {
        const double d = wakeDeficit(p);
        return 1 - (1 - (1 - d) * (1 - d)) * 1.15 - coreLoss(p);
    }
};

const FlowModel& model(double speed)
{
    // Rebuilt only when the inlet speed changes (overrides).
    static std::unique_ptr<FlowModel> m;
    if (!m || m->U != speed) m = std::make_unique<FlowModel>(speed);
    return *m;
}

// Smooth pseudo-random field for Fluent-vs-Kutta differences.
double smoothNoise(const QVector3D& p, int seed, double scale = 1.0)
{
    double v = 0;
    for (int i = 0; i < 6; ++i) {
        const double a = std::sin(seed * 12.9898 + i * 78.233) * 43758.5453;
        const double kx = (std::fmod(a, 1.0) - 0.5) * 6 / scale, ky = (std::fmod(a * 1.7, 1.0) - 0.5) * 6 / scale,
                     kz = (std::fmod(a * 2.3, 1.0) - 0.5) * 6 / scale;
        v += std::sin(kx * p.x() + ky * p.y() + kz * p.z() + i * 1.3 + seed);
    }
    return v / 6;
}

// The "disagreement" the comparison view should expose: rear-wing endplate
// vortex and the front-wing tip, plus small broadband noise.
double kuttaPerturbation(const QVector3D& p)
{
    double d = 0.04 * smoothNoise(p, 99, 0.7);
    d += -0.32 * std::exp(-std::pow(std::abs(p.y()) - 0.48, 2) / 0.003 - std::pow(p.z() - 1.0, 2) / 0.02) *
         smoothstep(1.5, 1.9, p.x());
    d += 0.18 * std::exp(-std::pow(std::abs(p.y()) - 0.70, 2) / 0.004 - std::pow(p.x() + 1.3, 2) / 0.05) *
         (p.z() < 0.3 ? 1 : 0);
    return d;
}

// ---- tessellation --------------------------------------------------------------------------------

struct Part {
    std::vector<QVector3D> pos, nrm;
    std::vector<std::uint32_t> tri;
    // Wing attributes per vertex (chord fraction, side +1 top/-1 bottom, span fraction).
    std::vector<float> chord, side, span;
    enum Kind { Body, Wing, Tyre, Plane } kind = Body;
    int wingIndex = -1;
};

using Surf = std::function<QVector3D(double u, double v)>;

// Parametric grid; u wraps when `wrapU`. Normals from central differences,
// oriented by `outward(p, n)`.
void grid(Part& part, int nu, int nv, bool wrapU, const Surf& f, const std::function<bool(QVector3D, QVector3D)>& outward)
{
    const std::uint32_t base = std::uint32_t(part.pos.size());
    const int cols = wrapU ? nu : nu + 1;
    for (int j = 0; j <= nv; ++j)
        for (int i = 0; i < cols; ++i) {
            const double u = double(i) / nu, v = double(j) / nv;
            const double eu = 1e-4, ev = 1e-4;
            const QVector3D p = f(u, v);
            const QVector3D du = f(u + eu, v) - f(u - eu, v);
            const QVector3D dv = f(u, std::min(1.0, v + ev)) - f(u, std::max(0.0, v - ev));
            QVector3D n = QVector3D::crossProduct(du, dv).normalized();
            if (n.isNull()) n = QVector3D(0, 0, 1);
            if (!outward(p, n)) n = -n;
            part.pos.push_back(p);
            part.nrm.push_back(n);
        }
    // Orient each triangle by its own face normal: vertex normals are unreliable
    // at degenerate grid points (nose tips, disc centres, superellipse corners).
    auto addTri = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
        const QVector3D pa = part.pos[a], pb = part.pos[b], pc = part.pos[c];
        const QVector3D fn = QVector3D::crossProduct(pb - pa, pc - pa);
        if (fn.lengthSquared() < 1e-24f) return;  // degenerate (collapsed row)
        if (outward((pa + pb + pc) / 3.0f, fn.normalized())) part.tri.insert(part.tri.end(), {a, b, c});
        else part.tri.insert(part.tri.end(), {a, c, b});
    };
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
            const std::uint32_t a = base + std::uint32_t(j * cols + i), b = base + std::uint32_t(j * cols + (i + 1) % cols);
            const std::uint32_t c = base + std::uint32_t((j + 1) * cols + i), d = base + std::uint32_t((j + 1) * cols + (i + 1) % cols);
            addTri(a, b, c);
            addTri(b, d, c);
        }
}

void loft(Part& part, double x0, double x1, double yc, Section (*section)(double), int nu, int nv)
{
    grid(part, nu, nv, true,
         [=](double u, double v) {
             const double x = x0 + (x1 - x0) * (0.5 - 0.5 * std::cos(kPi * std::clamp(v, 0.0, 1.0)));
             const Section q = section(x);
             const double t = 2 * kPi * u;
             const double sz = std::sin(t);
             const double hh = sz > 0 ? q.h + q.hump : q.h;
             return QVector3D(float(x), float(yc + q.w * spow(std::cos(t), 2.0 / 3.0)), float(q.c + hh * spow(sz, 2.0 / 3.0)));
         },
         [=](QVector3D p, QVector3D n) {
             const Section q = section(p.x());
             return QVector3D::dotProduct(n, QVector3D(0, float(p.y() - yc), float(p.z() - q.c))) >= 0 ||
                    (p.x() <= x0 + 1e-3 ? n.x() < 0 : p.x() >= x1 - 1e-3 ? n.x() > 0 : false);
         });
}

void wingPart(Part& part, const Wing& w, int nu, int nv)
{
    const std::size_t start = part.pos.size();
    auto loop = [&](double u, double* s, double* n, int* side) {
        // u in [0,1): top surface TE->LE, then bottom LE->TE (cosine spacing at LE).
        const double t = u < 0.5 ? 1 - 2 * u : 2 * u - 1;
        const double beta = 0.5 - 0.5 * std::cos(kPi * t);
        double yc, yt;
        airfoil(beta, w.camber, w.thick, &yc, &yt);
        *s = beta;
        *side = u < 0.5 ? 1 : -1;
        *n = yc + *side * yt;
    };
    grid(part, nu, nv, true,
         [&](double u, double v) {
             double s, n;
             int side;
             loop(std::fmod(u + 1.0, 1.0), &s, &n, &side);
             return fromChord(w, s, n, w.y0 + (w.y1 - w.y0) * std::clamp(v, 0.0, 1.0));
         },
         [&](QVector3D p, QVector3D nrm) {
             double s, n;
             toChord(w, p, &s, &n);
             double yc, yt;
             airfoil(std::clamp(s, 0.0, 1.0), w.camber, w.thick, &yc, &yt);
             const QVector3D mid = fromChord(w, std::clamp(s, 0.0, 1.0), yc, p.y());
             return QVector3D::dotProduct(nrm, p - mid) >= 0;
         });
    const int cols = nu;
    for (int j = 0; j <= nv; ++j)
        for (int i = 0; i < cols; ++i) {
            double s, n;
            int side;
            loop(double(i) / nu, &s, &n, &side);
            part.chord.push_back(float(s));
            part.side.push_back(float(side));
            const double yv = w.y0 + (w.y1 - w.y0) * double(j) / nv;
            const double mid = 0.5 * (w.y0 + w.y1), half = 0.5 * (w.y1 - w.y0);
            part.span.push_back(float((yv - mid) / half));
        }
    // End caps: fans over each tip section.
    for (int end = 0; end < 2; ++end) {
        const double y = end ? w.y1 : w.y0;
        const std::uint32_t centre = std::uint32_t(part.pos.size());
        part.pos.push_back(fromChord(w, 0.45, -w.camber * 0.99, y));
        part.nrm.push_back(QVector3D(0, end ? 1.f : -1.f, 0));
        part.chord.push_back(0.45f);
        part.side.push_back(0);
        part.span.push_back(end ? 1.f : -1.f);
        const std::uint32_t ring = std::uint32_t(part.pos.size());
        for (int i = 0; i < nu; ++i) {
            const std::size_t src = start + std::size_t((end ? nv : 0) * cols + i);
            part.pos.push_back(part.pos[src]);
            part.nrm.push_back(QVector3D(0, end ? 1.f : -1.f, 0));
            part.chord.push_back(part.chord[src]);
            part.side.push_back(part.side[src]);
            part.span.push_back(end ? 1.f : -1.f);
        }
        for (int i = 0; i < nu; ++i) {
            const std::uint32_t a = centre, b = ring + std::uint32_t(i), c = ring + std::uint32_t((i + 1) % nu);
            const QVector3D fn = QVector3D::crossProduct(part.pos[b] - part.pos[a], part.pos[c] - part.pos[a]);
            if ((fn.y() > 0) == (end == 1)) part.tri.insert(part.tri.end(), {a, b, c});
            else part.tri.insert(part.tri.end(), {a, c, b});
        }
    }
}

void boxPart(Part& part, const Box& b, int n)
{
    const QVector3D lo = b.lo, hi = b.hi;
    for (int axis = 0; axis < 3; ++axis)
        for (int sgn = -1; sgn <= 1; sgn += 2) {
            const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
            grid(part, n, n, false,
                 [=](double u, double v) {
                     QVector3D p;
                     p[axis] = sgn < 0 ? lo[axis] : hi[axis];
                     p[a1] = float(lo[a1] + (hi[a1] - lo[a1]) * std::clamp(u, 0.0, 1.0));
                     p[a2] = float(lo[a2] + (hi[a2] - lo[a2]) * std::clamp(v, 0.0, 1.0));
                     return p;
                 },
                 [=](QVector3D, QVector3D nn) { return nn[axis] * sgn > 0; });
        }
}

void tyrePart(Part& part, const Tyre& t, int nu, int nv)
{
    const QVector3D c = t.centre;
    // Tyre: superelliptic section revolved about the axle (y).
    grid(part, nu, nv, true,
         [=](double u, double v) {
             const double a = 2 * kPi * u, th = 2 * kPi * std::fmod(v + 1.0, 1.0);
             const double rr = 0.185 + 0.048 * spow(std::cos(a), 0.45);
             const double yy = kTyreHalfW * spow(std::sin(a), 0.45);
             return QVector3D(float(c.x() + rr * std::cos(th)), float(c.y() + yy), float(c.z() + rr * std::sin(th)));
         },
         [=](QVector3D p, QVector3D n) {
             const QVector3D radial(p.x() - c.x(), 0, p.z() - c.z());
             const QVector3D ring = radial.normalized() * 0.185f + QVector3D(c.x(), p.y(), c.z());
             return QVector3D::dotProduct(n, p - ring) >= 0;
         });
    // Rims: flat discs on both faces.
    for (int sgn = -1; sgn <= 1; sgn += 2)
        grid(part, nu / 2, 6, true,
             [=](double u, double v) {
                 const double th = 2 * kPi * u, r = 0.142 * std::clamp(v, 0.0, 1.0);
                 return QVector3D(float(c.x() + r * std::cos(th)), float(c.y() + sgn * 0.075), float(c.z() + r * std::sin(th)));
             },
             [=](QVector3D, QVector3D n) { return n.y() * sgn > 0; });
}

void planePart(Part& part, QVector3D origin, QVector3D du, QVector3D dv, int nu, int nv, QVector3D outward)
{
    grid(part, nu, nv, false,
         [=](double u, double v) { return origin + du * float(std::clamp(u, 0.0, 1.0)) + dv * float(std::clamp(v, 0.0, 1.0)); },
         [=](QVector3D, QVector3D n) { return QVector3D::dotProduct(n, outward) > 0; });
}

std::vector<Part> zoneParts(const QString& zone, int d)
{
    std::vector<Part> parts;
    auto add = [&](Part::Kind k) -> Part& { parts.emplace_back(); parts.back().kind = k; return parts.back(); };
    if (zone == "monocoque") {
        loft(add(Part::Body), kNoseX, kTailX, 0, bodySection, 140 * d, 220 * d);
    } else if (zone == "sidepods") {
        loft(add(Part::Body), kPodX0, kPodX1, kPodY, podSection, 80 * d, 110 * d);
        loft(add(Part::Body), kPodX0, kPodX1, -kPodY, podSection, 80 * d, 110 * d);
    } else if (zone == "front-wing" || zone == "rear-wing") {
        const bool front = zone == "front-wing";
        for (std::size_t i = 0; i < wings().size(); ++i) {
            const Wing& w = wings()[i];
            if (w.front != front) continue;
            Part& p = add(Part::Wing);
            p.wingIndex = int(i);
            const double span = w.y1 - w.y0;
            wingPart(p, w, 120 * d, std::max(20, int(110 * d * span / 1.2)));
        }
        Part& plates = add(Part::Body);
        for (const Box& b : endplates())
            if ((b.lo.x() < 0) == front) boxPart(plates, b, 18 * d);
    } else if (zone == "tyre-front" || zone == "tyre-rear") {
        for (const Tyre& t : tyres())
            if ((t.centre.x() < 0) == (zone == "tyre-front")) tyrePart(add(Part::Tyre), t, 56 * d, 150 * d);
    } else if (zone == "ground") {
        planePart(add(Part::Plane), {-3, -2.5, 0}, {9, 0, 0}, {0, 5, 0}, 180 * d, 100 * d, {0, 0, 1});
    } else if (zone == "inlet") {
        planePart(add(Part::Plane), {-6, -4, 0}, {0, 8, 0}, {0, 0, 4}, 24, 12, {1, 0, 0});
    } else if (zone == "outlet") {
        planePart(add(Part::Plane), {14, -4, 0}, {0, 8, 0}, {0, 0, 4}, 24, 12, {-1, 0, 0});
    } else if (zone == "tunnel") {
        planePart(add(Part::Plane), {-6, -4, 4}, {20, 0, 0}, {0, 8, 0}, 40, 16, {0, 0, -1});
        planePart(add(Part::Plane), {-6, -4, 0}, {20, 0, 0}, {0, 0, 4}, 40, 8, {0, 1, 0});
        planePart(add(Part::Plane), {-6, 4, 0}, {20, 0, 0}, {0, 0, 4}, 40, 8, {0, -1, 0});
    }
    return parts;
}

// Pressure coefficient on a surface vertex.
double surfaceCp(const Part& part, std::size_t i, const FlowModel& m, double convergence)
{
    const QVector3D p = part.pos[i], n = part.nrm[i];
    double cp;
    if (part.kind == Part::Wing) {
        const Wing& w = wings()[std::size_t(part.wingIndex)];
        const double s = std::clamp(double(part.chord[i]), 1e-3, 1.0), side = part.side[i];
        const double tip = 1 - std::pow(std::min(1.0, std::abs(double(part.span[i]))), 8);
        const double le = std::exp(-s / 0.03);
        const double load = std::abs(w.cl) * (w.front ? 1.25 : 1.0);
        if (side < 0) cp = -load * (1.15 * std::pow(1 - s, 0.7) + 0.9 * le) * (0.55 + 0.45 * tip) + 0.15 * (1 - tip);
        else if (side > 0) cp = 0.55 * (1 - s) * load * 0.5 + 0.8 * le * (1 - s);
        else cp = -0.2;
        if (w.front) cp *= 1 + 0.6 * std::exp(-p.z() / 0.06);  // ground effect on the main plane
    } else if (part.kind == Part::Tyre) {
        const double f = -n.x();
        cp = 0.95 * (f > 0 ? f * f : 0) - 0.75 * (1 - f * f) + 0.1 * (p.z() < 0.05 ? 1 : 0);
    } else {
        const QVector3D q = p + n * 0.006f;
        const QVector3D u = m.velocity(q);
        cp = 1 - double(u.lengthSquared()) / (m.U * m.U) - (1 - m.totalPressureCoefficient(q)) * 0.3;
        cp = std::clamp(cp, -3.5, 1.0);
    }
    // Live snapshots: early iterations have weaker loading plus transient noise.
    if (convergence < 1) cp = cp * convergence + (1 - convergence) * 0.25 * smoothNoise(p, 7, 0.5);
    return cp;
}

MeshData::Field makeField(const QString& name, const QString& unit, int comps, std::vector<float> v)
{
    MeshData::Field f;
    f.desc.name = name;
    f.desc.unit = unit;
    f.desc.components = comps;
    f.values = std::move(v);
    return f;
}

// Add a field and, with reference data, its @fluent and @delta companions.
void addWithReference(MeshData& m, const QString& name, const QString& unit, int comps, const std::vector<float>& kutta,
                      const std::vector<float>* fluent)
{
    m.fields.push_back(makeField(name, unit, comps, kutta));
    if (!fluent) return;
    m.fields.push_back(makeField(name + "@fluent", unit, comps, *fluent));
    std::vector<float> delta(kutta.size());
    for (std::size_t i = 0; i < kutta.size(); ++i) delta[i] = kutta[i] - (*fluent)[i];
    m.fields.push_back(makeField(name + "@delta", unit, comps, delta));
}

} // namespace

const std::vector<CarZone>& carZones()
{
    static const std::vector<CarZone> z{
        {"monocoque", "wall", 10, true},     {"sidepods", "wall", 11, true},   {"front-wing", "wall", 12, true},
        {"rear-wing", "wall", 13, true},     {"tyre-front", "wall", 14, true}, {"tyre-rear", "wall", 15, true},
        {"ground", "wall", 6, true},         {"inlet", "velocity-inlet", 2, true},
        {"outlet", "pressure-outlet", 3, true}, {"tunnel", "symmetry", 4, true},
    };
    return z;
}

MeshData carZoneMesh(const QString& zone, const CarParams& prm, bool reference)
{
    const FlowModel& m = model(prm.speed);
    const std::vector<Part> parts = zoneParts(zone, prm.detail);
    MeshData mesh;
    std::vector<float> cpK, cpF, tauK, tauF, ypK, ypF, velK, velF;
    const bool wall = zone != "inlet" && zone != "outlet" && zone != "tunnel";
    const double q = 0.5 * prm.density * prm.speed * prm.speed;
    const double qRef = 0.5 * prm.refDensity * prm.refVelocity * prm.refVelocity;
    for (const Part& part : parts) {
        const std::uint32_t base = std::uint32_t(mesh.positions.size() / 3);
        for (std::size_t i = 0; i < part.pos.size(); ++i) {
            const QVector3D p = part.pos[i];
            mesh.positions.insert(mesh.positions.end(), {p.x(), p.y(), p.z()});
            const double cp = surfaceCp(part, i, m, prm.convergence);
            const double fluentCp = cp + 0.025 * smoothNoise(p, 7);
            const double kuttaCp = fluentCp + kuttaPerturbation(p);
            cpK.push_back(float(kuttaCp));
            cpF.push_back(float(fluentCp));
            if (wall) {
                const double edge = prm.speed * std::sqrt(std::clamp(1 - fluentCp, 0.05, 6.0));
                // Skin friction falls slowly with Reynolds number (cf ~ Re^-0.2), so
                // a speed sweep moves the drag coefficient a little.
                const double cf = (zone == "ground" ? 0.0015 : 0.0042) * std::pow(prm.speed / 16.0, -0.2);
                const double tauF0 = 0.5 * prm.density * edge * edge * cf;
                const double tauK0 = tauF0 * (1 + 0.04 * smoothNoise(p, 31));
                // Shear along the free-stream direction projected on the surface.
                QVector3D t = QVector3D(1, 0, 0) - part.nrm[i] * QVector3D::dotProduct(QVector3D(1, 0, 0), part.nrm[i]);
                t = t.isNull() ? QVector3D(0, 0, 0) : t.normalized();
                for (int c = 0; c < 3; ++c) {
                    tauK.push_back(float(tauK0 * t[c]));
                    tauF.push_back(float(tauF0 * t[c]));
                }
                const double ypF0 = 48 * edge / prm.speed * (1 + 0.2 * smoothNoise(p, 5));
                ypF.push_back(float(ypF0));
                ypK.push_back(float(ypF0 * (1 + 0.03 * smoothNoise(p, 77))));
            } else {
                const QVector3D u = m.velocity(p) * float(prm.convergence < 1 ? 0.6 + 0.4 * prm.convergence : 1.0);
                for (int c = 0; c < 3; ++c) {
                    velF.push_back(u[c]);
                    velK.push_back(float(u[c] * (1 + 0.01 * smoothNoise(p, 3))));
                }
            }
        }
        for (std::uint32_t t : part.tri) mesh.indices.push_back(base + t);
    }
    auto scaled = [](const std::vector<float>& v, double k, double add = 0) {
        std::vector<float> o(v.size());
        for (std::size_t i = 0; i < v.size(); ++i) o[i] = float(v[i] * k + add);
        return o;
    };
    auto magnitude = [](const std::vector<float>& v) {
        std::vector<float> o(v.size() / 3);
        for (std::size_t i = 0; i < o.size(); ++i) o[i] = std::sqrt(v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]);
        return o;
    };
    const auto pK = scaled(cpK, q, 0), pF = scaled(cpF, q, 0);
    const auto pF2 = pF;
    addWithReference(mesh, "pressure", "Pa", 1, pK, reference ? &pF2 : nullptr);
    // Cp uses the reference values, exactly like Fluent's pressure coefficient.
    const auto cpRefK = scaled(pK, 1 / qRef, -prm.refPressure / qRef), cpRefF = scaled(pF, 1 / qRef, -prm.refPressure / qRef);
    addWithReference(mesh, "pressure_coefficient", "1", 1, cpRefK, reference ? &cpRefF : nullptr);
    if (wall) {
        addWithReference(mesh, "wall_shear", "Pa", 3, tauK, reference ? &tauF : nullptr);
        const auto mk = magnitude(tauK), mf = magnitude(tauF);
        addWithReference(mesh, "wall_shear_magnitude", "Pa", 1, mk, reference ? &mf : nullptr);
        addWithReference(mesh, "y_plus", "1", 1, ypK, reference ? &ypF : nullptr);
    } else {
        addWithReference(mesh, "velocity", "m/s", 3, velK, reference ? &velF : nullptr);
        const auto mk = magnitude(velK), mf = magnitude(velF);
        addWithReference(mesh, "velocity_magnitude", "m/s", 1, mk, reference ? &mf : nullptr);
    }
    return mesh;
}

MeshData carSlice(const QVector3D& origin, const QVector3D& normalIn, const CarParams& prm, bool reference)
{
    const FlowModel& m = model(prm.speed);
    const QVector3D n = normalIn.normalized();
    const QVector3D helper = std::abs(n.z()) < 0.9f ? QVector3D(0, 0, 1) : QVector3D(1, 0, 0);
    const QVector3D e1 = QVector3D::crossProduct(helper, n).normalized();
    const QVector3D e2 = QVector3D::crossProduct(n, e1);
    // Region of interest around the car; the slice spans its projection.
    const QVector3D lo(-2.4f, -1.5f, 0.0f), hi(5.5f, 1.5f, 1.7f);
    double a0 = 1e9, a1 = -1e9, b0 = 1e9, b1 = -1e9;
    for (int k = 0; k < 8; ++k) {
        const QVector3D c((k & 1) ? hi.x() : lo.x(), (k & 2) ? hi.y() : lo.y(), (k & 4) ? hi.z() : lo.z());
        const double a = QVector3D::dotProduct(c - origin, e1), b = QVector3D::dotProduct(c - origin, e2);
        a0 = std::min(a0, a); a1 = std::max(a1, a); b0 = std::min(b0, b); b1 = std::max(b1, b);
    }
    const double h = 0.016 / std::sqrt(double(prm.detail));
    const int na = std::max(2, int((a1 - a0) / h)), nb = std::max(2, int((b1 - b0) / h));
    MeshData mesh;
    std::vector<char> solid;
    std::vector<float> vel, velF, vmag, vmagF, cp, cpF, cp0, cp0F, kk, kkF, om, omF, vr, vrF;
    const double q = 0.5 * prm.density * prm.speed * prm.speed;
    const double conv = prm.convergence;
    for (int j = 0; j <= nb; ++j)
        for (int i = 0; i <= na; ++i) {
            QVector3D p = origin + e1 * float(a0 + (a1 - a0) * i / na) + e2 * float(b0 + (b1 - b0) * j / nb);
            const bool outside = p.x() < lo.x() - 1e-3f || p.x() > hi.x() + 1e-3f || p.y() < lo.y() - 1e-3f ||
                                 p.y() > hi.y() + 1e-3f || p.z() < -1e-3f || p.z() > hi.z() + 1e-3f;
            const bool inSolid = outside || insideSolid(p);
            solid.push_back(inSolid ? 1 : 0);
            mesh.positions.insert(mesh.positions.end(), {p.x(), p.y(), p.z()});
            QVector3D u = inSolid ? QVector3D() : m.velocity(p);
            if (conv < 1) u = QVector3D(float(prm.speed), 0, 0) * float(1 - conv) + u * float(conv);
            const double c0 = inSolid ? 1 : 1 - (1 - m.totalPressureCoefficient(p)) * conv;
            const double cpv = c0 - double(u.lengthSquared()) / (prm.speed * prm.speed);
            const double noiseF = 0.012 * smoothNoise(p, 3);
            const double pert = kuttaPerturbation(p) * 0.5;
            const double turb = std::max(0.0, 1 - c0);
            for (int c = 0; c < 3; ++c) {
                velF.push_back(float(u[c] * (1 + noiseF)));
                vel.push_back(float(u[c] * (1 + noiseF) * (1 + 0.6 * pert)));
            }
            vmagF.push_back(float(u.length() * (1 + noiseF)));
            vmag.push_back(float(u.length() * (1 + noiseF) * (1 + 0.6 * pert)));
            cpF.push_back(float(cpv + noiseF));
            cp.push_back(float(cpv + noiseF + pert));
            cp0F.push_back(float(c0));
            cp0.push_back(float(c0 + pert * 0.8));
            const double kv = 1.5 * std::pow(0.01 * prm.speed, 2) + 2.4 * turb * prm.speed * prm.speed * 0.02;
            kkF.push_back(float(kv));
            kk.push_back(float(kv * (1 + 0.08 * smoothNoise(p, 11))));
            const double ov = 40 + 900 * turb;
            omF.push_back(float(ov));
            om.push_back(float(ov * (1 + 0.06 * smoothNoise(p, 13))));
            const double vrv = 10 + 4000 * turb;
            vrF.push_back(float(vrv));
            vr.push_back(float(vrv * (1 + 0.1 * smoothNoise(p, 17))));
        }
    // Keep only fluid triangles; then compact unused vertices.
    const int cols = na + 1;
    std::vector<std::uint32_t> tris;
    for (int j = 0; j < nb; ++j)
        for (int i = 0; i < na; ++i) {
            const std::uint32_t a = std::uint32_t(j * cols + i), b = a + 1, c = a + std::uint32_t(cols), d = c + 1;
            if (!solid[a] && !solid[b] && !solid[c]) tris.insert(tris.end(), {a, b, c});
            if (!solid[b] && !solid[d] && !solid[c]) tris.insert(tris.end(), {b, d, c});
        }
    std::vector<std::int64_t> remap(solid.size(), -1);
    std::vector<std::uint32_t> keep;
    for (std::uint32_t& t : tris) {
        if (remap[t] < 0) { remap[t] = std::int64_t(keep.size()); keep.push_back(t); }
        t = std::uint32_t(remap[t]);
    }
    auto compact = [&](const std::vector<float>& v, int comps) {
        std::vector<float> o;
        o.reserve(keep.size() * std::size_t(comps));
        for (std::uint32_t k : keep)
            for (int c = 0; c < comps; ++c) o.push_back(v[std::size_t(k) * std::size_t(comps) + std::size_t(c)]);
        return o;
    };
    MeshData out;
    out.positions = compact(mesh.positions, 3);
    out.indices = std::move(tris);
    auto add = [&](const QString& name, const QString& unit, int comps, const std::vector<float>& k, const std::vector<float>& f) {
        const auto ck = compact(k, comps), cf = compact(f, comps);
        addWithReference(out, name, unit, comps, ck, reference ? &cf : nullptr);
    };
    add("velocity", "m/s", 3, vel, velF);
    add("velocity_magnitude", "m/s", 1, vmag, vmagF);
    std::vector<float> pK(cp.size()), pF(cpF.size());
    for (std::size_t i = 0; i < cp.size(); ++i) { pK[i] = float(cp[i] * q); pF[i] = float(cpF[i] * q); }
    add("pressure", "Pa", 1, pK, pF);
    add("pressure_coefficient", "1", 1, cp, cpF);
    add("total_pressure_coefficient", "1", 1, cp0, cp0F);
    add("k", "m^2/s^2", 1, kk, kkF);
    add("omega", "1/s", 1, om, omF);
    add("viscosity_ratio", "1", 1, vr, vrF);
    return out;
}

MeshData carStreamlines(const QVector3D& from, const QVector3D& to, int count, double maxLength, const CarParams& prm)
{
    const FlowModel& m = model(prm.speed);
    MeshData out;
    out.polylines = true;
    out.indices.push_back(0);
    std::vector<float> vmag, cp0;
    const double h = 0.012;
    const int maxSteps = int(maxLength / h);
    auto inDomain = [](const QVector3D& p) {
        return p.x() > -3 && p.x() < 7 && std::abs(p.y()) < 2.5f && p.z() > 0.002f && p.z() < 2.5f;
    };
    for (int s = 0; s < count; ++s) {
        QVector3D p = count == 1 ? from : from + (to - from) * (float(s) / float(count - 1));
        std::size_t n = 0;
        for (int step = 0; step < maxSteps && inDomain(p) && !insideSolid(p); ++step) {
            const QVector3D u = m.velocity(p);
            out.positions.insert(out.positions.end(), {p.x(), p.y(), p.z()});
            vmag.push_back(u.length());
            cp0.push_back(float(m.totalPressureCoefficient(p)));
            ++n;
            // RK2 with a step of fixed length.
            const float sp = std::max(u.length(), 1e-3f);
            const QVector3D mid = p + u / sp * float(h * 0.5);
            const QVector3D um = m.velocity(mid);
            p += um / std::max(um.length(), 1e-3f) * float(h);
        }
        if (n < 2) {
            out.positions.resize(out.positions.size() - 3 * n);
            vmag.resize(vmag.size() - n);
            cp0.resize(cp0.size() - n);
            continue;
        }
        out.indices.push_back(std::uint32_t(out.positions.size() / 3));
    }
    out.fields.push_back(makeField("velocity_magnitude", "m/s", 1, vmag));
    out.fields.push_back(makeField("total_pressure_coefficient", "1", 1, cp0));
    return out;
}

MeshData carIsoSurface(const QString& field, double value, const CarParams& prm)
{
    const FlowModel& m = model(prm.speed);
    const double h = 0.03 / std::sqrt(double(prm.detail));
    const QVector3D lo(-1.9f, -1.15f, 0.0f), hi(5.0f, 1.15f, 1.55f);
    const int nx = int((hi.x() - lo.x()) / h), ny = int((hi.y() - lo.y()) / h), nz = int((hi.z() - lo.z()) / h);
    auto idx = [&](int i, int j, int k) { return (std::size_t(k) * std::size_t(ny + 1) + std::size_t(j)) * std::size_t(nx + 1) + std::size_t(i); };
    auto pos = [&](int i, int j, int k) { return QVector3D(lo.x() + float(i * h), lo.y() + float(j * h), lo.z() + float(k * h)); };
    std::vector<float> f(std::size_t(nx + 1) * std::size_t(ny + 1) * std::size_t(nz + 1));
    const bool cp0 = field == "total_pressure_coefficient";
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j <= ny; ++j)
            for (int i = 0; i <= nx; ++i) {
                const QVector3D p = pos(i, j, k);
                double v;
                if (insideSolid(p)) v = cp0 ? 1.0 : prm.speed;
                else if (cp0) v = m.totalPressureCoefficient(p);
                else v = m.velocity(p).length();
                f[idx(i, j, k)] = float(v);
            }
    MeshData out;
    std::unordered_map<std::uint64_t, std::uint32_t> edgeVertex;
    std::vector<float> vmag;
    auto vertexOn = [&](std::size_t a, std::size_t b, const QVector3D& pa, const QVector3D& pb) {
        const std::uint64_t key = (std::uint64_t(std::min(a, b)) << 32) ^ std::uint64_t(std::max(a, b));
        if (auto it = edgeVertex.find(key); it != edgeVertex.end()) return it->second;
        const double fa = f[a], fb = f[b];
        const double t = std::clamp((value - fa) / (fb - fa), 0.0, 1.0);
        const QVector3D p = pa + (pb - pa) * float(t);
        const std::uint32_t id = std::uint32_t(out.positions.size() / 3);
        out.positions.insert(out.positions.end(), {p.x(), p.y(), p.z()});
        vmag.push_back(m.velocity(p).length());
        edgeVertex.emplace(key, id);
        return id;
    };
    // Six tetrahedra per cube sharing the main diagonal 0-6.
    static const int cube[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    static const int tets[6][4] = {{0,5,1,6},{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6}};
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i) {
                std::size_t ci[8];
                QVector3D cp[8];
                bool any = false, all = true;
                for (int c = 0; c < 8; ++c) {
                    ci[c] = idx(i + cube[c][0], j + cube[c][1], k + cube[c][2]);
                    cp[c] = pos(i + cube[c][0], j + cube[c][1], k + cube[c][2]);
                    const bool below = f[ci[c]] < value;
                    any |= below;
                    all &= below;
                }
                if (!any || all) continue;
                for (const auto& t : tets) {
                    int in[4], out4[4], ni = 0, no = 0;
                    for (int v = 0; v < 4; ++v) (f[ci[t[v]]] < value ? in[ni++] : out4[no++]) = t[v];
                    if (ni == 0 || ni == 4) continue;
                    auto V = [&](int a, int b) { return vertexOn(ci[a], ci[b], cp[a], cp[b]); };
                    if (ni == 1 || ni == 3) {
                        const int apex = ni == 1 ? in[0] : out4[0];
                        const int* others = ni == 1 ? out4 : in;
                        out.indices.insert(out.indices.end(), {V(apex, others[0]), V(apex, others[1]), V(apex, others[2])});
                    } else {
                        const std::uint32_t a = V(in[0], out4[0]), b = V(in[0], out4[1]), c = V(in[1], out4[0]), d = V(in[1], out4[1]);
                        out.indices.insert(out.indices.end(), {a, b, c, b, d, c});
                    }
                }
            }
    out.fields.push_back(makeField("velocity_magnitude", "m/s", 1, vmag));
    return out;
}

ProbeResult carLineProbe(const QVector3D& from, const QVector3D& to, int samples, const QStringList& fields, const CarParams& prm)
{
    const FlowModel& m = model(prm.speed);
    ProbeResult r;
    samples = std::max(2, samples);
    for (const QString& f : fields) {
        QString unit = f == "velocity_magnitude" ? "m/s" : f == "pressure" ? "Pa" : "1";
        r.series.push_back({f, unit, {}});
    }
    const double len = double((to - from).length());
    for (int i = 0; i < samples; ++i) {
        const double t = double(i) / (samples - 1);
        const QVector3D p = from + (to - from) * float(t);
        r.distance.push_back(t * len);
        r.points.push_back(p);
        const bool solid = insideSolid(p);
        const QVector3D u = solid ? QVector3D() : m.velocity(p);
        const double c0 = m.totalPressureCoefficient(p);
        const double cp = c0 - double(u.lengthSquared()) / (prm.speed * prm.speed);
        for (auto& s : r.series) {
            double v = std::nan("");
            if (!solid) {
                if (s.name == "velocity_magnitude") v = u.length();
                else if (s.name == "pressure_coefficient") v = cp;
                else if (s.name == "total_pressure_coefficient") v = c0;
                else if (s.name == "pressure") v = cp * 0.5 * prm.density * prm.speed * prm.speed;
            }
            s.y.push_back(v);
        }
    }
    return r;
}

std::vector<ZoneForce> carForces(const CarParams& prm)
{
    std::vector<ZoneForce> out;
    for (const CarZone& z : carZones()) {
        if (z.type != "wall" || z.name == "ground") continue;
        CarParams coarse = prm;
        coarse.detail = 1;
        const MeshData m = carZoneMesh(z.name, coarse, false);
        const auto* p = m.field("pressure");
        const auto* tau = m.field("wall_shear");
        ZoneForce f{z.name, {}, {}, {}, {}, 0};
        for (std::size_t t = 0; t + 2 < m.indices.size(); t += 3) {
            const std::uint32_t a = m.indices[t], b = m.indices[t + 1], c = m.indices[t + 2];
            const QVector3D pa(m.positions[3 * a], m.positions[3 * a + 1], m.positions[3 * a + 2]);
            const QVector3D pb(m.positions[3 * b], m.positions[3 * b + 1], m.positions[3 * b + 2]);
            const QVector3D pc(m.positions[3 * c], m.positions[3 * c + 1], m.positions[3 * c + 2]);
            const QVector3D areaVec = QVector3D::crossProduct(pb - pa, pc - pa) * 0.5f;  // outward (into fluid)
            const float pm = (p->values[a] + p->values[b] + p->values[c]) / 3.0f;
            const QVector3D centre = (pa + pb + pc) / 3.0f;
            const QVector3D dp = -areaVec * pm;
            f.pressure += dp;
            f.pressureMoment += QVector3D::crossProduct(centre, dp);
            QVector3D tm;
            for (int k = 0; k < 3; ++k)
                tm[k] = (tau->values[3 * a + std::size_t(k)] + tau->values[3 * b + std::size_t(k)] + tau->values[3 * c + std::size_t(k)]) / 3.0f;
            const QVector3D dv = tm * areaVec.length();
            f.viscous += dv;
            f.viscousMoment += QVector3D::crossProduct(centre, dv);
            f.area += areaVec.length();
        }
        out.push_back(f);
    }
    return out;
}

QJsonObject carInspect(const CarParams& prm, bool withData, const QJsonArray& rawSettings)
{
    auto num = [](double v) { return QJsonValue(v); };
    QJsonArray zones;
    for (const CarZone& z : carZones()) {
        QJsonObject cond{{"velocity", QJsonArray{QJsonValue(), QJsonValue(), QJsonValue()}},
                         {"velocity_magnitude", z.type == "velocity-inlet" ? num(prm.speed) : QJsonValue()},
                         {"direction", z.type == "velocity-inlet" ? QJsonArray{1, 0, 0} : QJsonArray{QJsonValue(), QJsonValue(), QJsonValue()}},
                         {"pressure", z.type == "pressure-outlet" || z.type == "velocity-inlet" ? num(0) : QJsonValue()},
                         {"roughness_height", z.type == "wall" ? num(0) : QJsonValue()},
                         {"wall_motion", z.name == "ground" ? QJsonValue("moving") : z.name.startsWith("tyre") ? QJsonValue("rotating") : z.type == "wall" ? QJsonValue("stationary") : QJsonValue()}};
        // Invented face counts in proportion to a ~18M-cell car mesh.
        const qint64 faces = z.name == "monocoque" ? 402118 : z.name == "sidepods" ? 211040 : z.name == "front-wing" ? 611904
                           : z.name == "rear-wing" ? 744230 : z.name.startsWith("tyre") ? 96412 : z.name == "ground" ? 210544 : 14400;
        zones.append(QJsonObject{{"id", z.id}, {"name", z.name}, {"type", z.type}, {"faces", faces}, {"condition", cond}});
    }
    auto setting = [](const QString& path, const QString& label, const QString& group, QJsonValue value, const QString& unit,
                      const QString& kind, QJsonValue min, QJsonValue max, bool editable, QJsonArray options = {}) {
        return QJsonObject{{"path", path}, {"label", label}, {"group", group}, {"value", value}, {"unit", unit}, {"kind", kind},
                           {"options", options}, {"min", min}, {"max", max}, {"editable", editable},
                           {"reason", editable ? "" : "not yet configurable in Kutta"}};
    };
    const double fmin = 1.1754943508222875e-38;
    QJsonArray s;
    // Canonical controls, mirroring the real CLI's paths and editability.
    s.append(setting("boundary/inlet/velocity-magnitude", "Velocity magnitude", "Boundary Conditions", prm.speed, "m/s", "number", 0, 1e15, true));
    s.append(setting("boundary/inlet/direction", "Flow direction", "Boundary Conditions", QJsonArray{1, 0, 0}, "1", "vec3", QJsonValue(), QJsonValue(), true));
    s.append(setting("boundary/inlet/turb-intensity", "Turbulent intensity", "Boundary Conditions", 0.01, "1", "number", 0, 1, true));
    s.append(setting("boundary/inlet/turb-viscosity-ratio", "Turbulent viscosity ratio", "Boundary Conditions", 10, "1", "number", fmin, 1e30, true));
    s.append(setting("boundary/inlet/gauge-pressure", "Gauge pressure", "Boundary Conditions", 0, "Pa", "number", QJsonValue(), QJsonValue(), false));
    s.append(setting("boundary/outlet/gauge-pressure", "Gauge pressure", "Boundary Conditions", 0, "Pa", "number", -1e30, 1e30, true));
    s.append(setting("boundary/outlet/turb-intensity", "Backflow turbulent intensity", "Boundary Conditions", 0.01, "1", "number", 0, 1, true));
    s.append(setting("boundary/outlet/turb-viscosity-ratio", "Backflow turbulent viscosity ratio", "Boundary Conditions", 10, "1", "number", fmin, 1e30, true));
    for (const CarZone& z : carZones()) {
        if (z.type != "wall") continue;
        const QString motion = z.name == "ground" ? "moving" : z.name.startsWith("tyre") ? "rotating" : "stationary";
        s.append(setting("boundary/" + z.name + "/wall-motion", "Wall motion", "Boundary Conditions", motion, "", "string", QJsonValue(), QJsonValue(), false));
        s.append(setting("boundary/" + z.name + "/shear-condition", "Shear condition", "Boundary Conditions", "no-slip", "", "string", QJsonValue(), QJsonValue(), false));
        s.append(setting("boundary/" + z.name + "/roughness-height", "Roughness height", "Boundary Conditions", 0, "m", "number", QJsonValue(), QJsonValue(), false));
        if (z.name == "ground")
            s.append(setting("boundary/ground/velocity", "Wall velocity", "Boundary Conditions", QJsonArray{prm.speed, 0, 0}, "m/s", "vec3", QJsonValue(), QJsonValue(), false));
        if (z.name.startsWith("tyre"))
            s.append(setting("boundary/" + z.name + "/rotation-speed", "Rotational speed", "Boundary Conditions", prm.speed / kTyreR, "rad/s", "number", QJsonValue(), QJsonValue(), false));
    }
    s.append(setting("boundary/tunnel/type", "Zone type", "Boundary Conditions", "symmetry", "", "string", QJsonValue(), QJsonValue(), false));
    s.append(setting("cell-zone/fluid/material", "Material", "Cell Zone Conditions", "air", "", "string", QJsonValue(), QJsonValue(), false));
    s.append(setting("cell-zone/fluid/frame-motion", "Frame motion", "Cell Zone Conditions", false, "", "bool", QJsonValue(), QJsonValue(), false));
    s.append(setting("general/solver-type", "Solver type", "General", "pressure-based", "", "string", QJsonValue(), QJsonValue(), false));
    s.append(setting("general/time", "Time", "General", "steady", "", "string", QJsonValue(), QJsonValue(), false));
    s.append(setting("general/gravity", "Gravity", "General", false, "", "bool", QJsonValue(), QJsonValue(), false));
    s.append(setting("general/operating-pressure", "Operating pressure", "General", 101325, "Pa", "number", QJsonValue(), QJsonValue(), false));
    s.append(setting("models/viscous", "Viscous model", "Models", "k-omega-sst", "", "string", QJsonValue(), QJsonValue(), false));
    s.append(setting("models/energy", "Energy", "Models", false, "", "bool", QJsonValue(), QJsonValue(), false));
    s.append(setting("models/curvature-correction", "Curvature correction", "Models", false, "", "bool", QJsonValue(), QJsonValue(), false));
    s.append(setting("models/production-limiter", "Production limiter", "Models", true, "", "bool", QJsonValue(), QJsonValue(), false));
    s.append(setting("materials/air/density", "Density", "Materials", prm.density, "kg/m^3", "number", fmin, 1e30, true));
    s.append(setting("materials/air/viscosity", "Dynamic viscosity", "Materials", prm.viscosity, "Pa s", "number", fmin, 1e30, true));
    s.append(setting("reference/area", "Reference area", "Reference Values", prm.refArea, "m^2", "number", fmin, 1e30, true));
    s.append(setting("reference/density", "Reference density", "Reference Values", prm.refDensity, "kg/m^3", "number", fmin, 1e30, true));
    s.append(setting("reference/velocity", "Reference velocity", "Reference Values", prm.refVelocity, "m/s", "number", fmin, 1e30, true));
    s.append(setting("reference/pressure", "Reference pressure", "Reference Values", prm.refPressure, "Pa", "number", -1e30, 1e30, true));
    s.append(setting("reference/length", "Reference length", "Reference Values", prm.refLength, "m", "number", fmin, 1e30, true));
    s.append(setting("reference/compute-from", "Compute from", "Reference Values", "inlet", "", "string", QJsonValue(), QJsonValue(), false));
    const QJsonArray schemes{"first-order-upwind", "second-order-upwind"};
    // Fluent's Methods panel, with the real CLI's paths, labels and activity.
    auto labelled = [](QJsonObject o, const QStringList& labels) {
        o["option_labels"] = QJsonArray::fromStringList(labels);
        const QJsonArray opts = o.value("options").toArray();
        for (qsizetype i = 0; i < opts.size() && i < labels.size(); ++i)
            if (opts[i] == o.value("value")) o["value_label"] = labels[i];
        return o;
    };
    const QStringList schemeLabels{"First Order Upwind", "Second Order Upwind"};
    s.append(labelled(setting("methods/pressure-velocity-coupling", "Pressure-velocity coupling", "Methods", "coupled", "", "enum", QJsonValue(), QJsonValue(), true,
                              QJsonArray{"simple", "simplec", "coupled"}), {"SIMPLE", "SIMPLEC", "Coupled"}));
    s.append(setting("methods/pseudo-transient", "Pseudo-transient formulation", "Methods", true, "", "bool", QJsonValue(), QJsonValue(), false));
    s.append(labelled(setting("methods/gradient", "Gradient", "Methods", "least-squares-cell-based", "", "enum", QJsonValue(), QJsonValue(), false,
                              QJsonArray{"least-squares-cell-based", "green-gauss-cell-based", "green-gauss-node-based"}),
                      {"Least Squares Cell Based", "Green-Gauss Cell Based", "Green-Gauss Node Based"}));
    s.append(labelled(setting("methods/pressure-scheme", "Pressure", "Methods", "second-order", "", "enum", QJsonValue(), QJsonValue(), false,
                              QJsonArray{"standard", "second-order"}), {"Standard", "Second Order"}));
    s.append(labelled(setting("methods/momentum-scheme", "Momentum", "Methods", "second-order-upwind", "", "enum", QJsonValue(), QJsonValue(), true, schemes), schemeLabels));
    s.append(labelled(setting("methods/k-scheme", "Turbulent kinetic energy", "Methods", "second-order-upwind", "", "enum", QJsonValue(), QJsonValue(), true, schemes), schemeLabels));
    s.append(labelled(setting("methods/omega-scheme", "Specific dissipation rate", "Methods", "second-order-upwind", "", "enum", QJsonValue(), QJsonValue(), true, schemes), schemeLabels));
    QJsonObject epsilonScheme = labelled(setting("methods/epsilon-scheme", "Turbulent dissipation rate", "Methods", QJsonValue(), "", "enum", QJsonValue(), QJsonValue(), false, schemes), schemeLabels);
    epsilonScheme["active"] = false;
    epsilonScheme["reason"] = "not solved by the case's viscous model";
    s.append(epsilonScheme);
    s.append(setting("methods/gradient-limiting", "Gradient limiting", "Methods", true, "", "bool", QJsonValue(), QJsonValue(), true));
    s.append(setting("controls/pseudo-time-scale-factor", "Pseudo-time scale factor", "Controls", 1, "1", "number", fmin, 1e15, true));
    struct Relax { const char* key; const char* label; double value; };
    for (const Relax& r : {Relax{"mom", "Momentum relaxation", 0.5}, Relax{"pressure", "Pressure relaxation", 0.5},
                           Relax{"k", "Turbulent kinetic energy relaxation", 0.75}, Relax{"omega", "Specific dissipation rate relaxation", 0.75}})
        s.append(setting(QString("controls/%1-relaxation").arg(r.key), r.label, "Controls", r.value, "1", "number", fmin, 1, true));
    // Convergence criteria, as the real CLI lists them. Residual values stand
    // for the case's saved Check Convergence and Absolute Criteria.
    s.append(setting("monitors/check-convergence", "Stop when converged", "Monitors", false, "", "bool", QJsonValue(), QJsonValue(), true));
    QJsonObject condition = setting("monitors/condition", "Convergence condition", "Monitors", "all", "", "enum", QJsonValue(), QJsonValue(), true,
                                    QJsonArray{"all", "any"});
    condition["option_labels"] = QJsonArray{"All conditions are met", "Any condition is met"};
    condition["value_label"] = "All conditions are met";
    s.append(condition);
    for (const char* eq : {"continuity", "x-velocity", "y-velocity", "z-velocity", "k", "omega"}) {
        s.append(setting(QStringLiteral("monitors/residual/%1/check").arg(eq), "Check convergence", "Monitors", true, "", "bool", QJsonValue(), QJsonValue(), true));
        s.append(setting(QStringLiteral("monitors/residual/%1/absolute-criterion").arg(eq), "Absolute criteria", "Monitors", 1e-3, "1", "number", 0, 1e30, true));
    }
    for (const char* r : {"downforce", "drag", "front-wing-downforce", "rear-wing-downforce"}) {
        s.append(setting(QStringLiteral("monitors/report/%1/active").arg(r), "Active", "Monitors", false, "", "bool", QJsonValue(), QJsonValue(), true));
        s.append(setting(QStringLiteral("monitors/report/%1/stop-criterion").arg(r), "Stop criterion", "Monitors", 1e-3, "1", "number", fmin, 1e30, true));
        s.append(setting(QStringLiteral("monitors/report/%1/previous-values").arg(r), "Use previous values", "Monitors", 5, "", "integer", 1, 100000, true));
        s.append(setting(QStringLiteral("monitors/report/%1/ignore-before").arg(r), "Ignore iterations before", "Monitors", 100, "", "integer", 0, 1e9, true));
    }
    s.append(setting("initialization/method", "Initialization method", "Initialization", "hybrid", "", "string", QJsonValue(), QJsonValue(), false));
    s.append(setting("initialization/velocity-x", "X velocity", "Initialization", prm.speed, "m/s", "number", -1e15, 1e15, !withData));
    s.append(setting("run/iterations", "Number of iterations", "Run Calculation", 1000, "", "integer", 1, 1000000000, true));
    s.append(setting("run/monitor-interval", "Monitor interval", "Run Calculation", 10, "", "integer", 1, 100000, true));
    s.append(setting("reports/downforce/direction", "Force vector", "Report Definitions", QJsonArray{0, 0, -1}, "1", "vec3", QJsonValue(), QJsonValue(), true));
    s.append(setting("reports/drag/direction", "Force vector", "Report Definitions", QJsonArray{1, 0, 0}, "1", "vec3", QJsonValue(), QJsonValue(), true));
    for (const QJsonValue& r : rawSettings) s.append(r);

    QJsonArray reports{
        QJsonObject{{"name", "downforce"}, {"kind", "lift"}, {"zones", QJsonArray{"monocoque", "sidepods", "front-wing", "rear-wing", "tyre-front", "tyre-rear"}}, {"direction", QJsonArray{0, 0, -1}}},
        QJsonObject{{"name", "drag"}, {"kind", "drag"}, {"zones", QJsonArray{"monocoque", "sidepods", "front-wing", "rear-wing", "tyre-front", "tyre-rear"}}, {"direction", QJsonArray{1, 0, 0}}},
        QJsonObject{{"name", "front-wing-downforce"}, {"kind", "lift"}, {"zones", QJsonArray{"front-wing"}}, {"direction", QJsonArray{0, 0, -1}}},
        QJsonObject{{"name", "rear-wing-downforce"}, {"kind", "lift"}, {"zones", QJsonArray{"rear-wing"}}, {"direction", QJsonArray{0, 0, -1}}},
    };
    QJsonArray issues{QJsonObject{{"path", "rp/recon/cell-lsf?"}, {"setting", "recon/cell-lsf?"}, {"label", "Warped-face gradient correction"},
                                  {"class", "approximable"}, {"approximation", "wfgc"},
                                  {"explanation", "Fluent corrects gradients on warped faces with an undocumented algorithm. Kutta can run this case only with the correction off, which changes the converged answer. SYNTHETIC FIXTURE: copied from the real CLI's wording for demonstration."},
                                  {"evidence", "src/kutta/fv/VERIFICATION.md#v6-initialization-comparison"}},
                      QJsonObject{{"path", "rp/ke-rough-wall-treatment-r14?"}, {"setting", "rp/ke-rough-wall-treatment-r14?"}, {"label", "Wall roughness"},
                                  {"class", "informational"},
                                  {"explanation", "Kutta implements the smooth-wall SST law. Every wall in this case has zero roughness, so this has no effect."}}};
    QJsonObject data{{"present", withData}, {"fields", withData ? QJsonArray{"SV_P", "SV_U", "SV_V", "SV_W", "SV_K", "SV_O", "SV_WALL_SHEAR", "SV_WALL_YPLUS"} : QJsonArray{}},
                     {"iterations_hint", withData ? QJsonValue(1200) : QJsonValue()}};
    return QJsonObject{
        {"v", 1}, {"synthetic", true},
        {"file", QJsonObject{{"path", "fsae_car.cas.h5"}, {"sha256", "synthetic-fixture"}, {"fluent_version", "24.2"}}},
        {"mesh", QJsonObject{{"cells", 18412903}, {"faces", 61220418}, {"nodes", 15904211},
                             {"cell_types", QJsonObject{{"tet", 0}, {"hex", 2204119}, {"poly", 16208784}}},
                             {"bounds", QJsonArray{QJsonArray{-6, -4, 0}, QJsonArray{14, 4, 4}}}}},
        {"zones", zones},
        {"models", QJsonObject{{"viscous", "k-omega-sst"}, {"energy", false}, {"energy_coupled", false}}},
        {"materials", QJsonArray{QJsonObject{{"name", "air"}, {"density", prm.density}, {"viscosity", prm.viscosity}}}},
        {"reports", reports},
        {"support", QJsonObject{{"status", "approximate"}, {"issues", issues}}},
        {"settings", s},
        {"data", data}};
}

} // namespace fixture
