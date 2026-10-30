#include "render/SceneView.h"

#include "core/ChartScale.h"

#include <QFontDatabase>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QTimer>
#include <QWheelEvent>

#include <algorithm>
#include <numeric>

namespace kapp {

namespace {
QColor kAccent(56, 204, 110);  // light-sheet green: focus, cursor, probe
}

class ViewOverlay : public QWidget {
public:
    explicit ViewOverlay(SceneView* parent) : QWidget(parent), m_view(parent)
    {
        setAttribute(Qt::WA_TransparentForMouseEvents);
        setAttribute(Qt::WA_NoSystemBackground);
        setAttribute(Qt::WA_TranslucentBackground);
    }
    QColor accent = kAccent;

protected:
    void paintEvent(QPaintEvent*) override
    {
        QPainter p(this);
        paintOverlay(p, rect(), m_view->overlayState(), accent, false);
    }

private:
    SceneView* m_view;
};

SceneView::SceneView(QWidget* parent) : QRhiWidget(parent)
{
    setSampleCount(4);
    setFocusPolicy(Qt::StrongFocus);
    setMouseTracking(false);
    setMinimumSize(160, 120);
    m_overlay = new ViewOverlay(this);
    m_overlay->setGeometry(rect());
    setAccessibleName(tr("Graphics viewport"));
}

SceneView::~SceneView() = default;

void SceneView::setScene(std::shared_ptr<const SceneData> scene)
{
    m_renderer.setScene(std::move(scene));
    ensureCamera();
    update();
}

void SceneView::setDisplay(const ViewDisplay& d)
{
    m_display = d;
    ensureCamera();
    m_renderer.setDisplay(m_display);
    update();
}

void SceneView::ensureCamera()
{
    if (m_display.cameraInitialized || !scene() || !scene()->bounds.valid) return;
    m_display.camera.fit(scene()->bounds, aspect());
    m_display.cameraInitialized = true;
    m_renderer.setDisplay(m_display);
}

void SceneView::setOverlay(const OverlayState& o)
{
    m_overlayState = o;
    m_overlayState.view = m_display.camera.view();
    m_overlay->update();
}

void SceneView::setTheme(bool dark, QColor clear, QColor accent)
{
    m_dark = dark;
    m_renderer.setClearColor(clear);
    m_renderer.setMissingGrey(dark ? 0.42f : 0.62f);
    m_overlay->accent = accent;
    m_overlayState.dark = dark;
    update();
    m_overlay->update();
}

float SceneView::aspect() const { return height() > 0 ? float(width()) / float(height()) : 1.0f; }

QPointF SceneView::toNdc(QPointF p) const
{
    return {2.0 * p.x() / std::max(1, width()) - 1.0, 1.0 - 2.0 * p.y() / std::max(1, height())};
}

void SceneView::cameraEdited()
{
    m_renderer.setDisplay(m_display);
    m_overlayState.view = m_display.camera.view();
    // The probe stays on its world point as the camera moves.
    if (m_overlayState.probe && m_overlayState.probe->hit) m_overlayState.probePixel = projectToLogical(m_overlayState.probe->world);
    m_overlay->update();
    update();
    emit cameraChanged(m_display.camera);
}

QPointF SceneView::projectToLogical(const QVector3D& world) const
{
    const qreal dpr = devicePixelRatio();
    const QVector3D p = m_renderer.projectToPixel(world, (QSizeF(size()) * dpr).toSize());
    return {p.x() / dpr, p.y() / dpr};
}

Bounds SceneView::visibleBounds() const
{
    Bounds b;
    if (!scene()) return b;
    for (const SceneItem& it : scene()->items) {
        const bool shown = m_display.visibleInitialized ? m_display.visible.contains(it.id) : it.defaultVisible;
        if (shown && it.mesh) b.expand(it.mesh->bounds());
    }
    return b.valid ? b : scene()->bounds;
}

void SceneView::fitView()
{
    if (!scene()) return;
    // Keep the title (top), legend (left) and triad (bottom) clear of the geometry.
    const double w = std::max(1, width()), h = std::max(1, height());
    const QMarginsF reserve(m_display.contours ? 118 / w : 16 / w, 58 / h, 16 / w, 30 / h);
    m_display.camera.fit(visibleBounds(), aspect(), reserve);
    if (qEnvironmentVariableIsSet("KUTTA_DEBUG_CAMERA")) {
        const Bounds b = visibleBounds();
        qDebug("fit: size %dx%d aspect %.3f bounds (%.3f %.3f %.3f)-(%.3f %.3f %.3f) dist %.4f yaw %.1f pitch %.1f fov %.1f", width(), height(), aspect(),
               b.min.x(), b.min.y(), b.min.z(), b.max.x(), b.max.y(), b.max.z(), m_display.camera.distance, m_display.camera.yaw, m_display.camera.pitch, m_display.camera.fovY);
    }
    m_display.cameraInitialized = true;
    cameraEdited();
}

void SceneView::applyPreset(Camera::Preset p)
{
    m_display.camera.preset(p);
    fitView();
}

void SceneView::setCamera(const Camera& c)
{
    m_display.camera = c;
    m_display.cameraInitialized = true;
    cameraEdited();
}

void SceneView::initialize(QRhiCommandBuffer*)
{
    m_renderer.initialize(rhi(), renderTarget()->renderPassDescriptor(), sampleCount());
}

void SceneView::releaseResources() { m_renderer.releaseResources(); }

void SceneView::render(QRhiCommandBuffer* cb)
{
    const bool benchmarking = m_benchRemaining > 0;
    if (benchmarking && m_benchClock.isValid()) m_benchTimes.push_back(m_benchClock.nsecsElapsed() / 1e6);
    if (benchmarking) m_benchClock.restart();
    m_renderer.render(cb, renderTarget(), renderTarget()->pixelSize());
    ++m_framesRendered;
    // Readbacks complete on a later frame; keep frames coming only while one is
    // outstanding. Otherwise nothing schedules another frame: no idle redraw.
    if (m_renderer.pickPending()) update();
    if (benchmarking) {
        if (--m_benchRemaining > 0) {
            update();
        } else {
            FrameTimes t;
            auto& v = m_benchTimes;
            if (!v.empty()) v.erase(v.begin());  // first interval includes warm-up
            t.frames = int(v.size());
            if (!v.empty()) {
                t.meanMs = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
                std::vector<double> s = v;
                std::sort(s.begin(), s.end());
                t.p95Ms = s[std::min(s.size() - 1, std::size_t(double(s.size()) * 0.95))];
                t.maxMs = s.back();
            }
            t.triangles = m_renderer.trianglesDrawn();
            m_benchClock.invalidate();
            if (auto done = std::move(m_benchDone)) QTimer::singleShot(0, this, [done, t] { done(t); });
        }
    }
}

void SceneView::measureFrames(int frames, std::function<void(FrameTimes)> done)
{
    m_benchRemaining = frames + 1;
    m_benchTimes.clear();
    m_benchClock.invalidate();
    m_benchDone = std::move(done);
    update();
}

QImage SceneView::renderImage(bool footer)
{
    QImage img = grabFramebuffer();
    img.setDevicePixelRatio(devicePixelRatio());
    if (img.isNull()) return img;
    QPainter p(&img);
    p.setRenderHint(QPainter::Antialiasing);
    paintOverlay(p, rect(), m_overlayState, m_overlay->accent, footer);
    return img;
}

void SceneView::resizeEvent(QResizeEvent* e)
{
    QRhiWidget::resizeEvent(e);
    m_overlay->setGeometry(rect());
}

void SceneView::focusInEvent(QFocusEvent* e)
{
    QRhiWidget::focusInEvent(e);
    emit activated(this);
}

void SceneView::mousePressEvent(QMouseEvent* e)
{
    emit activated(this);
    setFocus(Qt::MouseFocusReason);
    m_pressPos = m_lastPos = e->position().toPoint();
    m_dragButton = e->button();
    m_dragged = false;
}

void SceneView::mouseMoveEvent(QMouseEvent* e)
{
    if (m_dragButton == Qt::NoButton) return;
    const QPoint pos = e->position().toPoint();
    const QPoint d = pos - m_lastPos;
    m_lastPos = pos;
    if ((pos - m_pressPos).manhattanLength() > 3) m_dragged = true;
    if (!m_dragged) return;
    const bool pan = m_dragButton == Qt::MiddleButton || (m_dragButton == Qt::LeftButton && (e->modifiers() & Qt::ShiftModifier));
    if (pan) {
        m_display.camera.pan(QPointF(2.0 * d.x() / width(), -2.0 * d.y() / height()), aspect());
    } else if (m_dragButton == Qt::RightButton) {
        m_display.camera.zoomAt(toNdc(m_pressPos), aspect(), std::pow(1.006f, float(d.y())));
    } else {
        m_display.camera.orbit(-0.35f * d.x(), 0.35f * d.y());
    }
    cameraEdited();
}

void SceneView::mouseReleaseEvent(QMouseEvent* e)
{
    const bool click = !m_dragged && e->button() == Qt::LeftButton;
    m_dragButton = Qt::NoButton;
    if (!click || !scene()) return;
    // Probe: GPU pick at the clicked pixel.
    const QPointF logical = e->position();
    const QPointF fb = logical * devicePixelRatio();
    m_renderer.requestPick(fb, [this, logical](const PickResult& r) { emit picked(r, logical); });
    update();
}

void SceneView::wheelEvent(QWheelEvent* e)
{
    const double steps = e->angleDelta().y() / 120.0 + e->pixelDelta().y() / 240.0;
    if (steps == 0) return;
    m_display.camera.zoomAt(toNdc(e->position()), aspect(), float(std::pow(0.88, steps)));
    cameraEdited();
}

void SceneView::keyPressEvent(QKeyEvent* e)
{
    switch (e->key()) {
    case Qt::Key_F: fitView(); return;
    case Qt::Key_1: applyPreset(Camera::Preset::Front); return;
    case Qt::Key_2: applyPreset(Camera::Preset::Back); return;
    case Qt::Key_3: applyPreset(Camera::Preset::Left); return;
    case Qt::Key_4: applyPreset(Camera::Preset::Right); return;
    case Qt::Key_5: applyPreset(Camera::Preset::Top); return;
    case Qt::Key_6: applyPreset(Camera::Preset::Bottom); return;
    case Qt::Key_7:
    case Qt::Key_0: applyPreset(Camera::Preset::Isometric); return;
    case Qt::Key_Left: m_display.camera.orbit(10, 0); cameraEdited(); return;
    case Qt::Key_Right: m_display.camera.orbit(-10, 0); cameraEdited(); return;
    case Qt::Key_Up: m_display.camera.orbit(0, 10); cameraEdited(); return;
    case Qt::Key_Down: m_display.camera.orbit(0, -10); cameraEdited(); return;
    case Qt::Key_Plus:
    case Qt::Key_Equal: m_display.camera.zoomAt({0, 0}, aspect(), 0.85f); cameraEdited(); return;
    case Qt::Key_Minus: m_display.camera.zoomAt({0, 0}, aspect(), 1.0f / 0.85f); cameraEdited(); return;
    default: QRhiWidget::keyPressEvent(e);
    }
}

// ---- overlay painting ------------------------------------------------------------------

namespace {

void badge(QPainter& p, QRectF r, const QString& text, QColor fill, QColor fg, QColor outline)
{
    QPainterPath path;
    path.addRoundedRect(r, 3, 3);
    if (fill.alpha()) p.fillPath(path, fill);
    if (outline.alpha()) {
        p.setPen(QPen(outline, 1));
        p.drawPath(path);
    }
    p.setPen(fg);
    p.drawText(r, Qt::AlignCenter, text);
}

} // namespace

namespace {
// Text over geometry: a soft halo in the background colour keeps labels
// readable on any field colour (white deltas, yellow maxima).
void haloText(QPainter& p, const QRect& rect, int flags, const QString& text, const QColor& fg, bool dark)
{
    const QColor halo = dark ? QColor(16, 18, 20, 210) : QColor(255, 255, 255, 220);
    p.setPen(halo);
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            if (dx || dy) p.drawText(rect.translated(dx, dy), flags, text);
    p.setPen(fg);
    p.drawText(rect, flags, text);
}
} // namespace

void paintOverlay(QPainter& p, const QRect& r, const OverlayState& s, QColor accent, bool forExport)
{
    p.save();
    p.setRenderHint(QPainter::Antialiasing);
    p.setRenderHint(QPainter::TextAntialiasing);
    const QColor fg = s.dark ? QColor(226, 229, 233) : QColor(28, 31, 35);
    const QColor muted = s.dark ? QColor(150, 156, 164) : QColor(96, 102, 110);
    const QColor amber(232, 160, 32);
    QFont base = p.font();
    base.setPixelSize(12);
    QFont small = base;
    small.setPixelSize(11);
    QFont strong = base;
    strong.setWeight(QFont::DemiBold);
    strong.setPixelSize(13);
    QFont tabular = small;
    tabular.setFeature(QFont::Tag("tnum"), 1);

    const int pad = 12;
    // Badges, right-aligned. Approximate is filled amber and always first.
    p.setFont(QFont(small.family(), -1, QFont::DemiBold));
    QFont bf = small;
    bf.setWeight(QFont::DemiBold);
    bf.setPixelSize(10);
    bf.setLetterSpacing(QFont::AbsoluteSpacing, 0.6);
    p.setFont(bf);
    qreal x = r.right() - pad;
    for (const QString& b : s.badges) {
        const qreal w = p.fontMetrics().horizontalAdvance(b) + 14;
        const QRectF br(x - w, r.top() + 10, w, 18);
        if (b.startsWith("APPROXIMATE")) badge(p, br, b, amber, QColor(30, 22, 6), Qt::transparent);
        else if (b.startsWith("MODIFIED")) badge(p, br, b, Qt::transparent, amber, amber);
        else if (b.startsWith("LIVE")) badge(p, br, b, Qt::transparent, accent, accent);
        else if (b.startsWith("FAILED") || b.startsWith("BLOCKED")) badge(p, br, b, QColor(200, 64, 56), Qt::white, Qt::transparent);
        else badge(p, br, b, Qt::transparent, muted, muted);
        x -= w + 6;
    }

    const qreal titleRight = x;
    // Title and dataset.
    if (!s.title.isEmpty()) {
        p.setFont(strong);
        haloText(p, QRect(r.left() + pad, r.top() + 8, r.width() - 2 * pad, 18), Qt::AlignLeft | Qt::AlignVCenter,
                 p.fontMetrics().elidedText(s.title, Qt::ElideRight, int(titleRight) - r.left() - pad - 10), fg, s.dark);
        p.setFont(small);
        haloText(p, QRect(r.left() + pad, r.top() + 26, r.width() - 2 * pad, 16), Qt::AlignLeft | Qt::AlignVCenter,
                 p.fontMetrics().elidedText(s.subtitle, Qt::ElideRight, r.width() - 2 * pad - 200), muted, s.dark);
    }

    // Colour legend: vertical bar on the left, as in Fluent.
    if (s.legend && s.range.valid()) {
        const int top = r.top() + 56, barW = 12;
        const int h = std::clamp(int(r.height() * 0.42), 90, 320);
        const QRect bar(r.left() + pad, top, barW, h);
        const Colormap& cm = colormap(s.colormap);
        if (s.bands > 0) {
            for (int b = 0; b < s.bands; ++b) {
                const double t = (b + 0.5) / s.bands;
                const int y0 = bar.bottom() - int(std::round(double(b + 1) / s.bands * bar.height()));
                const int y1 = bar.bottom() - int(std::round(double(b) / s.bands * bar.height()));
                p.fillRect(QRect(bar.left(), y0 + 1, bar.width(), y1 - y0), cm.sample(t));
            }
        } else {
            for (int y = 0; y < bar.height(); ++y)
                p.fillRect(QRect(bar.left(), bar.bottom() - y, bar.width(), 1), cm.sample(double(y) / std::max(1, bar.height() - 1)));
        }
        p.setPen(QPen(muted, 1));
        p.setBrush(Qt::NoBrush);
        p.drawRect(QRectF(bar).adjusted(-0.5, -0.5, 0.5, 0.5));
        Axis ax;
        ax.setRange(s.range.min, s.range.max);
        ax.setPixels(bar.bottom(), bar.top());
        p.setFont(tabular);
        p.setPen(fg);
        std::vector<Tick> ticks = ax.ticks(std::max(3, h / 36));
        // Always label the ends so the range is explicit.
        const auto endLabel = [&](double v, int y) {
            haloText(p, QRect(bar.right() + 7, y - 8, 120, 16), Qt::AlignLeft | Qt::AlignVCenter, formatNumber(v, 4), fg, s.dark);
        };
        for (const Tick& t : ticks) {
            const int y = int(std::round(ax.toPixel(t.value)));
            if (std::abs(y - bar.top()) < 16 || std::abs(y - bar.bottom()) < 16) continue;
            p.setPen(muted);
            p.drawLine(bar.right() + 1, y, bar.right() + 4, y);
            haloText(p, QRect(bar.right() + 7, y - 8, 120, 16), Qt::AlignLeft | Qt::AlignVCenter, t.label, fg, s.dark);
        }
        endLabel(s.range.max, bar.top());
        endLabel(s.range.min, bar.bottom());
        if (!s.unit.isEmpty()) {
            p.setFont(small);
            p.setPen(muted);
            p.drawText(QRect(bar.left(), bar.bottom() + 6, 160, 16), Qt::AlignLeft | Qt::AlignVCenter, s.unit);
        }
    }

    // Axis triad, bottom-left (only with geometry to orient).
    if (s.triad) {
        const QPointF o(r.left() + 34, r.bottom() - 34);
        const float len = 20;
        const struct { QVector3D axis; QColor c; const char* label; } axes[] = {
            {{1, 0, 0}, QColor(222, 84, 74), "X"}, {{0, 1, 0}, QColor(92, 182, 96), "Y"}, {{0, 0, 1}, QColor(78, 140, 230), "Z"}};
        p.setFont(small);
        for (const auto& a : axes) {
            const QVector3D v = s.view.mapVector(a.axis);
            const QPointF end = o + QPointF(v.x(), -v.y()) * len;
            p.setPen(QPen(a.c, 2, Qt::SolidLine, Qt::RoundCap));
            p.drawLine(o, end);
            p.setPen(a.c);
            p.drawText(QRectF(end + QPointF(v.x(), -v.y()) * 7 - QPointF(6, 7), QSizeF(12, 14)), Qt::AlignCenter, a.label);
        }
    }

    // Probe crosshair and readout.
    if (s.probe && s.probe->hit) {
        const QPointF c = s.probePixel;
        p.setPen(QPen(accent, 1.5));
        p.drawLine(c + QPointF(-9, 0), c + QPointF(-3, 0));
        p.drawLine(c + QPointF(3, 0), c + QPointF(9, 0));
        p.drawLine(c + QPointF(0, -9), c + QPointF(0, -3));
        p.drawLine(c + QPointF(0, 3), c + QPointF(0, 9));
        p.setFont(tabular);
        const QStringList lines = s.probeText.split('\n');
        int w = 0;
        for (const QString& l : lines) w = std::max(w, p.fontMetrics().horizontalAdvance(l));
        QRectF box(c + QPointF(14, -10), QSizeF(w + 16, 16.0 * lines.size() + 8));
        if (box.right() > r.right() - 4) box.moveRight(c.x() - 14);
        if (box.bottom() > r.bottom() - 4) box.moveBottom(r.bottom() - 4);
        QPainterPath path;
        path.addRoundedRect(box, 3, 3);
        p.fillPath(path, s.dark ? QColor(16, 18, 21, 236) : QColor(255, 255, 255, 240));
        p.setPen(QPen(accent, 1));
        p.drawPath(path);
        p.setPen(fg);
        for (int i = 0; i < lines.size(); ++i)
            p.drawText(QRectF(box.left() + 8, box.top() + 4 + 16 * i, w + 4, 16), Qt::AlignLeft | Qt::AlignVCenter, lines[i]);
    }

    if (!s.emptyText.isEmpty()) {
        p.setFont(base);
        p.setPen(muted);
        p.drawText(r.adjusted(40, 40, -40, -40), Qt::AlignCenter | Qt::TextWordWrap, s.emptyText);
    }

    if (s.active) {
        p.setPen(QPen(accent, 2));
        p.setBrush(Qt::NoBrush);
        p.drawRect(QRectF(r).adjusted(1, 1, -1, -1));
    }

    // Exports carry the full provenance notice as a footer.
    if (forExport && !s.provenance.isEmpty()) {
        p.setFont(small);
        const QRect band(r.left(), r.bottom() - 44, r.width(), 44);
        p.fillRect(band, s.dark ? QColor(10, 11, 13, 230) : QColor(255, 248, 232, 240));
        p.setPen(s.dark ? amber : QColor(120, 78, 0));
        p.drawText(band.adjusted(pad, 4, -pad, -4), Qt::AlignLeft | Qt::AlignVCenter | Qt::TextWordWrap, s.provenance);
    }
    p.restore();
}

} // namespace kapp
