#include "ui/PlotWidget.h"

#include "ui/Theme.h"

#include <QApplication>
#include <QClipboard>
#include <QContextMenuEvent>
#include <QFile>
#include <QFileDialog>
#include <QMenu>
#include <QPainter>
#include <QPainterPath>
#include <QStandardPaths>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>

namespace kapp {

namespace {

constexpr double kHeader = 26;

QFont tabular(QFont f)
{
    f.setFeature(QFont::Tag("tnum"), 1);
    return f;
}

bool finite(double v) { return std::isfinite(v); }

// Index of the sample nearest to x in a sorted x array.
std::optional<std::size_t> nearestIndex(const std::vector<double>& xs, double x)
{
    if (xs.empty()) return std::nullopt;
    auto it = std::lower_bound(xs.begin(), xs.end(), x);
    if (it == xs.end()) return xs.size() - 1;
    if (it == xs.begin()) return 0;
    const std::size_t i = std::size_t(it - xs.begin());
    return std::abs(xs[i] - x) < std::abs(xs[i - 1] - x) ? i : i - 1;
}

} // namespace

PlotWidget::PlotWidget(QWidget* parent) : QWidget(parent)
{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    setMinimumSize(220, 140);
    setAttribute(Qt::WA_OpaquePaintEvent);
    connect(&Theme::instance(), &Theme::changed, this, qOverload<>(&QWidget::update));
}

void PlotWidget::setTitle(const QString& t) { m_title = t; setAccessibleName(t); update(); }
void PlotWidget::setAxisLabels(const QString& x, const QString& y) { m_xLabel = x; m_yLabel = y; update(); }
void PlotWidget::setBadges(const QStringList& b) { m_badges = b; update(); }
void PlotWidget::setEmptyText(const QString& t) { m_emptyText = t; update(); }

void PlotWidget::setYScale(Axis::Scale s)
{
    m_y.setScale(s);
    fitData();
    update();
}

void PlotWidget::setSeries(std::vector<PlotSeries> series)
{
    m_series = std::move(series);
    if (m_follow) fitData();
    update();
}

void PlotWidget::resetView()
{
    m_follow = true;
    fitData();
    update();
}

void PlotWidget::fitData()
{
    double x0 = INFINITY, x1 = -INFINITY, y0 = INFINITY, y1 = -INFINITY;
    const bool log = m_y.scale() == Axis::Scale::Log;
    for (const auto& s : m_series) {
        if (!s.visible) continue;
        for (std::size_t i = 0; i < s.x.size() && i < s.y.size(); ++i) {
            if (!finite(s.x[i]) || !finite(s.y[i]) || (log && s.y[i] <= 0)) continue;
            x0 = std::min(x0, s.x[i]);
            x1 = std::max(x1, s.x[i]);
            y0 = std::min(y0, s.y[i]);
            y1 = std::max(y1, s.y[i]);
        }
    }
    if (!(x1 >= x0)) { x0 = 0; x1 = 1; }
    if (!(y1 >= y0)) { y0 = log ? 1e-6 : 0; y1 = 1; }
    if (x1 == x0) x1 = x0 + 1;
    if (y1 == y0) {
        const double pad = log ? y0 : std::max(std::abs(y0) * 0.05, 1e-9);
        y0 -= log ? y0 * 0.5 : pad;
        y1 += pad;
    }
    Axis::niceRange(m_x.scale(), x0, x1);
    // Iterations start at the first sample, not at a rounded-down zero-ish value.
    Axis::niceRange(m_y.scale(), y0, y1);
    m_x.setRange(x0, x1);
    m_y.setRange(y0, y1);
}

PlotWidget::Layout PlotWidget::layout(const QRectF& r, const QFontMetricsF& fm) const
{
    Layout l;
    double yTickWidth = 0;
    for (const Tick& t : m_y.ticks(6))
        if (t.major) yTickWidth = std::max(yTickWidth, fm.horizontalAdvance(t.label));
    const double left = r.left() + yTickWidth + (m_yLabel.isEmpty() ? 12 : 26);
    // Legend chips wrap below the header.
    double x = r.left() + 10, y = r.top() + kHeader;
    const double h = fm.height() + 6;
    bool any = false;
    QStringList keys;
    for (std::size_t i = 0; i < m_series.size(); ++i) {
        const QString key = m_series[i].legendKey();
        if (keys.contains(key)) continue;
        keys << key;
        const double w = fm.horizontalAdvance(key) + 26;
        if (x + w > r.right() - 8 && x > r.left() + 10) { x = r.left() + 10; y += h; }
        l.legend.push_back({QRectF(x, y, w, h), i});
        x += w + 6;
        any = true;
    }
    const double top = any ? y + h + 6 : r.top() + kHeader + 4;
    const double bottom = r.bottom() - fm.height() - (m_xLabel.isEmpty() ? 10 : 26);
    l.plot = QRectF(QPointF(left, top), QPointF(r.right() - 14, bottom));
    return l;
}

void PlotWidget::configureAxes(const QRectF& plot) const
{
    m_x.setPixels(plot.left(), plot.right());
    m_y.setPixels(plot.bottom(), plot.top());
}

void PlotWidget::paintEvent(QPaintEvent*)
{
    QPainter p(this);
    paintTo(p, rect(), false);
}

void PlotWidget::paintTo(QPainter& p, const QRectF& full, bool forExport) const
{
    const auto& t = tokens();
    p.fillRect(full, t.panel);
    p.setRenderHint(QPainter::Antialiasing);
    const QFont base = tabular(font());
    QFont small = base;
    small.setPointSizeF(base.pointSizeF() * 0.9);
    p.setFont(small);
    const QFontMetricsF fm(small);
    // Exports reserve a footer line for the provenance notice.
    const bool footer = forExport && !m_provenance.isEmpty();
    const QRectF r = footer ? full.adjusted(0, 0, 0, -(fm.height() + 8)) : full;
    const Layout l = layout(r, fm);
    configureAxes(l.plot);

    // Header: title and badges.
    QFont titleFont = base;
    titleFont.setWeight(QFont::DemiBold);
    p.setFont(titleFont);
    p.setPen(t.text);
    const QFontMetricsF tfm(titleFont);
    p.drawText(QPointF(r.left() + 10, r.top() + 8 + tfm.ascent()), m_title);
    double bx = r.left() + 10 + tfm.horizontalAdvance(m_title) + 10;
    QFont badgeFont = base;
    badgeFont.setPointSizeF(base.pointSizeF() * 0.78);
    badgeFont.setWeight(QFont::Bold);
    badgeFont.setLetterSpacing(QFont::AbsoluteSpacing, 0.6);
    const QFontMetricsF bfm(badgeFont);
    p.setFont(badgeFont);
    for (const QString& b : m_badges) {
        const QRectF br(bx, r.top() + 8, bfm.horizontalAdvance(b) + 12, bfm.height() + 3);
        const bool amber = b.startsWith("APPROX") || b.startsWith("MODIFIED") || b.startsWith("SYNTHETIC") || b.startsWith("FIXTURE");
        const bool live = b.startsWith("LIVE");
        const bool failed = b.startsWith("FAILED") || b.startsWith("ERROR");
        QColor fill = Qt::transparent, edge = t.borderStrong, ink = t.muted;
        if (b.startsWith("APPROX")) { fill = t.amber; edge = t.amber; ink = QColor(0x2a1d05); }
        else if (amber) { edge = t.amber; ink = t.dark ? t.amber : t.amberText; }
        else if (live) { edge = t.accent; ink = t.dark ? t.accent : t.accent.darker(115); }
        else if (failed) { fill = t.red; edge = t.red; ink = Qt::white; }
        p.setPen(QPen(edge, 1));
        p.setBrush(fill);
        p.drawRoundedRect(br.adjusted(0.5, 0.5, -0.5, -0.5), 3, 3);
        p.setPen(ink);
        p.drawText(br, Qt::AlignCenter, b);
        bx += br.width() + 6;
    }

    // Legend chips (click to hide a series).
    p.setFont(small);
    for (const Chip& chip : l.legend) {
        const auto& s = m_series[chip.series];
        const QRectF c = chip.rect;
        const double cy = c.center().y();
        QPen pen(s.visible ? s.color : t.faint, 2);
        if (s.dashed) pen.setStyle(Qt::DashLine);
        p.setPen(pen);
        p.drawLine(QPointF(c.left() + 2, cy), QPointF(c.left() + 16, cy));
        p.setPen(s.visible ? t.muted : t.faint);
        p.drawText(QRectF(c.left() + 20, c.top(), c.width() - 20, c.height()), Qt::AlignVCenter | Qt::AlignLeft, s.legendKey());
    }

    const QRectF plot = l.plot;
    if (plot.width() < 20 || plot.height() < 20) return;

    // Grid and ticks.
    p.setFont(small);
    const QColor grid = t.dark ? QColor(255, 255, 255, 16) : QColor(0, 0, 0, 16);
    const QColor gridMinor = t.dark ? QColor(255, 255, 255, 7) : QColor(0, 0, 0, 7);
    for (const Tick& tk : m_y.ticks(6)) {
        const double py = m_y.toPixel(tk.value);
        if (py < plot.top() - 0.5 || py > plot.bottom() + 0.5) continue;
        p.setPen(QPen(tk.major ? grid : gridMinor, 1));
        p.drawLine(QPointF(plot.left(), py), QPointF(plot.right(), py));
        if (tk.major) {
            p.setPen(t.muted);
            p.drawText(QRectF(r.left(), py - fm.height() / 2, plot.left() - r.left() - 6, fm.height()),
                       Qt::AlignRight | Qt::AlignVCenter, tk.label);
        }
    }
    for (const Tick& tk : m_x.ticks(std::max(2, int(plot.width() / 90)))) {
        const double px = m_x.toPixel(tk.value);
        if (px < plot.left() - 0.5 || px > plot.right() + 0.5) continue;
        p.setPen(QPen(tk.major ? grid : gridMinor, 1));
        p.drawLine(QPointF(px, plot.top()), QPointF(px, plot.bottom()));
        if (tk.major) {
            p.setPen(t.muted);
            p.drawText(QRectF(px - 50, plot.bottom() + 4, 100, fm.height()), Qt::AlignHCenter | Qt::AlignTop, tk.label);
        }
    }
    p.setPen(QPen(t.border, 1));
    p.setBrush(Qt::NoBrush);
    p.drawRect(plot);
    if (!m_xLabel.isEmpty()) {
        p.setPen(t.muted);
        p.drawText(QRectF(plot.left(), r.bottom() - fm.height() - 6, plot.width(), fm.height()), Qt::AlignCenter, m_xLabel);
    }
    if (!m_yLabel.isEmpty()) {
        // Centred on the plot when it fits; otherwise on the whole widget height, elided.
        const bool fits = fm.horizontalAdvance(m_yLabel) <= plot.height();
        const double span = fits ? plot.height() : r.height() - 8;
        p.save();
        p.translate(r.left() + 6 + fm.ascent(), fits ? plot.center().y() : r.center().y());
        p.rotate(-90);
        p.setPen(t.muted);
        p.drawText(QRectF(-span / 2, -fm.ascent(), span, fm.height()), Qt::AlignCenter, fm.elidedText(m_yLabel, Qt::ElideRight, int(span)));
        p.restore();
    }

    // Series.
    bool anyData = false;
    p.save();
    p.setClipRect(plot.adjusted(-1, -1, 1, 1));
    const bool log = m_y.scale() == Axis::Scale::Log;
    for (const auto& s : m_series) {
        if (!s.visible) continue;
        QPainterPath path;
        bool pen = false;
        const std::size_t n = std::min(s.x.size(), s.y.size());
        // Decimate to about two points per pixel column, keeping the min and max.
        const double colW = 0.5;
        double colX = -INFINITY, colMin = 0, colMax = 0;
        auto flush = [&](double px) {
            if (!finite(colX)) return;
            if (!pen) { path.moveTo(px, colMin); pen = true; }
            else path.lineTo(px, colMin);
            if (colMax != colMin) path.lineTo(px, colMax);
        };
        for (std::size_t i = 0; i < n; ++i) {
            const double yv = s.y[i];
            if (!finite(s.x[i]) || !finite(yv) || (log && yv <= 0)) {
                flush(colX);
                colX = -INFINITY;
                pen = false;
                continue;
            }
            anyData = true;
            const double px = m_x.toPixel(s.x[i]), py = m_y.toPixel(yv);
            if (finite(colX) && std::abs(px - colX) < colW) {
                colMin = std::min(colMin, py);
                colMax = std::max(colMax, py);
                continue;
            }
            flush(colX);
            colX = px;
            colMin = colMax = py;
        }
        flush(colX);
        QPen lp(s.color, 1.6);
        lp.setCapStyle(Qt::RoundCap);
        lp.setJoinStyle(Qt::RoundJoin);
        if (s.dashed) lp.setDashPattern({4, 3});
        p.setPen(lp);
        p.setBrush(Qt::NoBrush);
        p.drawPath(path);
        if (s.markers && n < 400) {
            p.setBrush(s.color);
            p.setPen(Qt::NoPen);
            for (std::size_t i = 0; i < n; ++i)
                if (finite(s.x[i]) && finite(s.y[i]) && !(log && s.y[i] <= 0))
                    p.drawEllipse(QPointF(m_x.toPixel(s.x[i]), m_y.toPixel(s.y[i])), 2.2, 2.2);
        }
    }
    p.restore();

    if (!anyData) {
        p.setPen(t.faint);
        p.setFont(base);
        p.drawText(plot, Qt::AlignCenter | Qt::TextWordWrap, m_emptyText.isEmpty() ? tr("No data yet") : m_emptyText);
    }

    // Linked cursor: a vertical line and the values of every series at that x.
    if (m_cursorX && anyData && !forExport) {
        const double px = m_x.toPixel(*m_cursorX);
        if (px >= plot.left() && px <= plot.right()) {
            p.setPen(QPen(t.accent, 1));
            p.drawLine(QPointF(px, plot.top()), QPointF(px, plot.bottom()));
            QStringList lines;
            std::vector<QColor> colors;
            double xAt = *m_cursorX;
            for (const auto& s : m_series) {
                if (!s.visible) continue;
                auto i = nearestIndex(s.x, *m_cursorX);
                if (!i || *i >= s.y.size()) continue;
                xAt = s.x[*i];
                lines << QStringLiteral("%1  %2").arg(s.name, formatNumber(s.y[*i], 4));
                colors.push_back(s.color);
                const double py = m_y.toPixel(s.y[*i]);
                if (finite(py) && py >= plot.top() && py <= plot.bottom()) {
                    p.setBrush(s.color);
                    p.setPen(QPen(t.panel, 1.5));
                    p.drawEllipse(QPointF(m_x.toPixel(s.x[*i]), py), 3.2, 3.2);
                }
            }
            if (!lines.isEmpty()) {
                const QString head = (m_xLabel.isEmpty() ? QString("x") : m_xLabel) + "  " + formatNumber(xAt, 6);
                double w = fm.horizontalAdvance(head);
                for (const QString& s : lines) w = std::max(w, fm.horizontalAdvance(s) + 14);
                const double h = (lines.size() + 1) * fm.height() + 10;
                QRectF box(px + 10, plot.top() + 8, w + 16, h);
                if (box.right() > plot.right()) box.moveRight(px - 10);
                p.setPen(QPen(t.borderStrong, 1));
                p.setBrush(t.raised);
                p.drawRoundedRect(box, 4, 4);
                p.setPen(t.muted);
                p.drawText(QPointF(box.left() + 8, box.top() + 5 + fm.ascent()), head);
                for (int i = 0; i < lines.size(); ++i) {
                    const double y = box.top() + 5 + (i + 1) * fm.height();
                    p.setPen(Qt::NoPen);
                    p.setBrush(colors[std::size_t(i)]);
                    p.drawRect(QRectF(box.left() + 8, y + fm.height() / 2 - 1, 8, 2));
                    p.setPen(t.text);
                    p.drawText(QPointF(box.left() + 22, y + fm.ascent()), lines[i]);
                }
            }
        }
    }

    // Box zoom rubber band.
    if (m_zoomBox && m_dragging && !forExport) {
        p.setPen(QPen(t.accent, 1, Qt::DashLine));
        p.setBrush(QColor(t.accent.red(), t.accent.green(), t.accent.blue(), 24));
        p.drawRect(QRectF(m_pressPos, m_lastPos).normalized().intersected(plot));
    }

    if (footer) {
        p.setPen(t.muted);
        p.setFont(small);
        p.drawText(QRectF(full.left() + 10, full.bottom() - fm.height() - 4, full.width() - 20, fm.height()),
                   Qt::AlignLeft | Qt::AlignVCenter, fm.elidedText(m_provenance, Qt::ElideRight, full.width() - 20));
    }
}

void PlotWidget::setCursorX(std::optional<double> x)
{
    if (x == m_cursorX) return;
    m_cursorX = x;
    update();
}

void PlotWidget::mouseMoveEvent(QMouseEvent* e)
{
    const QPointF pos = e->position();
    const Layout l = layout(rect(), QFontMetricsF(font()));
    configureAxes(l.plot);
    if (m_dragging) {
        if (m_zoomBox) {
            m_lastPos = pos;
        } else {
            const QPointF d = pos - m_lastPos;
            m_x.pan(-d.x());
            m_y.pan(-d.y());
            m_lastPos = pos;
            m_follow = false;
        }
        update();
        return;
    }
    if (l.plot.contains(pos)) {
        m_cursorX = m_x.fromPixel(pos.x());
        emit cursorMoved(m_cursorX);
    } else if (m_cursorX) {
        m_cursorX.reset();
        emit cursorMoved(std::nullopt);
    }
    update();
}

void PlotWidget::mousePressEvent(QMouseEvent* e)
{
    const Layout l = layout(rect(), QFontMetricsF(font()));
    for (const Chip& chip : l.legend)
        if (chip.rect.contains(e->position()) && e->button() == Qt::LeftButton) {
            const QString key = m_series[chip.series].legendKey();
            const bool visible = !m_series[chip.series].visible;
            for (auto& s : m_series)
                if (s.legendKey() == key) s.visible = visible;
            emit seriesToggled(key, visible);
            if (m_follow) fitData();
            update();
            return;
        }
    if (!l.plot.contains(e->position())) return;
    m_dragging = e->button() == Qt::LeftButton;
    m_zoomBox = m_dragging && (e->modifiers() & Qt::ShiftModifier);
    m_pressPos = m_lastPos = e->position();
}

void PlotWidget::mouseReleaseEvent(QMouseEvent* e)
{
    if (m_dragging && m_zoomBox) {
        const QRectF box = QRectF(m_pressPos, e->position()).normalized();
        if (box.width() > 6 && box.height() > 6) {
            const double x0 = m_x.fromPixel(box.left()), x1 = m_x.fromPixel(box.right());
            const double y0 = m_y.fromPixel(box.bottom()), y1 = m_y.fromPixel(box.top());
            m_x.setRange(x0, x1);
            m_y.setRange(y0, y1);
            m_follow = false;
        }
    }
    m_dragging = m_zoomBox = false;
    update();
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent*) { resetView(); }

void PlotWidget::wheelEvent(QWheelEvent* e)
{
    const Layout l = layout(rect(), QFontMetricsF(font()));
    configureAxes(l.plot);
    const double steps = e->angleDelta().y() / 120.0;
    if (steps == 0) return;
    const double f = std::pow(0.85, steps);
    const QPointF pos = e->position();
    // Wheel zooms x (the iteration axis); Shift zooms y; Ctrl zooms both.
    if (!(e->modifiers() & Qt::ShiftModifier) || (e->modifiers() & Qt::ControlModifier)) m_x.zoom(m_x.fromPixel(pos.x()), f);
    if (e->modifiers() & (Qt::ShiftModifier | Qt::ControlModifier)) m_y.zoom(m_y.fromPixel(pos.y()), f);
    m_follow = false;
    update();
}

void PlotWidget::leaveEvent(QEvent*)
{
    if (m_cursorX) {
        m_cursorX.reset();
        emit cursorMoved(std::nullopt);
        update();
    }
}

void PlotWidget::keyPressEvent(QKeyEvent* e)
{
    if (e->key() == Qt::Key_0 || e->key() == Qt::Key_Home) return resetView();
    if (e->key() == Qt::Key_L) return setYScale(m_y.scale() == Axis::Scale::Log ? Axis::Scale::Linear : Axis::Scale::Log);
    QWidget::keyPressEvent(e);
}

void PlotWidget::contextMenuEvent(QContextMenuEvent* e)
{
    QMenu m(this);
    m.addAction(tr("Reset view"), this, &PlotWidget::resetView);
    QAction* log = m.addAction(tr("Logarithmic y axis"));
    log->setCheckable(true);
    log->setChecked(m_y.scale() == Axis::Scale::Log);
    connect(log, &QAction::toggled, this, [this](bool on) { setYScale(on ? Axis::Scale::Log : Axis::Scale::Linear); });
    m.addSeparator();
    m.addAction(tr("Copy data as CSV"), this, [this] { QApplication::clipboard()->setText(toCsv()); });
    m.addAction(tr("Export CSV…"), this, [this] {
        const QString f = QFileDialog::getSaveFileName(this, tr("Export CSV"), m_title + ".csv", tr("CSV (*.csv)"));
        if (!f.isEmpty()) exportCsv(f);
    });
    m.addAction(tr("Export PNG…"), this, [this] {
        const QString f = QFileDialog::getSaveFileName(this, tr("Export PNG"), m_title + ".png", tr("PNG image (*.png)"));
        if (!f.isEmpty()) exportPng(f);
    });
    m.exec(e->globalPos());
}

QString PlotWidget::toCsv() const
{
    // One x column per series keeps series with different samples honest.
    QString out;
    if (!m_provenance.isEmpty()) out += "# " + m_provenance + "\n";
    if (!m_badges.isEmpty()) out += "# " + m_badges.join(" | ") + "\n";
    QStringList head;
    std::size_t rows = 0;
    for (const auto& s : m_series) {
        head << QStringLiteral("\"%1 %2\"").arg(s.name, m_xLabel.isEmpty() ? "x" : m_xLabel) << QStringLiteral("\"%1\"").arg(s.name);
        rows = std::max(rows, std::min(s.x.size(), s.y.size()));
    }
    out += head.join(',') + "\n";
    for (std::size_t i = 0; i < rows; ++i) {
        QStringList row;
        for (const auto& s : m_series) {
            const bool has = i < s.x.size() && i < s.y.size();
            row << (has ? QString::number(s.x[i], 'g', 10) : QString()) << (has && finite(s.y[i]) ? QString::number(s.y[i], 'g', 10) : QString());
        }
        out += row.join(',') + "\n";
    }
    return out;
}

bool PlotWidget::exportCsv(const QString& path) const
{
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    f.write(toCsv().toUtf8());
    return true;
}

QImage PlotWidget::toImage(QSize logicalSize, qreal dpr, bool withProvenance) const
{
    QImage img(logicalSize * dpr, QImage::Format_ARGB32_Premultiplied);
    img.setDevicePixelRatio(dpr);
    QPainter p(&img);
    p.setFont(font());
    paintTo(p, QRectF(QPointF(0, 0), QSizeF(logicalSize)), withProvenance);
    return img;
}

bool PlotWidget::exportPng(const QString& path) const
{
    return toImage(size(), 2.0, true).save(path, "PNG");
}

} // namespace kapp
