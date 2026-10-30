#pragma once

// A 2D line chart for residuals, report histories and XY plots. Painted with
// QPainter on demand: nothing repaints unless data, size or the cursor changes.

#include "core/ChartScale.h"

#include <QColor>
#include <QWidget>

#include <optional>
#include <vector>

namespace kapp {

struct PlotSeries {
    QString name;
    std::vector<double> x, y;  // NaN breaks the line
    QColor color;
    bool dashed = false;
    bool visible = true;
    bool markers = false;
    // Series sharing a legend key get one chip, which toggles them together
    // (one per residual when several runs overlay). Empty: the name.
    QString legend;
    QString legendKey() const { return legend.isEmpty() ? name : legend; }
};

class PlotWidget : public QWidget {
    Q_OBJECT
public:
    explicit PlotWidget(QWidget* parent = nullptr);

    void setTitle(const QString& t);
    void setAxisLabels(const QString& x, const QString& y);
    void setYScale(Axis::Scale s);
    Axis::Scale yScale() const { return m_y.scale(); }
    // Badges drawn in the header and into exports ("APPROXIMATE", "LIVE").
    void setBadges(const QStringList& badges);
    void setProvenance(const QString& p) { m_provenance = p; }
    void setEmptyText(const QString& t);

    void setSeries(std::vector<PlotSeries> series);
    const std::vector<PlotSeries>& series() const { return m_series; }
    // Data changed and the view follows it unless the user zoomed or panned.
    void resetView();
    bool following() const { return m_follow; }

    // The shared iteration cursor: views that plot against the same x link it.
    void setCursorX(std::optional<double> x);
    std::optional<double> cursorX() const { return m_cursorX; }

    QString toCsv() const;
    QImage toImage(QSize logicalSize, qreal dpr, bool withProvenance) const;
    bool exportCsv(const QString& path) const;
    bool exportPng(const QString& path) const;

signals:
    void cursorMoved(std::optional<double> x);
    void seriesToggled(const QString& name, bool visible);

protected:
    void paintEvent(QPaintEvent*) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void leaveEvent(QEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;

private:
    // series: the first with the chip's key; npos marks the "+N more" chip.
    struct Chip { QRectF rect; std::size_t series; QString text; };
    // Every label owns its own rectangle: legend rows, the y title (rotated
    // beside the axis when it fits, else on one line above the plot), tick
    // labels and the x title. Nothing is drawn where another element is.
    struct Layout {
        QRectF plot;
        std::vector<Chip> legend;
        QStringList hiddenLegend;  // keys folded into the "+N more" chip
        QRectF yTitle;             // empty when there is no y title
        bool yRotated = false;
        QRectF xTitle;
        int xTicks = 6, yTicks = 6;
    };
    Layout layout(const QRectF& r) const;
    QFont chartFont() const;
    void paintTo(QPainter& p, const QRectF& r, bool forExport) const;
    void fitData();
    void configureAxes(const QRectF& plot) const;

    QString m_title, m_xLabel, m_yLabel, m_emptyText, m_provenance;
    QStringList m_badges;
    std::vector<PlotSeries> m_series;
    mutable Axis m_x, m_y;
    bool m_follow = true;
    std::optional<double> m_cursorX;
    QPointF m_pressPos, m_lastPos;
    bool m_dragging = false;
    bool m_zoomBox = false;
};

} // namespace kapp
