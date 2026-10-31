#pragma once

// Residuals, report histories and XY plots. Follows the live run while one is
// going; otherwise shows the runs chosen in history (several overlay, the
// later ones dashed). Residual and report charts share the iteration cursor.

#include "core/Session.h"
#include "ui/PlotWidget.h"

#include <QTimer>
#include <QWidget>

class QComboBox;
class QTabWidget;

namespace kapp {

// convergenceMarks() as chart guides, coloured by outcome.
std::vector<PlotGuide> convergenceGuides(const RunProgress& p, bool planar);

class ChartsPanel : public QWidget {
    Q_OBJECT
public:
    explicit ChartsPanel(Session* session, QWidget* parent = nullptr);

    enum class Tab { Residuals, Reports, XY };
    void showTab(Tab t);
    void showRuns(const QStringList& runIds);  // empty: latest run
    QStringList shownRuns() const { return m_runs; }
    void setXY(const QString& title, const QString& xLabel, const QString& yLabel, std::vector<PlotSeries> series,
               const QStringList& badges, const QString& provenance);
    // Plot an export-view probe: one series per field and curve, against arc
    // length (xAxis -1) or a coordinate (0..2). Fluent fields are dashed.
    void plotProbe(const PostOutput& probe, const QString& title, int xAxis);
    PlotWidget* residuals() const { return m_residuals; }
    PlotWidget* reports() const { return m_reports; }
    PlotWidget* xy() const { return m_xy; }

private:
    void rebuild();
    void scheduleLive();
    struct Source { QString title; RunProgress progress; QStringList badges; QString provenance; };
    std::vector<Source> sources() const;

    Session* m_session;
    QTabWidget* m_tabs;
    PlotWidget* m_residuals;
    PlotWidget* m_reports;
    PlotWidget* m_xy;
    QComboBox* m_reportQuantity;
    QStringList m_runs;
    QTimer m_liveTimer;
};

} // namespace kapp
