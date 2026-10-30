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
