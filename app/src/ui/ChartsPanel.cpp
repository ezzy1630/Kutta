#include "ui/ChartsPanel.h"

#include "ui/Theme.h"

#include <QComboBox>
#include <QHBoxLayout>
#include <QTabWidget>
#include <QVBoxLayout>

namespace kapp {

namespace {

QStringList runBadges(const RunRecord& r)
{
    QStringList b;
    if (r.approximate()) b << "APPROXIMATE";
    if (r.modified()) b << "MODIFIED";
    if (r.synthetic) b << "SYNTHETIC";
    else if (r.fixture) b << "FIXTURE REPLAY";
    if (r.status == RunStatus::Failed) b << "FAILED";
    if (r.status == RunStatus::Running) b << "LIVE";
    return b;
}

QString runProvenance(const RunRecord& r)
{
    QStringList p;
    if (r.approximate()) p << "APPROXIMATE: " + r.approximations.join(", ");
    if (r.modified()) p << QStringLiteral("MODIFIED: %1 settings differ from the case").arg(r.overrides.size());
    if (r.synthetic) p << "SYNTHETIC fixture data, not a Kutta solution";
    else if (r.fixture) p << "FIXTURE REPLAY of a recorded Kutta run";
    return p.join(" | ");
}

} // namespace

ChartsPanel::ChartsPanel(Session* session, QWidget* parent) : QWidget(parent), m_session(session)
{
    auto* v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);
    m_tabs = new QTabWidget(this);
    m_tabs->setDocumentMode(true);
    v->addWidget(m_tabs);
    m_residuals = new PlotWidget;
    m_residuals->setTitle(tr("Scaled Residuals"));
    m_residuals->setAxisLabels(tr("Iteration"), {});
    m_residuals->setYScale(Axis::Scale::Log);
    m_residuals->setEmptyText(tr("Residuals appear here when a calculation runs."));
    m_reports = new PlotWidget;
    m_reports->setTitle(tr("Report Plots"));
    m_reports->setAxisLabels(tr("Iteration"), {});
    m_reports->setEmptyText(tr("Report definitions are plotted per iteration.\nThis case may define none; see Solution ▸ Report Definitions."));
    m_xy = new PlotWidget;
    m_xy->setTitle(tr("XY Plot"));
    m_xy->setEmptyText(tr("Create a line probe in Results ▸ Plots ▸ XY Plot."));
    m_tabs->addTab(m_residuals, tr("Residuals"));
    m_tabs->addTab(m_reports, tr("Reports"));
    m_tabs->addTab(m_xy, tr("XY Plot"));

    auto* corner = new QWidget(this);
    auto* ch = new QHBoxLayout(corner);
    ch->setContentsMargins(0, 0, 6, 0);
    m_reportQuantity = new QComboBox(corner);
    m_reportQuantity->addItem(tr("Force (N)"), "force");
    m_reportQuantity->addItem(tr("Coefficient"), "coefficient");
    m_reportQuantity->setAccessibleName(tr("Report quantity"));
    ch->addWidget(m_reportQuantity);
    m_tabs->setCornerWidget(corner);
    corner->setVisible(false);
    connect(m_tabs, &QTabWidget::currentChanged, this, [corner](int i) { corner->setVisible(i == 1); });
    connect(m_reportQuantity, &QComboBox::activated, this, &ChartsPanel::rebuild);

    // The linked cursor: residuals and reports share the iteration under the mouse.
    connect(m_residuals, &PlotWidget::cursorMoved, m_reports, &PlotWidget::setCursorX);
    connect(m_reports, &PlotWidget::cursorMoved, m_residuals, &PlotWidget::setCursorX);

    m_liveTimer.setSingleShot(true);
    m_liveTimer.setInterval(120);
    connect(&m_liveTimer, &QTimer::timeout, this, &ChartsPanel::rebuild);
    connect(session, &Session::iteration, this, &ChartsPanel::scheduleLive);
    connect(session, &Session::runStarted, this, [this] {
        m_runs.clear();
        m_residuals->resetView();
        m_reports->resetView();
        rebuild();
    });
    connect(session, &Session::runFinished, this, &ChartsPanel::rebuild);
    connect(session, &Session::caseOpened, this, [this] { m_runs.clear(); rebuild(); });
    connect(&Theme::instance(), &Theme::changed, this, &ChartsPanel::rebuild);
}

void ChartsPanel::showTab(Tab t) { m_tabs->setCurrentIndex(int(t)); }

void ChartsPanel::showRuns(const QStringList& runIds)
{
    m_runs = runIds;
    m_residuals->resetView();
    m_reports->resetView();
    rebuild();
}

void ChartsPanel::scheduleLive()
{
    if (!m_liveTimer.isActive()) m_liveTimer.start();
}

std::vector<ChartsPanel::Source> ChartsPanel::sources() const
{
    std::vector<Source> out;
    const auto& h = m_session->history();
    if (m_session->running() && m_runs.isEmpty()) {
        const RunRecord* r = h.find(m_session->live().runId);
        Source s{r ? r->title() : tr("Live"), m_session->live().progress, r ? runBadges(*r) : QStringList{"LIVE"}, r ? runProvenance(*r) : QString()};
        if (!s.badges.contains("LIVE")) s.badges << "LIVE";
        out.push_back(std::move(s));
        return out;
    }
    QStringList ids = m_runs;
    if (ids.isEmpty())
        if (const RunRecord* r = h.runs().empty() ? nullptr : &h.runs().back()) ids << r->id;
    for (const QString& id : ids) {
        const RunRecord* r = h.find(id);
        if (!r) continue;
        RunProgress p = (m_session->running() && m_session->live().runId == id) ? m_session->live().progress : h.progress(id);
        out.push_back({r->title(), std::move(p), runBadges(*r), runProvenance(*r)});
    }
    return out;
}

void ChartsPanel::rebuild()
{
    const auto& t = tokens();
    const auto src = sources();
    const bool multi = src.size() > 1;
    std::vector<PlotSeries> res, rep;
    QStringList badges;
    QStringList provenance;
    QStringList hidden;
    for (const auto& s : m_residuals->series())
        if (!s.visible) hidden << s.name;
    for (std::size_t k = 0; k < src.size(); ++k) {
        const auto& s = src[k];
        for (const QString& b : s.badges)
            if (!badges.contains(b)) badges << b;
        if (!s.provenance.isEmpty()) provenance << (multi ? s.title + ": " : QString()) + s.provenance;
        const auto& r = s.progress.residuals;
        const QStringList order = residualOrder(r.names);
        std::vector<double> xs(r.iterations.begin(), r.iterations.end());
        for (int i = 0; i < order.size(); ++i) {
            const int n = int(r.names.indexOf(order[i]));
            PlotSeries ps;
            ps.name = (multi ? s.title + " · " : QString()) + residualLabel(order[i]);
            ps.x = xs;
            ps.y = r.values[std::size_t(n)];
            ps.color = t.series[std::size_t(i) % t.series.size()];
            ps.dashed = k > 0;
            ps.visible = !hidden.contains(ps.name);
            res.push_back(std::move(ps));
        }
        const bool coeff = m_reportQuantity->currentData().toString() == "coefficient";
        for (std::size_t i = 0; i < s.progress.reports.size(); ++i) {
            const auto& rs = s.progress.reports[i];
            PlotSeries ps;
            ps.name = (multi ? s.title + " · " : QString()) + rs.name;
            ps.x.assign(rs.iter.begin(), rs.iter.end());
            ps.y = coeff ? rs.coefficient : rs.force;
            ps.color = t.series[(i + 2) % t.series.size()];
            ps.dashed = k > 0;
            rep.push_back(std::move(ps));
        }
    }
    QString title = src.empty() ? QString() : src.size() == 1 ? src.front().title : tr("%1 runs").arg(src.size());
    m_residuals->setTitle(title.isEmpty() ? tr("Scaled Residuals") : tr("Scaled Residuals · %1").arg(title));
    m_reports->setTitle(title.isEmpty() ? tr("Report Plots") : tr("Report Plots · %1").arg(title));
    m_residuals->setBadges(badges);
    m_reports->setBadges(badges);
    m_residuals->setProvenance(provenance.join(" | "));
    m_reports->setProvenance(provenance.join(" | "));
    m_reports->setAxisLabels(tr("Iteration"), m_reportQuantity->currentData().toString() == "coefficient" ? tr("Coefficient") : tr("Force [N]"));
    m_residuals->setSeries(std::move(res));
    m_reports->setSeries(std::move(rep));
}

void ChartsPanel::setXY(const QString& title, const QString& xLabel, const QString& yLabel, std::vector<PlotSeries> series,
                        const QStringList& badges, const QString& provenance)
{
    m_xy->setTitle(title);
    m_xy->setAxisLabels(xLabel, yLabel);
    m_xy->setBadges(badges);
    m_xy->setProvenance(provenance);
    m_xy->setSeries(std::move(series));
    m_xy->resetView();
    showTab(Tab::XY);
}

} // namespace kapp
