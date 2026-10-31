#include "ui/ChartsPanel.h"

#include "ui/Theme.h"

#include <QComboBox>
#include <QJsonArray>

#include <cmath>
#include <QHBoxLayout>
#include <QTabWidget>
#include <QVBoxLayout>

namespace kapp {

namespace {

QStringList runBadges(const RunRecord& r)
{
    QStringList b;
    if (r.approximate()) b << r.approximateBadge();
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
    if (r.approximate()) p << r.approximateBadge();
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

std::vector<PlotGuide> convergenceGuides(const RunProgress& p, bool planar)
{
    const auto& t = tokens();
    std::vector<PlotGuide> out;
    for (const ConvergenceMark& m : convergenceMarks(p, planar)) {
        using K = ConvergenceMark::Kind;
        const QColor c = m.kind == K::Converged ? (t.dark ? t.accent : t.accent.darker(130))
                       : m.kind == K::Stopped   ? (t.dark ? t.amber : t.amberText)
                                                : QColor();
        out.push_back({m.kind == K::Criterion ? PlotGuide::Along::Y : PlotGuide::Along::X, m.value, m.label, c});
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
        if (!s.visible) hidden << s.legendKey();
    for (std::size_t k = 0; k < src.size(); ++k) {
        const auto& s = src[k];
        for (const QString& b : s.badges)
            if (!badges.contains(b)) badges << b;
        if (!s.provenance.isEmpty()) provenance << (multi ? s.title + ": " : QString()) + s.provenance;
        const auto& r = s.progress.residuals;
        const bool planar = m_session->hasCase() && m_session->caseSummary().planar();
        QStringList order;
        for (const QString& name : residualOrder(r.names))
            if (residualShown(name, planar)) order << name;
        std::vector<double> xs(r.iterations.begin(), r.iterations.end());
        for (int i = 0; i < order.size(); ++i) {
            const int n = int(r.names.indexOf(order[i]));
            PlotSeries ps;
            ps.name = (multi ? s.title + " · " : QString()) + residualLabel(order[i]);
            ps.legend = residualLabel(order[i]);  // runs differ by dash; one chip per equation
            ps.x = xs;
            ps.y = r.values[std::size_t(n)];
            ps.color = t.series[std::size_t(i) % t.series.size()];
            ps.dashed = k > 0;
            ps.visible = !hidden.contains(ps.legendKey());
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
    QString title = src.empty() ? QString() : src.front().title;
    if (multi) {
        QStringList dashed;
        for (std::size_t k = 1; k < src.size(); ++k) dashed << src[k].title;
        title = tr("%1 solid, %2 dashed").arg(src.front().title, dashed.join(", "));
    }
    m_residuals->setTitle(title.isEmpty() ? tr("Scaled Residuals") : tr("Scaled Residuals · %1").arg(title));
    m_reports->setTitle(title.isEmpty() ? tr("Report Plots") : tr("Report Plots · %1").arg(title));
    m_residuals->setBadges(badges);
    m_reports->setBadges(badges);
    m_residuals->setProvenance(provenance.join(" | "));
    m_reports->setProvenance(provenance.join(" | "));
    m_reports->setAxisLabels(tr("Iteration"), m_reportQuantity->currentData().toString() == "coefficient" ? tr("Coefficient") : tr("Force [N]"));
    // One run: its residual criteria and where (and why) it stopped.
    std::vector<PlotGuide> guides;
    if (src.size() == 1) guides = convergenceGuides(src.front().progress, m_session->hasCase() && m_session->caseSummary().planar());
    m_residuals->setGuides(std::move(guides));
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

namespace kapp {

void ChartsPanel::plotProbe(const PostOutput& probe, const QString& title, int xAxis)
{
    const auto* item = probe.scene ? probe.scene->item("post:" + probe.name) : nullptr;
    if (!item || !item->mesh) return;
    const MeshData& m = item->mesh->data();
    const MeshData::Field* arc = m.field("arc_length");
    const MeshData::Field* valid = m.field("valid");
    const QJsonArray lines = probe.part.meta.value("lines").toArray();
    QStringList quantities;
    QString unit;
    std::vector<PlotSeries> series;
    int colour = -1;
    for (const auto& f : m.fields) {
        if (f.desc.name == "arc_length" || f.desc.name == "valid") continue;
        const FieldRef ref = FieldRef::parse(f.desc.name);
        // Kutta and Fluent values of one quantity share a colour.
        const qsizetype q = quantities.indexOf(ref.quantity);
        if (q < 0) {
            quantities << ref.quantity;
            unit = f.desc.unit;
            ++colour;
        }
        const int c = f.desc.components;
        for (std::size_t j = 0; j + 1 < m.indices.size(); ++j) {
            PlotSeries ps;
            ps.name = displayLabel(ref);
            const QString zone = lines.at(qsizetype(j)).toObject().value("zone").toString();
            if (m.indices.size() > 2) ps.name += zone.isEmpty() ? QStringLiteral(" · curve %1").arg(j + 1) : QStringLiteral(" · %1 (%2)").arg(zone).arg(j + 1);
            ps.color = tokens().series[std::size_t(q < 0 ? colour : q) % tokens().series.size()];
            ps.dashed = ref.role == FieldRole::Fluent;
            for (std::uint32_t v = m.indices[j]; v < m.indices[j + 1]; ++v) {
                ps.x.push_back(xAxis >= 0 ? m.positions[3 * v + std::size_t(xAxis)] : arc ? arc->values[v] : double(v));
                double y = f.values[std::size_t(c) * v];
                if (c == 3) y = std::hypot(y, double(f.values[3 * v + 1]), double(f.values[3 * v + 2]));
                if (valid && valid->values[v] == 0) y = NAN;
                ps.y.push_back(y);
            }
            series.push_back(std::move(ps));
        }
    }
    QStringList badges;
    if (probe.synthetic) badges << "SYNTHETIC";
    if (probe.provenance.contains("APPROXIMATE")) badges << "APPROXIMATE";
    if (probe.provenance.contains("MODIFIED")) badges << "MODIFIED";
    const QString xLabel = xAxis >= 0 ? QStringLiteral("%1 [m]").arg(QChar('x' + xAxis)) : tr("Distance along curve [m]");
    QString yLabel = tr("Value");
    if (quantities.size() == 1) {
        const QString u = prettyUnit(unit);
        yLabel = u.isEmpty() ? fluentLabel(quantities[0]) : QStringLiteral("%1 [%2]").arg(fluentLabel(quantities[0]), u);
    }
    setXY(title, xLabel, yLabel, std::move(series), badges, probe.provenance);
}

} // namespace kapp
