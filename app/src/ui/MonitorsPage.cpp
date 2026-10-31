// Solution ▸ Monitors: Fluent's residual monitor criteria and report
// convergence conditions, and where the live or last run stands against them.

#include "ui/ChartsPanel.h"
#include "ui/Inspector.h"
#include "ui/PageKit.h"
#include "ui/SettingWidgets.h"

#include <QGridLayout>
#include <QLocale>

namespace kapp::pages {

using namespace kit;

namespace {

// An inventory editor that follows the model (undo, restore, reset).
QWidget* editor(const InspectorContext& c, const QString& path, QWidget* parent, int maxWidth = 0)
{
    SettingsModel* m = c.session->settings();
    const Setting* s = m->find(path);
    if (!s) return label(QStringLiteral("—"), "faint", parent);
    auto* e = new SettingEditor(m, *s, parent);
    QObject::connect(m, &SettingsModel::valueChanged, e, [e, path](const QString& p) { if (p == path) e->refresh(); });
    e->refresh();
    if (maxWidth > 0) e->setMaximumWidth(maxWidth);
    return e;
}

QString pct(double fraction) { return QLocale().toString(fraction * 100, 'g', 3) + QStringLiteral(" %"); }

} // namespace

namespace {

// Checked residuals that were above their criterion at the sample before the
// last: the ones whose drop ended the run.
QStringList lastToConverge(const RunProgress& p, const ConvergenceSetup& setup)
{
    const auto& r = p.residuals;
    QStringList out;
    if (r.iterations.size() < 2) return out;
    const std::size_t prev = r.iterations.size() - 2;
    for (auto it = setup.residuals.begin(); it != setup.residuals.end(); ++it) {
        const qsizetype n = r.names.indexOf(it.key());
        if (n >= 0 && r.values[std::size_t(n)][prev] > it.value()) out << it.key();
    }
    return residualOrder(out);
}

} // namespace

QWidget* convergenceStatus(const InspectorContext& c)
{
    Session* s = c.session;
    auto* box = new QFrame;
    box->setProperty("callout", "neutral");
    box->setAccessibleName(QObject::tr("Convergence status"));
    auto* v = new QVBoxLayout(box);
    v->setContentsMargins(10, 8, 10, 10);
    v->setSpacing(3);
    auto* title = label({});
    QFont tf = title->font();
    tf.setWeight(QFont::DemiBold);
    title->setFont(tf);
    v->addWidget(title);
    auto* note = label({}, "muted");
    v->addWidget(note);
    auto* grid = new QGridLayout;
    grid->setHorizontalSpacing(12);
    grid->setVerticalSpacing(2);
    v->addLayout(grid);

    auto refresh = [s, box, title, grid, note] {
        while (QLayoutItem* it = grid->takeAt(0)) {
            delete it->widget();
            delete it;
        }
        const auto& t = tokens();
        const bool planar = s->hasCase() && s->caseSummary().planar();
        auto row = [&](int r, const QString& name, const QString& value, const QString& criterion, bool met) {
            grid->addWidget(label(name), r, 0);
            auto* val = label(value);
            val->setAlignment(Qt::AlignRight);
            grid->addWidget(val, r, 1);
            grid->addWidget(label(criterion, "faint"), r, 2);
            auto* m = label(met ? QObject::tr("met") : QObject::tr("not met"), met ? QString() : "muted");
            if (met) m->setStyleSheet(QStringLiteral("color: %1").arg((t.dark ? t.accent : t.accent.darker(125)).name()));
            grid->addWidget(m, r, 3);
        };
        // Each condition at one monitor sample: residuals first, then reports.
        auto conditions = [&](const ConvergenceState& state) {
            int r = 0;
            for (const ConditionState& cs : state.conditions) {
                if (cs.kind == "residuals") {
                    for (const QString& eq : residualOrder(cs.equations.keys())) {
                        if (!residualShown(eq, planar)) continue;
                        const auto& e = cs.equations[eq];
                        row(r++, residualLabel(eq), QString::number(e.value, 'e', 2),
                            QObject::tr("≤ %1").arg(QString::number(e.criterion, 'g', 3)), e.met);
                    }
                } else {
                    row(r++, cs.report, cs.value ? pct(*cs.value) : QObject::tr("too few values"),
                        QObject::tr("< %1").arg(pct(cs.criterion.value_or(0))), cs.met);
                }
            }
            grid->setColumnStretch(4, 1);
        };
        const auto& live = s->live();
        if (s->running()) {
            const std::optional<ConvergenceState> state = live.last ? live.last->convergence : std::nullopt;
            const ConvergenceSetup setup = live.setup ? live.setup->convergence : ConvergenceSetup{};
            if (!setup.enabled) {
                title->setText(QObject::tr("Running to the iteration limit"));
                note->setText(QObject::tr("Convergence checking was off when this run started."));
            } else {
                title->setText(state ? QObject::tr("Iteration %1: %2").arg(live.last->iter).arg(state->met ? QObject::tr("converged") : QObject::tr("not converged yet"))
                                     : QObject::tr("Waiting for the first monitor sample"));
                note->setText(setup.condition == "any" ? QObject::tr("Stops when any condition is met.") : QObject::tr("Stops when all conditions are met at the same sample."));
            }
            if (state) conditions(*state);
            box->setVisible(true);
            return;
        }
        const RunRecord* last = s->history().runs().empty() ? nullptr : &s->history().runs().back();
        if (!last || last->status == RunStatus::Running) {
            box->setVisible(false);
            return;
        }
        box->setVisible(true);
        const RunProgress progress = s->history().progress(last->id);
        const auto summary = s->history().summary(last->id);
        const ConvergenceSetup conv = summary ? summary->convergence : progress.setup ? progress.setup->convergence : ConvergenceSetup{};
        if (last->converged()) {
            title->setText(QObject::tr("%1 converged at iteration %2").arg(last->title()).arg(last->iterationsDone));
            QStringList met;
            for (const QString& id : last->convergedBy) met << describeCondition(id, conv);
            QString text = QObject::tr("Stopped because %1.").arg(met.join(QObject::tr("; and ")));
            if (last->convergedBy.contains("residuals")) {
                QStringList lastEq;
                for (const QString& eq : lastToConverge(progress, conv))
                    if (residualShown(eq, planar)) lastEq << eq;
                if (!lastEq.isEmpty()) text += ' ' + QObject::tr("Last to reach its criterion: %1.").arg(lastEq.join(", "));
            }
            note->setText(text);
        } else if (last->status == RunStatus::Completed) {
            title->setText(QObject::tr("%1 reached the iteration limit (%2) without converging").arg(last->title()).arg(last->iterationsDone));
            note->setText(conv.enabled ? QObject::tr("Not every condition was met at the last sample; check the residual and report plots before trusting the numbers.")
                                       : QObject::tr("Convergence checking was off for that run, so it ran to the iteration limit."));
        } else if (last->status == RunStatus::Stopped) {
            title->setText(QObject::tr("%1 stopped by you at iteration %2").arg(last->title()).arg(last->iterationsDone));
            note->setText(conv.enabled ? QObject::tr("Not converged: you stopped it before its conditions were met. The results so far are kept and can be post-processed.")
                                       : QObject::tr("Not converged (checking was off). The results so far are kept and can be post-processed."));
        } else {
            title->setText(QObject::tr("%1 %2 at iteration %3").arg(last->title(), runStatusText(*last).toLower()).arg(last->iterationsDone));
            note->setText(QObject::tr("Not converged."));
        }
        if (progress.convergence) conditions(*progress.convergence);
    };
    refresh();
    QObject::connect(s, &Session::iteration, box, refresh);
    QObject::connect(s, &Session::runStateChanged, box, refresh);
    QObject::connect(s, &Session::historyChanged, box, refresh);
    return box;
}

QWidget* monitors(const InspectorContext& c)
{
    QWidget* page = new QWidget;
    auto* outer = new QVBoxLayout(page);
    outer->setContentsMargins(0, 0, 0, 0);
    outer->setSpacing(0);
    outer->addWidget(header(QObject::tr("Monitors"), QObject::tr("Solution"), page));
    auto* content = new QWidget(page);
    QVBoxLayout* v = body(content);
    outer->addWidget(content);
    outer->addStretch(1);
    SettingsModel* m = c.session->settings();

    if (!m->find("monitors/check-convergence")) {
        // An older CLI: no criteria, runs stop at the iteration limit.
        v->addWidget(label(QObject::tr("This kutta CLI has no convergence criteria: runs stop at the iteration count. Scaled residuals and report values are still plotted live."), "muted"));
        auto* f = new SettingsForm(m);
        f->addSetting("run/monitor-interval");
        f->finish();
        v->addWidget(f);
        return page;
    }

    auto* top = new QGridLayout;
    top->setHorizontalSpacing(10);
    top->addWidget(label(QObject::tr("Stop when converged")), 0, 0);
    top->addWidget(editor(c, "monitors/check-convergence", content), 0, 1);
    top->addWidget(label(QObject::tr("Converged when")), 1, 0);
    top->addWidget(editor(c, "monitors/condition", content), 1, 1);
    top->addWidget(label(QObject::tr("Monitor interval")), 2, 0);
    top->addWidget(editor(c, "run/monitor-interval", content), 2, 1);
    top->setColumnStretch(1, 1);
    v->addLayout(top);
    v->addWidget(label(QObject::tr("Unchecked, runs go to the iteration count. These settings say when to stop, not what is solved, so they never make a run MODIFIED."), "faint"));
    v->addWidget(convergenceStatus(c));

    // Residuals, as Fluent's Residual Monitors dialog lists them.
    v->addWidget(section(QObject::tr("Residuals")));
    auto* rg = new QGridLayout;
    rg->setHorizontalSpacing(12);
    rg->setVerticalSpacing(4);
    rg->addWidget(label(QObject::tr("Equation"), "muted"), 0, 0);
    rg->addWidget(label(QObject::tr("Check"), "muted"), 0, 1);
    rg->addWidget(label(QObject::tr("Absolute criteria"), "muted"), 0, 2);
    QStringList eqs;
    for (const Setting& s : m->settings())
        if (s.path.startsWith("monitors/residual/") && s.path.endsWith("/check")) eqs << s.path.section('/', 2, 2);
    eqs = residualOrder(eqs);
    const bool planar = c.session->caseSummary().planar();
    eqs.removeIf([planar](const QString& eq) { return !residualShown(eq, planar); });
    int r = 1;
    for (const QString& eq : eqs) {
        rg->addWidget(label(residualLabel(eq)), r, 0);
        rg->addWidget(editor(c, "monitors/residual/" + eq + "/check", content), r, 1);
        rg->addWidget(editor(c, "monitors/residual/" + eq + "/absolute-criterion", content, 120), r, 2);
        ++r;
    }
    rg->setColumnStretch(3, 1);
    v->addLayout(rg);
    v->addWidget(label(QObject::tr("Met when every checked equation's scaled residual is at or below its criterion. Values are read from the case."), "faint"));
    // Fluent's definition, which Kutta reports: users compare these numbers with Fluent's.
    auto* def = label(QObject::tr("Scaled residuals, as Fluent defines them: each momentum and turbulence equation's imbalance summed over the cells, "
                                  "divided by the sum of |a_P φ_P|; continuity's mass imbalance divided by its largest value in the first five iterations."), "faint");
    def->setAccessibleName(QObject::tr("Scaled residual definition"));
    v->addWidget(def);
    if (planar) v->addWidget(label(QObject::tr("2D case: no z-velocity. Kutta solves w on the one-cell extrusion, where it stays at round-off."), "faint"));

    // Report definitions, as Fluent's Convergence Conditions dialog lists them.
    v->addWidget(section(QObject::tr("Report definitions")));
    QStringList reports;
    for (const Setting& s : m->settings())
        if (s.path.startsWith("monitors/report/") && s.path.endsWith("/active")) reports << s.path.section('/', 2, 2);
    if (reports.isEmpty()) {
        v->addWidget(label(QObject::tr("This case defines no force reports, so only residuals can decide convergence."), "muted"));
    } else {
        auto* g = new QGridLayout;
        g->setHorizontalSpacing(8);
        g->setVerticalSpacing(4);
        const QStringList heads{QObject::tr("Report"), QObject::tr("Active"), QObject::tr("Stop at"), QObject::tr("Values"), QObject::tr("Ignore <")};
        for (int i = 0; i < heads.size(); ++i) g->addWidget(label(heads[i], "muted"), 0, i);
        int row = 1;
        for (const QString& rep : reports) {
            const QString base = "monitors/report/" + rep + "/";
            g->addWidget(label(decodePathSegment(rep)), row, 0);
            g->addWidget(editor(c, base + "active", content), row, 1);
            g->addWidget(editor(c, base + "stop-criterion", content, 76), row, 2);
            g->addWidget(editor(c, base + "previous-values", content, 56), row, 3);
            g->addWidget(editor(c, base + "ignore-before", content, 64), row, 4);
            ++row;
        }
        g->setColumnStretch(0, 1);
        v->addLayout(g);
        // Reports the CLI cannot use as a condition say why, once.
        QStringList fixed;
        QString why;
        for (const QString& rep : reports)
            if (const Setting* a = m->find("monitors/report/" + rep + "/active"); a && !a->editable) {
                fixed << decodePathSegment(rep);
                why = a->reason;
            }
        if (!fixed.isEmpty())
            v->addWidget(label(QObject::tr("%1: read only (%2).").arg(fixed.join(", "), why), "faint"));
        v->addWidget(label(QObject::tr("Stop at: Fluent's stop criterion. Values: how many previous values to compare (Fluent's Use Iterations). Ignore <: iterations to skip at the start."), "faint"));
        v->addWidget(label(QObject::tr("A report is converged when the change between each of its last N values, divided by the current value, is below the stop criterion (Fluent's rule). Kutta takes report values at the monitor interval, so with an interval above 1 this is stricter than Fluent's per-iteration check."), "faint"));
    }

    auto* row = new QHBoxLayout;
    auto* res = button(QObject::tr("Residual plot"), "chart-line");
    auto* rep = button(QObject::tr("Report plot"), "chart-spline");
    row->addWidget(res);
    row->addWidget(rep);
    row->addStretch(1);
    v->addLayout(row);
    QObject::connect(res, &QPushButton::clicked, res, [c] { c.charts->showTab(ChartsPanel::Tab::Residuals); });
    QObject::connect(rep, &QPushButton::clicked, rep, [c] { c.charts->showTab(ChartsPanel::Tab::Reports); });
    return page;
}

} // namespace kapp::pages
