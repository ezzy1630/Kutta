#include "ui/Inspector.h"
#include "core/ChartScale.h"

#include "ui/ChartsPanel.h"
#include "ui/GraphicsArea.h"
#include "ui/Outline.h"
#include "ui/PageKit.h"
#include "ui/ScreenChecks.h"
#include "ui/SettingWidgets.h"

#include <QCheckBox>
#include <QCoreApplication>
#include <QComboBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QGridLayout>
#include <QHeaderView>
#include <QKeySequence>
#include <QJsonArray>
#include <QLocale>
#include <QScrollBar>
#include <QSettings>
#include <QTableWidget>
#include <QTimer>

namespace kapp {

using namespace kit;

Inspector::Inspector(Session* session, GraphicsArea* graphics, ChartsPanel* charts, QWidget* parent) : QScrollArea(parent)
{
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    setAccessibleName(tr("Inspector"));
    m_ctx = {session, graphics, charts, [this](const QString& id) { emit navigate(id); },
             [this](const QString& f) { showSettingsFiltered(f); }};
    m_ctx.refresh = [this] { refresh(); };
    m_ctx.current = [this] { return m_id; };
    connect(session, &Session::caseOpened, this, &Inspector::refresh);
    connect(session, &Session::caseFailed, this, &Inspector::refresh);
    connect(&Theme::instance(), &Theme::changed, this, &Inspector::refresh);
    connect(graphics, &GraphicsArea::activeChanged, this, [this] {
        if (m_id.startsWith("results/") && m_id != "results/residuals" && m_id != "results/report-plots") refresh();
    });
    connect(session, &Session::newCaseChanged, this, [this] {
        if (m_id == "new-case") refresh();
    });
    connect(session, &Session::approvalsChanged, this, [this] {
        if (m_id == "setup/general" || m_id == "solution/run") refresh();
    });
    connect(session, &Session::historyChanged, this, [this] {
        // Pages listing runs follow the history; editing pages keep their state.
        if (m_id.startsWith("runs/") || m_id == "results/forces" || m_id == "results/aero" || m_id == "setup/general") refresh();
    });
    show("setup/general");
}

// Pages are exactly as wide as the panel: wrapped text and tables adapt, and
// nothing hides past the right edge (there is no horizontal scrolling).
void Inspector::fitPageWidth()
{
    if (QWidget* w = widget()) w->setFixedWidth(viewport()->width());
}

void Inspector::resizeEvent(QResizeEvent* e)
{
    QScrollArea::resizeEvent(e);
    fitPageWidth();
}

QStringList Inspector::layoutProblemsAtWidths(const QList<int>& widths)
{
    QWidget* page = widget();
    if (!page) return {};
    QStringList problems;
    auto settle = [this, page](int width) {
        page->setFixedWidth(width);
        // As the scroll area sizes a resizable page: its height for that width.
        for (int pass = 0; pass < 2; ++pass) {
            QCoreApplication::sendPostedEvents(nullptr, QEvent::LayoutRequest);
            const int h = page->hasHeightForWidth() ? page->heightForWidth(width) : page->sizeHint().height();
            page->resize(width, std::max(h, viewport()->height()));
        }
        QCoreApplication::sendPostedEvents(nullptr, QEvent::LayoutRequest);
    };
    for (int width : widths) {
        settle(width);
        for (const QString& p : layoutProblems(page)) problems << QStringLiteral("%1 at %2 px: %3").arg(m_id).arg(width).arg(p);
    }
    fitPageWidth();
    page->updateGeometry();
    QCoreApplication::sendPostedEvents(nullptr, QEvent::LayoutRequest);
    return problems;
}

void Inspector::showSettingsFiltered(const QString& filter)
{
    m_settingsFilter = filter;
    emit navigate("settings");
    show("settings");
}

void Inspector::show(const QString& rawId)
{
    QString id = rawId;
    static const QHash<QString, QString> groups = {
        {"setup", "setup/general"},           {"solution", "solution/run"},        {"results", "results/contours"},
        {"results/graphics", "results/contours"}, {"results/plots", "results/residuals"}, {"results/reports", "results/forces"},
        {"runs", "runs/history"}};
    id = groups.value(id, id);
    const int scroll = id == m_id ? verticalScrollBar()->value() : 0;
    m_id = id;
    QWidget* old = takeWidget();
    if (old) old->deleteLater();
    QWidget* page = build(id);
    page->setObjectName("Panel");
    // Pickers size to the panel, not to their longest entry, so no page is
    // wider than the inspector.
    for (auto* cb : page->findChildren<QComboBox*>()) {
        cb->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
        cb->setMinimumContentsLength(8);
    }
    setWidget(page);
    fitPageWidth();
    if (scroll) QTimer::singleShot(0, this, [this, scroll] { verticalScrollBar()->setValue(scroll); });
    if (id != "settings") m_settingsFilter.clear();
}

QWidget* Inspector::build(const QString& id)
{
    const InspectorContext& c = m_ctx;
    if (id == "new-case") return pages::newCase(c);
    const bool needsCase = id.startsWith("setup/") || id.startsWith("solution/");
    if (needsCase && !c.session->hasCase()) return pages::noCase(c, Outline::titleOf(id));
    if (id == "setup/general") return pages::general(c);
    if (id == "setup/models") return pages::models(c);
    if (id == "setup/materials") return pages::materials(c, {});
    if (id.startsWith("setup/material/")) return pages::materials(c, id.mid(15));
    if (id == "setup/cell-zones") return pages::cellZones(c);
    if (id == "setup/boundaries") return pages::boundaries(c);
    if (id.startsWith("setup/boundary/")) return pages::boundary(c, id.mid(15));
    if (id == "setup/reference") return pages::reference(c);
    if (id == "solution/methods") return pages::methods(c);
    if (id == "solution/controls") return pages::controls(c);
    if (id == "solution/reports") return pages::reportDefinitions(c);
    if (id == "solution/monitors") return pages::monitors(c);
    if (id == "solution/initialization") return pages::initialization(c);
    if (id == "solution/run") return pages::run(c);
    if (id == "results/contours") return pages::contours(c);
    if (id == "results/vectors") return pages::vectors(c);
    if (id == "results/pathlines") return pages::pathlines(c);
    if (id == "results/iso") return pages::isoSurfaces(c);
    if (id == "results/planes") return pages::planes(c);
    if (id == "results/mesh") return pages::mesh(c);
    if (id == "results/xy") return pages::xyPlot(c);
    if (id == "results/residuals") return pages::residuals(c);
    if (id == "results/report-plots") return pages::reportPlots(c);
    if (id == "results/forces") return pages::forces(c);
    if (id == "results/integrals") return pages::integrals(c);
    if (id == "runs/history") return pages::runHistory(c);
    if (id == "runs/compare") return pages::compare(c);
    if (id == "runs/studies") return pages::studies(c);
    if (id == "results/aero") return pages::aero(c);
    if (id == "settings") return c.session->hasCase() ? pages::allSettings(c, m_settingsFilter) : pages::noCase(c, tr("All Case Settings"));
    return pages::noCase(c, Outline::titleOf(id));
}

namespace pages {

namespace {

QWidget* makePage(const QString& title, const QString& crumb, QVBoxLayout** bodyOut)
{
    auto* page = new QWidget;
    auto* v = new QVBoxLayout(page);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    v->addWidget(header(title, crumb, page));
    auto* content = new QWidget(page);
    *bodyOut = body(content);
    v->addWidget(content);
    v->addStretch(1);
    return page;
}

// A settings form page: header plus one SettingsForm.
SettingsForm* formPage(const InspectorContext& c, const QString& title, const QString& crumb, QWidget** pageOut)
{
    auto* page = new QWidget;
    auto* v = new QVBoxLayout(page);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    v->addWidget(header(title, crumb, page));
    auto* form = new SettingsForm(c.session->settings(), page);
    v->addWidget(form);
    v->addStretch(1);
    *pageOut = page;
    return form;
}

std::vector<const Setting*> settingsUnder(const InspectorContext& c, const QString& prefix, bool includeRaw = false)
{
    std::vector<const Setting*> out;
    for (const auto& s : c.session->settings()->settings())
        if (s.path.startsWith(prefix) && (includeRaw || !s.raw)) out.push_back(&s);
    return out;
}

int rawCount(const InspectorContext& c, const QString& prefix)
{
    int n = 0;
    for (const auto& s : c.session->settings()->settings()) n += s.path.startsWith(prefix) && s.raw;
    return n;
}

QPushButton* rawLink(const InspectorContext& c, const QString& prefix, const QString& what)
{
    const int n = rawCount(c, prefix);
    if (n == 0) return nullptr;
    auto* b = button(QObject::tr("%1 retained Fluent settings").arg(QLocale().toString(n)), "list-filter");
    b->setProperty("flat", true);
    b->setToolTip(QObject::tr("Retained settings for %1: saved in the case and read by Kutta, but not exposed as controls.").arg(what));
    QObject::connect(b, &QPushButton::clicked, b, [c, prefix] { c.showSettings(prefix); });
    return b;
}

QString count(qint64 n) { return QLocale().toString(n); }

} // namespace

QWidget* header(const QString& title, const QString& crumb, QWidget* parent)
{
    auto* w = new QWidget(parent);
    auto* v = new QVBoxLayout(w);
    v->setContentsMargins(14, 12, 14, 6);
    v->setSpacing(1);
    if (!crumb.isEmpty()) v->addWidget(label(crumb, "faint", w));
    auto* t = label(title, "title", w);
    t->setAccessibleName(title);
    v->addWidget(t);
    return w;
}


namespace {

// The failed open, in place of the page's usual content: no modal dialog.
QWidget* openErrorCallout(const InspectorContext& c)
{
    const Session::OpenError& e = c.session->openError();
    if (e.message.isEmpty()) return nullptr;
    auto* f = callout("red");
    auto* top = new QHBoxLayout;
    auto* t = label(QObject::tr("Could not open %1").arg(QFileInfo(e.casePath).fileName()), "title");
    top->addWidget(t, 1);
    auto* dismiss = button(QObject::tr("Dismiss"), "x");
    dismiss->setProperty("flat", true);
    QObject::connect(dismiss, &QPushButton::clicked, dismiss, [c] {
        c.session->clearOpenError();
        QTimer::singleShot(0, c.session, [c] { if (c.refresh) c.refresh(); });
    });
    top->addWidget(dismiss, 0, Qt::AlignTop);
    static_cast<QVBoxLayout*>(f->layout())->addLayout(top);
    f->layout()->addWidget(label(e.message));
    f->layout()->addWidget(label(c.session->hasCase()
                                     ? QObject::tr("%1 is still open.").arg(c.session->caseName())
                                     : QObject::tr("Kutta reads Fluent cases (.cas.h5, .cas) and Kutta cases (.kutta.json). Check the file is one, and that it is readable."),
                                 "muted"));
    return f;
}

// Open, approve, run, read: the first-time path, one compact card. Each step
// says what is true now and offers the one action that moves it on; a power
// user can hide it (View ▸ Show Getting Started brings it back).
QWidget* gettingStarted(const InspectorContext& c, bool hideable)
{
    if (hideable && QSettings().value("ui/hideGettingStarted").toBool()) return nullptr;
    Session* s = c.session;
    auto* w = new QFrame;
    w->setProperty("callout", "neutral");
    w->setAccessibleName(QObject::tr("Getting started"));
    auto* v = new QVBoxLayout(w);
    v->setContentsMargins(10, 8, 10, 10);
    v->setSpacing(6);
    auto* head = new QHBoxLayout;
    head->addWidget(section(QObject::tr("Getting started")), 1);
    if (hideable) {
        auto* hide = button(QObject::tr("Hide"));
        hide->setProperty("flat", true);
        hide->setToolTip(QObject::tr("View ▸ Show Getting Started brings it back."));
        QObject::connect(hide, &QPushButton::clicked, hide, [c] {
            QSettings().setValue("ui/hideGettingStarted", true);
            QTimer::singleShot(0, c.session, [c] { if (c.refresh) c.refresh(); });
        });
        head->addWidget(hide);
    }
    v->addLayout(head);
    auto* grid = new QGridLayout;
    grid->setHorizontalSpacing(8);
    grid->setVerticalSpacing(6);
    grid->setColumnStretch(1, 1);
    v->addLayout(grid);

    enum class Mark { Done, Current, Blocked, Later };
    struct Link {
        QString text, page;
    };
    int row = 0;
    // Everything a step offers sits under its text, never beside it, so the
    // card is as narrow as the panel allows and wraps instead of clipping.
    auto step = [&](Mark m, const QString& title, const QString& detail, QPushButton* action = nullptr, const QList<Link>& links = {}) {
        const auto& t = tokens();
        auto* glyph = new QLabel;
        const QString icon = m == Mark::Done ? "circle-check" : m == Mark::Blocked ? "circle-alert" : "circle";
        const QColor col = m == Mark::Done ? t.accent : m == Mark::Blocked ? t.red : m == Mark::Current ? t.text : t.faint;
        glyph->setPixmap(Theme::instance().icon(icon, col).pixmap(16, 16));
        glyph->setFixedWidth(18);
        grid->addWidget(glyph, row, 0, Qt::AlignTop);
        auto* text = new QWidget;
        auto* tv = new QVBoxLayout(text);
        tv->setContentsMargins(0, 0, 0, 0);
        tv->setSpacing(1);
        auto* tl = label(QStringLiteral("%1. %2").arg(row + 1).arg(title), m == Mark::Later ? "faint" : QString());
        if (m == Mark::Current || m == Mark::Blocked) {
            QFont f = tl->font();
            f.setWeight(QFont::DemiBold);
            tl->setFont(f);
        }
        tv->addWidget(tl);
        if (!detail.isEmpty()) tv->addWidget(label(detail, m == Mark::Blocked ? "error" : "muted"));
        if (!links.isEmpty()) {
            // One wrapping line of page links rather than a row of buttons.
            QStringList anchors;
            for (const Link& l : links)
                anchors << QStringLiteral("<a href=\"%1\" style=\"text-decoration:none\">%2</a>").arg(l.page.toHtmlEscaped(), l.text.toHtmlEscaped());
            auto* ll = new QLabel(anchors.join(QStringLiteral(" &nbsp;·&nbsp; ")));
            ll->setTextFormat(Qt::RichText);
            ll->setWordWrap(true);
            ll->setTextInteractionFlags(Qt::TextBrowserInteraction);
            ll->setAccessibleName(QObject::tr("Go to"));
            QObject::connect(ll, &QLabel::linkActivated, ll, [c](const QString& page) { c.navigate(page); });
            tv->addWidget(ll);
        }
        if (action) {
            tv->addSpacing(3);
            tv->addWidget(action, 0, Qt::AlignLeft);
        }
        grid->addWidget(text, row, 1);
        ++row;
    };
    auto primary = [&](const QString& text, const QString& icon, std::function<void()> fn) {
        auto* b = button(text, icon, "primary");
        QObject::connect(b, &QPushButton::clicked, b, std::move(fn));
        return b;
    };

    // 1. Open
    if (!s->hasCase()) {
        step(Mark::Current, QObject::tr("Open your Fluent case"), QObject::tr("The .cas.h5; add its .dat.h5 to compare with Fluent."),
             primary(QObject::tr("Open Case…"), "folder-open", [c] { if (c.openDialog) c.openDialog(); }));
        step(Mark::Later, QObject::tr("Review what Kutta needs from you"), {});
        step(Mark::Later, QObject::tr("Run the calculation"), {});
        step(Mark::Later, QObject::tr("Read the results"), {});
        return w;
    }
    const bool kutta = s->caseSummary().kuttaCase();
    step(Mark::Done, kutta ? QObject::tr("Open your case") : QObject::tr("Open your Fluent case"),
         kutta                    ? QObject::tr("%1, a Kutta case on %2").arg(breakable(s->caseName()), breakable(QFileInfo(s->caseSummary().meshFile).fileName()))
         : s->dataPath().isEmpty() ? QObject::tr("%1, without Fluent data (no comparison)").arg(breakable(s->caseName()))
                                 : QObject::tr("%1 with Fluent's data, for comparison").arg(breakable(s->caseName())));

    // 2. Approvals
    const CaseSummary& cs = s->caseSummary();
    const QStringList approved = s->approvedApproximations();
    const auto blocking = cs.blocking(approved);
    const bool unsupported = std::any_of(blocking.begin(), blocking.end(), [](const SupportIssue* i) { return i->cls == IssueClass::Unsupported; });
    const bool onGeneral = c.current && c.current() == "setup/general";
    const QList<Link> review = onGeneral ? QList<Link>{} : QList<Link>{{QObject::tr("Review on General"), "setup/general"}};
    if (unsupported) {
        step(Mark::Blocked, QObject::tr("Review what Kutta needs from you"),
             QObject::tr("Kutta cannot solve this case as it stands; the reasons are listed %1.").arg(onGeneral ? QObject::tr("below") : QObject::tr("on General")),
             nullptr, review);
    } else if (!blocking.empty()) {
        step(Mark::Current, QObject::tr("Review what Kutta needs from you"),
             blocking.size() == 1 ? QObject::tr("1 approximation needs your decision%1.").arg(onGeneral ? QObject::tr(" below") : QString())
                                  : QObject::tr("%1 approximations need your decision%2.").arg(blocking.size()).arg(onGeneral ? QObject::tr(" below") : QString()),
             nullptr, review);
    } else {
        step(Mark::Done, QObject::tr("Review what Kutta needs from you"),
             approved.isEmpty() ? (kutta ? QObject::tr("Nothing to approve: Kutta solves this case as its setup defines it.")
                                         : QObject::tr("Nothing to approve: Kutta solves this case as Fluent defines it."))
                                : QObject::tr("Approved; results will be marked APPROXIMATE."));
    }

    // 3. Run
    const auto& runs = s->history().runs();
    const bool finished = std::any_of(runs.begin(), runs.end(), [](const RunRecord& r) {
        return r.status == RunStatus::Completed || r.status == RunStatus::Stopped;
    });
    const QString keys = QKeySequence("Ctrl+R").toString(QKeySequence::NativeText);
    QString device = QObject::tr("the device kutta picks");
    for (const auto& d : s->devices())
        if (d.backend == s->backend()) device = d.label();
    if (s->running()) {
        step(Mark::Current, QObject::tr("Run the calculation"), QObject::tr("Running; residuals are in Charts below."));
    } else if (!blocking.empty()) {
        step(Mark::Later, QObject::tr("Run the calculation"), {});
    } else {
        step(finished ? Mark::Done : Mark::Current, QObject::tr("Run the calculation"),
             finished ? QObject::tr("%1 finished. Run again after changing settings (%2).").arg(runs.size() == 1 ? QObject::tr("1 run") : QObject::tr("%1 runs").arg(runs.size()), keys)
                      : QObject::tr("%1 iterations on %2 (%3).").arg(s->requestedIterations()).arg(device, keys),
             finished ? nullptr : primary(QObject::tr("Run"), "play", [s] { s->run(); }));
    }

    // 4. Results
    if (!finished) {
        step(Mark::Later, QObject::tr("Read the results"), QObject::tr("Contours, forces and the comparison with Fluent."));
    } else {
        step(Mark::Current, QObject::tr("Read the results"),
             s->dataPath().isEmpty() ? QObject::tr("The viewports show the latest run; Forces lists the report definitions.")
                                     : QObject::tr("The viewports show the latest run; pick a Fluent or Δ field to compare."),
             nullptr,
             {{QObject::tr("Contours"), "results/contours"}, {QObject::tr("Forces"), "results/forces"}, {QObject::tr("Run History"), "runs/history"}});
    }
    return w;
}

} // namespace

QWidget* noCase(const InspectorContext& c, const QString& title)
{
    QVBoxLayout* v;
    if (c.session->hasCase()) {
        QWidget* page = makePage(title, {}, &v);
        v->addWidget(label(QObject::tr("Nothing to show here yet."), "muted"));
        return page;
    }
    // Welcome: what Kutta does, and the ways to start.
    QWidget* page = makePage(QObject::tr("Open a case"), {}, &v);
    if (QWidget* e = openErrorCallout(c)) v->addWidget(e);
    v->addWidget(label(QObject::tr("Kutta reads a Fluent case (.cas.h5), solves it on this computer's GPU and shows the results next to Fluent's. Open a case, or drop one on the window."), "muted"));
    v->addWidget(gettingStarted(c, false));
    // No Fluent case: start from the mesh alone.
    auto* fromMesh = new QHBoxLayout;
    auto* nb = button(QObject::tr("New Case from Mesh…"), "plus");
    nb->setProperty("flat", true);
    QObject::connect(nb, &QPushButton::clicked, nb, [c] { if (c.newCaseDialog) c.newCaseDialog(); });
    fromMesh->addWidget(nb);
    fromMesh->addWidget(label(QObject::tr("No Fluent case? Start from your mesh (.msh) and set it up here."), "faint"), 1);
    v->addLayout(fromMesh);
    if (!c.samples.isEmpty()) {
        v->addSpacing(6);
        v->addWidget(section(QObject::tr("Samples")));
        auto* grid = new QGridLayout;
        grid->setHorizontalSpacing(10);
        grid->setVerticalSpacing(4);
        grid->setColumnStretch(1, 1);
        int row = 0;
        for (const QString& p : c.samples) {
            const QString name = QFileInfo(p).fileName();
            const bool aero = name.startsWith("fsae_car") || name.startsWith("naca");
            auto* b = button(name, aero ? "wind" : "waves");
            b->setProperty("flat", true);
            QObject::connect(b, &QPushButton::clicked, b, [c, p] { c.openPath(p); });
            grid->addWidget(b, row, 0, Qt::AlignLeft | Qt::AlignTop);
            const QString replays = QObject::tr("Replays a recorded Kutta run.");
            const QString about = name.startsWith("fsae_car")     ? QObject::tr("Synthetic FSAE car for trying the workbench. Not a CFD result.")
                                  : name.startsWith("naca")       ? QObject::tr("2D ground-effect wing: SST with intermittency transition, moving ground. %1").arg(replays)
                                  : name.startsWith("elbow_param") ? QObject::tr("Ansys elbow: standard k-epsilon, SIMPLE. %1").arg(replays)
                                  : name.startsWith("cylinder")   ? QObject::tr("Ansys cylinder: laminar, 2D. %1").arg(replays)
                                                                  : QObject::tr("Ansys mixing elbow: SST, Coupled. %1").arg(replays);
            grid->addWidget(label(about, "faint"), row++, 1);
        }
        v->addLayout(grid);
    }
    const QStringList recent = c.session->recentCases();
    if (!recent.isEmpty()) {
        v->addSpacing(6);
        v->addWidget(section(QObject::tr("Recent")));
        for (const QString& p : recent.mid(0, 6)) {
            auto* b = button(QFileInfo(p).fileName(), "file-text");
            b->setProperty("flat", true);
            b->setToolTip(QDir::toNativeSeparators(p));
            b->setEnabled(QFileInfo::exists(p));
            QObject::connect(b, &QPushButton::clicked, b, [c, p] { c.openPath(p); });
            v->addWidget(b, 0, Qt::AlignLeft);
        }
    }
    return page;
}

// ---- Setup ----------------------------------------------------------------------------

QWidget* general(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("General"), QObject::tr("Setup"), &v);
    const CaseSummary& cs = c.session->caseSummary();

    // Support status first: it decides whether and how this case can run.
    auto* statusRow = new QHBoxLayout;
    statusRow->setSpacing(8);
    const QStringList approved = c.session->approvedApproximations();
    const auto blocking = cs.blocking(approved);
    const bool unsupported = std::any_of(blocking.begin(), blocking.end(), [](const SupportIssue* i) { return i->cls == IssueClass::Unsupported; });
    if (cs.status == SupportStatus::Supported) statusRow->addWidget(badge("SUPPORTED", "supported"));
    else if (unsupported) statusRow->addWidget(badge("NOT SUPPORTED", "unsupported"));
    else statusRow->addWidget(badge("APPROXIMATE", "approximate"));
    if (cs.synthetic) statusRow->addWidget(badge("SYNTHETIC", "modified"));
    if (c.session->cliInfo().fixture) statusRow->addWidget(badge("FIXTURE CLI", "neutral"));
    QString statusText = cs.status == SupportStatus::Supported ? (cs.kuttaCase() ? QObject::tr("Kutta solves this case as its setup defines it.")
                                                                                 : QObject::tr("Kutta solves this case as Fluent defines it."))
                         : unsupported ? QObject::tr("Kutta cannot solve this case; see below.")
                                       : QObject::tr("Kutta can solve this case only with the approximations below. Results are marked APPROXIMATE everywhere.");
    statusRow->addStretch(1);
    v->addLayout(statusRow);
    v->addWidget(label(statusText, "muted"));
    if (QWidget* e = openErrorCallout(c)) v->addWidget(e);
    if (QWidget* g = gettingStarted(c, true)) v->addWidget(g);

    std::vector<const SupportIssue*> approx, blockers, notes;
    for (const auto& i : cs.issues) {
        if (i.cls == IssueClass::Unsupported) blockers.push_back(&i);
        else if (i.cls == IssueClass::Approximable) approx.push_back(&i);
        else notes.push_back(&i);
    }
    for (const auto* i : blockers) {
        auto* f = callout("red");
        f->layout()->addWidget(label(i->label, "title"));
        f->layout()->addWidget(label(i->explanation));
        if (!i->evidence.isEmpty()) f->layout()->addWidget(label(i->evidence, "muted"));
        f->layout()->addWidget(label(i->setting, "faint"));
        v->addWidget(f);
    }
    for (const auto* i : approx) {
        const bool on = approved.contains(i->approximation);
        auto* f = callout(on ? "amber" : "neutral");
        auto* top = new QHBoxLayout;
        auto* t = label(i->label);
        QFont bf = t->font();
        bf.setWeight(QFont::DemiBold);
        t->setFont(bf);
        top->addWidget(t, 1);
        if (on) top->addWidget(badge("APPROXIMATE", "approximate"));
        static_cast<QVBoxLayout*>(f->layout())->addLayout(top);
        f->layout()->addWidget(label(i->explanationOnly()));
        // The measured effect is what the user weighs before opting in; older
        // CLIs fold it into the explanation above.
        if (!i->measuredEffect.isEmpty()) {
            QString effect = i->measuredEffect;
            effect[0] = effect[0].toUpper();
            auto* head = label(QObject::tr("Measured effect"), "faint");
            QFont hf = head->font();
            hf.setWeight(QFont::DemiBold);
            head->setFont(hf);
            f->layout()->addWidget(head);
            f->layout()->addWidget(label(effect + '.', "muted"));
        }
        if (!i->evidence.isEmpty()) {
            auto* ev = label(QObject::tr("Evidence: %1").arg(i->evidence), "faint");
            ev->setTextInteractionFlags(Qt::TextSelectableByMouse);
            f->layout()->addWidget(ev);
        }
        auto* cb = new QCheckBox(QObject::tr("Run with this approximation"), f);
        cb->setChecked(on);
        cb->setAccessibleDescription(QObject::tr("Results will be marked APPROXIMATE in every view, chart and export."));
        cb->setToolTip(QObject::tr("%1\nkutta solve --approximate %2").arg(cb->accessibleDescription(), i->approximation));
        const QString name = i->approximation;
        QObject::connect(cb, &QCheckBox::toggled, cb, [c, name](bool b) {
            // Defer: the page is rebuilt on approval changes.
            QTimer::singleShot(0, c.session, [c, name, b] { c.session->setApproximationApproved(name, b); });
        });
        f->layout()->addWidget(cb);
        v->addWidget(f);
    }

    v->addWidget(section(QObject::tr("Case")));
    auto* grid = new QGridLayout;
    grid->setHorizontalSpacing(12);
    grid->setVerticalSpacing(5);
    int row = 0;
    auto addLabel = [&](const QString& k, QLabel* l) {
        grid->addWidget(label(k, "muted"), row, 0, Qt::AlignTop);
        l->setTextInteractionFlags(Qt::TextSelectableByMouse);
        grid->addWidget(l, row++, 1);
    };
    auto addRow = [&](const QString& k, const QString& val) { addLabel(k, label(val)); };
    grid->setColumnStretch(1, 1);
    addLabel(QObject::tr("File"), pathLabel(c.session->casePath()));
    if (c.session->dataPath().isEmpty()) addRow(QObject::tr("Data"), QObject::tr("none (initialized by Kutta)"));
    else addLabel(QObject::tr("Data"), pathLabel(c.session->dataPath()));
    if (cs.kuttaCase()) {
        addRow(QObject::tr("Format"), QObject::tr("Kutta case (mesh and JSON setup)"));
        addLabel(QObject::tr("Mesh"), pathLabel(cs.meshFile));
    } else {
        addRow(QObject::tr("Written by"), cs.fluentVersion.isEmpty() ? QObject::tr("unknown Fluent version") : "Fluent " + cs.fluentVersion);
    }
    addRow(QObject::tr("SHA-256"), cs.sha256.left(16) + "…");
    v->addLayout(grid);

    v->addWidget(section(QObject::tr("Mesh")));
    auto* mg = new QGridLayout;
    mg->setHorizontalSpacing(12);
    mg->setVerticalSpacing(5);
    row = 0;
    grid = mg;
    addRow(QObject::tr("Cells"), count(cs.cells));
    QStringList types;
    for (auto it = cs.cellTypes.begin(); it != cs.cellTypes.end(); ++it)
        if (it.value() > 0) types << QStringLiteral("%1 %2").arg(count(it.value()), it.key());
    if (!types.isEmpty()) addRow(QObject::tr("Cell types"), types.join(", "));
    addRow(QObject::tr("Faces"), count(cs.faces));
    addRow(QObject::tr("Nodes"), count(cs.nodes));
    if (cs.bounds.valid) {
        const QVector3D e = cs.bounds.max - cs.bounds.min;
        addRow(QObject::tr("Extent"), QStringLiteral("%1 × %2 × %3 m").arg(formatNumber(e.x(), 4), formatNumber(e.y(), 4), formatNumber(e.z(), 4)));
    }
    mg->setColumnStretch(1, 1);
    v->addLayout(mg);

    v->addWidget(section(QObject::tr("Solver")));
    auto* form = new SettingsForm(c.session->settings());
    form->layout()->setContentsMargins(0, 0, 0, 0);
    form->addRow(QObject::tr("Type"), label(QObject::tr("Pressure-based, steady (as read from the case)")));
    form->addSetting("general/operating-pressure");
    if (cs.dimension) {
        form->addRow(QObject::tr("Dimension"), label(cs.planar() ? QObject::tr("2D (planar)") : QObject::tr("%1D").arg(*cs.dimension)));
        if (cs.planar() && cs.bounds.valid)
            form->addNote(QObject::tr("Kutta solves a 2D case as one cell through its %1 m reference depth, between two symmetry planes: "
                                      "the same planar equations. Forces are for that depth, as in Fluent; views show the plane face on.")
                              .arg(formatNumber(cs.bounds.max.z() - cs.bounds.min.z(), 4)));
    } else {
        form->addSetting("general/dimension");
    }
    v->addWidget(form);

    if (!notes.empty()) {
        v->addWidget(section(QObject::tr("Notes")));
        // The count and what the notes are about wrap as text; the button
        // stays short enough for the narrowest panel.
        v->addWidget(label(notes.size() == 1 ? QObject::tr("1 note on a setting that does not affect this case.")
                                             : QObject::tr("%1 notes on settings that do not affect this case.").arg(notes.size()),
                           "muted"));
        const QString showText = notes.size() == 1 ? QObject::tr("Show Note") : QObject::tr("Show Notes");
        auto* toggle = button(showText, "info");
        toggle->setProperty("flat", true);
        toggle->setCheckable(true);
        QObject::connect(toggle, &QPushButton::toggled, toggle, [toggle, showText](bool on) { toggle->setText(on ? QObject::tr("Hide") : showText); });
        auto* list = new QWidget;
        auto* lv = new QVBoxLayout(list);
        lv->setContentsMargins(0, 0, 0, 0);
        for (const auto* i : notes) {
            auto* f = callout("neutral");
            f->layout()->addWidget(label(i->label));
            f->layout()->addWidget(label(i->explanation, "muted"));
            lv->addWidget(f);
        }
        list->setVisible(false);
        QObject::connect(toggle, &QPushButton::toggled, list, &QWidget::setVisible);
        v->addWidget(toggle, 0, Qt::AlignLeft);
        v->addWidget(list);
    }
    if (auto* b = rawLink(c, "case/", QObject::tr("the case"))) v->addWidget(b, 0, Qt::AlignLeft);
    return page;
}

QWidget* models(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Models"), QObject::tr("Setup"), &page);
    const CaseSummary& cs = c.session->caseSummary();
    f->addSection(QObject::tr("Viscous"));
    if (!f->addSetting("models/viscous")) f->addRow(QObject::tr("Viscous model"), label(cs.viscousModel));
    // Fluent's Viscous Model dialog: the transition model belongs to the closure.
    if (f->addSetting("models/transition", QObject::tr("Transition"))) {
        const QJsonValue t = c.session->settings()->value("models/transition");
        bool blocked = false;
        for (const auto& i : cs.issues)
            if (i.cls == IssueClass::Unsupported && (i.setting.contains("rp-trans") || i.setting.contains("rp-transition-model"))) {
                f->addNote(i.explanation, "error");
                blocked = true;
            }
        if (t == QJsonValue("gamma-transport-eqn")) {
            f->addNote(QObject::tr("Fluent's Intermittency Transition Model: one transport equation for intermittency γ on top of "
                                   "k-omega SST, which damps turbulence production in laminar boundary layers until they transition. "
                                   "Kutta solves it with the case's constants below and adds the intermit residual and the Intermittency field."));
            // The model constants as saved; Kutta reads and checks them.
            const std::pair<const char*, const char*> constants[] = {
                {"trans-one-eqn-ctu1", "CTU1"}, {"trans-one-eqn-ctu2", "CTU2"}, {"trans-one-eqn-ctu3", "CTU3"},
                {"trans-one-eqn-cpg1", "CPG1"}, {"trans-one-eqn-cpg2", "CPG2"}, {"trans-one-eqn-cpg3", "CPG3"},
                {"trans-one-eqn-lim-cpg1", "CPG1 lim"}, {"trans-one-eqn-lim-cpg2", "CPG2 lim"},
                {"trans-one-eqn-csep", "CSEP"}, {"intermit-sig", "σ (intermittency Prandtl number)"}};
            for (const auto& [key, name] : constants) f->addSetting(QStringLiteral("case/%1").arg(key), QString::fromUtf8(name));
            f->addNote(QObject::tr("Its discretization is under Methods, its relaxation under Controls and the inflow intermittency on each "
                                   "inlet and outlet. Changing the model or its constants is read-only in Kutta."), "faint");
        } else if (!blocked && cs.viscousModel.contains("sst")) {
            f->addNote(QObject::tr("Fully turbulent. Kutta also runs cases saved with Fluent's Intermittency Transition Model; "
                                   "Transition SST and the algebraic model are not implemented."), "faint");
        }
    }
    f->addSection(QObject::tr("Energy"));
    f->addSetting("models/energy");
    for (const auto& i : cs.issues)
        if (i.setting.contains("energy")) f->addNote(i.explanation);
    for (const Setting* s : settingsUnder(c, "models/"))
        if (s->path != "models/viscous" && s->path != "models/energy" && s->path != "models/transition") f->addSetting(s->path);
    f->finish();
    return page;
}

QWidget* materials(const InspectorContext& c, const QString& only)
{
    QWidget* page;
    SettingsForm* f = formPage(c, only.isEmpty() ? QObject::tr("Materials") : only,
                               only.isEmpty() ? QObject::tr("Setup") : QObject::tr("Setup ▸ Materials"), &page);
    for (const auto& m : c.session->caseSummary().materials) {
        if (!only.isEmpty() && m.name != only) continue;
        const QString prefix = "materials/" + m.name + "/";
        const Setting* type = c.session->settings()->find(prefix + "type");
        f->addSection(type ? QStringLiteral("%1 · %2").arg(m.name, type->value.toString()) : m.name);
        f->addSetting(prefix + "density");
        f->addSetting(prefix + "viscosity");
        for (const Setting* s : settingsUnder(c, prefix))
            if (!s->path.endsWith("/density") && !s->path.endsWith("/viscosity") && !s->path.endsWith("/type")) f->addSetting(s->path);
    }
    f->finish();
    return page;
}

QWidget* cellZones(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Cell Zone Conditions"), QObject::tr("Setup"), &page);
    QString owner;
    for (const Setting* s : settingsUnder(c, "cell-zone/")) {
        if (s->owner != owner) {
            owner = s->owner;
            f->addSection(owner);
        }
        f->addSetting(s->path);
    }
    if (auto* b = rawLink(c, "cell-zone/", QObject::tr("cell zones"))) f->addWidget(b);
    f->finish();
    return page;
}

QWidget* boundaries(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Boundary Conditions"), QObject::tr("Setup"), &v);
    const CaseSummary& cs = c.session->caseSummary();
    auto* table = new QTableWidget(0, 4);
    table->setHorizontalHeaderLabels({QObject::tr("Zone"), QObject::tr("Type"), QObject::tr("ID"), QObject::tr("Faces")});
    table->verticalHeader()->hide();
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setShowGrid(false);
    table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    table->setAccessibleName(QObject::tr("Boundary zones"));
    for (const auto& z : cs.zones) {
        if (z.type == "interior") continue;
        const int r = table->rowCount();
        table->insertRow(r);
        auto* name = new QTableWidgetItem(z.name);
        const int mods = [&] {
            int n = 0;
            for (const QString& p : c.session->settings()->modifiedPaths()) n += p.startsWith("boundary/" + z.name + "/");
            return n;
        }();
        if (mods) name->setForeground(tokens().dark ? tokens().amber : tokens().amberText);
        table->setItem(r, 0, name);
        table->setItem(r, 1, new QTableWidgetItem(z.type));
        auto* id = new QTableWidgetItem(QString::number(z.id));
        id->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        table->setItem(r, 2, id);
        auto* faces = new QTableWidgetItem(count(z.faces));
        faces->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        table->setItem(r, 3, faces);
    }
    table->resizeColumnsToContents();
    table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    table->setMinimumHeight(std::min(420, 30 + table->rowCount() * 26));
    QObject::connect(table, &QTableWidget::cellActivated, table, [c, table](int r) { c.navigate("setup/boundary/" + table->item(r, 0)->text()); });
    QObject::connect(table, &QTableWidget::cellClicked, table, [c, table](int r) { c.navigate("setup/boundary/" + table->item(r, 0)->text()); });
    v->addWidget(table);
    v->addWidget(label(QObject::tr("Choose a zone to edit its conditions. Interior zones have no conditions."), "muted"));
    return page;
}

QWidget* boundary(const InspectorContext& c, const QString& zoneName)
{
    QWidget* page;
    SettingsForm* f = formPage(c, zoneName, QObject::tr("Setup ▸ Boundary Conditions"), &page);
    const Zone* z = c.session->caseSummary().zone(zoneName);
    const QString prefix = "boundary/" + zoneName + "/";
    QStringList placed;
    auto add = [&](const QString& leaf, const QString& label = {}) {
        if (f->addSetting(prefix + leaf, label)) placed << prefix + leaf;
    };
    const QString type = z ? z->type : QString();
    f->addRow(QObject::tr("Type"), label(type.isEmpty() ? QObject::tr("unknown") : type));
    if (z) f->addRow(QObject::tr("Zone ID · faces"), label(QStringLiteral("%1 · %2").arg(z->id).arg(count(z->faces))));
    placed << prefix + "type";
    if (type.contains("inlet") || type.contains("outlet") || type.contains("far-field")) {
        f->addSection(QObject::tr("Momentum"));
        // Fluent's panel words; the section heading says which method.
        add("velocity-specification", QObject::tr("Velocity specification method"));
        add("velocity-magnitude");
        add("velocity");
        add("direction");
        add("gauge-pressure");
        f->addSection(type.contains("outlet") ? QObject::tr("Turbulence (backflow)") : QObject::tr("Turbulence"));
        add("turbulence-specification", QObject::tr("Specification method"));
        add("turb-intensity");
        add("turb-hydraulic-diam");
        add("turb-length-scale");
        add("turb-viscosity-ratio");
        add("k");
        add("omega");
        add("epsilon");
    } else if (type == "wall" && z) {
        f->addSection(QObject::tr("Wall"));
        const QString motion = z->condition.value("wall_motion").toString();
        f->addRow(QObject::tr("Wall motion"), label(motion.isEmpty() ? QStringLiteral("—") : motion));
        placed << prefix + "wall-motion" << prefix + "roughness-height";
        if (motion == "rotating") {
            placed << prefix + "rotation-speed";
            const auto omega = z->condition.value("angular_speed");
            f->addRow(QObject::tr("Angular speed"), label(omega.isDouble() ? formatNumber(omega.toDouble(), 6) : QStringLiteral("—")), "rad/s");
            for (const auto& item : {std::pair{"axis_origin", QObject::tr("Axis origin")}, std::pair{"axis_direction", QObject::tr("Axis direction")}}) {
                const auto a = z->condition.value(item.first).toArray();
                if (a.size() == 3)
                    f->addRow(item.second, label(QStringLiteral("(%1, %2, %3)").arg(formatNumber(a[0].toDouble(), 6), formatNumber(a[1].toDouble(), 6), formatNumber(a[2].toDouble(), 6))),
                              QString(item.first) == "axis_origin" ? "m" : "");
            }
        }
        const QJsonValue rough = z->condition.value("roughness_height");
        f->addRow(QObject::tr("Roughness height"), label(rough.isDouble() ? formatNumber(rough.toDouble(), 4) : QStringLiteral("—")), "m");
        const QJsonValue vel = z->condition.value("velocity");
        if (vel.isArray() && vel.toArray().at(0).isDouble()) {
            const QJsonArray a = vel.toArray();
            f->addRow(QObject::tr("Wall velocity"), label(QStringLiteral("(%1, %2, %3)").arg(formatNumber(a[0].toDouble(), 4), formatNumber(a[1].toDouble(), 4), formatNumber(a[2].toDouble(), 4))), "m/s");
        }
        for (const QString& l : {"velocity-specification", "velocity-magnitude", "velocity", "direction", "gauge-pressure", "turbulence-specification"}) placed << prefix + l;
    } else if (type == "symmetry") {
        f->addNote(QObject::tr("Symmetry boundaries have no inputs: zero normal velocity and zero normal gradients."));
        for (const QString& l : {"velocity-specification", "velocity-magnitude", "velocity", "direction", "gauge-pressure", "turbulence-specification"}) placed << prefix + l;
    }
    bool otherHeader = false;
    for (const Setting* s : settingsUnder(c, prefix)) {
        if (placed.contains(s->path)) continue;
        if (!otherHeader) { f->addSection(QObject::tr("Other")); otherHeader = true; }
        f->addSetting(s->path);
    }
    if (auto* b = rawLink(c, prefix, zoneName)) f->addWidget(b);
    f->finish();
    return page;
}

QWidget* reference(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Reference Values"), QObject::tr("Setup"), &page);
    f->addNote(QObject::tr("Used for force and pressure coefficients (Cl, Cd, Cp). They do not change the flow solution."));
    for (const QString& p : {"area", "density", "length", "pressure", "velocity", "viscosity"}) f->addSetting("reference/" + p);
    for (const Setting* s : settingsUnder(c, "reference/"))
        if (!QStringList{"area", "density", "length", "pressure", "velocity", "viscosity"}.contains(s->leaf)) f->addSetting(s->path);
    f->finish();
    return page;
}

// ---- Solution ------------------------------------------------------------------------

// What the selected pressure-velocity coupling does for this case, in one
// sentence: Kutta's Coupled solver needs pseudo time, so a case saved with
// Coupled but without it runs SIMPLEC, and such a case is not offered Coupled.
static QString couplingNote(const SettingsModel& m)
{
    const Setting* s = m.find("methods/pressure-velocity-coupling");
    if (!s) return {};
    const QString v = m.value(s->path).toString();
    const QJsonValue pt = m.value("methods/pseudo-transient");
    const bool pseudo = pt.isBool() && pt.toBool();
    const bool offered = s->options.contains(QJsonValue("coupled"));
    if (v == "coupled" && pt.isBool() && !pseudo)
        return QObject::tr("Saved with Coupled but without the pseudo-transient formulation, which Kutta's Coupled solver needs: "
                           "this case runs SIMPLEC. Same discrete equations; only the convergence path differs.");
    if (v == "coupled")
        return QObject::tr("Pressure and momentum are solved together, with pseudo-transient steps and the flow Courant number "
                           "as implicit relaxation (Kutta's native Coupled solver; see Controls).");
    QString note = v == "simplec" ? QObject::tr("Segregated, with SIMPLEC's consistent pressure correction.")
                                  : QObject::tr("Segregated, with SIMPLE's pressure correction.");
    if (!offered) note += ' ' + QObject::tr("Coupled is not offered: Kutta's Coupled solver needs the pseudo-transient formulation, which this case does not use.");
    return note;
}

QWidget* methods(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Methods"), QObject::tr("Solution"), &page);
    SettingsModel* m = c.session->settings();
    QStringList placed;
    auto add = [&](const QString& leaf, const QString& label = {}) {
        if (f->addSetting("methods/" + leaf, label)) placed << leaf;
    };
    // Fluent's Methods panel, top to bottom.
    f->addSection(QObject::tr("Pressure-Velocity Coupling"));
    add("pressure-velocity-coupling", QObject::tr("Scheme"));
    add("pseudo-transient", QObject::tr("Pseudo-transient formulation"));
    if (m->find("methods/pressure-velocity-coupling")) {
        auto* note = new QLabel(couplingNote(*m));
        note->setProperty("role", "note");
        note->setWordWrap(true);
        f->addWidget(note);
        QObject::connect(m, &SettingsModel::valueChanged, note, [m, note](const QString& p) {
            if (p == "methods/pressure-velocity-coupling") note->setText(couplingNote(*m));
        });
        // What the last run actually ran, from its setup event.
        const auto& runs = c.session->history().runs();
        auto last = std::find_if(runs.rbegin(), runs.rend(), [](const RunRecord& r) { return !r.coupling.isEmpty(); });
        if (last != runs.rend()) f->addNote(QObject::tr("%1 ran %2.").arg(last->title(), couplingLabel(last->coupling)), "faint");
    }
    f->addSection(QObject::tr("Spatial Discretization"));
    add("gradient");
    add("pressure-scheme");
    add("momentum-scheme");
    add("k-scheme");
    add("omega-scheme");
    add("epsilon-scheme");
    add("intermittency-scheme");
    add("gradient-limiting");
    bool other = false;
    for (const Setting* s : settingsUnder(c, "methods/")) {
        if (placed.contains(s->leaf)) continue;
        if (!other) { f->addSection(QObject::tr("Other")); other = true; }
        f->addSetting(s->path);
    }
    f->finish();
    return page;
}

QWidget* controls(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Controls"), QObject::tr("Solution"), &page);
    f->addSection(QObject::tr("Under-relaxation factors"));
    f->addSetting("controls/pressure-relaxation");
    f->addSetting("controls/mom-relaxation");
    f->addSetting("controls/k-relaxation");
    f->addSetting("controls/omega-relaxation");
    f->addSetting("controls/epsilon-relaxation");
    f->addSetting("controls/intermit-relaxation");
    f->addSection(QObject::tr("Pseudo time"));
    f->addSetting("controls/pseudo-time-scale-factor");
    // Fluent's Flow Courant Number, the Coupled solver's implicit relaxation.
    SettingsModel* m = c.session->settings();
    if (f->addSetting("controls/courant-number")) {
        auto* note = new QLabel(QObject::tr("The Coupled solver applies the flow Courant number as implicit relaxation on top of the "
                                            "pseudo time step. Lower values are more robust and slower; Fluent's default is 200."));
        note->setProperty("role", "note");
        note->setWordWrap(true);
        f->addWidget(note);
        note->setVisible(m->isActive("controls/courant-number"));
        QObject::connect(m, &SettingsModel::valueChanged, note, [m, note](const QString&) { note->setVisible(m->isActive("controls/courant-number")); });
    }
    for (const Setting* s : settingsUnder(c, "controls/"))
        if (!s->leaf.endsWith("relaxation") && s->leaf != "pseudo-time-scale-factor" && s->leaf != "courant-number") f->addSetting(s->path);
    f->finish();
    return page;
}

QWidget* reportDefinitions(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Report Definitions"), QObject::tr("Solution"), &v);
    const auto& reps = c.session->caseSummary().reports;
    if (reps.empty()) {
        v->addWidget(label(QObject::tr("This case defines no force reports. Kutta still writes the force on every wall zone to forces.csv after each run; see Results ▸ Reports ▸ Forces."), "muted"));
        return page;
    }
    auto* t = new QTableWidget(int(reps.size()), 4);
    t->setHorizontalHeaderLabels({QObject::tr("Name"), QObject::tr("Kind"), QObject::tr("Zones"), QObject::tr("Direction")});
    t->verticalHeader()->hide();
    t->setEditTriggers(QAbstractItemView::NoEditTriggers);
    t->setShowGrid(false);
    for (int i = 0; i < int(reps.size()); ++i) {
        const auto& r = reps[std::size_t(i)];
        t->setItem(i, 0, new QTableWidgetItem(r.name));
        t->setItem(i, 1, new QTableWidgetItem(r.kind));
        t->setItem(i, 2, new QTableWidgetItem(r.zones.join(", ")));
        t->setItem(i, 3, new QTableWidgetItem(r.direction ? QStringLiteral("(%1, %2, %3)").arg(r.direction->x()).arg(r.direction->y()).arg(r.direction->z()) : QStringLiteral("—")));
    }
    t->resizeColumnsToContents();
    t->horizontalHeader()->setStretchLastSection(true);
    t->setMinimumHeight(40 + int(reps.size()) * 26);
    t->setAccessibleName(QObject::tr("Report definitions"));
    v->addWidget(t);
    v->addWidget(label(QObject::tr("Report definitions come from the case. Values are plotted each iteration under Results ▸ Plots ▸ Report Plots."), "muted"));
    return page;
}

QWidget* initialization(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Initialization"), QObject::tr("Solution"), &page);
    f->addSetting("initialization/method");
    f->addSection(QObject::tr("Start from Fluent data"));
    const QString data = c.session->dataPath();
    f->addRow(QObject::tr("Data file"), label(data.isEmpty() ? QObject::tr("none") : QFileInfo(data).fileName()));
    f->addNote(data.isEmpty() ? QObject::tr("Without a data file Kutta initializes the flow itself. With one, it starts from Fluent's solution and the viewports can show Fluent's fields and the difference.")
                              : QObject::tr("Kutta starts from this Fluent solution. Its fields are also the reference for comparisons (Fluent and Δ fields in the viewports)."));
    auto* pick = button(data.isEmpty() ? QObject::tr("Use a data file…") : QObject::tr("Use another data file…"), "file-input");
    QObject::connect(pick, &QPushButton::clicked, pick, [c, pick] {
        const QString start = QFileInfo(c.session->casePath()).absolutePath();
        const QString file = QFileDialog::getOpenFileName(pick->window(), QObject::tr("Fluent data file"), start, QObject::tr("Fluent data (*.dat.h5 *.dat)"));
        if (!file.isEmpty()) c.session->openCase(c.session->casePath(), file);
    });
    f->addWidget(pick);
    if (!data.isEmpty()) {
        auto* clear = button(QObject::tr("Start without data"), "x");
        clear->setProperty("flat", true);
        QObject::connect(clear, &QPushButton::clicked, clear, [c] { c.session->openCase(c.session->casePath()); });
        f->addWidget(clear);
    }
    f->finish();
    return page;
}

QWidget* run(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Run Calculation"), QObject::tr("Solution"), &page);
    Session* s = c.session;
    f->addSection(QObject::tr("Calculation"));
    f->addSetting("run/iterations");

    auto* device = new QComboBox;
    device->setAccessibleName(QObject::tr("Device"));
    device->addItem(QObject::tr("Automatic"), "auto");
    for (const auto& d : s->devices()) device->addItem(d.label(), d.backend);
    device->setCurrentIndex(std::max(0, device->findData(s->backend())));
    QObject::connect(device, &QComboBox::activated, device, [s, device](int i) { s->setBackend(device->itemData(i).toString()); });
    f->addRow(QObject::tr("Device"), device);

    auto* live = new QComboBox;
    live->setAccessibleName(QObject::tr("Live contours"));
    live->addItem(QObject::tr("Off"), 0);
    for (int n : {20, 50, 100, 200}) live->addItem(QObject::tr("Every %1 iterations").arg(n), n);
    live->setCurrentIndex(std::max(0, live->findData(s->liveViewInterval())));
    live->setToolTip(QObject::tr("Kutta writes a field snapshot during the run so viewports update while it solves."));
    QObject::connect(live, &QComboBox::activated, live, [s, live](int i) { s->setLiveViewInterval(live->itemData(i).toInt()); });
    f->addRow(QObject::tr("Live contours"), live);

    const QStringList approved = s->approvedApproximations();
    if (!approved.isEmpty()) {
        auto* box = callout("amber");
        QStringList names;
        for (const auto& i : s->caseSummary().issues)
            if (approved.contains(i.approximation) && !names.contains(i.label)) names << i.label;
        box->layout()->addWidget(label(QObject::tr("This run will be APPROXIMATE (%1). Its results are marked wherever they appear.").arg(names.join(", "))));
        f->addWidget(box);
    }
    const int mods = s->settings()->caseModifiedCount();
    if (mods) {
        auto* box = callout("neutral");
        box->layout()->addWidget(label(mods == 1 ? QObject::tr("1 setting differs from the case. The run and its results are marked MODIFIED.")
                                                 : QObject::tr("%1 settings differ from the case. The run and its results are marked MODIFIED.").arg(mods)));
        f->addWidget(box);
    }

    auto* actions = new QWidget;
    auto* ah = new QHBoxLayout(actions);
    ah->setContentsMargins(0, 6, 0, 0);
    auto* runBtn = button(QObject::tr("Run Calculation"), "play", "primary");
    auto* stopBtn = button(QObject::tr("Stop"), "circle-stop", "danger");
    ah->addWidget(runBtn);
    ah->addWidget(stopBtn);
    ah->addStretch(1);
    f->addWidget(actions);
    auto* status = label({}, "muted");
    f->addWidget(status);
    auto* why = label({}, "error");
    f->addWidget(why);
    // A run that did not complete says why and what to do, until the next one.
    auto* outcome = callout("red");
    auto* outcomeTitle = label({}, "title");
    auto* outcomeError = label({});
    outcomeError->setTextInteractionFlags(Qt::TextSelectableByMouse);
    auto* outcomeAdvice = label({}, "muted");
    outcome->layout()->addWidget(outcomeTitle);
    outcome->layout()->addWidget(outcomeError);
    outcome->layout()->addWidget(outcomeAdvice);
    f->addWidget(outcome);
    // With convergence criteria, the status card says why the last run stopped;
    // the outcome callout is then only for runs that failed.
    const bool hasCriteria = s->settings()->find("monitors/check-convergence") != nullptr;
    if (hasCriteria) f->addWidget(convergenceStatus(c));

    auto refresh = [s, hasCriteria, runBtn, stopBtn, status, why, outcome, outcomeTitle, outcomeError, outcomeAdvice] {
        QString reason;
        const bool can = s->canRun(&reason);
        runBtn->setEnabled(can);
        stopBtn->setVisible(s->running());
        stopBtn->setEnabled(s->state() == RunState::Running || s->state() == RunState::Starting);
        why->setText(!can && !s->running() ? reason : QString());
        why->setVisible(!why->text().isEmpty());
        const auto& l = s->live();
        if (s->running()) {
            QString t;
            if (s->state() == RunState::Starting) t = QObject::tr("Starting kutta…");
            else if (s->state() == RunState::Stopping) t = QObject::tr("Stopping after the current iteration…");
            else if (l.last) {
                t = QObject::tr("Iteration %1 of %2").arg(l.last->iter).arg(l.requestedIterations);
                if (l.last->secPerIter) t += QObject::tr(" · %1 s per iteration").arg(*l.last->secPerIter, 0, 'f', 3);
                const double el = l.clock.elapsed() / 1000.0;
                t += QObject::tr(" · %1 s elapsed").arg(el, 0, 'f', 0);
                if (l.last->secPerIter && l.requestedIterations > l.last->iter)
                    t += QObject::tr(" · about %1 s left").arg(*l.last->secPerIter * (l.requestedIterations - l.last->iter), 0, 'f', 0);
            } else t = QObject::tr("Setting up on %1…").arg(l.setup ? l.setup->device : QObject::tr("the device"));
            status->setText(t);
        } else if (const RunRecord* r = s->history().runs().empty() ? nullptr : &s->history().runs().back()) {
            QString t = QObject::tr("Last run: %1 · %2 · %3 iterations").arg(r->title(), runStatusText(*r)).arg(r->iterationsDone);
            if (!r->coupling.isEmpty()) t += QObject::tr(" · %1").arg(couplingLabel(r->coupling));
            status->setText(t);
        } else {
            status->setText(QObject::tr("No runs yet for this case."));
        }
        const RunRecord* last = s->running() || s->history().runs().empty() ? nullptr : &s->history().runs().back();
        const QString advice = last && !(hasCriteria && last->status == RunStatus::Stopped) ? runOutcomeAdvice(*last) : QString();
        outcome->setVisible(!advice.isEmpty());
        if (!advice.isEmpty()) {
            const bool failed = last->status == RunStatus::Failed || last->status == RunStatus::Interrupted;
            outcome->setProperty("callout", failed ? "red" : "neutral");
            repolish(outcome);
            outcomeTitle->setText(failed ? QObject::tr("%1 %2 at iteration %3").arg(last->title(), runStatusLabel(last->status).toLower()).arg(last->iterationsDone)
                                         : QObject::tr("%1 stopped").arg(last->title()));
            outcomeError->setText(last->error);
            outcomeError->setVisible(failed && !last->error.isEmpty());
            outcomeAdvice->setText(advice);
        }
    };
    refresh();
    QObject::connect(runBtn, &QPushButton::clicked, page, [s] { s->run(); });
    QObject::connect(stopBtn, &QPushButton::clicked, page, [s] { s->stop(); });
    QObject::connect(s, &Session::runStateChanged, page, refresh);
    QObject::connect(s, &Session::historyChanged, page, refresh);
    QObject::connect(s, &Session::iteration, page, refresh);
    QObject::connect(s, &Session::approvalsChanged, page, refresh);
    QObject::connect(s->settings(), &SettingsModel::overridesChanged, page, refresh);
    f->finish();
    return page;
}

} // namespace pages

} // namespace kapp
