#include "ui/Inspector.h"
#include "core/ChartScale.h"

#include "ui/ChartsPanel.h"
#include "ui/GraphicsArea.h"
#include "ui/Outline.h"
#include "ui/PageKit.h"
#include "ui/SettingWidgets.h"

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QHeaderView>
#include <QJsonArray>
#include <QLocale>
#include <QScrollBar>
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
    connect(session, &Session::caseOpened, this, &Inspector::refresh);
    connect(&Theme::instance(), &Theme::changed, this, &Inspector::refresh);
    connect(graphics, &GraphicsArea::activeChanged, this, [this] {
        if (m_id.startsWith("results/") && m_id != "results/residuals" && m_id != "results/report-plots") refresh();
    });
    connect(session, &Session::approvalsChanged, this, [this] {
        if (m_id == "setup/general" || m_id == "solution/run") refresh();
    });
    connect(session, &Session::historyChanged, this, [this] {
        // Pages listing runs follow the history; editing pages keep their state.
        if (m_id.startsWith("runs/") || m_id == "results/forces") refresh();
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
    auto* b = button(QObject::tr("%1 retained Fluent settings for %2").arg(QLocale().toString(n), what), "list-filter");
    b->setProperty("flat", true);
    b->setToolTip(QObject::tr("Settings saved in the case that Kutta reads but does not expose as controls."));
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

QWidget* noCase(const InspectorContext& c, const QString& title)
{
    QVBoxLayout* v;
    if (c.session->hasCase()) {
        QWidget* page = makePage(title, {}, &v);
        v->addWidget(label(QObject::tr("Nothing to show here yet."), "muted"));
        return page;
    }
    // Welcome: what Kutta does, and the ways to start.
    QWidget* page = makePage(QObject::tr("Open a Fluent case"), {}, &v);
    v->addWidget(label(QObject::tr("Kutta reads a Fluent case (.cas.h5), solves it on this computer's GPU and shows the results next to Fluent's. Open a case, or drop one on the window."), "muted"));
    auto* open = button(QObject::tr("Open Case…"), "folder-open", "primary");
    QObject::connect(open, &QPushButton::clicked, open, [c] { if (c.openDialog) c.openDialog(); });
    v->addWidget(open, 0, Qt::AlignLeft);
    if (!c.samples.isEmpty()) {
        v->addSpacing(6);
        v->addWidget(section(QObject::tr("Samples")));
        auto* grid = new QGridLayout;
        grid->setHorizontalSpacing(10);
        grid->setVerticalSpacing(4);
        grid->setColumnStretch(1, 1);
        int row = 0;
        for (const QString& p : c.samples) {
            const bool synthetic = QFileInfo(p).fileName().startsWith("fsae_car");
            auto* b = button(QFileInfo(p).fileName(), synthetic ? "wind" : "waves");
            b->setProperty("flat", true);
            QObject::connect(b, &QPushButton::clicked, b, [c, p] { c.openPath(p); });
            grid->addWidget(b, row, 0, Qt::AlignLeft | Qt::AlignTop);
            grid->addWidget(label(synthetic ? QObject::tr("Synthetic FSAE car for trying the workbench. Not a CFD result.")
                                            : QObject::tr("Ansys mixing elbow. Replays a recorded Kutta run."), "faint"), row++, 1);
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
    QString statusText = cs.status == SupportStatus::Supported ? QObject::tr("Kutta solves this case as Fluent defines it.")
                         : unsupported ? QObject::tr("Kutta cannot solve this case; see below.")
                                       : QObject::tr("Kutta can solve this case only with the approximations below. Results are marked APPROXIMATE everywhere.");
    statusRow->addStretch(1);
    v->addLayout(statusRow);
    v->addWidget(label(statusText, "muted"));

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
        f->layout()->addWidget(label(i->explanation));
        if (!i->evidence.isEmpty()) f->layout()->addWidget(label(QObject::tr("Measured: %1").arg(i->evidence), "muted"));
        auto* cb = new QCheckBox(QObject::tr("Run with this approximation"), f);
        cb->setChecked(on);
        cb->setAccessibleDescription(QObject::tr("Results will be marked APPROXIMATE in every view, chart and export."));
        cb->setToolTip(cb->accessibleDescription());
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
    auto addRow = [&](const QString& k, const QString& val) {
        grid->addWidget(label(k, "muted"), row, 0, Qt::AlignTop);
        auto* l = label(val);
        l->setTextInteractionFlags(Qt::TextSelectableByMouse);
        grid->addWidget(l, row++, 1);
    };
    grid->setColumnStretch(1, 1);
    addRow(QObject::tr("File"), QDir::toNativeSeparators(c.session->casePath()));
    addRow(QObject::tr("Data"), c.session->dataPath().isEmpty() ? QObject::tr("none (initialized by Kutta)") : QDir::toNativeSeparators(c.session->dataPath()));
    addRow(QObject::tr("Written by"), cs.fluentVersion.isEmpty() ? QObject::tr("unknown Fluent version") : "Fluent " + cs.fluentVersion);
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
    form->addSetting("general/dimension");
    v->addWidget(form);

    if (!notes.empty()) {
        v->addWidget(section(QObject::tr("Notes")));
        auto* toggle = button(notes.size() == 1 ? QObject::tr("1 note on a setting that does not affect this case")
                                                : QObject::tr("%1 notes on settings that do not affect this case").arg(notes.size()), "info");
        toggle->setProperty("flat", true);
        toggle->setCheckable(true);
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
    f->addSection(QObject::tr("Energy"));
    f->addSetting("models/energy");
    for (const auto& i : cs.issues)
        if (i.setting.contains("energy")) f->addNote(i.explanation);
    for (const Setting* s : settingsUnder(c, "models/"))
        if (s->path != "models/viscous" && s->path != "models/energy") f->addSetting(s->path);
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

QWidget* methods(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Methods"), QObject::tr("Solution"), &page);
    f->addSection(QObject::tr("Spatial discretization"));
    f->addSetting("methods/momentum-scheme");
    f->addSetting("methods/k-scheme");
    f->addSetting("methods/omega-scheme");
    f->addSection(QObject::tr("Gradients"));
    f->addSetting("methods/gradient-limiting");
    for (const Setting* s : settingsUnder(c, "methods/"))
        if (!QStringList{"momentum-scheme", "k-scheme", "omega-scheme", "gradient-limiting"}.contains(s->leaf)) f->addSetting(s->path);
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
    f->addSection(QObject::tr("Pseudo time"));
    f->addSetting("controls/pseudo-time-scale-factor");
    for (const Setting* s : settingsUnder(c, "controls/"))
        if (!s->leaf.endsWith("relaxation") && s->leaf != "pseudo-time-scale-factor") f->addSetting(s->path);
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

QWidget* monitors(const InspectorContext& c)
{
    QWidget* page;
    SettingsForm* f = formPage(c, QObject::tr("Monitors"), QObject::tr("Solution"), &page);
    f->addSection(QObject::tr("Residuals"));
    f->addNote(QObject::tr("Scaled residuals for continuity, momentum and the turbulence equations are reported every iteration and plotted live."));
    f->addSetting("monitors/convergence-criteria");
    f->addSetting("run/monitor-interval");
    for (const Setting* s : settingsUnder(c, "monitors/"))
        if (s->path != "monitors/convergence-criteria") f->addSetting(s->path);
    auto* b = button(QObject::tr("Show residual plot"), "chart-line");
    QObject::connect(b, &QPushButton::clicked, b, [c] { c.charts->showTab(ChartsPanel::Tab::Residuals); });
    f->addWidget(b);
    f->finish();
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

    auto refresh = [s, runBtn, stopBtn, status, why] {
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
            status->setText(QObject::tr("Last run: %1 · %2 · %3 iterations").arg(r->title(), runStatusLabel(r->status)).arg(r->iterationsDone));
        } else {
            status->setText(QObject::tr("No runs yet for this case."));
        }
    };
    refresh();
    QObject::connect(runBtn, &QPushButton::clicked, page, [s] { s->run(); });
    QObject::connect(stopBtn, &QPushButton::clicked, page, [s] { s->stop(); });
    QObject::connect(s, &Session::runStateChanged, page, refresh);
    QObject::connect(s, &Session::iteration, page, refresh);
    QObject::connect(s, &Session::approvalsChanged, page, refresh);
    QObject::connect(s->settings(), &SettingsModel::overridesChanged, page, refresh);
    f->finish();
    return page;
}

} // namespace pages

} // namespace kapp
