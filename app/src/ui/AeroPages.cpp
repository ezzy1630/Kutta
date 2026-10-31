// Results ▸ Reports ▸ Aerodynamics and Runs ▸ Design Studies: the numbers an
// FSAE aero engineer designs with (CL, CD, L/D, balance, centre of pressure),
// per run and across a sweep of one setting. Every number comes from
// `kutta aero` on the run's loads.json; nothing here re-derives the physics.

#include "ui/ChartsPanel.h"
#include "ui/GraphicsArea.h"
#include "ui/Inspector.h"
#include "ui/PageKit.h"
#include "ui/SettingWidgets.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QGridLayout>
#include <QHeaderView>
#include <QLineEdit>
#include <QLocale>
#include <QMessageBox>
#include <QPointer>
#include <QSaveFile>
#include <QSettings>
#include <QRegularExpression>
#include <QTableWidget>
#include <QTimer>

#include <cmath>

namespace kapp::pages {

using namespace kit;

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
    // Content takes the whole height and ends in its own stretch (see
    // endPage): sections that grow after results arrive are never squeezed
    // by a height measured before the panel had its width.
    v->addWidget(content, 1);
    return page;
}

QWidget* endPage(QWidget* page, QVBoxLayout* body)
{
    body->addStretch(1);
    return page;
}

void clearLayout(QLayout* l)
{
    while (QLayoutItem* it = l->takeAt(0)) {
        if (it->layout()) clearLayout(it->layout());
        delete it->widget();
        delete it;
    }
}

QString fixed(double v, int decimals) { return QLocale().toString(v, 'f', decimals); }
QString fixed(std::optional<double> v, int decimals) { return v ? fixed(*v, decimals) : QStringLiteral("—"); }

bool finished(const RunRecord& r) { return r.status == RunStatus::Completed || r.status == RunStatus::Stopped; }

// The run a page acts on: the active viewport's, else the latest finished.
QString defaultRun(const InspectorContext& c)
{
    if (auto* p = c.graphics->active())
        if (const QString k = p->datasetKey(); k.startsWith("run:")) return k.mid(4);
    if (const RunRecord* r = c.session->history().latestFinished()) return r->id;
    return {};
}

void addBadges(QHBoxLayout* row, const RunRecord& r)
{
    if (r.approximate()) row->addWidget(badge(r.approximateBadge(), "approximate"));
    if (r.modified()) row->addWidget(badge("MODIFIED", "modified"));
    if (r.synthetic) row->addWidget(badge("SYNTHETIC", "modified"));
    else if (r.fixture) row->addWidget(badge("FIXTURE REPLAY", "modified"));
}

QTableWidget* table(int columns, const QStringList& heads, const QString& name)
{
    auto* t = new QTableWidget(0, columns);
    t->setHorizontalHeaderLabels(heads);
    t->verticalHeader()->hide();
    t->setEditTriggers(QAbstractItemView::NoEditTriggers);
    t->setShowGrid(false);
    t->setAccessibleName(name);
    t->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    QFont tab = t->font();
    tab.setFeature(QFont::Tag("tnum"), 1);
    t->setFont(tab);
    return t;
}

QTableWidgetItem* numberItem(const QString& text)
{
    auto* it = new QTableWidgetItem(text);
    it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
    return it;
}

// Axis choices as Fluent users set them: a signed coordinate axis.
struct AxisChoice { QString label; QVector3D v; };
const std::vector<AxisChoice>& axisChoices()
{
    static const std::vector<AxisChoice> a{{"+X", {1, 0, 0}}, {"−X", {-1, 0, 0}}, {"+Y", {0, 1, 0}},
                                           {"−Y", {0, -1, 0}}, {"+Z", {0, 0, 1}}, {"−Z", {0, 0, -1}}};
    return a;
}

QComboBox* axisCombo(const QVector3D& current, const QString& name)
{
    auto* box = new QComboBox;
    box->setAccessibleName(name);
    for (const AxisChoice& a : axisChoices()) box->addItem(a.label, QVariant::fromValue(a.v));
    int at = -1;
    for (int i = 0; i < box->count(); ++i)
        if ((box->itemData(i).value<QVector3D>() - current.normalized()).length() < 1e-5) at = i;
    if (at < 0) {
        box->addItem(QStringLiteral("(%1, %2, %3)").arg(current.x()).arg(current.y()).arg(current.z()), QVariant::fromValue(current));
        at = box->count() - 1;
    }
    box->setCurrentIndex(at);
    return box;
}

// Three coordinates in metres; commits when editing finishes.
struct PointEdit {
    QWidget* widget;
    std::array<QLineEdit*, 3> edits;
    std::optional<QVector3D> value() const
    {
        QVector3D v;
        for (int i = 0; i < 3; ++i) {
            bool ok = false;
            const double x = QLocale().toDouble(edits[std::size_t(i)]->text().trimmed(), &ok);
            if (!ok || !std::isfinite(x)) return std::nullopt;
            v[i] = float(x);
        }
        return v;
    }
};
PointEdit pointEdit(const std::optional<QVector3D>& v, const QString& name)
{
    PointEdit p{new QWidget, {}};
    auto* h = new QHBoxLayout(p.widget);
    h->setContentsMargins(0, 0, 0, 0);
    h->setSpacing(4);
    for (int i = 0; i < 3; ++i) {
        auto* e = new QLineEdit(v ? QLocale().toString(double((*v)[i]), 'g', 6) : QString());
        e->setPlaceholderText(QString("xyz"[i]));
        e->setAccessibleName(QStringLiteral("%1 %2").arg(name, QString("XYZ"[i])));
        h->addWidget(e, 1);
        p.edits[std::size_t(i)] = e;
    }
    return p;
}

// One aero quantity of a study variant, for the table and the plot.
struct Quantity { QString key, label, unit; int decimals; };
const std::vector<Quantity>& quantities()
{
    static const std::vector<Quantity> q{
        {"cl", "CL", "", 3},           {"cd", "CD", "", 3},          {"ld", "L/D", "", 2},
        {"downforce", "Downforce", "N", 1}, {"drag", "Drag", "N", 1}, {"front", "Front", "%", 1},
        {"cop", "CoP from front axle", "m", 3}, {"cla", "ClA", "m²", 3}, {"cda", "CdA", "m²", 3}};
    return q;
}
std::optional<double> quantity(const AeroResult& a, const QString& key)
{
    const AeroRow* t = a.total();
    if (!t) return std::nullopt;
    const auto& k = t->coefficients;
    if (key == "cl") return k ? std::optional<double>(k->cl) : std::nullopt;
    if (key == "cd") return k ? std::optional<double>(k->cd) : std::nullopt;
    if (key == "cla") return t->cla;
    if (key == "cda") return t->cda;
    if (key == "ld") return a.liftToDrag;
    if (key == "downforce") return t->downforce;
    if (key == "drag") return t->total.drag;
    if (key == "front") return a.balance.frontPercent;
    if (key == "cop") return a.balance.copFromFrontAxle;
    return std::nullopt;
}

} // namespace

// ---- Aerodynamics ---------------------------------------------------------------------

QWidget* aero(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Aerodynamics"), QObject::tr("Results ▸ Reports"), &v);
    Session* s = c.session;
    if (!s->hasCase()) {
        v->addWidget(label(QObject::tr("Open a case to see its aerodynamic coefficients."), "muted"));
        return endPage(page, v);
    }
    std::vector<const RunRecord*> runs;
    for (auto it = s->history().runs().rbegin(); it != s->history().runs().rend(); ++it)
        if (finished(*it)) runs.push_back(&*it);
    // A 2D case is a section: no vehicle, no balance, forces for the reference depth.
    const bool planar = s->caseSummary().planar();
    if (runs.empty()) {
        v->addWidget(label(planar ? QObject::tr("Downforce, drag and their coefficients appear here when a run finishes.")
                                  : QObject::tr("Downforce, drag, balance and centre of pressure appear here when a run finishes. Set the wheelbase below now if you like; it is kept with the case."), "muted"));
    }

    // ---- result
    auto* picker = new QComboBox;
    picker->setAccessibleName(QObject::tr("Run"));
    for (const RunRecord* r : runs) picker->addItem(r->title(), r->id);
    picker->setCurrentIndex(std::max(0, picker->findData(defaultRun(c))));
    auto* top = new QHBoxLayout;
    top->addWidget(picker, 1);
    auto* badges = new QHBoxLayout;
    top->addLayout(badges);
    if (!runs.empty()) v->addLayout(top);
    auto* status = label({}, "muted");
    v->addWidget(status);
    auto* headline = new QGridLayout;
    headline->setHorizontalSpacing(16);
    headline->setVerticalSpacing(0);
    v->addLayout(headline);
    auto* summary = new QGridLayout;
    summary->setHorizontalSpacing(12);
    summary->setVerticalSpacing(3);
    v->addLayout(summary);
    auto* balanceNote = callout("amber");
    balanceNote->layout()->addWidget(label(QObject::tr("Enter the front axle position and wheelbase below to see the aero balance (% front) and the centre of pressure.")));
    balanceNote->hide();
    v->addWidget(balanceNote);

    auto* breakdownHead = new QHBoxLayout;
    breakdownHead->addWidget(section(QObject::tr("Breakdown")), 1);
    auto* component = new QComboBox;
    component->setAccessibleName(QObject::tr("Force component"));
    component->addItems({QObject::tr("Total"), QObject::tr("Pressure"), QObject::tr("Viscous")});
    breakdownHead->addWidget(component);
    v->addLayout(breakdownHead);
    auto* breakdown = table(6, {QObject::tr("Group"), QObject::tr("Down [N]"), QObject::tr("Drag [N]"), "CL", "CD", QObject::tr("Share")},
                            QObject::tr("Aero breakdown"));
    breakdown->horizontalHeaderItem(1)->setToolTip(QObject::tr("Downforce, N"));
    breakdown->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    breakdown->verticalHeader()->setDefaultSectionSize(26);
    v->addWidget(breakdown);
    auto* refs = label({}, "faint");
    v->addWidget(refs);
    auto* actions = new QHBoxLayout;
    auto* save = button(QObject::tr("Export CSV…"), "file-down");
    auto* toRefs = button(QObject::tr("Reference Values"), "chevron-right");
    toRefs->setProperty("flat", true);
    actions->addWidget(save);
    actions->addStretch(1);
    actions->addWidget(toRefs);
    v->addLayout(actions);
    QObject::connect(toRefs, &QPushButton::clicked, toRefs, [c] { c.navigate("setup/reference"); });

    // ---- vehicle and axes (kept with the case, applied to every run)
    const AeroSpec spec = s->aeroSpec();
    v->addWidget(section(planar ? QObject::tr("Axes") : QObject::tr("Vehicle")));
    auto* vg = new QGridLayout;
    vg->setHorizontalSpacing(10);
    vg->setVerticalSpacing(4);
    PointEdit axle = pointEdit(spec.frontAxle, QObject::tr("Front axle"));
    auto* wheelbase = new QLineEdit(spec.wheelbase ? QLocale().toString(*spec.wheelbase, 'g', 6) : QString());
    wheelbase->setAccessibleName(QObject::tr("Wheelbase"));
    wheelbase->setPlaceholderText(QObject::tr("e.g. 1.53"));
    int row = 0;
    if (!planar) {
        vg->addWidget(label(QObject::tr("Front axle (tyre contact)")), row, 0);
        vg->addWidget(axle.widget, row, 1);
        vg->addWidget(label("m", "faint"), row++, 2);
        vg->addWidget(label(QObject::tr("Wheelbase")), row, 0);
        vg->addWidget(wheelbase, row, 1);
        vg->addWidget(label("m", "faint"), row++, 2);
    } else {
        axle.widget->setParent(page);
        axle.widget->hide();
        wheelbase->setParent(page);
        wheelbase->hide();
    }
    auto* flow = axisCombo(spec.drag, QObject::tr("Free-stream direction"));
    auto* up = axisCombo(spec.lift, QObject::tr("Up direction"));
    vg->addWidget(label(QObject::tr("Free stream (drag)")), row, 0);
    vg->addWidget(flow, row++, 1);
    vg->addWidget(label(QObject::tr("Up (lift)")), row, 0);
    vg->addWidget(up, row++, 1);
    PointEdit centre = pointEdit(spec.momentCentre, QObject::tr("Moment centre"));
    auto* centreToAxle = button(QObject::tr("Front axle"), {});
    centreToAxle->setProperty("flat", true);
    centreToAxle->setToolTip(QObject::tr("Take moments about the front tyre contact point"));
    centreToAxle->setText(QObject::tr("Use the front axle"));
    vg->addWidget(label(QObject::tr("Moment centre")), row, 0);
    vg->addWidget(centre.widget, row, 1);
    vg->addWidget(label("m", "faint"), row++, 2);
    if (!planar) vg->addWidget(centreToAxle, row++, 1, Qt::AlignLeft);
    else {
        centreToAxle->setParent(page);
        centreToAxle->hide();
    }
    vg->setColumnStretch(1, 1);
    v->addLayout(vg);
    auto* specError = label({}, "error");
    specError->hide();
    v->addWidget(specError);
    v->addWidget(label(planar ? QObject::tr("A 2D case is a section: forces are for the case's reference depth, so CL and CD are section coefficients. Kept with the case; changing it recomputes every run without solving again.")
                              : QObject::tr("Balance splits the downforce between the axles using the pitching moment about the front contact point, so drag acting above the ground counts. CoP is where the total downforce alone would give that split. Kept with the case; changing it recomputes every run without solving again."), "faint"));

    // ---- groups
    v->addWidget(section(QObject::tr("Groups")));
    const QStringList car = s->carZones();
    QStringList groupNames;
    for (const auto& g : suggestAeroGroups({"front-wing", "rear-wing", "undertray", "wheel", "body"})) groupNames << g.name;
    for (const auto& g : spec.groups)
        if (!groupNames.contains(g.name)) groupNames << g.name;
    auto* groups = table(2, {QObject::tr("Zone"), QObject::tr("Group")}, QObject::tr("Aero groups"));
    groups->setRowCount(int(car.size()));
    std::vector<QComboBox*> groupBoxes;
    for (int i = 0; i < car.size(); ++i) {
        groups->setItem(i, 0, new QTableWidgetItem(car[i]));
        auto* box = new QComboBox;
        box->setAccessibleName(QObject::tr("Group of %1").arg(car[i]));
        box->addItem(QObject::tr("(none)"), QString());
        for (const QString& g : groupNames) box->addItem(g, g);
        QString current;
        for (const auto& g : spec.groups)
            if (g.zones.contains(car[i])) current = g.name;
        box->setCurrentIndex(std::max(0, box->findData(current)));
        groups->setCellWidget(i, 1, box);
        groupBoxes.push_back(box);
    }
    groups->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    groups->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
    groups->verticalHeader()->setDefaultSectionSize(30);
    groups->setFixedHeight(34 + int(car.size()) * 30);
    v->addWidget(groups);
    auto* suggest = button(QObject::tr("Suggest from zone names"), "wand-sparkles");
    suggest->setProperty("flat", true);
    v->addWidget(suggest, 0, Qt::AlignLeft);
    v->addWidget(label(planar ? QObject::tr("Zones are the walls of the case's force reports, so the total matches them; the ground is not included.")
                              : QObject::tr("Car zones are the walls of the case's force reports, so the total matches them; the ground is not part of the car."), "faint"));

    // Build the spec from the inputs; invalid input says why and keeps the saved spec.
    auto commit = [s, axle, wheelbase, flow, up, centre, groupBoxes, car, groupNames, specError] {
        AeroSpec n = s->aeroSpec();
        QString error;
        n.drag = flow->currentData().value<QVector3D>();
        n.lift = up->currentData().value<QVector3D>();
        if (std::abs(QVector3D::dotProduct(n.drag.normalized(), n.lift.normalized())) > 1e-6)
            error = QObject::tr("The free stream and up directions must be perpendicular.");
        const bool axleBlank = std::all_of(axle.edits.begin(), axle.edits.end(), [](QLineEdit* e) { return e->text().trimmed().isEmpty(); });
        n.frontAxle = axleBlank ? std::nullopt : axle.value();
        if (!axleBlank && !n.frontAxle) error = QObject::tr("The front axle needs three numbers (m).");
        bool ok = false;
        const double wb = QLocale().toDouble(wheelbase->text().trimmed(), &ok);
        n.wheelbase = ok ? std::optional<double>(wb) : std::nullopt;
        if (!wheelbase->text().trimmed().isEmpty() && (!ok || wb <= 0)) error = QObject::tr("The wheelbase must be a positive length (m).");
        if (const auto mc = centre.value()) n.momentCentre = *mc;
        else error = QObject::tr("The moment centre needs three numbers (m).");
        n.groups.clear();
        for (const QString& g : groupNames) {
            AeroSpec::Group grp{g, {}};
            for (std::size_t i = 0; i < groupBoxes.size(); ++i)
                if (groupBoxes[i]->currentData().toString() == g) grp.zones << car[int(i)];
            if (!grp.zones.isEmpty()) n.groups.push_back(grp);
        }
        specError->setText(error);
        specError->setVisible(!error.isEmpty());
        if (error.isEmpty()) s->setAeroSpec(n);
    };
    for (QLineEdit* e : axle.edits) QObject::connect(e, &QLineEdit::editingFinished, page, commit);
    for (QLineEdit* e : centre.edits) QObject::connect(e, &QLineEdit::editingFinished, page, commit);
    QObject::connect(wheelbase, &QLineEdit::editingFinished, page, commit);
    QObject::connect(flow, &QComboBox::activated, page, commit);
    QObject::connect(up, &QComboBox::activated, page, commit);
    for (QComboBox* b : groupBoxes) QObject::connect(b, &QComboBox::activated, page, commit);
    QObject::connect(centreToAxle, &QPushButton::clicked, page, [axle, centre, commit] {
        for (int i = 0; i < 3; ++i) centre.edits[std::size_t(i)]->setText(axle.edits[std::size_t(i)]->text());
        commit();
    });
    QObject::connect(suggest, &QPushButton::clicked, page, [groupBoxes, car, commit] {
        const auto suggested = suggestAeroGroups(car);
        for (std::size_t i = 0; i < groupBoxes.size(); ++i) {
            QString g;
            for (const auto& x : suggested)
                if (x.zones.contains(car[int(i)])) g = x.name;
            groupBoxes[i]->setCurrentIndex(std::max(0, groupBoxes[i]->findData(g)));
        }
        commit();
    });

    // ---- fill from kutta aero
    QPointer<QWidget> guard(page);
    auto fill = [=](const AeroResult& a) {
        clearLayout(headline);
        clearLayout(summary);
        const AeroRow* t = a.total();
        const auto& k = t->coefficients;
        const auto& b = a.balance;
        // The four numbers a design review asks for first.
        const std::vector<std::pair<QString, QString>> tiles{
            {"CL", k ? fixed(k->cl, 3) : QStringLiteral("—")},
            {"CD", k ? fixed(k->cd, 3) : QStringLiteral("—")},
            {QObject::tr("L/D"), fixed(a.liftToDrag, 2)},
            planar ? std::pair{QObject::tr("Downforce"), fixed(t->downforce, 1) + " N"}
                   : std::pair{QObject::tr("Balance"), b.frontPercent ? QObject::tr("%1 % F").arg(fixed(b.frontPercent, 1)) : QStringLiteral("—")}};
        for (int i = 0; i < int(tiles.size()); ++i) {
            headline->addWidget(label(tiles[std::size_t(i)].first, "faint"), 0, i);
            auto* big = label(tiles[std::size_t(i)].second);
            QFont bf = big->font();
            bf.setPointSizeF(bf.pointSizeF() * 1.45);
            bf.setWeight(QFont::DemiBold);
            bf.setFeature(QFont::Tag("tnum"), 1);
            big->setFont(bf);
            big->setTextInteractionFlags(Qt::TextSelectableByMouse);
            headline->addWidget(big, 1, i);
        }
        int row = 0;
        auto line = [&](const QString& name, const QString& value) {
            summary->addWidget(label(name, "muted"), row, 0, Qt::AlignTop);
            auto* l = label(value);
            QFont f = l->font();
            f.setFeature(QFont::Tag("tnum"), 1);
            l->setFont(f);
            l->setTextInteractionFlags(Qt::TextSelectableByMouse);
            summary->addWidget(l, row++, 1);
        };
        const QString sep = QStringLiteral("  ·  ");
        line(QObject::tr("Downforce"), fixed(t->downforce, 1) + " N" + (t->cla ? sep + "ClA " + fixed(*t->cla, 3) + " m²" : QString()));
        line(QObject::tr("Drag"), fixed(t->total.drag, 1) + " N" + (t->cda ? sep + "CdA " + fixed(*t->cda, 3) + " m²" : QString()) +
                                      sep + QObject::tr("%1 % viscous").arg(fixed(t->total.drag != 0 ? 100 * t->viscous.drag / t->total.drag : 0, 0)));
        if (!planar) line(QObject::tr("Side force"), fixed(t->total.side, 1) + " N" + (k ? sep + "CS " + fixed(k->cs, 3) : QString()));
        line(QObject::tr("Pitching moment"), fixed(t->pitch, 1) + " N m" + (k && k->cmPitch ? sep + "CM " + fixed(*k->cmPitch, 3) : QString()) + sep +
                                                 QObject::tr("about (%1, %2, %3) m").arg(a.momentCentre.x()).arg(a.momentCentre.y()).arg(a.momentCentre.z()));
        if (b.status == "defined") {
            line(QObject::tr("Balance"), QObject::tr("%1 % front, %2 % rear").arg(fixed(b.frontPercent, 1), fixed(100 - b.frontPercent.value_or(0), 1)) + sep +
                                             QObject::tr("front %1 N, rear %2 N").arg(fixed(b.frontLoad, 1), fixed(b.rearLoad, 1)));
            line(QObject::tr("Centre of pressure"), QObject::tr("%1 m behind the front axle").arg(fixed(b.copFromFrontAxle, 3)) +
                                                        (b.wheelbase ? sep + QObject::tr("%1 % of wheelbase").arg(fixed(100 * b.copFromFrontAxle.value_or(0) / *b.wheelbase, 1)) : QString()));
        } else if (b.status == "no_downforce") {
            line(QObject::tr("Balance"), QObject::tr("no net downforce") + sep + QObject::tr("front %1 N, rear %2 N").arg(fixed(b.frontLoad, 1), fixed(b.rearLoad, 1)));
        }
        summary->setColumnStretch(1, 1);
        balanceNote->setVisible(!planar && b.status == "needs_vehicle");

        const int comp = component->currentIndex();
        breakdown->setRowCount(0);
        auto pick = [comp](const AeroRow& r) { return comp == 1 ? r.pressure : comp == 2 ? r.viscous : r.total; };
        const double qa = a.dynamicPressure && a.area ? *a.dynamicPressure * *a.area : 0;
        std::vector<const AeroRow*> order;
        for (const auto& r : a.rows)
            if (r.kind == "group") order.push_back(&r);
        for (const auto& r : a.rows)
            if (r.kind == "zone" && a.ungrouped.contains(r.name)) order.push_back(&r);
        if (order.empty() || std::none_of(a.rows.begin(), a.rows.end(), [](const AeroRow& r) { return r.kind == "group"; }))
            for (const auto& r : a.rows)
                if (r.kind == "zone" && std::find(order.begin(), order.end(), &r) == order.end()) order.push_back(&r);
        order.push_back(t);
        for (const AeroRow* r : order) {
            const int i = breakdown->rowCount();
            breakdown->insertRow(i);
            const AeroComponents f = pick(*r);
            auto* name = new QTableWidgetItem(r->kind == "zone" && !a.ungrouped.isEmpty() && r->kind != "total" &&
                                                      std::any_of(a.rows.begin(), a.rows.end(), [](const AeroRow& x) { return x.kind == "group"; })
                                                  ? r->name + QObject::tr("  (ungrouped)")
                                                  : r->name);
            name->setToolTip(r->zones.join(", "));
            if (r->kind == "total") {
                QFont fb = name->font();
                fb.setBold(true);
                name->setFont(fb);
            }
            breakdown->setItem(i, 0, name);
            breakdown->setItem(i, 1, numberItem(fixed(-f.lift, 1)));
            breakdown->setItem(i, 2, numberItem(fixed(f.drag, 1)));
            breakdown->setItem(i, 3, numberItem(qa > 0 ? fixed(f.lift / qa, 3) : QStringLiteral("—")));
            breakdown->setItem(i, 4, numberItem(qa > 0 ? fixed(f.drag / qa, 3) : QStringLiteral("—")));
            breakdown->setItem(i, 5, numberItem(r->downforceShare ? fixed(100 * *r->downforceShare, 0) + " %" : QStringLiteral("—")));
        }
        // Names take the room the numbers leave; numbers never elide.
        for (int col = 1; col < breakdown->columnCount(); ++col) breakdown->horizontalHeader()->setSectionResizeMode(col, QHeaderView::ResizeToContents);
        breakdown->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
        breakdown->setFixedHeight(std::max(breakdown->horizontalHeader()->sizeHint().height(), breakdown->fontMetrics().height() + 14) + breakdown->rowCount() * 26 + 4);
        refs->setText(a.dynamicPressure
                          ? QObject::tr("Coefficients use the run's Reference Values: A %1 m², ρ %2 kg/m³, V %3 m/s, L %4 m (q = %5 Pa). Share is of the total downforce. Positive pitch raises the nose.")
                                .arg(fixed(a.area, 3), QLocale().toString(a.density.value_or(0), 'g', 5), QLocale().toString(a.velocity.value_or(0), 'g', 5),
                                     a.length ? QLocale().toString(*a.length, 'g', 5) : QStringLiteral("—"), fixed(a.dynamicPressure, 1))
                          : QObject::tr("The run has no valid reference density, velocity or area, so coefficients are not available; forces are."));
        refs->setToolTip(a.provenance);
    };
    auto compute = [=] {
        const QString id = picker->currentData().toString();
        const RunRecord* r = s->history().find(id);
        if (!r) return;
        clearLayout(badges);
        addBadges(badges, *r);
        status->setText(QObject::tr("Computing from %1…").arg(r->title()));
        s->aero(id, [=](Result<AeroResult> res) {
            if (!guard) return;
            if (!res) {
                status->setText(res.error());
                status->setProperty("role", "error");
                repolish(status);
                return;
            }
            status->setText(QObject::tr("%1 · %2 iterations · %3").arg(r->title()).arg(res->iterations).arg(runStatusText(*r)));
            fill(res.value());
        });
    };
    if (!runs.empty()) {
        compute();
        QObject::connect(picker, &QComboBox::activated, page, compute);
        QObject::connect(component, &QComboBox::activated, page, compute);
        QObject::connect(s, &Session::aeroSpecChanged, page, compute);
    } else {
        component->hide();
        breakdown->hide();
        save->hide();
        status->hide();
    }

    QObject::connect(save, &QPushButton::clicked, page, [s, picker, save] {
        const auto a = s->aeroIfComputed(picker->currentData().toString());
        if (!a) return;
        const QString f = QFileDialog::getSaveFileName(save->window(), QObject::tr("Export aero breakdown"), "aero.csv", "CSV (*.csv)");
        if (f.isEmpty()) return;
        QString csv;
        for (const QString& part : a->provenance.split(" | ", Qt::SkipEmptyParts)) csv += "# " + part + "\n";
        csv += "name,kind,downforce_N,drag_N,side_N,pressure_drag_N,viscous_drag_N,pitch_Nm,cl,cd,cs,cla_m2,cda_m2,cm_pitch,downforce_share\n";
        const QLocale c = QLocale::c();
        auto n = [&](std::optional<double> x) { return x ? c.toString(*x, 'g', 10) : QString(); };
        for (const auto& r : a->rows) {
            const auto& k = r.coefficients;
            csv += QStringList{"\"" + r.name + "\"", r.kind, n(r.downforce), n(r.total.drag), n(r.total.side), n(r.pressure.drag), n(r.viscous.drag), n(r.pitch),
                               n(k ? std::optional(k->cl) : std::nullopt), n(k ? std::optional(k->cd) : std::nullopt), n(k ? std::optional(k->cs) : std::nullopt),
                               n(r.cla), n(r.cda), n(k ? k->cmPitch : std::nullopt), n(r.downforceShare)}
                       .join(',') + "\n";
        }
        QSaveFile out(f);
        if (out.open(QIODevice::WriteOnly)) {
            out.write(csv.toUtf8());
            out.commit();
        }
    });
    return endPage(page, v);
}

// ---- Design studies -------------------------------------------------------------------

QWidget* studies(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Design Studies"), QObject::tr("Runs"), &v);
    Session* s = c.session;
    if (!s->hasCase()) {
        v->addWidget(label(QObject::tr("Open a case to run a study of it."), "muted"));
        return endPage(page, v);
    }
    SettingsModel* m = s->settings();

    // ---- new study
    v->addWidget(section(QObject::tr("New study")));
    auto* f = new QGridLayout;
    f->setHorizontalSpacing(10);
    f->setVerticalSpacing(4);
    auto* param = new QComboBox;
    param->setAccessibleName(QObject::tr("Setting to vary"));
    int allowed = 0, editable = 0;
    for (const Setting& st : m->settings()) {
        editable += st.editable && !st.raw;
        if (!studyParameterAllowed(st)) continue;
        ++allowed;
        const DisplayUnit du = displayUnit(st);
        QString text = settingLabel(st);
        if (!st.owner.isEmpty()) text += QStringLiteral(" · %1").arg(st.owner);
        else if (!st.section.isEmpty()) text += QStringLiteral(" · %1").arg(st.group);
        if (!du.unit.isEmpty() && du.unit != "1") text += QStringLiteral(" [%1]").arg(du.unit);
        param->addItem(text, st.path);
    }
    // An inlet speed sweep is the common first study.
    for (int i = 0; i < param->count(); ++i)
        if (param->itemData(i).toString().endsWith("/velocity-magnitude")) {
            param->setCurrentIndex(i);
            break;
        }
    auto* values = new QLineEdit;
    values->setAccessibleName(QObject::tr("Values"));
    values->setPlaceholderText(QObject::tr("16, 20, 24   or   12 to 28 step 4"));
    auto* name = new QLineEdit;
    name->setAccessibleName(QObject::tr("Study name"));
    name->setPlaceholderText(QObject::tr("optional, e.g. Speed sweep"));
    auto* follow = new QCheckBox;
    auto* current = label({}, "faint");
    f->addWidget(label(QObject::tr("Vary")), 0, 0);
    f->addWidget(param, 0, 1);
    f->addWidget(current, 1, 1);
    f->addWidget(label(QObject::tr("Values")), 2, 0);
    f->addWidget(values, 2, 1);
    f->addWidget(follow, 3, 1);
    f->addWidget(label(QObject::tr("Name")), 4, 0);
    f->addWidget(name, 4, 1);
    f->setColumnStretch(1, 1);
    v->addLayout(f);
    auto* createError = label({}, "error");
    createError->hide();
    v->addWidget(createError);
    auto* create = button(QObject::tr("Create study"), "plus");
    v->addWidget(create, 0, Qt::AlignLeft);
    v->addWidget(label(QObject::tr("Only the %1 settings Kutta can change in this case and that change what is solved can vary (%2 editable in all). Every variant is an ordinary run with your other settings, iteration count and convergence criteria; variants that change the case are marked MODIFIED.")
                           .arg(allowed).arg(editable), "faint"));
    auto updateParam = [m, param, follow, current] {
        const Setting* st = m->find(param->currentData().toString());
        if (!st) return;
        current->setText(QObject::tr("Case value: %1").arg(displayValue(*st, st->value)));
        const QStringList followers = referenceFollowers(*m, st->path);
        follow->setVisible(!followers.isEmpty());
        follow->setChecked(!followers.isEmpty());
        follow->setText(QObject::tr("Reference velocity follows it"));
        follow->setToolTip(QObject::tr("The case's reference velocity equals this inlet's, so coefficients stay referenced to the speed actually run."));
    };
    updateParam();
    QObject::connect(param, &QComboBox::currentIndexChanged, page, updateParam);

    // ---- existing studies
    const auto& all = s->studies();
    if (all.empty()) {
        v->addWidget(label(QObject::tr("No studies yet. A study runs one variant after another on the GPU and tabulates CL, CD, L/D and balance for each."), "muted"));
    }
    QComboBox* pick = nullptr;
    if (!all.empty()) {
        v->addWidget(section(QObject::tr("Studies")));
        pick = new QComboBox;
        pick->setAccessibleName(QObject::tr("Study"));
        for (const Study& st : all) pick->addItem(QStringLiteral("%1 · %2, %3 values").arg(st.title(), st.parameterLabel).arg(st.variants.size()), st.id);
        const QString remembered = QSettings().value("studies/selected").toString();
        pick->setCurrentIndex(pick->findData(remembered) >= 0 ? pick->findData(remembered) : pick->count() - 1);
        v->addWidget(pick);
    }
    // A nested layout, not a wrapper widget: wrapped labels then get their
    // height for the panel's width.
    auto* dv = new QVBoxLayout;
    dv->setSpacing(6);
    v->addLayout(dv);

    QPointer<QWidget> guard(page);
    auto plotKey = std::make_shared<QString>(QSettings().value("studies/plot", "cl").toString());
    // The table rebuilds itself as results arrive; the page's connections own it.
    auto build = std::make_shared<std::function<void()>>();
    std::weak_ptr<std::function<void()>> weak = build;
    auto rebuild = [weak] { if (auto b = weak.lock()) (*b)(); };
    *build = [=] {
        clearLayout(dv);
        if (!pick) return;
        const Study* st = s->study(pick->currentData().toString());
        if (!st) return;
        QSettings().setValue("studies/selected", st->id);
        const bool active = s->activeStudy() == st->id;
        int done = 0;
        for (const auto& var : st->variants) done += var.status == VariantStatus::Done;
        auto* head = new QHBoxLayout;
        head->addWidget(label(active ? QObject::tr("Running: %1 of %2 done").arg(done).arg(st->variants.size())
                                     : QObject::tr("%1 of %2 variants done").arg(done).arg(st->variants.size()), "muted"), 1);
        auto* runBtn = button(done == 0 ? QObject::tr("Run study") : QObject::tr("Run remaining"), "play", "primary");
        auto* stopBtn = button(QObject::tr("Stop study"), "circle-stop", "danger");
        runBtn->setVisible(!active && st->variants.size() > std::size_t(done));
        runBtn->setEnabled(!s->running() && st->variants.size() > std::size_t(done));
        stopBtn->setVisible(active);
        head->addWidget(runBtn);
        head->addWidget(stopBtn);
        dv->addLayout(head);
        auto* runError = label({}, "error");
        runError->hide();
        dv->addWidget(runError);
        QObject::connect(runBtn, &QPushButton::clicked, runBtn, [s, st, runError] {
            QString why;
            if (!s->startStudy(st->id, &why)) {
                runError->setText(why);
                runError->show();
            }
        });
        QObject::connect(stopBtn, &QPushButton::clicked, stopBtn, [s] { s->stopStudy(); });

        const QString unit = prettyUnit(st->unit);
        const QString axisLabel = unit.isEmpty() || unit == "1" ? st->parameterLabel : QStringLiteral("%1 [%2]").arg(st->parameterLabel, unit);
        QStringList heads{unit.isEmpty() || unit == "1" ? QObject::tr("Value") : QStringLiteral("[%1]").arg(unit),
                          QObject::tr("Run"), "CL", "CD", "L/D", QObject::tr("Front %"), QObject::tr("CoP [m]")};
        // CoP is measured behind the front axle.
        auto* caption = label(QObject::tr("%1 of %2").arg(st->parameterLabel, decodePathSegment(st->parameter.section('/', 1, 1))), "faint");
        caption->setWordWrap(false);
        dv->addWidget(caption);
        auto* t = table(int(heads.size()), heads, QObject::tr("Study results"));
        t->horizontalHeaderItem(0)->setToolTip(axisLabel);
        t->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        t->verticalHeader()->setDefaultSectionSize(26);
        t->setSelectionBehavior(QAbstractItemView::SelectRows);
        t->setSelectionMode(QAbstractItemView::ExtendedSelection);
        t->setRowCount(int(st->variants.size()));
        const auto& tk = tokens();
        std::vector<PlotSeries> series;
        const Quantity* q = &quantities()[0];
        for (const Quantity& x : quantities())
            if (x.key == *plotKey) q = &x;
        PlotSeries line{q->label, {}, {}, tk.series.empty() ? tk.accent : tk.series[0], false, true, true, {}};
        QStringList marks;
        QString provenance;
        bool missing = false;
        for (int i = 0; i < int(st->variants.size()); ++i) {
            const StudyVariant& var = st->variants[std::size_t(i)];
            t->setItem(i, 0, numberItem(QLocale().toString(var.value, 'g', 8)));
            const RunRecord* r = s->history().find(var.runId);
            QString runText = variantStatusLabel(var.status);
            if (r) {
                runText = QStringLiteral("Run %1").arg(r->number);
                if (var.status != VariantStatus::Done) runText += " · " + variantStatusLabel(var.status).toLower();
                else if (r->converged()) runText += QObject::tr(" · converged");
                QStringList m;
                if (r->approximate()) m << "APPROXIMATE";
                if (r->modified()) m << "MODIFIED";
                if (r->synthetic) m << "SYNTHETIC";
                for (const QString& x : m)
                    if (!marks.contains(x)) marks << x;
            }
            auto* runItem = new QTableWidgetItem(runText);
            runItem->setData(Qt::UserRole, var.runId);
            if (var.status == VariantStatus::Failed) runItem->setForeground(tk.red);
            // MODIFIED variants differ from the case; the one at the case value does not.
            if (r && r->modified()) {
                runItem->setText(runText + " *");
                runItem->setForeground(tk.dark ? tk.amber : tk.amberText);
            }
            if (r) runItem->setToolTip(QStringLiteral("%1\n%2").arg(r->title(), r->modified() ? QObject::tr("MODIFIED: %1").arg(r->overrides.keys().join(", ")) : QObject::tr("as in the case")));
            t->setItem(i, 1, runItem);
            const auto a = var.status == VariantStatus::Done ? s->aeroIfComputed(var.runId) : std::nullopt;
            if (var.status == VariantStatus::Done && !a) missing = true;
            auto cell = [&](int col, const QString& key, int decimals) {
                const auto x = a ? quantity(*a, key) : std::nullopt;
                t->setItem(i, col, numberItem(x ? fixed(*x, decimals) : QStringLiteral("—")));
            };
            cell(2, "cl", 3);
            cell(3, "cd", 3);
            cell(4, "ld", 2);
            cell(5, "front", 1);
            cell(6, "cop", 3);
            if (a) {
                if (provenance.isEmpty()) provenance = a->provenance;
                if (const auto y = quantity(*a, q->key)) {
                    line.x.push_back(var.value);
                    line.y.push_back(*y);
                }
            }
        }
        for (int col = 0; col < t->columnCount(); ++col) t->horizontalHeader()->setSectionResizeMode(col, QHeaderView::ResizeToContents);
        t->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
        // Before the table is shown its header has no size hint yet.
        t->setFixedHeight(std::max(t->horizontalHeader()->sizeHint().height(), t->fontMetrics().height() + 14) + t->rowCount() * 26 + 4);
        dv->addWidget(t);

        auto* plotRow = new QHBoxLayout;
        plotRow->addWidget(label(QObject::tr("Plot"), "muted"));
        auto* what = new QComboBox;
        what->setAccessibleName(QObject::tr("Quantity to plot"));
        for (const Quantity& x : quantities()) what->addItem(x.unit.isEmpty() ? x.label : QStringLiteral("%1 [%2]").arg(x.label, x.unit), x.key);
        what->setCurrentIndex(std::max(0, what->findData(*plotKey)));
        plotRow->addWidget(what, 1);
        auto* compare = button(QObject::tr("Compare"), "git-compare-arrows");
        compare->setToolTip(QObject::tr("Select two rows to compare their settings, forces and fields"));
        auto* exportCsv = button(QObject::tr("Export CSV…"), "file-down");
        plotRow->addWidget(compare);
        plotRow->addWidget(exportCsv);
        dv->addLayout(plotRow);
        auto* foot = new QHBoxLayout;
        auto* del = button(QObject::tr("Delete study…"), "trash-2");
        del->setProperty("flat", true);
        del->setEnabled(!active);
        foot->addWidget(label(QObject::tr("* MODIFIED: the run changes the case. Balance and CoP use the vehicle on Results ▸ Reports ▸ Aerodynamics. Runs stay in Run History; deleting the study keeps them."), "faint"), 1);
        foot->addWidget(del, 0, Qt::AlignTop);
        dv->addLayout(foot);

        if (!line.x.empty()) {
            const QString yLabel = q->unit.isEmpty() ? q->label : QStringLiteral("%1 [%2]").arg(q->label, q->unit);
            c.charts->setXY(QStringLiteral("%1 vs %2 · %3").arg(q->label, st->parameterLabel, st->title()), axisLabel, yLabel, {line}, marks, provenance);
        }
        // Results computed after the table was built fill it in.
        if (missing)
            for (const auto& var : st->variants)
                if (var.status == VariantStatus::Done && !s->aeroIfComputed(var.runId))
                    s->aero(var.runId, [guard, rebuild](Result<AeroResult>) {
                        if (guard) QTimer::singleShot(0, guard, rebuild);
                    });

        QObject::connect(what, &QComboBox::activated, what, [what, plotKey, rebuild, guard] {
            *plotKey = what->currentData().toString();
            QSettings().setValue("studies/plot", *plotKey);
            QTimer::singleShot(0, guard, rebuild);
        });
        QObject::connect(compare, &QPushButton::clicked, compare, [c, t] {
            QStringList ids;
            for (const auto& idx : t->selectionModel()->selectedRows())
                if (const QString id = t->item(idx.row(), 1)->data(Qt::UserRole).toString(); !id.isEmpty()) ids << id;
            if (ids.size() != 2) return;
            QSettings().setValue("compare/a", ids[1]);
            QSettings().setValue("compare/b", ids[0]);
            c.navigate("runs/compare");
        });
        QObject::connect(exportCsv, &QPushButton::clicked, exportCsv, [s, st, exportCsv] {
            const QString file = QFileDialog::getSaveFileName(exportCsv->window(), QObject::tr("Export study"), st->title() + ".csv", "CSV (*.csv)");
            if (file.isEmpty()) return;
            QString csv;
            const QLocale lc = QLocale::c();
            QStringList head{st->parameter, "run", "status"};
            for (const Quantity& x : quantities()) head << x.key;
            head << "marks";
            csv += head.join(',') + "\n";
            for (const auto& var : st->variants) {
                const RunRecord* r = s->history().find(var.runId);
                const auto a = s->aeroIfComputed(var.runId);
                QStringList row{lc.toString(var.value, 'g', 10), r ? r->id : QString(), variantStatusLabel(var.status)};
                for (const Quantity& x : quantities()) {
                    const auto y = a ? quantity(*a, x.key) : std::nullopt;
                    row << (y ? lc.toString(*y, 'g', 10) : QString());
                }
                QStringList m;
                if (r && r->approximate()) m << "APPROXIMATE";
                if (r && r->modified()) m << "MODIFIED";
                if (r && r->synthetic) m << "SYNTHETIC";
                row << m.join(' ');
                csv += row.join(',') + "\n";
            }
            QSaveFile out(file);
            if (out.open(QIODevice::WriteOnly)) {
                out.write(csv.toUtf8());
                out.commit();
            }
        });
        QObject::connect(del, &QPushButton::clicked, del, [s, st, del] {
            if (QMessageBox::question(del->window(), QObject::tr("Delete study"), QObject::tr("Delete %1? Its runs stay in Run History.").arg(st->title()),
                                      QMessageBox::Cancel | QMessageBox::Yes, QMessageBox::Cancel) == QMessageBox::Yes)
                s->deleteStudy(st->id);
        });
    };
    (*build)();
    if (pick) QObject::connect(pick, &QComboBox::activated, page, [build] { (*build)(); });
    QObject::connect(s, &Session::studiesChanged, page, [c] { QTimer::singleShot(0, c.session, [c] { if (c.current && c.current() == "runs/studies") c.refresh(); }); });
    QObject::connect(s, &Session::aeroSpecChanged, page, [build] { (*build)(); });

    QObject::connect(create, &QPushButton::clicked, page, [s, m, param, values, name, follow, createError, c] {
        QString error;
        const auto vals = parseSweepValues(values->text(), &error);
        if (vals.empty()) {
            createError->setText(error);
            createError->show();
            return;
        }
        const QString path = param->currentData().toString();
        const QStringList followers = follow->isVisible() && follow->isChecked() ? referenceFollowers(*m, path) : QStringList{};
        auto study = makeSweep(*m, path, vals, followers);
        if (!study) {
            createError->setText(study.error());
            createError->show();
            return;
        }
        study->name = name->text().trimmed();
        const QString id = s->addStudy(study.value(), &error);
        if (id.isEmpty()) {
            createError->setText(error);
            createError->show();
            return;
        }
        QSettings().setValue("studies/selected", id);
        c.refresh();
    });
    return endPage(page, v);
}

} // namespace kapp::pages
