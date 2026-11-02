// File ▸ New Case from Mesh: a bare mesh, its zones and the role each takes
// in an FSAE external-aero case. Roles start from `kutta new-case`'s
// suggestions (zone names, then mesh zone types); Create writes the Kutta
// case and opens it, after which the usual Setup and Solution pages edit it.

#include "ui/Inspector.h"
#include "ui/PageKit.h"

#include <QComboBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileInfo>
#include <QGridLayout>
#include <QHash>
#include <QLocale>
#include <QTimer>

#include <algorithm>

namespace kapp {

using namespace kit;

namespace pages {

namespace {

QString roleLabel(const QString& role)
{
    static const QHash<QString, QString> labels = {
        {"velocity-inlet", QObject::tr("Velocity inlet")}, {"pressure-outlet", QObject::tr("Pressure outlet")},
        {"moving-ground", QObject::tr("Moving ground")},   {"wall", QObject::tr("Wall")},
        {"symmetry", QObject::tr("Symmetry")},             {"interior", QObject::tr("Interior")}};
    return labels.value(role, role);
}

QWidget* newCasePage(QVBoxLayout** bodyOut)
{
    auto* page = new QWidget;
    auto* v = new QVBoxLayout(page);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    v->addWidget(header(QObject::tr("New Case from Mesh"), QObject::tr("Setup"), page));
    auto* content = new QWidget(page);
    *bodyOut = body(content);
    v->addWidget(content);
    v->addStretch(1);
    return page;
}

} // namespace

QWidget* newCase(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = newCasePage(&v);
    Session* s = c.session;
    const auto& draft = s->newCaseDraft();
    const auto chooseMesh = [c] { if (c.newCaseDialog) c.newCaseDialog(); };

    if (!s->newCaseError().isEmpty()) {
        auto* f = callout("red");
        f->layout()->addWidget(label(draft ? QObject::tr("Could not create the case") : QObject::tr("Could not read the mesh"), "title"));
        f->layout()->addWidget(label(s->newCaseError()));
        v->addWidget(f);
    }
    if (s->newCaseBusy() && !draft) {
        v->addWidget(label(QObject::tr("Reading the mesh…"), "muted"));
        return page;
    }
    if (!draft) {
        v->addWidget(label(QObject::tr("Start from a mesh exported by your mesher (.msh or .msh.h5; a Fluent .cas also works, only its mesh is used). "
                                       "Kutta lists its zones, suggests what each one is from its name, and writes a Kutta case you can edit, run and save."),
                           "muted"));
        auto* b = button(QObject::tr("Choose Mesh…"), "folder-open", "primary");
        QObject::connect(b, &QPushButton::clicked, b, chooseMesh);
        v->addWidget(b, 0, Qt::AlignLeft);
        return page;
    }
    const NewCaseDraft& d = *draft;

    v->addWidget(section(QObject::tr("Mesh")));
    auto* mesh = new QGridLayout;
    mesh->setHorizontalSpacing(12);
    mesh->setVerticalSpacing(5);
    mesh->setColumnStretch(1, 1);
    int row = 0;
    auto factLabel = [&](const QString& k, QLabel* l) {
        mesh->addWidget(label(k, "muted"), row, 0, Qt::AlignTop);
        l->setTextInteractionFlags(Qt::TextSelectableByMouse);
        mesh->addWidget(l, row++, 1);
    };
    auto fact = [&](const QString& k, const QString& value) { factLabel(k, label(value)); };
    factLabel(QObject::tr("File"), pathLabel(d.meshPath));
    fact(QObject::tr("Size"), QObject::tr("%1, %2 cells, %3 faces%4")
                                  .arg(d.dimension == 2 ? QObject::tr("2D") : QObject::tr("3D"), QLocale().toString(d.cells),
                                       QLocale().toString(d.faces), d.polyhedral ? QObject::tr(", polyhedral") : QString()));
    v->addLayout(mesh);
    auto* other = button(QObject::tr("Choose Another Mesh…"), "folder-open");
    other->setProperty("flat", true);
    QObject::connect(other, &QPushButton::clicked, other, chooseMesh);
    v->addWidget(other, 0, Qt::AlignLeft);

    // One template for now: stated as a fact, not a picker with one entry.
    v->addWidget(section(QObject::tr("Template")));
    auto* templ = label(QObject::tr("FSAE external aero"));
    templ->setAccessibleName(QObject::tr("Template: FSAE external aero"));
    QFont tf = templ->font();
    tf.setWeight(QFont::DemiBold);
    templ->setFont(tf);
    v->addWidget(templ);
    v->addWidget(label(d.dimension == 2
                           ? QObject::tr("Flow along +x, up +y. k-ω SST, air, a moving ground at the freestream speed, drag and downforce coefficients on every car wall.")
                           : QObject::tr("Flow along +x, up +z. k-ω SST, air, a moving ground at the freestream speed, drag and downforce coefficients on every car wall."),
                       "muted"));
    auto* speedRow = new QHBoxLayout;
    speedRow->addWidget(label(QObject::tr("Freestream speed"), "muted"));
    auto* speed = new QDoubleSpinBox;
    speed->setAccessibleName(QObject::tr("Freestream speed"));
    speed->setRange(0.1, 400);
    speed->setDecimals(2);
    speed->setSuffix(QObject::tr(" m/s"));
    speed->setValue(d.speed);
    speed->setKeyboardTracking(false);
    QObject::connect(speed, &QDoubleSpinBox::valueChanged, s, [s](double value) { s->setNewCaseSpeed(value); });
    speedRow->addWidget(speed);
    speedRow->addStretch(1);
    v->addLayout(speedRow);

    v->addWidget(section(QObject::tr("Zones")));
    v->addWidget(label(QObject::tr("Each boundary zone needs a role. Suggestions come from the zone names; check them before creating the case."), "muted"));
    auto* zones = new QGridLayout;
    zones->setHorizontalSpacing(10);
    zones->setVerticalSpacing(6);
    zones->setColumnStretch(0, 1);
    row = 0;
    int changed = 0;
    QStringList walls;
    for (const NewCaseZone* z : d.boundaries()) {
        auto* name = new QWidget;
        auto* nv = new QVBoxLayout(name);
        nv->setContentsMargins(0, 0, 0, 0);
        nv->setSpacing(0);
        nv->addWidget(label(z->name));
        const bool overridden = z->role != z->suggestedRole;
        changed += overridden;
        const QString why = overridden ? QObject::tr("Your choice; suggested: %1 (%2)").arg(roleLabel(z->suggestedRole), z->reason)
                                       : z->reason;
        nv->addWidget(label(QObject::tr("%1 faces · %2").arg(QLocale().toString(z->faces), why), "faint"));
        // The template's walls do not move: a wheel kept as a wall stands
        // still under a moving road, which changes the car's loads.
        if (z->role == "wall" && z->group == "Wheels")
            nv->addWidget(label(QObject::tr("Stationary wall: this wheel will not turn (see Before you run, below)"), "warning"));
        zones->addWidget(name, row, 0);
        auto* role = new QComboBox;
        role->setAccessibleName(QObject::tr("Role of %1").arg(z->name));
        for (const QString& r : d.roles)
            if (r != "interior") role->addItem(roleLabel(r), r);
        role->setCurrentIndex(role->findData(z->role));
        const QString zoneName = z->name;
        QObject::connect(role, &QComboBox::activated, role, [s, role, zoneName](int i) {
            const QString r = role->itemData(i).toString();
            // The page rebuilds on the change; let the popup close first.
            QTimer::singleShot(0, s, [s, zoneName, r] { s->setNewCaseRole(zoneName, r); });
        });
        zones->addWidget(role, row, 1, Qt::AlignTop);
        if (z->role == "wall") walls << (z->group.isEmpty() ? z->name : QStringLiteral("%1 (%2)").arg(z->name, z->group));
        ++row;
    }
    v->addLayout(zones);
    const auto interior = std::count_if(d.zones.begin(), d.zones.end(), [](const NewCaseZone& z) { return z.kind == "interior"; });
    if (interior) v->addWidget(label(QObject::tr("%1 interior face zone(s) need no role; the cells are fluid.").arg(interior), "faint"));
    if (walls.isEmpty())
        v->addWidget(label(QObject::tr("No zone is a wall: the case will have no car surfaces and no drag or downforce reports."), "error"));
    else
        v->addWidget(label(QObject::tr("Drag and downforce are reported on: %1.").arg(walls.join(", ")), "muted"));
    // What kutta new-case says to check: easy to miss in faint text, and each
    // one changes the results (reference area, wheels, approximations).
    if (!d.notes.isEmpty()) {
        auto* f = callout("neutral");
        auto* head = label(QObject::tr("Before you run"));
        QFont hf = head->font();
        hf.setWeight(QFont::DemiBold);
        head->setFont(hf);
        f->layout()->addWidget(head);
        for (const QString& note : d.notes) f->layout()->addWidget(label(note));
        // Where the note's rotating walls are set: the app does not edit wall motion yet.
        const bool wheels = std::any_of(d.zones.begin(), d.zones.end(), [](const NewCaseZone& z) { return z.role == "wall" && z.group == "Wheels"; });
        if (wheels)
            f->layout()->addWidget(label(QObject::tr("Wall motion is set in the case file (a boundary's motion, docs/kutta-case.md); "
                                                     "the app shows it but does not edit it yet."),
                                         "muted"));
        v->addWidget(f);
    }

    auto* actions = new QHBoxLayout;
    auto* create = button(QObject::tr("Create Case…"), "plus", "primary");
    create->setEnabled(!s->newCaseBusy() && !s->running());
    create->setToolTip(QObject::tr("Writes a Kutta case (.kutta.json) beside the mesh, or where you choose, and opens it."));
    QObject::connect(create, &QPushButton::clicked, create, [c] { if (c.createCaseDialog) c.createCaseDialog(); });
    actions->addWidget(create);
    auto* cancel = button(QObject::tr("Cancel"));
    cancel->setProperty("flat", true);
    QObject::connect(cancel, &QPushButton::clicked, cancel, [c, s] {
        QTimer::singleShot(0, s, [c, s] {
            s->discardNewCase();
            c.navigate("setup/general");
        });
    });
    actions->addWidget(cancel);
    actions->addStretch(1);
    v->addLayout(actions);
    if (s->newCaseBusy()) v->addWidget(label(QObject::tr("Writing the case…"), "muted"));
    if (changed) v->addWidget(label(QObject::tr("%1 role(s) changed from the suggestions.").arg(changed), "faint"));
    return page;
}

} // namespace pages

} // namespace kapp
