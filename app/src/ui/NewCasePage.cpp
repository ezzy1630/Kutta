// File ▸ New Case from Mesh: a bare mesh, its zones and the role each takes
// in an FSAE external-aero case. Roles start from `kutta new-case`'s
// suggestions (zone names, then mesh zone types); Create writes the Kutta
// case and opens it, after which the usual Setup and Solution pages edit it.

#include "ui/Inspector.h"
#include "ui/PageKit.h"

#include <QComboBox>
#include <QCheckBox>
#include <QDir>
#include <QDoubleSpinBox>
#include <QFileInfo>
#include <QGridLayout>
#include <QHash>
#include <QLocale>
#include <QLineEdit>
#include <QTimer>

#include <algorithm>
#include <cmath>
#include <limits>

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
        QObject::connect(s, &Session::newCaseValuesChanged, f, [s, f] { f->setVisible(!s->newCaseError().isEmpty()); });
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
    speed->setEnabled(!s->newCaseBusy());
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
        if (z->role == "wall") {
            auto* rotate = new QCheckBox(QObject::tr("Rotate as a wheel"));
            rotate->setAccessibleName(QObject::tr("Rotate %1 as a wheel").arg(z->name));
            rotate->setChecked(z->wheel.has_value());
            rotate->setEnabled(!s->newCaseBusy());
            nv->addWidget(rotate);
            const QString name = z->name;
            const int dimension = d.dimension;
            QObject::connect(rotate, &QCheckBox::toggled, rotate, [s, name, dimension](bool on) {
                QTimer::singleShot(0, s, [s, name, dimension, on] {
                    WheelGeometry wheel;
                    wheel.axis = dimension == 2 ? std::array<double, 3>{0, 0, 1} : std::array<double, 3>{0, -1, 0};
                    s->setNewCaseWheel(name, on ? std::optional(wheel) : std::nullopt);
                });
            });
            if (!z->wheel && z->group == "Wheels")
                nv->addWidget(label(QObject::tr("Stationary: enable rotation if this wheel turns."), "warning"));
        }
        zones->addWidget(name, row, 0);
        auto* role = new QComboBox;
        role->setEnabled(!s->newCaseBusy());
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
        if (z->role == "wall" && z->wheel) {
            auto* editor = new QWidget;
            editor->setEnabled(!s->newCaseBusy());
            auto* fields = new QVBoxLayout(editor);
            fields->setContentsMargins(0, 2, 0, 10);
            fields->setSpacing(5);
            auto number = [](double value, const QString& name) {
                auto* spin = new QLineEdit(QLocale().toString(value, 'g', 10));
                spin->setAccessibleName(name);
                spin->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
                return spin;
            };
            auto* radiusRow = new QHBoxLayout;
            radiusRow->addWidget(label(QObject::tr("Radius (m)"), "muted"));
            auto* radius = number(z->wheel->radius, QObject::tr("Radius of %1 (m)").arg(z->name));
            radiusRow->addWidget(radius, 1);
            fields->addLayout(radiusRow);
            auto vector = [&](const QString& title, const QString& accessible, const std::array<double, 3>& value) {
                fields->addWidget(label(title, "muted"));
                auto* grid = new QGridLayout;
                grid->setHorizontalSpacing(6);
                grid->setVerticalSpacing(2);
                std::array<QLineEdit*, 3> spins;
                for (int i = 0; i < 3; ++i) {
                    const QString component = QString("XYZ").mid(i, 1);
                    grid->addWidget(label(component, "faint"), 0, i);
                    spins[i] = number(value[i], accessible.arg(component, z->name));
                    grid->addWidget(spins[i], 1, i);
                    grid->setColumnStretch(i, 1);
                }
                fields->addLayout(grid);
                return spins;
            };
            const auto centre = vector(QObject::tr("Wheel centre (m)"), QObject::tr("Centre %1 of %2 (m)"), z->wheel->centre);
            const auto axis = vector(QObject::tr("Axle direction"), QObject::tr("Axis %1 of %2"), z->wheel->axis);
            fields->addWidget(label(d.dimension == 2
                                       ? QObject::tr("Enter geometry from your mesh in metres; check all three centre coordinates. The axle points out of the 2D plane (X = Y = 0).")
                                       : QObject::tr("Enter geometry from your mesh in metres; check all three centre coordinates. The axle runs across the flow (X = 0)."), "muted"));
            auto* status = label({});
            fields->addWidget(status);
            const QString name = z->name;
            auto update = [s, name, status] {
                if (!s->newCaseDraft()) return;
                const auto* zone = s->newCaseDraft()->zone(name);
                if (!zone || !zone->wheel) return;
                const auto& d = *s->newCaseDraft();
                QString error = zone->wheel->error(d.dimension);
                if (error.isEmpty() && !std::isfinite(d.speed / zone->wheel->radius)) error = QObject::tr("Radius is too small for a finite angular speed.");
                const auto axis = zone->wheel->rollingAxis(d.dimension);
                status->setText(error.isEmpty() ? QObject::tr("Rotating wall · %1 rad/s\nRolling axis: (%2, %3, %4)")
                                   .arg(QLocale().toString(d.speed / zone->wheel->radius, 'g', 6), QString::number(axis[0]), QString::number(axis[1]), QString::number(axis[2])) : error);
                status->setProperty("role", error.isEmpty() ? "muted" : "error");
                status->style()->unpolish(status);
                status->style()->polish(status);
            };
            QObject::connect(s, &Session::newCaseValuesChanged, editor, update);
            update();
            auto save = [s, name, radius, centre, axis] {
                auto value = [](QLineEdit* field) {
                    bool ok = false;
                    const double x = QLocale().toDouble(field->text(), &ok);
                    return ok ? x : std::numeric_limits<double>::quiet_NaN();
                };
                WheelGeometry w;
                w.radius = value(radius);
                for (int i = 0; i < 3; ++i) { w.centre[i] = value(centre[i]); w.axis[i] = value(axis[i]); }
                s->setNewCaseWheel(name, w);
            };
            QObject::connect(radius, &QLineEdit::textEdited, editor, save);
            for (auto* spin : centre) QObject::connect(spin, &QLineEdit::textEdited, editor, save);
            for (auto* spin : axis) QObject::connect(spin, &QLineEdit::textEdited, editor, save);
            zones->addWidget(editor, row++, 0, 1, 2);
        }
    }
    v->addLayout(zones);
    v->addWidget(label(QObject::tr("Wheel angular speed = freestream speed ÷ radius (rad/s). Kutta chooses the axle's sign so the wheel's bottom moves along +x with the moving ground."), "muted"));
    const auto interior = std::count_if(d.zones.begin(), d.zones.end(), [](const NewCaseZone& z) { return z.kind == "interior"; });
    if (interior) v->addWidget(label(QObject::tr("%1 interior face zone(s) need no role; the cells are fluid.").arg(interior), "faint"));
    if (walls.isEmpty())
        v->addWidget(label(QObject::tr("No zone is a wall: the case will have no car surfaces and no drag or downforce reports."), "error"));
    else
        v->addWidget(label(QObject::tr("Drag and downforce are reported on: %1.").arg(walls.join(", ")), "muted"));
    // What kutta new-case says to check: easy to miss in faint text, and each
    // one changes the results (reference area, wheels, approximations).
    {
        auto* f = callout("neutral");
        auto* head = label(QObject::tr("Before you run"));
        QFont hf = head->font();
        hf.setWeight(QFont::DemiBold);
        head->setFont(hf);
        f->layout()->addWidget(head);
        auto* notes = label({});
        f->layout()->addWidget(notes);
        auto updateNotes = [s, notes] {
            if (!s->newCaseDraft()) return;
            if (!s->newCaseDraft()->wheelError().isEmpty()) notes->setText(QObject::tr("Complete the wheel geometry to check this setup."));
            else if (!s->newCaseNotesError().isEmpty()) notes->setText(withHomeAsTilde(s->newCaseNotesError()));
            else if (!s->newCaseNotesCurrent()) notes->setText(QObject::tr("Checking the current choices…"));
            else {
                QStringList text;
                for (const auto& note : s->newCaseDraft()->notes)
                    text << (note.startsWith("Wheels without --wheel") || note.startsWith("Wheels are stationary walls")
                                 ? QObject::tr("Some wheel zones are stationary. Enable ‘Rotate as a wheel’ and enter their geometry if they should turn.") : note);
                notes->setText(withHomeAsTilde(text.join("\n\n")));
            }
        };
        QObject::connect(s, &Session::newCaseNotesChanged, notes, updateNotes);
        updateNotes();
        v->addWidget(f);
    }

    auto* actions = new QHBoxLayout;
    auto* create = button(QObject::tr("Create Case…"), "plus", "primary");
    auto updateCreate = [s, create] {
        create->setEnabled(s->newCaseDraft() && s->newCaseDraft()->wheelError().isEmpty() && !s->newCaseBusy() && !s->running());
    };
    QObject::connect(s, &Session::newCaseValuesChanged, create, updateCreate);
    updateCreate();
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
