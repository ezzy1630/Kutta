// Meshing: New Case from Geometry in Fluent Meshing's order. Import the STL
// solids, set the wind tunnel, refinement regions, surface sizes and boundary
// layers, generate the volume mesh with `kutta mesh`, read its report, then
// set up the case. Design Variants (Runs) compares a design's variants.

#include "ui/Inspector.h"
#include "ui/ChartsPanel.h"
#include "ui/FlowLayout.h"
#include "ui/GraphicsArea.h"
#include "ui/PageKit.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFileInfo>
#include <QGridLayout>
#include <QHeaderView>
#include <QLineEdit>
#include <QLocale>
#include <QPainter>
#include <QTableWidget>
#include <QTimer>

#include <cmath>

namespace kapp {

using namespace kit;

namespace {

QWidget* makePage(const QString& title, QVBoxLayout** bodyOut)
{
    auto* page = new QWidget;
    auto* v = new QVBoxLayout(page);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    v->addWidget(pages::header(title, QObject::tr("Meshing"), page));
    auto* content = new QWidget(page);
    *bodyOut = body(content);
    v->addWidget(content);
    v->addStretch(1);
    return page;
}

QString count(double n)
{
    if (n >= 1e6) return QObject::tr("%1 M").arg(QLocale().toString(n / 1e6, 'f', n >= 1e7 ? 1 : 2));
    if (n >= 1e4) return QObject::tr("%1 k").arg(QLocale().toString(n / 1e3, 'f', 0));
    return QLocale().toString(qint64(std::llround(n)));
}

QString bytes(double b)
{
    return b >= 1e9 ? QObject::tr("%1 GB").arg(QLocale().toString(b / 1e9, 'f', 1)) : QObject::tr("%1 MB").arg(QLocale().toString(b / 1e6, 'f', 0));
}

// A length shown and typed in the project's unit (mm for most CAD), kept in metres.
// Lengths to 0.1 mm (or 1 µm in metres), without trailing zeros.
QString lengthText(double metres, LengthUnit u)
{
    const double x = metres / unitScale(u);
    const int decimals = u == LengthUnit::Millimetre ? 1 : 6;
    QString t = QLocale().toString(x, 'f', decimals);
    const QChar point = QLocale().decimalPoint().at(0);
    if (t.contains(point)) {
        while (t.endsWith('0')) t.chop(1);
        if (t.endsWith(point)) t.chop(1);
    }
    return t == "-0" ? QStringLiteral("0") : t;
}

std::optional<double> parseNumber(const QString& text)
{
    bool ok = false;
    const double x = QLocale().toDouble(text.trimmed(), &ok);
    if (ok && std::isfinite(x)) return x;
    const double y = QLocale::c().toDouble(text.trimmed(), &ok);
    return ok && std::isfinite(y) ? std::optional(y) : std::nullopt;
}

// A numeric field that applies on Enter or leaving it, and shows its unit.
QWidget* numberField(const QString& text, const QString& unit, const QString& accessible, std::function<bool(double)> apply)
{
    auto* w = new QWidget;
    auto* h = new QHBoxLayout(w);
    h->setContentsMargins(0, 0, 0, 0);
    h->setSpacing(6);
    auto* e = new QLineEdit(text);
    e->setAccessibleName(accessible);
    e->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    e->setMinimumWidth(70);
    h->addWidget(e, 1);
    if (!unit.isEmpty()) {
        auto* u = label(unit, "unit");
        u->setWordWrap(false);
        h->addWidget(u);
    }
    QObject::connect(e, &QLineEdit::editingFinished, e, [e, apply, text] {
        const auto x = parseNumber(e->text());
        const bool ok = x && apply(*x);
        e->setProperty("invalid", !ok);
        repolish(e);
        if (!ok) e->setToolTip(QObject::tr("Not applied: enter a valid number."));
        else e->setToolTip({});
    });
    return w;
}

// Two columns: what, then the field.
QGridLayout* form()
{
    auto* g = new QGridLayout;
    g->setHorizontalSpacing(12);
    g->setVerticalSpacing(6);
    g->setColumnStretch(0, 1);
    g->setColumnStretch(1, 1);
    return g;
}

void formRow(QGridLayout* g, const QString& what, QWidget* field, const QString& hint = {})
{
    const int row = g->rowCount();
    g->addWidget(label(what, "muted"), row, 0, Qt::AlignVCenter);
    g->addWidget(field, row, 1);
    if (!hint.isEmpty()) g->addWidget(label(hint, "faint"), row + 1, 0, 1, 2);
}

QPushButton* nextButton(const InspectorContext& c, const QString& id, const QString& text)
{
    auto* b = button(text, "chevron-right");
    b->setProperty("flat", true);
    QObject::connect(b, &QPushButton::clicked, b, [c, id] { c.navigate(id); });
    return b;
}

// Edits the open project's settings through a copy, as one change.
void editSetup(Session* s, const std::function<void(MeshSetup&)>& fn)
{
    MeshSetup m = s->meshing()->project().setup;
    fn(m);
    s->meshing()->setSetup(m);
}

QWidget* noProject(const InspectorContext& c, const QString& title)
{
    QVBoxLayout* v;
    QWidget* page = makePage(title, &v);
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->error().isEmpty()) {
        auto* f = callout("red");
        f->layout()->addWidget(label(QObject::tr("Could not import the geometry"), "title"));
        f->layout()->addWidget(label(mw->error()));
        v->addWidget(f);
    }
    if (mw->loading()) {
        v->addWidget(label(QObject::tr("Reading the geometry…"), "muted"));
        return page;
    }
    v->addWidget(label(QObject::tr("Start from your car's CAD exported as STL: one or more files, each solid a part (front wing, rear wing, wheels, body). "
                                   "Kutta shows the surfaces at once, puts a wind tunnel around them, meshes it with kutta mesh and sets up the case."),
                       "muted"));
    auto* b = button(QObject::tr("Import Geometry…"), "file-input", "primary");
    QObject::connect(b, &QPushButton::clicked, b, [c] { if (c.importGeometryDialog) c.importGeometryDialog(); });
    v->addWidget(b, 0, Qt::AlignLeft);
    if (c.session->cliInfo().fixture && c.sampleGeometry) {
        v->addSpacing(6);
        v->addWidget(section(QObject::tr("Sample")));
        v->addWidget(label(QObject::tr("The fixture CLI's synthetic car as STL in millimetres, with a second, larger rear wing to try a design variant. SYNTHETIC: its mesh and results are invented."), "faint"));
        auto* sb = button(QObject::tr("Import the Synthetic Car"), "box");
        sb->setProperty("flat", true);
        QObject::connect(sb, &QPushButton::clicked, sb, [c] { c.sampleGeometry(); });
        v->addWidget(sb, 0, Qt::AlignLeft);
    }
    return page;
}

// Thin bar: the mesher's progress through its stages.
class ProgressStrip : public QWidget {
public:
    explicit ProgressStrip(QWidget* parent = nullptr) : QWidget(parent) { setFixedHeight(6); }
    void setValue(double f) { m_value = std::clamp(f, 0.0, 1.0); update(); }
    QSize sizeHint() const override { return {200, 6}; }

protected:
    void paintEvent(QPaintEvent*) override
    {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing);
        p.setPen(Qt::NoPen);
        p.setBrush(tokens().border);
        p.drawRoundedRect(rect(), 3, 3);
        p.setBrush(tokens().accent);
        p.drawRoundedRect(QRectF(0, 0, width() * m_value, height()), 3, 3);
    }

private:
    double m_value = 0;
};

// A quality metric's distribution as bars; the worse end is tinted.
class Histogram : public QWidget {
public:
    Histogram(const MeshQuality& q, QWidget* parent = nullptr) : QWidget(parent), m_q(q)
    {
        setMinimumHeight(56);
        setAccessibleName(QObject::tr("%1 histogram").arg(q.label));
    }
    QSize sizeHint() const override { return {240, 56}; }

protected:
    void paintEvent(QPaintEvent*) override
    {
        QPainter p(this);
        const auto& t = tokens();
        const int n = int(m_q.counts.size());
        if (n == 0) return;
        qint64 peak = 1;
        for (qint64 c : m_q.counts) peak = std::max(peak, c);
        const double w = double(width()) / n;
        for (int i = 0; i < n; ++i) {
            // Log scale: the few poor cells matter more than the many good ones.
            const double f = m_q.counts[std::size_t(i)] > 0 ? std::log10(1.0 + double(m_q.counts[std::size_t(i)])) / std::log10(1.0 + double(peak)) : 0;
            const double h = std::max(f > 0 ? 2.0 : 0.0, f * (height() - 2));
            const bool poor = m_q.higherIsBetter ? i < n / 5 : i >= n - n / 4;
            p.fillRect(QRectF(i * w + 1, height() - h, w - 2, h), poor ? t.amber : t.accent);
        }
        p.setPen(t.border);
        p.drawLine(0, height() - 1, width(), height() - 1);
    }

private:
    MeshQuality m_q;
};

std::optional<Device> selectedDevice(const Session* s)
{
    for (const auto& d : s->devices())
        if (s->backend() == "auto" ? d.isDefault : d.backend == s->backend()) return d;
    if (!s->devices().empty()) return s->devices().front();
    return std::nullopt;
}

} // namespace

namespace pages {

// ---- Import Geometry ------------------------------------------------------------------

QWidget* geometry(const InspectorContext& c)
{
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Import Geometry"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Import Geometry"), &v);
    const GeometryProject& p = mw->project();
    if (!mw->error().isEmpty()) {
        auto* f = callout("red");
        f->layout()->addWidget(label(mw->error()));
        v->addWidget(f);
    }
    if (mw->loading()) v->addWidget(label(QObject::tr("Reading the geometry…"), "muted"));

    auto* facts = form();
    facts->addWidget(label(QObject::tr("Project"), "muted"), 0, 0);
    facts->addWidget(pathLabel(p.dir), 0, 1);
    auto* units = new QComboBox;
    units->setAccessibleName(QObject::tr("Geometry units"));
    units->addItem(QObject::tr("Millimetres (mm)"), "mm");
    units->addItem(QObject::tr("Metres (m)"), "m");
    units->setCurrentIndex(units->findData(unitName(p.setup.units)));
    units->setEnabled(!mw->meshing());
    QObject::connect(units, &QComboBox::activated, units, [mw, units](int i) {
        const LengthUnit u = unitFromName(units->itemData(i).toString());
        QTimer::singleShot(0, mw, [mw, u] { mw->setUnits(u); });
    });
    formRow(facts, QObject::tr("Units"), units);
    const Bounds b = mw->geometryBounds();
    if (b.valid) {
        const QVector3D e = b.max - b.min;
        formRow(facts, QObject::tr("Size (L × W × H)"),
                label(QObject::tr("%1 × %2 × %3 m").arg(QLocale().toString(e.x(), 'f', 3), QLocale().toString(e.y(), 'f', 3), QLocale().toString(e.z(), 'f', 3))));
        formRow(facts, QObject::tr("Surface"), label(QObject::tr("%1 m², %2 triangles").arg(QLocale().toString(mw->surfaceArea(), 'f', 2), count(double(mw->triangles())))));
    }
    v->addLayout(facts);
    if (b.valid && std::max({b.max.x() - b.min.x(), b.max.y() - b.min.y(), b.max.z() - b.min.z()}) > 20)
        v->addWidget(label(QObject::tr("The geometry is over 20 m in metres: if it was exported in millimetres, choose Millimetres."), "warning"));
    if (!p.variantOf.isEmpty())
        v->addWidget(label(QObject::tr("Design variant of %1. Its case will reuse the settings of %2.")
                               .arg(QFileInfo(p.variantOf).fileName(), p.baseCase.isEmpty() ? QObject::tr("the base case once it exists") : QFileInfo(p.baseCase).fileName()),
                           "muted"));

    v->addWidget(section(QObject::tr("Files")));
    for (int i = 0; i < int(mw->files().size()); ++i) {
        const auto& f = mw->files()[std::size_t(i)];
        auto* row = new QHBoxLayout;
        auto* name = new QVBoxLayout;
        name->setSpacing(0);
        name->addWidget(pathLabel(f.path));
        name->addWidget(label(f.error.isEmpty() ? QObject::tr("%1 triangles · %2").arg(count(double(f.triangles)), f.binary ? QObject::tr("binary STL") : QObject::tr("ASCII STL"))
                                                : f.error,
                              f.error.isEmpty() ? "faint" : "error"));
        row->addLayout(name, 1);
        auto* replace = button(QObject::tr("Replace…"), "refresh-cw");
        replace->setProperty("flat", true);
        replace->setAccessibleName(QObject::tr("Replace %1").arg(QFileInfo(f.path).fileName()));
        replace->setToolTip(QObject::tr("Swap this file for another (a new rear wing, say); solids with the same name keep their zone name and role."));
        replace->setEnabled(!mw->meshing());
        QObject::connect(replace, &QPushButton::clicked, replace, [c, i] { if (c.replaceGeometryDialog) c.replaceGeometryDialog(i); });
        row->addWidget(replace, 0, Qt::AlignTop);
        auto* remove = button({}, "trash-2");
        remove->setProperty("flat", true);
        remove->setAccessibleName(QObject::tr("Remove %1").arg(QFileInfo(f.path).fileName()));
        remove->setToolTip(QObject::tr("Remove this file from the geometry"));
        remove->setEnabled(!mw->meshing() && mw->files().size() > 1);
        QObject::connect(remove, &QPushButton::clicked, remove, [mw, i] { QTimer::singleShot(0, mw, [mw, i] { mw->removeGeometry(i); }); });
        row->addWidget(remove, 0, Qt::AlignTop);
        v->addLayout(row);
    }
    auto* add = button(QObject::tr("Add Geometry…"), "plus");
    add->setProperty("flat", true);
    add->setEnabled(!mw->meshing());
    QObject::connect(add, &QPushButton::clicked, add, [c] { if (c.addGeometryDialog) c.addGeometryDialog(); });
    v->addWidget(add, 0, Qt::AlignLeft);

    v->addWidget(section(QObject::tr("Solids → wall zones")));
    v->addWidget(label(QObject::tr("Each solid becomes a wall zone of the car, named here. Wheels rotate with the moving ground once the case is set up (Kutta proposes each axle from the tyre surface for you to check)."), "faint"));
    auto* grid = new QGridLayout;
    grid->setHorizontalSpacing(8);
    grid->setVerticalSpacing(5);
    grid->setColumnStretch(0, 3);
    grid->setColumnStretch(1, 2);
    for (int i = 0; i < int(p.solids.size()); ++i) {
        const auto& s = p.solids[std::size_t(i)];
        auto* name = new QLineEdit(s.name);
        name->setAccessibleName(QObject::tr("Zone name of solid %1").arg(s.source));
        name->setToolTip(s.source == s.name ? QObject::tr("Solid %1 in %2").arg(s.source, QFileInfo(s.file).fileName())
                                            : QObject::tr("Solid %1 in %2, renamed").arg(s.source, QFileInfo(s.file).fileName()));
        name->setEnabled(!mw->meshing());
        QObject::connect(name, &QLineEdit::editingFinished, name, [mw, name, i] {
            const QString n = name->text();
            QTimer::singleShot(0, mw, [mw, i, n] { mw->renameSolid(i, n); });
        });
        grid->addWidget(name, i, 0);
        auto* role = new QComboBox;
        role->setAccessibleName(QObject::tr("Role of %1").arg(s.name));
        role->addItem(QObject::tr("Wall"), false);
        role->addItem(QObject::tr("Rotating wheel"), true);
        role->setCurrentIndex(s.wheel ? 1 : 0);
        role->setEnabled(!mw->meshing());
        QObject::connect(role, &QComboBox::activated, role, [mw, role, i](int k) {
            const bool wheel = role->itemData(k).toBool();
            QTimer::singleShot(0, mw, [mw, i, wheel] { mw->setSolidWheel(i, wheel); });
        });
        grid->addWidget(role, i, 1);
    }
    v->addLayout(grid);
    v->addWidget(label(QObject::tr("The wind tunnel's zones are named by the mesher: inlet (velocity inlet), outlet (pressure outlet), ground (moving ground), sides and top (symmetry)."), "faint"));

    auto* actions = new FlowLayout;
    actions->addWidget(nextButton(c, "mesh/tunnel", QObject::tr("Next: Wind Tunnel")));
    auto* dup = button(QObject::tr("Duplicate as Design Variant"), "copy");
    dup->setProperty("flat", true);
    dup->setToolTip(QObject::tr("A copy of this project beside it, with the same geometry, mesh settings and case settings: swap or add a part and compare."));
    dup->setEnabled(!mw->meshing() && !mw->loading());
    QObject::connect(dup, &QPushButton::clicked, dup, [c] { if (c.duplicateVariant) c.duplicateVariant(); });
    actions->addWidget(dup);
    v->addLayout(actions);
    return page;
}

// ---- Wind Tunnel ----------------------------------------------------------------------

QWidget* tunnel(const InspectorContext& c)
{
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Wind Tunnel"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Wind Tunnel"), &v);
    Session* s = c.session;
    const MeshSetup& m = mw->project().setup;
    v->addWidget(label(QObject::tr("The box of air around the car: distances in car lengths from the geometry. The outlines show in the viewport as you type."), "muted"));
    const Bounds g = mw->geometryBounds();
    const double L = g.valid ? double(g.max[parseAxis(m.flow).first] - g.min[parseAxis(m.flow).first]) : 0;
    if (L > 0) v->addWidget(label(QObject::tr("Car length %1 m along the flow").arg(QLocale().toString(L, 'f', 3)), "faint"));

    v->addWidget(section(QObject::tr("Directions")));
    auto* dirs = form();
    auto axisBox = [&](const QString& value, const QString& accessible, bool flow) {
        auto* b = new QComboBox;
        b->setAccessibleName(accessible);
        for (const char* a : {"+x", "-x", "+y", "-y", "+z", "-z"}) b->addItem(a, a);
        b->setCurrentIndex(b->findData(value));
        QObject::connect(b, &QComboBox::activated, b, [s, b, flow](int i) {
            const QString a = b->itemData(i).toString();
            QTimer::singleShot(0, s, [s, a, flow] { editSetup(s, [&](MeshSetup& x) { (flow ? x.flow : x.up) = a; }); });
        });
        return b;
    };
    formRow(dirs, QObject::tr("Flow direction"), axisBox(m.flow, QObject::tr("Flow direction"), true));
    formRow(dirs, QObject::tr("Up"), axisBox(m.up, QObject::tr("Up direction"), false));
    v->addLayout(dirs);

    v->addWidget(section(QObject::tr("Extent (car lengths)")));
    auto* ext = form();
    auto lengths = [&](const QString& what, double value, double TunnelSetup::* member) {
        formRow(ext, what, numberField(QLocale().toString(value, 'g', 4), QObject::tr("L"), what, [s, member](double x) {
            if (!(x > 0) || x > 100) return false;
            editSetup(s, [&](MeshSetup& m) { m.tunnel.*member = x; });
            return true;
        }));
    };
    lengths(QObject::tr("Upstream of the nose"), m.tunnel.upstream, &TunnelSetup::upstream);
    lengths(QObject::tr("Downstream of the tail"), m.tunnel.downstream, &TunnelSetup::downstream);
    lengths(QObject::tr("Each side"), m.tunnel.side, &TunnelSetup::side);
    lengths(QObject::tr("Above the car"), m.tunnel.top, &TunnelSetup::top);
    v->addLayout(ext);

    v->addWidget(section(QObject::tr("Ground and symmetry")));
    auto* ground = new QCheckBox(QObject::tr("Moving ground at the car's lowest point"));
    ground->setAccessibleName(QObject::tr("Moving ground"));
    ground->setChecked(m.tunnel.ground);
    QObject::connect(ground, &QCheckBox::toggled, ground, [s](bool on) { editSetup(s, [&](MeshSetup& m) { m.tunnel.ground = on; }); });
    v->addWidget(ground);
    auto* half = new QCheckBox(QObject::tr("Half car: symmetry plane at y = 0"));
    half->setAccessibleName(QObject::tr("Symmetry plane at y = 0"));
    half->setChecked(m.tunnel.symmetry == "y");
    QObject::connect(half, &QCheckBox::toggled, half, [s](bool on) { editSetup(s, [&](MeshSetup& m) { m.tunnel.symmetry = on ? "y" : QString(); }); });
    v->addWidget(half);
    v->addWidget(label(QObject::tr("Half a car halves the cells; forces are then for half the car."), "faint"));
    auto* tb = label({}, "faint");
    auto update = [mw, tb] {
        if (!mw->hasProject()) return;
        const Bounds t = tunnelBounds(mw->project().setup, mw->geometryBounds());
        if (!t.valid) return tb->setText({});
        const QVector3D e = t.max - t.min;
        tb->setText(QObject::tr("Tunnel %1 × %2 × %3 m").arg(QLocale().toString(e.x(), 'f', 1), QLocale().toString(e.y(), 'f', 1), QLocale().toString(e.z(), 'f', 1)));
    };
    QObject::connect(mw, &MeshWorkflow::setupChanged, tb, update);
    update();
    v->addWidget(tb);
    v->addWidget(nextButton(c, "mesh/refinement", QObject::tr("Next: Local Refinement Regions")), 0, Qt::AlignLeft);
    return page;
}

// ---- Local Refinement Regions --------------------------------------------------------

QWidget* refinement(const InspectorContext& c)
{
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Local Refinement Regions"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Local Refinement Regions"), &v);
    Session* s = c.session;
    const MeshSetup& m = mw->project().setup;
    const LengthUnit u = m.units;
    const QString un = unitName(u);
    v->addWidget(label(QObject::tr("Boxes where cells are no larger than their size: around the car and through its wake. Coordinates in %1.").arg(un), "muted"));
    for (int i = 0; i < int(m.boxes.size()); ++i) {
        const RefinementBox& b = m.boxes[std::size_t(i)];
        auto* head = new QHBoxLayout;
        auto* name = new QLineEdit(b.name);
        name->setAccessibleName(QObject::tr("Name of refinement region %1").arg(i + 1));
        QObject::connect(name, &QLineEdit::editingFinished, name, [s, name, i] {
            const QString n = name->text().trimmed();
            if (n.isEmpty()) return;
            editSetup(s, [&](MeshSetup& m) { if (i < int(m.boxes.size())) m.boxes[std::size_t(i)].name = n; });
        });
        head->addWidget(name, 1);
        auto* del = button({}, "trash-2");
        del->setProperty("flat", true);
        del->setAccessibleName(QObject::tr("Remove refinement region %1").arg(b.name));
        QObject::connect(del, &QPushButton::clicked, del, [c, s, i] {
            QTimer::singleShot(0, s, [c, s, i] {
                editSetup(s, [&](MeshSetup& m) { if (i < int(m.boxes.size())) m.boxes.erase(m.boxes.begin() + i); });
                c.refresh();
            });
        });
        head->addWidget(del);
        v->addLayout(head);
        auto* grid = new QGridLayout;
        grid->setHorizontalSpacing(6);
        grid->setVerticalSpacing(3);
        for (int k = 0; k < 3; ++k) {
            grid->addWidget(label(QString("XYZ").mid(k, 1), "faint"), 0, k + 1);
            grid->setColumnStretch(k + 1, 1);
        }
        auto coord = [&](int row, const QString& what, const QVector3D& value, bool isMax) {
            grid->addWidget(label(what, "muted"), row, 0);
            for (int k = 0; k < 3; ++k) {
                auto* e = new QLineEdit(lengthText(value[k], u));
                e->setAccessibleName(QObject::tr("%1 %2 of %3 (%4)").arg(what, QString("XYZ").mid(k, 1), b.name, un));
                e->setMinimumWidth(48);
                QObject::connect(e, &QLineEdit::editingFinished, e, [s, e, i, k, isMax, u] {
                    const auto x = parseNumber(e->text());
                    if (!x) return;
                    editSetup(s, [&](MeshSetup& m) {
                        if (i >= int(m.boxes.size())) return;
                        auto& box = m.boxes[std::size_t(i)];
                        (isMax ? box.max : box.min)[k] = float(*x * unitScale(u));
                    });
                });
                grid->addWidget(e, row, k + 1);
            }
        };
        coord(1, QObject::tr("Min"), b.min, false);
        coord(2, QObject::tr("Max"), b.max, true);
        v->addLayout(grid);
        auto* sizeRow = form();
        formRow(sizeRow, QObject::tr("Cell size"), numberField(lengthText(b.size, u), un, QObject::tr("Cell size of %1 (%2)").arg(b.name, un), [s, i, u](double x) {
            if (!(x > 0)) return false;
            editSetup(s, [&](MeshSetup& m) { if (i < int(m.boxes.size())) m.boxes[std::size_t(i)].size = x * unitScale(u); });
            return true;
        }));
        v->addLayout(sizeRow);
        v->addWidget(rule());
    }
    if (m.boxes.empty()) v->addWidget(label(QObject::tr("No refinement regions: cells grow from the surface size straight to the base size."), "faint"));
    auto* row = new FlowLayout;
    auto* add = button(QObject::tr("Add Region"), "plus");
    add->setProperty("flat", true);
    QObject::connect(add, &QPushButton::clicked, add, [c, s, mw] {
        QTimer::singleShot(0, s, [c, s, mw] {
            const Bounds g = mw->geometryBounds();
            editSetup(s, [&](MeshSetup& m) {
                RefinementBox b;
                b.name = QObject::tr("region-%1").arg(m.boxes.size() + 1);
                b.min = g.valid ? g.min : QVector3D();
                b.max = g.valid ? g.max : QVector3D(1, 1, 1);
                b.size = m.maxSize * 2;
                m.boxes.push_back(b);
            });
            c.refresh();
        });
    });
    row->addWidget(add);
    auto* reset = button(QObject::tr("Suggested Regions"), "wand-sparkles");
    reset->setProperty("flat", true);
    reset->setToolTip(QObject::tr("Replace the regions with a box around the car and one through its near wake"));
    QObject::connect(reset, &QPushButton::clicked, reset, [c, s, mw] {
        QTimer::singleShot(0, s, [c, s, mw] {
            const MeshSetup d = defaultMeshSetup(mw->geometryBounds(), mw->project().setup.units);
            editSetup(s, [&](MeshSetup& m) { m.boxes = d.boxes; });
            c.refresh();
        });
    });
    row->addWidget(reset);
    row->addWidget(nextButton(c, "mesh/sizing", QObject::tr("Next: Sizing and Boundary Layers")));
    v->addLayout(row);
    return page;
}

// ---- Sizing and Boundary Layers ---------------------------------------------------------

QWidget* sizing(const InspectorContext& c)
{
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Sizing and Boundary Layers"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Sizing and Boundary Layers"), &v);
    Session* s = c.session;
    const MeshSetup& m = mw->project().setup;
    const LengthUnit u = m.units;
    const QString un = unitName(u);
    auto length = [&](QGridLayout* f, const QString& what, double value, std::function<void(MeshSetup&, double)> set, const QString& hint = {}) {
        formRow(f, what, numberField(lengthText(value, u), un, QObject::tr("%1 (%2)").arg(what, un), [s, set, u](double x) {
            if (!(x > 0)) return false;
            editSetup(s, [&](MeshSetup& m) { set(m, x * unitScale(u)); });
            return true;
        }), hint);
    };

    v->addWidget(section(QObject::tr("Surface mesh")));
    auto* sf = form();
    length(sf, QObject::tr("Minimum size"), m.minSize, [](MeshSetup& m, double x) { m.minSize = x; });
    length(sf, QObject::tr("Maximum size"), m.maxSize, [](MeshSetup& m, double x) { m.maxSize = x; });
    formRow(sf, QObject::tr("Curvature normal angle"), numberField(QLocale().toString(m.curvatureDegrees, 'g', 4), QObject::tr("°"), QObject::tr("Curvature normal angle"), [s](double x) {
        if (!(x > 0 && x < 90)) return false;
        editSetup(s, [&](MeshSetup& m) { m.curvatureDegrees = x; });
        return true;
    }));
    formRow(sf, QObject::tr("Cells per gap"), numberField(QString::number(m.proximityCells), {}, QObject::tr("Cells per gap"), [s](double x) {
        if (x < 0 || x > 20 || x != std::floor(x)) return false;
        editSetup(s, [&](MeshSetup& m) { m.proximityCells = int(x); });
        return true;
    }));
    v->addLayout(sf);
    v->addWidget(label(QObject::tr("Surface cells shrink from the maximum towards the minimum size where the surface curves more than the normal angle per cell, and to fit the cells per gap between close surfaces (wing slots, wheel and floor)."), "faint"));

    v->addWidget(section(QObject::tr("Boundary layers")));
    auto* bl = form();
    formRow(bl, QObject::tr("Number of layers"), numberField(QString::number(m.layers), {}, QObject::tr("Number of layers"), [s](double x) {
        if (x < 0 || x > 60 || x != std::floor(x)) return false;
        editSetup(s, [&](MeshSetup& m) { m.layers = int(x); });
        return true;
    }));
    length(bl, QObject::tr("First layer height"), m.firstHeight, [](MeshSetup& m, double x) { m.firstHeight = x; });
    formRow(bl, QObject::tr("Growth rate"), numberField(QLocale().toString(m.growth, 'g', 4), {}, QObject::tr("Growth rate"), [s](double x) {
        if (!(x >= 1 && x <= 3)) return false;
        editSetup(s, [&](MeshSetup& m) { m.growth = x; });
        return true;
    }));
    v->addLayout(bl);
    auto* yp = label({}, "faint");
    auto updateYPlus = [mw, yp] {
        if (!mw->hasProject()) return;
        const Bounds g = mw->geometryBounds();
        const double L = g.valid ? double(g.max.x() - g.min.x()) : 1.0;
        const auto& m = mw->project().setup;
        const double total = m.layers > 0 ? m.firstHeight * (std::abs(m.growth - 1) < 1e-9 ? m.layers : (std::pow(m.growth, m.layers) - 1) / (m.growth - 1)) : 0;
        yp->setText(QObject::tr("First-cell y+ ≈ %1 at 20 m/s over the car length (flat-plate estimate). Layers %2 %3 thick on every wall.")
                        .arg(QLocale().toString(yPlusForFirstHeight(m.firstHeight, 20, L), 'f', 0))
                        .arg(lengthText(total, m.units), unitName(m.units)));
    };
    QObject::connect(mw, &MeshWorkflow::setupChanged, yp, updateYPlus);
    updateYPlus();
    v->addWidget(yp);
    v->addWidget(label(QObject::tr("Kutta's SST walls work across y+ (wall functions above about 30, resolved below 1). Layers grow on every solid (zones \"*\")."), "faint"));
    v->addWidget(nextButton(c, "mesh/generate", QObject::tr("Next: Generate the Volume Mesh")), 0, Qt::AlignLeft);
    return page;
}

// ---- Generate the Volume Mesh -------------------------------------------------------------

QWidget* generateMesh(const InspectorContext& c)
{
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Generate the Volume Mesh"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Generate the Volume Mesh"), &v);
    Session* s = c.session;
    const MeshSetup& m = mw->project().setup;
    const LengthUnit u = m.units;
    auto* vf = form();
    formRow(vf, QObject::tr("Base size"), numberField(lengthText(m.baseSize, u), unitName(u), QObject::tr("Base size (%1)").arg(unitName(u)), [s, u](double x) {
        if (!(x > 0)) return false;
        editSetup(s, [&](MeshSetup& m) { m.baseSize = x * unitScale(u); });
        return true;
    }), QObject::tr("The largest cells, in the far field. Octree cells halve from it down to the surface and refinement sizes."));
    formRow(vf, QObject::tr("Maximum cells"), numberField(QLocale().toString(m.maxCells), {}, QObject::tr("Maximum cells"), [s](double x) {
        if (x < 1000 || x > 2e9) return false;
        editSetup(s, [&](MeshSetup& m) { m.maxCells = qint64(x); });
        return true;
    }), QObject::tr("kutta mesh stops with an error rather than build more cells than this."));
    vf->setEnabled(!mw->meshing());
    v->addLayout(vf);

    // The estimate against the device the case will run on.
    v->addWidget(section(QObject::tr("Estimate")));
    auto* est = label({});
    auto* mem = label({}, "muted");
    auto updateEstimate = [s, mw, est, mem] {
        if (!mw->hasProject()) return;
        const CellEstimate e = mw->estimate();
        est->setText(QObject::tr("≈ %1 cells (surface band %2, boundary layers %3, refinement regions %4, far field %5)")
                         .arg(count(e.total()), count(e.octreeCells), count(e.layerCells), count(e.boxCells), count(e.farFieldCells)));
        const auto dev = selectedDevice(s);
        QString text = dev ? QObject::tr("Solving needs about %1 of %2 on %3.").arg(bytes(e.solverBytes()), bytes(double(dev->memoryBytes)), dev->name)
                           : QObject::tr("Solving needs about %1 of device memory.").arg(bytes(e.solverBytes()));
        QString role = "muted";
        if (e.overLimit) {
            text += QObject::tr(" Over the maximum cells: raise it or coarsen the sizes.");
            role = "error";
        } else if (dev && e.solverBytes() > 0.85 * double(dev->memoryBytes)) {
            text += QObject::tr(" That is more than this device holds comfortably: coarsen the sizes or use a half car.");
            role = "warning";
        }
        mem->setText(text);
        mem->setProperty("role", role);
        repolish(mem);
    };
    QObject::connect(mw, &MeshWorkflow::setupChanged, est, updateEstimate);
    QObject::connect(mw, &MeshWorkflow::projectChanged, est, updateEstimate);
    QObject::connect(s, &Session::devicesChanged, est, updateEstimate);
    updateEstimate();
    v->addWidget(est);
    v->addWidget(mem);
    v->addWidget(label(QObject::tr("A planning figure (±50 %) from the surface area and sizes, at %1 kB of device memory per cell. The mesh report gives the real count.")
                           .arg(QLocale().toString(kPlanningBytesPerCell / 1000, 'f', 1)),
                       "faint"));

    const QStringList problems = mw->problems();
    if (!problems.isEmpty() && !mw->loading()) {
        auto* f = callout("red");
        f->layout()->addWidget(label(QObject::tr("Before generating"), "title"));
        f->layout()->addWidget(label(problems.join("\n")));
        v->addWidget(f);
    }
    if (mw->meshStale() && !mw->meshing())
        v->addWidget(label(QObject::tr("The geometry or settings changed since the mesh was generated: generate it again before setting up the case."), "warning"));

    // Progress, stage by stage, as the console and status bar also show it.
    auto* progress = new QWidget;
    auto* pv = new QVBoxLayout(progress);
    pv->setContentsMargins(0, 4, 0, 4);
    pv->setSpacing(4);
    auto* stageLabel = label({});
    auto* strip = new ProgressStrip;
    auto* stages = label({}, "faint");
    pv->addWidget(stageLabel);
    pv->addWidget(strip);
    pv->addWidget(stages);
    v->addWidget(progress);

    auto* outcome = callout("neutral");
    auto* outcomeTitle = label({}, "title");
    auto* outcomeText = label({});
    outcome->layout()->addWidget(outcomeTitle);
    outcome->layout()->addWidget(outcomeText);
    v->addWidget(outcome);

    auto* actions = new FlowLayout;
    auto* gen = button(mw->hasMesh() ? QObject::tr("Generate Again") : QObject::tr("Generate Mesh"), "play", "primary");
    gen->setAccessibleName(QObject::tr("Generate Mesh"));
    auto* cancel = button(QObject::tr("Cancel"), "circle-stop");
    auto* report = button(QObject::tr("Mesh Report"), "file-text");
    report->setProperty("flat", true);
    actions->addWidget(gen);
    actions->addWidget(cancel);
    actions->addWidget(report);
    v->addLayout(actions);
    QObject::connect(gen, &QPushButton::clicked, gen, [mw] { mw->generate(); });
    QObject::connect(cancel, &QPushButton::clicked, cancel, [mw] { mw->cancel(); });
    QObject::connect(report, &QPushButton::clicked, report, [c] { c.navigate("mesh/report"); });

    auto refresh = [mw, s, gen, cancel, report, progress, stageLabel, strip, stages, outcome, outcomeTitle, outcomeText, vf] {
        const bool busy = mw->meshing();
        gen->setVisible(!busy);
        gen->setEnabled(mw->problems().isEmpty() && !s->running());
        gen->setToolTip(s->running() ? QObject::tr("The GPU is solving; generate the mesh after the run.") : mw->problems().join("\n"));
        cancel->setVisible(busy);
        report->setVisible(!busy && mw->hasMesh());
        vf->setEnabled(!busy);
        const auto& p = mw->progress();
        progress->setVisible(busy);
        if (busy) {
            stageLabel->setText(QObject::tr("%1 · %2 % · %3 s").arg(p.stage.isEmpty() ? QObject::tr("Starting kutta mesh") : meshStageLabel(p.stage))
                                    .arg(int(std::round(p.fraction * 100))).arg(QLocale().toString(p.seconds, 'f', 0)));
            strip->setValue(p.fraction);
            QStringList lines;
            const QStringList all = meshStages();
            const qsizetype at = all.indexOf(p.stage);
            for (qsizetype i = 0; i < all.size(); ++i)
                lines << QStringLiteral("%1  %2").arg(i < at ? QStringLiteral("✓") : i == at ? QStringLiteral("▸") : QStringLiteral("·"), meshStageLabel(all[i]));
            stages->setText(lines.join("\n"));
        }
        using O = MeshWorkflow::Outcome;
        const O o = mw->outcome();
        outcome->setVisible(!busy && (o != O::None || mw->hasMesh()));
        if (o == O::Failed) {
            outcome->setProperty("callout", "red");
            outcomeTitle->setText(QObject::tr("Mesh generation failed"));
            outcomeText->setText(mw->failure() + (mw->advice().isEmpty() ? QString() : "\n\n" + mw->advice()));
        } else if (o == O::Cancelled) {
            outcome->setProperty("callout", "neutral");
            outcomeTitle->setText(QObject::tr("Cancelled"));
            outcomeText->setText(QObject::tr("Nothing was written; the previous mesh, if any, is unchanged."));
        } else if (const auto& r = mw->report()) {
            outcome->setProperty("callout", "neutral");
            outcomeTitle->setText(r->synthetic ? QObject::tr("Mesh ready (SYNTHETIC)") : QObject::tr("Mesh ready"));
            outcomeText->setText(QObject::tr("%1 cells, %2 faces. Read the report, then set up the case.").arg(count(double(r->cells)), count(double(r->faces))));
        }
        repolish(outcome);
    };
    QObject::connect(mw, &MeshWorkflow::meshStateChanged, gen, refresh);
    QObject::connect(mw, &MeshWorkflow::setupChanged, gen, refresh);
    QObject::connect(s, &Session::runStateChanged, gen, refresh);
    refresh();
    return page;
}

// ---- Mesh Report -------------------------------------------------------------------------

QWidget* meshReport(const InspectorContext& c)
{
    MeshWorkflow* mw = c.session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Mesh Report"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Mesh Report"), &v);
    const auto& r = mw->report();
    if (!r) {
        v->addWidget(label(mw->reportError().isEmpty() ? QObject::tr("No mesh yet. Generate it on Generate the Volume Mesh.") : mw->reportError(), "muted"));
        v->addWidget(nextButton(c, "mesh/generate", QObject::tr("Generate the Volume Mesh")), 0, Qt::AlignLeft);
        return page;
    }
    if (r->synthetic) {
        auto* row = new FlowLayout;
        row->addWidget(badge("SYNTHETIC", "modified"));
        v->addLayout(row);
        v->addWidget(label(r->provenance.isEmpty() ? QObject::tr("Fixture data, not a mesh.") : r->provenance, "warning"));
    }
    if (mw->meshStale()) v->addWidget(label(QObject::tr("The geometry or settings changed since this mesh was generated."), "warning"));

    v->addWidget(section(QObject::tr("Size")));
    auto* f = form();
    formRow(f, QObject::tr("Cells"), label(QLocale().toString(r->cells)));
    formRow(f, QObject::tr("Faces"), label(QLocale().toString(r->faces)));
    formRow(f, QObject::tr("Nodes"), label(QLocale().toString(r->nodes)));
    if (r->wallSeconds) formRow(f, QObject::tr("Time"), label(QObject::tr("%1 s").arg(QLocale().toString(*r->wallSeconds, 'f', 1))));
    if (r->peakMemoryBytes) formRow(f, QObject::tr("Peak memory"), label(bytes(double(*r->peakMemoryBytes))));
    if (r->yPlus)
        formRow(f, QObject::tr("First-cell y+"), label(r->yPlusSpeed ? QObject::tr("≈ %1 at %2 m/s").arg(QLocale().toString(*r->yPlus, 'f', 1), QLocale().toString(*r->yPlusSpeed, 'g', 4))
                                                                   : QObject::tr("≈ %1").arg(QLocale().toString(*r->yPlus, 'f', 1))));
    v->addLayout(f);
    if (!r->cellTypes.isEmpty()) {
        QStringList types;
        for (auto it = r->cellTypes.begin(); it != r->cellTypes.end(); ++it)
            if (it.value() > 0) types << QObject::tr("%1 %2 (%3 %)").arg(count(double(it.value())), it.key())
                                             .arg(QLocale().toString(r->cells ? 100.0 * double(it.value()) / double(r->cells) : 0, 'f', 1));
        v->addWidget(label(QObject::tr("Cells by type: %1").arg(types.join(", ")), "muted"));
    }

    v->addWidget(section(QObject::tr("Quality")));
    if (r->quality.empty()) v->addWidget(label(QObject::tr("The report has no quality figures."), "muted"));
    for (const auto& q : r->quality) {
        auto fmt = [](const std::optional<double>& x) { return x ? QLocale().toString(*x, 'g', 3) : QStringLiteral("—"); };
        const QString worst = q.higherIsBetter ? QObject::tr("min %1").arg(fmt(q.min)) : QObject::tr("max %1").arg(fmt(q.max));
        auto* head = label(QObject::tr("%1 · %2 · mean %3").arg(q.label, worst, fmt(q.mean)));
        v->addWidget(head);
        if (!q.counts.empty()) {
            v->addWidget(new Histogram(q));
            v->addWidget(label(QObject::tr("%1 to %2, log scale; %3").arg(QLocale().toString(q.edges.front(), 'g', 3), QLocale().toString(q.edges.back(), 'g', 3),
                                                                          q.higherIsBetter ? QObject::tr("low is poor") : QObject::tr("high is poor")),
                               "faint"));
        }
    }
    if (const auto* oq = r->metric("orthogonal_quality"); oq && oq->min && *oq->min < 0.1)
        v->addWidget(label(QObject::tr("Minimum orthogonal quality is below 0.1; Fluent's guidance is to keep it above 0.1. Kutta treats such cells with low-order face pressures."), "warning"));

    if (!r->worst.empty()) {
        v->addWidget(section(QObject::tr("Worst cells")));
        auto* t = new QTableWidget(int(r->worst.size()), 3);
        t->setHorizontalHeaderLabels({QObject::tr("Metric"), QObject::tr("Value"), QObject::tr("Near (m)")});
        t->verticalHeader()->hide();
        t->setEditTriggers(QAbstractItemView::NoEditTriggers);
        t->setShowGrid(false);
        t->setAccessibleName(QObject::tr("Worst cells"));
        for (int i = 0; i < int(r->worst.size()); ++i) {
            const auto& w = r->worst[std::size_t(i)];
            t->setItem(i, 0, new QTableWidgetItem(QObject::tr("%1%2").arg(QString(w.metric).replace('_', ' '), w.zone.isEmpty() ? QString() : " · " + w.zone)));
            t->setItem(i, 1, new QTableWidgetItem(QLocale().toString(w.value, 'g', 3)));
            t->setItem(i, 2, new QTableWidgetItem(w.position ? QStringLiteral("%1, %2, %3").arg(QLocale().toString(w.position->x(), 'f', 2), QLocale().toString(w.position->y(), 'f', 2),
                                                                                                QLocale().toString(w.position->z(), 'f', 2))
                                                             : QStringLiteral("—")));
        }
        t->resizeColumnsToContents();
        t->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
        t->setMaximumHeight(t->horizontalHeader()->height() + 4 + t->rowHeight(0) * int(r->worst.size()));
        v->addWidget(t);
    }
    v->addWidget(section(QObject::tr("Repairs")));
    v->addWidget(label(r->repairs.isEmpty() ? QObject::tr("None needed.") : r->repairs.join("\n"), r->repairs.isEmpty() ? "muted" : QString()));

    if (!r->zones.empty()) {
        v->addWidget(section(QObject::tr("Zones")));
        QStringList z;
        for (const auto& zone : r->zones) z << QObject::tr("%1 (%2, %3 faces)").arg(zone.name, zone.type, count(double(zone.faces)));
        v->addWidget(label(z.join("\n"), "muted"));
    }

    v->addWidget(section(QObject::tr("Display")));
    auto* show = new FlowLayout;
    auto* surf = button(QObject::tr("Show Surface Mesh"), "grid-3x3");
    surf->setProperty("flat", true);
    QObject::connect(surf, &QPushButton::clicked, surf, [c] {
        const QString k = c.session->generatedMeshSceneKey();
        if (!k.isEmpty()) c.graphics->show(k, 0);
    });
    show->addWidget(surf);
    auto* geo = button(QObject::tr("Show Geometry"), "box");
    geo->setProperty("flat", true);
    QObject::connect(geo, &QPushButton::clicked, geo, [c] { c.graphics->show("geometry", 0); });
    show->addWidget(geo);
    v->addLayout(show);
    v->addWidget(label(QObject::tr("Slices through the volume mesh need a section export of the mesh from the CLI (export-geometry draws boundary zones only); once the case has run, Results ▸ Planes cuts the solution."), "faint"));

    v->addWidget(section(QObject::tr("Next")));
    const bool variant = !mw->project().baseCase.isEmpty();
    v->addWidget(label(variant ? QObject::tr("This design variant's case takes every setting of %1 on the new mesh; new solids become walls in the whole-car reports. "
                                             "Unsaved changes to that case are not included: save it first.").arg(QFileInfo(mw->project().baseCase).fileName())
                               : QObject::tr("Set Up Case reads the mesh with kutta new-case: roles from the zone names, solids flagged as wheels rotating. Check them and the speed, then Create Case."),
                       "muted"));
    auto* setUp = button(variant ? QObject::tr("Set Up Variant Case") : QObject::tr("Set Up Case"), "arrow-right-left", "primary");
    setUp->setEnabled(!mw->meshStale() && !c.session->running());
    setUp->setToolTip(mw->meshStale() ? QObject::tr("Generate the mesh again first: the geometry or settings changed.") : QString());
    QObject::connect(setUp, &QPushButton::clicked, setUp, [c] { QTimer::singleShot(0, c.session, [c] { c.session->setUpCaseFromMesh(); }); });
    v->addWidget(setUp, 0, Qt::AlignLeft);
    if (!mw->project().caseFile.isEmpty()) {
        auto* open = button(QObject::tr("Open %1").arg(QFileInfo(mw->project().caseFile).fileName()), "folder-open");
        open->setProperty("flat", true);
        const QString path = mw->project().caseFile;
        QObject::connect(open, &QPushButton::clicked, open, [c, path] { if (c.openPath) c.openPath(path); });
        v->addWidget(open, 0, Qt::AlignLeft);
    }
    return page;
}

// ---- Design Variants ---------------------------------------------------------------------

QWidget* variants(const InspectorContext& c)
{
    QVBoxLayout* v;
    auto* page = new QWidget;
    {
        auto* pv = new QVBoxLayout(page);
        pv->setContentsMargins(0, 0, 0, 0);
        pv->setSpacing(0);
        pv->addWidget(header(QObject::tr("Design Variants"), QObject::tr("Runs"), page));
        auto* content = new QWidget(page);
        v = body(content);
        pv->addWidget(content);
        pv->addStretch(1);
    }
    Session* s = c.session;
    MeshWorkflow* mw = s->meshing();
    if (!mw->hasProject()) {
        v->addWidget(label(QObject::tr("Variants are copies of a geometry project with a part swapped or added (Meshing ▸ Import Geometry ▸ Duplicate as Design Variant). "
                                       "Open a case set up from geometry to compare its variants here."), "muted"));
        return page;
    }
    const auto runs = s->variantRuns();
    if (runs.size() < 2) {
        v->addWidget(label(QObject::tr("%1 has no variants yet. Duplicate it as a design variant on Import Geometry, swap a part, mesh and run it.").arg(mw->project().name()), "muted"));
        return page;
    }
    v->addWidget(label(QObject::tr("The latest finished run of each variant, with coefficients from the open case's aero settings so they compare like for like. Δ is against the base design."), "muted"));
    if (!s->hasCase()) {
        v->addWidget(label(QObject::tr("Open one of the variants' cases to compute the coefficients."), "warning"));
        return page;
    }
    // Variants across, quantities down: the deltas stay in view in a narrow panel.
    const QStringList rows{QObject::tr("Run"), QObject::tr("CL"), QObject::tr("CD"), QObject::tr("L/D"), QObject::tr("Front %"),
                           QObject::tr("ΔCL"), QObject::tr("ΔCD"), QObject::tr("Δ front % (pts)")};
    auto* t = new QTableWidget(int(rows.size()), int(runs.size()));
    t->setVerticalHeaderLabels(rows);
    QStringList heads;
    for (const auto& vr : runs) heads << vr.project + (vr.current ? QObject::tr(" (open)") : QString());
    t->setHorizontalHeaderLabels(heads);
    t->setEditTriggers(QAbstractItemView::NoEditTriggers);
    t->setSelectionBehavior(QAbstractItemView::SelectColumns);
    t->setSelectionMode(QAbstractItemView::SingleSelection);
    t->setShowGrid(false);
    t->setWordWrap(true);
    t->setAccessibleName(QObject::tr("Design variants"));
    QFont tab = t->font();
    tab.setFeature(QFont::Tag("tnum"), 1);
    t->setFont(tab);
    auto num = [](const QString& text) {
        auto* x = new QTableWidgetItem(text);
        x->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        return x;
    };
    struct Row { double cl = NAN, cd = NAN; std::optional<double> front; };
    auto results = std::make_shared<std::vector<Row>>(runs.size());
    int baseRow = 0;
    for (int i = 0; i < int(runs.size()); ++i) {
        const auto& vr = runs[std::size_t(i)];
        if (vr.base) baseRow = i;
        t->setItem(0, i, num(vr.run ? vr.run->title() + (vr.run->synthetic ? QObject::tr(" · SYNTHETIC") : QString())
                                    : vr.caseFile.isEmpty() ? QObject::tr("no case yet") : QObject::tr("no finished run")));
        for (int k = 1; k < int(rows.size()); ++k) t->setItem(k, i, num(vr.run ? QObject::tr("…") : QStringLiteral("—")));
    }
    auto fit = [t] {
        t->resizeColumnsToContents();
        t->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        t->resizeRowsToContents();
        int h = t->horizontalHeader()->height() + 6;
        for (int r = 0; r < t->rowCount(); ++r) h += t->rowHeight(r);
        t->setMinimumHeight(h);
        t->setMaximumHeight(h);
    };
    auto fill = [t, results, baseRow, num, fit] {
        const Row& b = (*results)[std::size_t(baseRow)];
        for (int i = 0; i < int(results->size()); ++i) {
            const Row& r = (*results)[std::size_t(i)];
            if (std::isnan(r.cl)) continue;
            t->setItem(1, i, num(QLocale().toString(r.cl, 'f', 3)));
            t->setItem(2, i, num(QLocale().toString(r.cd, 'f', 3)));
            t->setItem(3, i, num(r.cd != 0 ? QLocale().toString(r.cl / r.cd, 'f', 2) : QStringLiteral("—")));
            t->setItem(4, i, num(r.front ? QLocale().toString(*r.front, 'f', 1) : QStringLiteral("—")));
            auto delta = [&](double x, double base) {
                if (i == baseRow) return QObject::tr("base");
                if (std::isnan(base)) return QStringLiteral("—");
                const double pct = base != 0 ? (x - base) / std::abs(base) * 100 : 0;
                return QStringLiteral("%1%2 (%3%4 %)").arg(x - base >= 0 ? "+" : "", QLocale().toString(x - base, 'f', 3), pct >= 0 ? "+" : "", QLocale().toString(pct, 'f', 1));
            };
            t->setItem(5, i, num(delta(r.cl, b.cl)));
            t->setItem(6, i, num(delta(r.cd, b.cd)));
            t->setItem(7, i, num(i == baseRow ? QObject::tr("base")
                                 : r.front && b.front ? QStringLiteral("%1%2").arg(*r.front - *b.front >= 0 ? "+" : "", QLocale().toString(*r.front - *b.front, 'f', 1))
                                                      : QStringLiteral("—")));
        }
        fit();
    };
    for (int i = 0; i < int(runs.size()); ++i) {
        if (!runs[std::size_t(i)].run) continue;
        QPointer<QTableWidget> guard(t);
        s->aeroForRunDir(runs[std::size_t(i)].run->dir, [guard, results, i, fill, num](Result<AeroResult> a) {
            if (!guard) return;
            if (!a) {
                guard->setItem(1, i, num(QObject::tr("aero failed")));
                guard->item(1, i)->setToolTip(a.error());
                return;
            }
            if (const AeroRow* tot = a->total(); tot && tot->coefficients) {
                (*results)[std::size_t(i)].cl = tot->coefficients->cl;
                (*results)[std::size_t(i)].cd = tot->coefficients->cd;
            }
            (*results)[std::size_t(i)].front = a->balance.frontPercent;
            fill();
        });
    }
    fit();
    v->addWidget(t);
    v->addWidget(label(QObject::tr("Front % is the share of downforce on the front axle (needs the vehicle on the Aerodynamics page). Coefficients use the open case's reference values."), "faint"));

    auto* actions = new FlowLayout;
    auto* side = button(QObject::tr("Side by Side"), "columns-2", "primary");
    side->setToolTip(QObject::tr("The base design's run beside the selected variant's (click its column), cameras linked"));
    auto* openCase = button(QObject::tr("Open Variant Case"), "folder-open");
    openCase->setProperty("flat", true);
    actions->addWidget(side);
    actions->addWidget(openCase);
    v->addLayout(actions);
    t->selectColumn(baseRow == 0 && runs.size() > 1 ? 1 : 0);
    auto runsCopy = std::make_shared<std::vector<Session::VariantRun>>(runs);
    QObject::connect(side, &QPushButton::clicked, side, [c, t, runsCopy, baseRow] {
        const int row = t->currentColumn();
        if (row < 0) return;
        const auto& a = (*runsCopy)[std::size_t(baseRow)];
        const auto& b = (*runsCopy)[std::size_t(row)];
        if (!a.run || !b.run) return;
        auto key = [&](const Session::VariantRun& r) { return r.current ? "run:" + r.run->id : "ext:" + r.run->dir; };
        c.graphics->addDataset(key(a), QStringLiteral("%1 · %2").arg(a.project, a.run->title()));
        c.graphics->addDataset(key(b), QStringLiteral("%1 · %2").arg(b.project, b.run->title()));
        c.graphics->setViewLayout(GraphicsArea::Layout::SideBySide);
        c.graphics->setLinkedCameras(true);
        c.graphics->show(key(a), 0);
        c.graphics->show(key(b), 1);
    });
    QObject::connect(openCase, &QPushButton::clicked, openCase, [c, t, runsCopy] {
        const int row = t->currentColumn();
        if (row >= 0 && !(*runsCopy)[std::size_t(row)].caseFile.isEmpty() && c.openPath) c.openPath((*runsCopy)[std::size_t(row)].caseFile);
    });
    return page;
}

} // namespace pages

} // namespace kapp
