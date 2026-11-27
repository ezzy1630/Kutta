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
#include <QDesktopServices>
#include <QDir>
#include <QJsonArray>
#include <QUrl>
#include <QFileInfo>
#include <QGridLayout>
#include <QHeaderView>
#include <QLineEdit>
#include <QLocale>
#include <QFileDialog>
#include <QPainter>
#include <QPainterPath>
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

// A reported distribution as bars, one per bin (the bins need not be equal),
// on a log scale; bins on the poor side are tinted.
class Histogram : public QWidget {
public:
    Histogram(const QString& name, std::vector<MeshBin> bins, std::function<bool(const MeshBin&)> poor, QWidget* parent = nullptr)
        : QWidget(parent), m_bins(std::move(bins)), m_poor(std::move(poor))
    {
        setMinimumHeight(56);
        setAccessibleName(QObject::tr("%1 histogram").arg(name));
        QStringList tip;
        for (const auto& b : m_bins) tip << QObject::tr("%1 to %2: %3").arg(edge(b.from), edge(b.to), QLocale().toString(b.count));
        setToolTip(tip.join("\n"));
    }
    QSize sizeHint() const override { return {240, 56}; }
    static QString edge(double x) { return x >= 1e100 ? QStringLiteral("∞") : QLocale().toString(x, 'g', 3); }

protected:
    void paintEvent(QPaintEvent*) override
    {
        QPainter p(this);
        const auto& t = tokens();
        const int n = int(m_bins.size());
        if (n == 0) return;
        qint64 peak = 1;
        for (const auto& b : m_bins) peak = std::max(peak, b.count);
        const double w = double(width()) / n;
        for (int i = 0; i < n; ++i) {
            const auto& b = m_bins[std::size_t(i)];
            const double f = b.count > 0 ? std::log10(1.0 + double(b.count)) / std::log10(1.0 + double(peak)) : 0;
            const double h = std::max(f > 0 ? 2.0 : 0.0, f * (height() - 2));
            p.fillRect(QRectF(i * w + 1, height() - h, w - 2, h), m_poor(b) ? t.amber : t.accent);
        }
        p.setPen(t.border);
        p.drawLine(0, height() - 1, width(), height() - 1);
    }

private:
    std::vector<MeshBin> m_bins;
    std::function<bool(const MeshBin&)> m_poor;
};

// Shares of a whole as one bar (layer coverage, y+ bands), each with its tone.
class ShareBar : public QWidget {
public:
    struct Part { double share; QColor color; };
    ShareBar(const QString& name, std::vector<Part> parts, QWidget* parent = nullptr) : QWidget(parent), m_parts(std::move(parts))
    {
        setFixedHeight(10);
        setAccessibleName(name);
    }
    QSize sizeHint() const override { return {240, 10}; }

protected:
    void paintEvent(QPaintEvent*) override
    {
        QPainter p(this);
        p.setRenderHint(QPainter::Antialiasing);
        QPainterPath clip;
        clip.addRoundedRect(QRectF(rect()), 3, 3);
        p.setClipPath(clip);
        p.fillRect(rect(), tokens().border);
        double x = 0;
        for (const auto& part : m_parts) {
            const double w = std::clamp(part.share, 0.0, 1.0) * width();
            p.fillRect(QRectF(x, 0, w, height()), part.color);
            x += w;
        }
    }

private:
    std::vector<Part> m_parts;
};

QString percent(const std::optional<double>& f)
{
    if (!f) return QObject::tr("not reported");
    const double x = *f * 100;
    return QObject::tr("%1 %").arg(QLocale().toString(x, 'f', x > 0 && x < 0.1 ? 2 : 1));
}

QString millimetres(const std::optional<double>& m)
{
    if (!m) return QObject::tr("not reported");
    const double x = *m * 1000;
    return QObject::tr("%1 mm").arg(QLocale().toString(x, 'f', x < 0.1 ? 3 : x < 10 ? 2 : 1));
}

QString number3(const std::optional<double>& x) { return x ? QLocale().toString(*x, 'g', 3) : QObject::tr("not reported"); }
QString whole(const std::optional<qint64>& n) { return n ? QLocale().toString(*n) : QObject::tr("not reported"); }

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
    v->addWidget(label(QObject::tr("The wind tunnel's zones are named by the mesher: inlet (velocity inlet), outlet (pressure outlet), ground (moving or stationary, chosen on Wind Tunnel), side-low, side-high and top (symmetry). A solid cannot use those names, or bottom or symmetry."), "faint"));

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
    formRow(dirs, QObject::tr("Flow direction"), label(m.flow));
    formRow(dirs, QObject::tr("Up"), label(m.up));
    v->addLayout(dirs);
    v->addWidget(label(QObject::tr("kutta mesh meshes with the flow along +x and +z up: export the car from CAD in those axes."), "faint"));

    v->addWidget(section(QObject::tr("Extent (car lengths)")));
    auto* ext = form();
    auto lengths = [&](const QString& what, double value, double TunnelSetup::* member) {
        formRow(ext, what, numberField(QLocale().toString(value, 'g', 4), QObject::tr("L"), what, [s, member](double x) {
            if (!(x >= 0.01) || x > 100) return false;
            editSetup(s, [&](MeshSetup& m) { m.tunnel.*member = x; });
            return true;
        }));
    };
    lengths(QObject::tr("Upstream of the nose"), m.tunnel.upstream, &TunnelSetup::upstream);
    lengths(QObject::tr("Downstream of the tail"), m.tunnel.downstream, &TunnelSetup::downstream);
    lengths(QObject::tr("Each side (and below, without ground)"), m.tunnel.side, &TunnelSetup::side);
    lengths(QObject::tr("Above the car"), m.tunnel.top, &TunnelSetup::top);
    v->addLayout(ext);

    v->addWidget(section(QObject::tr("Ground and symmetry")));
    auto* ground = new QCheckBox(QObject::tr("Ground: the plane z = 0"));
    ground->setAccessibleName(QObject::tr("Ground plane at z = 0"));
    ground->setChecked(m.tunnel.ground);
    QObject::connect(ground, &QCheckBox::toggled, ground, [s](bool on) { editSetup(s, [&](MeshSetup& m) { m.tunnel.ground = on; }); });
    v->addWidget(ground);
    auto* motionForm = form();
    auto* motion = new QComboBox;
    motion->setAccessibleName(QObject::tr("Ground motion"));
    motion->addItem(QObject::tr("Moving (cars)"), groundMotionKey(GroundMotion::Moving));
    motion->addItem(QObject::tr("Stationary (tunnel)"), groundMotionKey(GroundMotion::Stationary));
    motion->setCurrentIndex(motion->findData(groundMotionKey(mw->project().ground)));
    motion->setEnabled(m.tunnel.ground && !mw->meshing());
    QObject::connect(motion, &QComboBox::activated, motion, [mw, motion](int i) {
        const GroundMotion g = groundMotionFrom(motion->itemData(i).toString());
        QTimer::singleShot(0, mw, [mw, g] { mw->setGroundMotion(g); });
    });
    QObject::connect(ground, &QCheckBox::toggled, motion, [motion, mw](bool on) { motion->setEnabled(on && !mw->meshing()); });
    formRow(motionForm, QObject::tr("Ground motion"), motion);
    v->addLayout(motionForm);
    v->addWidget(label(QObject::tr("A car on the road sees the ground move with the air; a wind-tunnel model such as the Ahmed body sits on a fixed floor. Either way the ground is not in the forces."), "faint"));
    auto* half = new QCheckBox(QObject::tr("Half car: symmetry plane at y = 0"));
    half->setAccessibleName(QObject::tr("Symmetry plane at y = 0"));
    half->setChecked(m.tunnel.symmetry == "y");
    QObject::connect(half, &QCheckBox::toggled, half, [s](bool on) { editSetup(s, [&](MeshSetup& m) { m.tunnel.symmetry = on ? "y" : QString(); }); });
    v->addWidget(half);
    v->addWidget(label(QObject::tr("The tyres stand on z = 0; geometry below it is refused. Half a car halves the cells; forces are then for half the car. Upper limits round out to whole base cells."), "faint"));
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
        if (!(x >= 0.1 && x <= 180)) return false;
        editSetup(s, [&](MeshSetup& m) { m.curvatureDegrees = x; });
        return true;
    }));
    formRow(sf, QObject::tr("Cells per gap"), numberField(QString::number(m.proximityCells), {}, QObject::tr("Cells per gap"), [s](double x) {
        if (x < 1 || x > 20 || x != std::floor(x)) return false;
        editSetup(s, [&](MeshSetup& m) { m.proximityCells = int(x); });
        return true;
    }));
    formRow(sf, QObject::tr("Smallest gap cell"), numberField(m.proximityMinSize > 0 ? lengthText(m.proximityMinSize, u) : QString(), un, QObject::tr("Smallest gap cell (%1)").arg(un), [s, u](double x) {
        if (x < 0) return false;
        editSetup(s, [&](MeshSetup& m) { m.proximityMinSize = x * unitScale(u); });
        return true;
    }), QObject::tr("Gaps may refine below the minimum size down to this; empty or 0: a quarter of the minimum size."));
    v->addLayout(sf);
    v->addWidget(label(QObject::tr("Surface cells shrink from the maximum towards the minimum size where the surface curves more than the normal angle per cell, and to fit the cells per gap between close surfaces (wing slots, wheel and floor)."), "faint"));

    v->addWidget(section(QObject::tr("Boundary layers")));
    auto* bl = form();
    formRow(bl, QObject::tr("Number of layers"), numberField(QString::number(m.layers), {}, QObject::tr("Number of layers"), [s](double x) {
        if (x < 0 || x > 100 || x != std::floor(x)) return false;
        editSetup(s, [&](MeshSetup& m) { m.layers = int(x); });
        return true;
    }));
    length(bl, QObject::tr("First layer height"), m.firstHeight, [](MeshSetup& m, double x) { m.firstHeight = x; });
    formRow(bl, QObject::tr("Growth rate"), numberField(QLocale().toString(m.growth, 'g', 4), {}, QObject::tr("Growth rate"), [s](double x) {
        if (!(x > 1 && x <= 3)) return false;
        editSetup(s, [&](MeshSetup& m) { m.growth = x; });
        return true;
    }));
    formRow(bl, QObject::tr("Freestream speed"), numberField(QLocale().toString(mw->project().speed, 'g', 4), QObject::tr("m/s"), QObject::tr("Freestream speed (m/s)"), [mw](double x) {
        if (!(x > 0) || x > 400) return false;
        QTimer::singleShot(0, mw, [mw, x] { mw->setSpeed(x); });
        return true;
    }), QObject::tr("For the mesher's y+ estimate, and the case's starting speed."));
    v->addLayout(bl);
    auto* yp = label({}, "faint");
    auto updateYPlus = [mw, yp] {
        if (!mw->hasProject()) return;
        const Bounds g = mw->geometryBounds();
        const double L = g.valid ? double(g.max.x() - g.min.x()) : 1.0;
        const auto& m = mw->project().setup;
        const double total = m.layers > 0 ? m.firstHeight * (std::abs(m.growth - 1) < 1e-9 ? m.layers : (std::pow(m.growth, m.layers) - 1) / (m.growth - 1)) : 0;
        // At the first cell's centroid, half the first height, as kutta mesh estimates it.
        const QString layers = m.layers > 0 ? QObject::tr("Layers %1 thick in all.").arg(millimetres(total)) : QObject::tr("No boundary layers.");
        yp->setText(QObject::tr("First-cell y+ ≈ %1 at %2 m/s over the car length (flat-plate estimate at the cell centre). %3")
                        .arg(QLocale().toString(yPlusForFirstHeight(m.firstHeight / 2, mw->project().speed, L), 'f', 0))
                        .arg(QLocale().toString(mw->project().speed, 'g', 4))
                        .arg(layers));
    };
    QObject::connect(mw, &MeshWorkflow::setupChanged, yp, updateYPlus);
    updateYPlus();
    v->addWidget(yp);
    v->addWidget(label(QObject::tr("Kutta's SST walls work across y+ (wall functions above about 30, resolved below 1). Layers grow on every solid (zones \"*\")."), "faint"));
    auto* actions = new FlowLayout;
    actions->addWidget(nextButton(c, "mesh/generate", QObject::tr("Next: Generate the Volume Mesh")));
    auto* load = button(QObject::tr("Load mesh.json…"), "file-input");
    load->setProperty("flat", true);
    load->setToolTip(QObject::tr("Take the settings of a kutta mesh setup file (lengths in metres); the geometry's units stay."));
    load->setEnabled(!mw->meshing());
    QObject::connect(load, &QPushButton::clicked, load, [mw, load] {
        const QString file = QFileDialog::getOpenFileName(load->window(), QObject::tr("Load Mesh Settings"), QFileInfo(mw->project().dir).absolutePath(), QObject::tr("Mesh settings (*.json)"));
        if (!file.isEmpty()) QTimer::singleShot(0, mw, [mw, file] { mw->loadSetupFile(file); });
    });
    actions->addWidget(load);
    v->addLayout(actions);
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
        } else if (e.solverBytes() > kTwelveGbUsableBytes) {
            text += QObject::tr(" Probably too much for a 12 GB GPU.");
            role = "warning";
        }
        if (mw->hasMesh() && mw->report())
            text += QObject::tr(" The last mesh had %1 cells.").arg(QLocale().toString(mw->report()->cells));
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
    v->addWidget(label(QObject::tr("A planning figure (±50 %) from the surface area and sizes, at about %1 kB of device memory per cell (measured on Metal). The mesh report gives the real count.")
                           .arg(QLocale().toString(kSolverBytesPerCell / 1000, 'f', 1)),
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
            outcomeTitle->setText(mw->meshSynthetic() ? QObject::tr("Mesh ready (SYNTHETIC)") : QObject::tr("Mesh ready"));
            QString text = QObject::tr("%1 cells, %2 faces. Read the report, then set up the case.").arg(count(double(r->cells)), count(double(r->faces)));
            if (!r->closedGaps().empty())
                text += QObject::tr(" A slot has an unresolved section. See Slots and gaps.");
            outcomeText->setText(text);
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

namespace {
class GapSection : public QWidget {
public:
    explicit GapSection(QWidget* parent) : QWidget(parent)
    {
        setMinimumHeight(184);
        setAccessibleName(tr("Fluid mesh cross-section through the gap"));
    }

    void setSection(const QJsonObject& data)
    {
        m_radius = data.value("radius_m").toDouble();
        m_width = data.value("width_m").toDouble();
        m_cells.clear();
        for (const auto& entry : data.value("polygons").toArray()) {
            QPolygonF polygon;
            for (const auto& pair : entry.toArray()) {
                const auto point = pair.toArray();
                if (point.size() == 2) polygon << QPointF(point[0].toDouble(), point[1].toDouble());
            }
            if (polygon.size() >= 3) m_cells.push_back(polygon);
        }
        update();
    }

protected:
    void paintEvent(QPaintEvent*) override
    {
        if (m_radius <= 0 || m_width <= 0) return;
        const auto& color = tokens();
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing);
        const QRectF area = QRectF(rect()).adjusted(1, 1, -1, -26);
        const double span = std::min(area.width(), area.height());
        const QRectF plot(area.center() - QPointF(span / 2, span / 2), QSizeF(span, span));
        painter.fillRect(plot, color.viewport);
        painter.save();
        painter.setClipRect(plot);
        const double scale = span / (2 * m_radius);
        QTransform toPixel;
        toPixel.translate(plot.center().x(), plot.center().y());
        toPixel.scale(scale, -scale);
        painter.setPen(QPen(color.borderStrong, 1));
        painter.setBrush(color.accentSoft);
        for (const QPolygonF& cell : m_cells) painter.drawPolygon(toPixel.map(cell));

        painter.setPen(QPen(color.amber, 2, Qt::DashLine));
        const auto a = toPixel.map(QPointF(0, -m_width / 2));
        const auto b = toPixel.map(QPointF(0, m_width / 2));
        painter.drawLine(a, b);
        painter.setPen(Qt::NoPen);
        painter.setBrush(color.amber);
        painter.drawEllipse(a, 3, 3);
        painter.drawEllipse(b, 3, 3);
        painter.restore();
        painter.setPen(color.borderStrong);
        painter.setBrush(Qt::NoBrush);
        painter.drawRect(plot);
        painter.setPen(color.amberText);
        painter.drawText(QRectF(0, plot.bottom() + 3, width(), 22), Qt::AlignCenter,
                         tr("Sampled gap · %1 mm").arg(QLocale().toString(m_width * 1000, 'g', 4)));
    }

private:
    double m_radius = 0, m_width = 0;
    std::vector<QPolygonF> m_cells;
};

void gapInspector(const InspectorContext& c, QVBoxLayout* layout, const std::vector<MeshGap>& gaps)
{
    auto* mw = c.session->meshing();
    std::vector<MeshGap> candidates;
    for (const auto& gap : gaps)
        if (!gap.contact && gap.cellsAcross < 3 && gap.width > 0 && gap.sample)
            candidates.push_back(gap);
    std::sort(candidates.begin(), candidates.end(), [](const MeshGap& a, const MeshGap& b) {
        return a.cellsAcross != b.cellsAcross ? a.cellsAcross < b.cellsAcross : a.width < b.width;
    });
    if (candidates.empty()) {
        if (std::any_of(gaps.begin(), gaps.end(), [](const MeshGap& gap) { return !gap.contact && gap.cellsAcross < 3; }))
            layout->addWidget(label(QObject::tr("Regenerate the mesh to record gap locations for inspection."), "muted"));
        return;
    }

    auto* gapList = new QComboBox;
    gapList->setAccessibleName(QObject::tr("Under-resolved gap"));
    gapList->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    gapList->setMinimumContentsLength(12);
    for (const auto& gap : candidates)
        gapList->addItem(gap.cellsAcross == 1 ? QObject::tr("%1 / %2 · 1 cell").arg(gap.a, gap.b)
                                            : QObject::tr("%1 / %2 · %3 cells").arg(gap.a, gap.b).arg(gap.cellsAcross));
    layout->addWidget(gapList);
    auto* actions = new FlowLayout;
    auto* inspect = button(QObject::tr("Inspect cross-section"), "scan-line");
    auto* refine = button(QObject::tr("Add local refinement"), "plus");
    refine->setEnabled(!mw->meshing() && mw->meshGeometryCurrent());
    refine->setToolTip(QObject::tr("Add an editable box around the sampled gap, with a target of four cells across its width."));
    actions->addWidget(inspect);
    actions->addWidget(refine);
    layout->addLayout(actions);
    auto* sectionPlot = new GapSection(gapList->parentWidget());
    auto* status = label({}, "muted");
    sectionPlot->hide();
    layout->addWidget(sectionPlot);
    layout->addWidget(status);
    status->hide();
    const auto refresh = [sectionPlot, status, inspect] {
        sectionPlot->hide();
        status->hide();
        inspect->setEnabled(true);
    };
    QObject::connect(gapList, &QComboBox::currentIndexChanged, sectionPlot, refresh);
    QObject::connect(gapList, &QComboBox::activated, inspect, &QPushButton::click);
    QObject::connect(inspect, &QPushButton::clicked, sectionPlot, [=] {
        const int selected = gapList->currentIndex();
        const auto& gap = candidates.at(std::size_t(selected));
        sectionPlot->hide();
        status->setText(QObject::tr("Cutting the volume mesh…"));
        status->show();
        inspect->setEnabled(false);
        gapList->setEnabled(false);
        QPointer<GapSection> guard(sectionPlot);
        c.session->cli()->meshSection(mw->meshFile(), gap.sample->first, gap.sample->second,
            [=](Result<QJsonObject> result) {
                if (!guard) return;
                gapList->setEnabled(true);
                inspect->setEnabled(true);
                if (!result) {
                    status->setText(result.error());
                    return;
                }
                sectionPlot->setSection(result.value());
                sectionPlot->show();
                const int skipped = result->value("skipped_cells").toInt();
                status->setText(skipped ? QObject::tr("Fluid-cell outlines. %n invalid cell cut(s) omitted.", nullptr, skipped)
                                        : QObject::tr("Fluid-cell outlines from the generated mesh."));
            });
    });
    QObject::connect(refine, &QPushButton::clicked, refine, [=] {
        if (mw->meshing() || !mw->meshGeometryCurrent()) return;
        const auto& gap = candidates.at(std::size_t(gapList->currentIndex()));
        auto setup = mw->project().setup;
        const auto [a, b] = *gap.sample;
        const float pad = float(2 * gap.width);
        QVector3D lo, hi;
        for (int axis = 0; axis < 3; ++axis) {
            lo[axis] = std::min(a[axis], b[axis]) - pad;
            hi[axis] = std::max(a[axis], b[axis]) + pad;
        }
        if (setup.tunnel.ground) lo.setZ(std::max(0.0f, lo.z()));
        if (setup.tunnel.symmetry == "y") lo.setY(std::max(0.0f, lo.y()));
        const QString base = QObject::tr("Gap %1 / %2").arg(gap.a, gap.b);
        QString name = base;
        for (int n = 2; std::any_of(setup.boxes.begin(), setup.boxes.end(),
                    [&](const RefinementBox& box) { return box.name == name; }); ++n)
            name = QStringLiteral("%1 (%2)").arg(base).arg(n);
        const double size = std::max(1e-7, gap.width / 4);
        setup.boxes.push_back({name, lo, hi, std::min(size, setup.baseSize)});
        mw->setSetup(setup);
        c.graphics->show("geometry", 0);
        c.navigate("mesh/refinement");
    });
}
} // namespace

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
    if (mw->meshSynthetic()) {
        auto* row = new FlowLayout;
        row->addWidget(badge("SYNTHETIC", "modified"));
        v->addLayout(row);
        v->addWidget(label(QObject::tr("The fixture CLI made this report for app development: no mesh was generated and every figure is invented."), "warning"));
    }
    if (mw->meshStale()) v->addWidget(label(QObject::tr("The geometry or settings changed since this mesh was generated."), "warning"));
    const QString notReported = QObject::tr("not reported");

    // ---- Size and what it costs to solve
    v->addWidget(section(QObject::tr("Size")));
    auto* f = form();
    formRow(f, QObject::tr("Cells"), label(QLocale().toString(r->cells)));
    formRow(f, QObject::tr("Faces"), label(QLocale().toString(r->faces)));
    formRow(f, QObject::tr("Nodes"), label(QLocale().toString(r->nodes)));
    formRow(f, QObject::tr("Mesh time"), label(r->wallSeconds ? QObject::tr("%1 s").arg(QLocale().toString(*r->wallSeconds, 'f', 1)) : notReported));
    formRow(f, QObject::tr("Mesher peak memory"), label(r->peakMemoryBytes ? bytes(double(*r->peakMemoryBytes)) : notReported));
    v->addLayout(f);
    {
        const auto& t = r->cellTypes;
        QStringList types;
        auto add = [&](const std::optional<qint64>& n, const QString& what) {
            if (n) types << QObject::tr("%1 %2").arg(QLocale().toString(*n), what);
        };
        add(t.hex, QObject::tr("hex"));
        add(t.prism, QObject::tr("prism"));
        add(t.poly, QObject::tr("polyhedra"));
        QString text = types.isEmpty() ? QObject::tr("Cells by type: not reported.") : QObject::tr("Cells by type: %1.").arg(types.join(", "));
        if (t.hexFractionOfCore) text += QObject::tr(" Hex share of the core (cells outside the layers): %1.").arg(percent(t.hexFractionOfCore));
        if (t.hexWithoutHangingNodes && t.hex) text += QObject::tr(" Hex counts octree cubes with 2:1 hanging faces; %1 have none.").arg(QLocale().toString(*t.hexWithoutHangingNodes));
        v->addWidget(label(text, "muted"));
    }
    {
        const double need = solverBytesEstimate(r->cells);
        auto* box = callout(need > kTwelveGbUsableBytes ? "amber" : "neutral");
        box->layout()->addWidget(label(QObject::tr("Solver memory ≈ %1 (estimate)").arg(bytes(need)), "title"));
        box->layout()->addWidget(label(QObject::tr("At about %1 kB per cell, measured on Metal for steady SST on kutta mesh meshes; not a measurement of this mesh.")
                                           .arg(QLocale().toString(kSolverBytesPerCell / 1000, 'f', 1)), "muted"));
        if (need > kTwelveGbUsableBytes)
            box->layout()->addWidget(label(QObject::tr("Probably will not fit a 12 GB GPU: about %1 cells is the limit at this rate once the CUDA context is allowed for. Coarsen the sizes, use fewer layers or a half car.")
                                               .arg(count(kTwelveGbUsableBytes / kSolverBytesPerCell)), "warning"));
        else
            box->layout()->addWidget(label(QObject::tr("Fits a 12 GB GPU at this rate (up to about %1 cells).").arg(count(kTwelveGbUsableBytes / kSolverBytesPerCell)), "muted"));
        if (const auto dev = selectedDevice(c.session); dev && dev->memoryBytes > 0) {
            const bool fits = need < 0.85 * double(dev->memoryBytes);
            box->layout()->addWidget(label(fits ? QObject::tr("%1 (%2): fits.").arg(dev->name, bytes(double(dev->memoryBytes)))
                                                : QObject::tr("%1 (%2): more than it holds comfortably.").arg(dev->name, bytes(double(dev->memoryBytes))),
                                           fits ? "muted" : "warning"));
        }
        v->addWidget(box);
    }

    // ---- Slots and gaps the octree resolved, or closed
    v->addWidget(section(QObject::tr("Slots and gaps")));
    if (!r->gaps) {
        v->addWidget(label(QObject::tr("Gaps: not reported."), "muted"));
    } else {
        std::vector<MeshGap> closed, contact, open;
        for (const auto& g : *r->gaps) {
            if (g.contact) contact.push_back(g);
            else if (g.closed()) closed.push_back(g);
            else open.push_back(g);
        }
        auto pair = [](const MeshGap& g) { return QObject::tr("%1 – %2").arg(g.a, g.b); };
        if (!closed.empty()) {
            auto* box = callout("amber");
            box->layout()->addWidget(label(closed.size() == 1 ? QObject::tr("One slot has an unresolved section") : QObject::tr("%1 slots have unresolved sections").arg(closed.size()), "title"));
            QStringList lines;
            for (const auto& g : closed) lines << QObject::tr("%1, %2 wide").arg(pair(g), millimetres(g.width));
            box->layout()->addWidget(label(lines.join("\n")));
            box->layout()->addWidget(label(QObject::tr("The mesh has zero fluid cells across the narrowest sampled section. Inspect the sampled location and refine the gap before comparing loads."), "muted"));
            v->addWidget(box);
        } else {
            v->addWidget(label(r->gaps->empty() ? QObject::tr("No narrow gaps between parts.") : QObject::tr("Every sampled open gap contains fluid cells."), "muted"));
        }
        std::sort(open.begin(), open.end(), [](const MeshGap& x, const MeshGap& y) { return x.cellsAcross != y.cellsAcross ? x.cellsAcross < y.cellsAcross : x.width < y.width; });
        if (!open.empty()) {
            QStringList lines;
            for (std::size_t i = 0; i < std::min<std::size_t>(open.size(), 5); ++i)
                lines << (open[i].cellsAcross == 1 ? QObject::tr("%1: 1 cell across %2").arg(pair(open[i]), millimetres(open[i].width))
                                                   : QObject::tr("%1: %2 cells across %3").arg(pair(open[i])).arg(open[i].cellsAcross).arg(millimetres(open[i].width)));
            v->addWidget(label(QObject::tr("Narrowest open gaps (of %1):\n%2").arg(open.size()).arg(lines.join("\n")), open.front().cellsAcross < 2 ? "warning" : "muted"));
        }
        if (!mw->meshSynthetic()) gapInspector(c, v, *r->gaps);
        if (!contact.empty()) {
            QStringList names;
            for (const auto& g : contact) names << pair(g);
            v->addWidget(label(QObject::tr("Touching, no gap to resolve: %1.").arg(names.join(", ")), "faint"));
        }
    }

    // ---- Prism layers
    v->addWidget(section(QObject::tr("Boundary layers")));
    if (!r->layers) {
        v->addWidget(label(QObject::tr("Layers: not reported."), "muted"));
    } else if (r->layers->requested && *r->layers->requested == 0) {
        v->addWidget(label(QObject::tr("None: no layers were requested (count 0)."), "muted"));
    } else {
        const auto& l = *r->layers;
        auto* lf = form();
        formRow(lf, QObject::tr("Layers"), label(QObject::tr("%1 generated of %2 requested").arg(whole(l.generated), whole(l.requested))));
        formRow(lf, QObject::tr("Prism cells"), label(whole(l.cells)));
        formRow(lf, QObject::tr("First height achieved"), label(l.firstHeightMean ? QObject::tr("mean %1, min %2").arg(millimetres(l.firstHeightMean), millimetres(l.firstHeightMin)) : notReported));
        formRow(lf, QObject::tr("Terminations"), label(l.splitPoints || l.liftedPoints ? QObject::tr("%1 split, %2 lifted (a step)").arg(whole(l.splitPoints), whole(l.liftedPoints)) : notReported));
        v->addLayout(lf);
        if (l.coverageFull || l.coverageThinned || l.coverageNone) {
            const auto& t = tokens();
            v->addWidget(new ShareBar(QObject::tr("Layer coverage"), {{l.coverageFull.value_or(0), t.accent}, {l.coverageThinned.value_or(0), t.amber}, {l.coverageNone.value_or(0), t.red}}));
            v->addWidget(label(QObject::tr("Coverage of the layered walls%1: full %2, thinned %3, none %4.")
                                   .arg(l.wallArea ? QObject::tr(" (%1 m²)").arg(QLocale().toString(*l.wallArea, 'f', 2)) : QString(),
                                        percent(l.coverageFull), percent(l.coverageThinned), percent(l.coverageNone)),
                               l.coverageNone.value_or(0) > 0.1 ? "warning" : "muted"));
            if (l.coverageNone.value_or(0) > 0)
                v->addWidget(label(QObject::tr("Where a stack would fall below orthogonal quality 0.15 it collapses: those walls have no layers, and their first cell is a core cell."), "faint"));
        } else {
            v->addWidget(label(QObject::tr("Coverage: not reported."), "muted"));
        }
    }

    // ---- Snap: how far the wall stays from the STL
    v->addWidget(section(QObject::tr("Distance from the geometry")));
    if (!r->snap || !r->snap->distanceMax) {
        v->addWidget(label(QObject::tr("Snap error: not reported."), "muted"));
    } else {
        const auto& sn = *r->snap;
        auto* sf = form();
        formRow(sf, QObject::tr("Wall points to the STL"), label(QObject::tr("p99 %1, max %2").arg(millimetres(sn.distanceP99), millimetres(sn.distanceMax))));
        formRow(sf, QObject::tr("Mean"), label(millimetres(sn.distanceMean)));
        if (sn.relativeMax) formRow(sf, QObject::tr("Largest, of the local cell size"), label(percent(sn.relativeMax)));
        v->addLayout(sf);
        v->addWidget(label(QObject::tr("Where snapping fully would push a cell below orthogonal quality 0.15, the point stays off the surface: this is the geometric error left, at sharp trailing edges and concave junctions."), "faint"));
    }

    // ---- y+
    v->addWidget(section(QObject::tr("First-cell y+ (estimate)")));
    if (!r->yPlus) {
        v->addWidget(label(QObject::tr("y+: not reported."), "muted"));
    } else {
        const auto& y = *r->yPlus;
        auto* yf = form();
        formRow(yf, QObject::tr("Median (p50)"), label(y.speed ? QObject::tr("%1 at %2 m/s").arg(number3(y.p50), QLocale().toString(*y.speed, 'g', 4)) : number3(y.p50)));
        formRow(yf, QObject::tr("p05 to p95"), label(QObject::tr("%1 to %2").arg(number3(y.p05), number3(y.p95))));
        formRow(yf, QObject::tr("Min to max"), label(QObject::tr("%1 to %2").arg(number3(y.min), number3(y.max))));
        v->addLayout(yf);
        if (y.below5 || y.from5to30 || y.from30to300 || y.above300) {
            const auto& t = tokens();
            v->addWidget(new ShareBar(QObject::tr("y+ bands"), {{y.below5.value_or(0), t.link}, {y.from5to30.value_or(0), t.amber}, {y.from30to300.value_or(0), t.accent}, {y.above300.value_or(0), t.red}}));
            v->addWidget(label(QObject::tr("Wall area: below 5: %1 · 5–30: %2 · 30–300: %3 · above 300: %4")
                                   .arg(percent(y.below5), percent(y.from5to30), percent(y.from30to300), percent(y.above300)), "muted"));
        } else {
            v->addWidget(label(QObject::tr("y+ bands: not reported."), "muted"));
        }
        v->addWidget(label(QObject::tr("A flat-plate estimate from the first cells' distances, not solved y+. 5–30 is the buffer layer, where neither wall-resolved nor wall-function modelling is accurate; Kutta's SST walls blend across it."), "faint"));
    }

    // ---- Quality
    v->addWidget(section(QObject::tr("Quality")));
    if (!r->quality) {
        v->addWidget(label(QObject::tr("Quality: not reported."), "muted"));
    } else {
        const auto& q = *r->quality;
        auto* qf = form();
        formRow(qf, QObject::tr("Orthogonal quality"), label(QObject::tr("min %1, mean %2").arg(number3(q.orthogonalMin), number3(q.orthogonalMean))));
        formRow(qf, QObject::tr("Cells below 0.15 / 0.05"), label(QObject::tr("%1 / %2").arg(whole(q.cellsBelow015), whole(q.cellsBelow005))));
        formRow(qf, QObject::tr("Worst cell at"), label(q.worstCell ? QStringLiteral("%1, %2, %3 m").arg(QLocale().toString(q.worstCell->x(), 'f', 3), QLocale().toString(q.worstCell->y(), 'f', 3),
                                                                                                         QLocale().toString(q.worstCell->z(), 'f', 3))
                                                                   : notReported));
        formRow(qf, QObject::tr("Maximum skewness"), label(number3(q.skewnessMax)));
        formRow(qf, QObject::tr("Minimum cell volume"), label(q.volumeMin ? QObject::tr("%1 m³").arg(QLocale().toString(*q.volumeMin, 'g', 3)) : notReported));
        formRow(qf, QObject::tr("Negative volumes"), label(whole(q.negativeCells)));
        formRow(qf, QObject::tr("Nonconvex cells"), label(whole(q.nonconvexCells)));
        v->addLayout(qf);
        if (!q.orthogonalHistogram.empty()) {
            v->addWidget(new Histogram(QObject::tr("Orthogonal quality"), q.orthogonalHistogram, [](const MeshBin& b) { return b.to <= 0.2 + 1e-9; }));
            v->addWidget(label(QObject::tr("Orthogonal quality, %1 to %2 by the report's bins, log scale; below 0.2 tinted.")
                                   .arg(Histogram::edge(q.orthogonalHistogram.front().from), Histogram::edge(q.orthogonalHistogram.back().to)), "faint"));
        }
        if (q.orthogonalMin && *q.orthogonalMin < 0.1)
            v->addWidget(label(QObject::tr("Minimum orthogonal quality is below 0.1; Fluent's guidance is to keep it above 0.1. Kutta treats such cells with low-order face pressures."), "warning"));
        auto ratio = [&](const std::optional<MeshDistribution>& d, const QString& where) {
            if (!d) return v->addWidget(label(QObject::tr("Volume ratio, %1: not reported.").arg(where), "muted"));
            qint64 over20 = 0;
            for (const auto& b : d->histogram) if (b.from >= 20) over20 += b.count;
            v->addWidget(label(QObject::tr("Adjacent volume ratio, %1%2: p50 %3, p99 %4, p999 %5, max %6%7")
                                   .arg(where, d->faces ? QObject::tr(" (%1 faces)").arg(QLocale().toString(*d->faces)) : QString(),
                                        number3(d->p50), number3(d->p99), number3(d->p999), number3(d->max),
                                        d->histogram.empty() ? QString() : QObject::tr("; %1 above 20").arg(QLocale().toString(over20)))));
            if (!d->histogram.empty()) v->addWidget(new Histogram(QObject::tr("Volume ratio, %1").arg(where), d->histogram, [](const MeshBin& b) { return b.from >= 20; }));
        };
        ratio(q.volumeRatioCore, QObject::tr("core"));
        ratio(q.volumeRatioLayers, QObject::tr("layers"));
        v->addWidget(label(QObject::tr("Ratio of neighbouring cells' volumes across each face; above 20 tinted. The largest sit beside 2:1 octree steps where snapping squeezed a fine cell."), "faint"));
    }

    // ---- Surface repairs, missing parts, the mesher's own warnings
    v->addWidget(section(QObject::tr("Repairs")));
    {
        static const QList<std::pair<const char*, const char*>> known{
            {"open_edges", QT_TR_NOOP("open edges")}, {"nonmanifold_edges", QT_TR_NOOP("non-manifold edges")},
            {"orientation_conflicts", QT_TR_NOOP("orientation conflicts")}, {"reoriented_triangles", QT_TR_NOOP("triangles reoriented")},
            {"removed_degenerate_triangles", QT_TR_NOOP("degenerate triangles removed")}, {"removed_duplicate_triangles", QT_TR_NOOP("duplicate triangles removed")},
            {"welded_vertex_occurrences", QT_TR_NOOP("vertices welded")}};
        QStringList done;
        for (const auto& [key, what] : known)
            if (r->repairs.value(key) > 0) done << QObject::tr("%1 %2").arg(QLocale().toString(r->repairs.value(key)), QObject::tr(what));
        QString text = r->repairs.isEmpty() ? QObject::tr("Repairs: not reported.") : done.isEmpty() ? QObject::tr("None needed.") : done.join(", ") + ".";
        if (r->repairs.contains("components")) text += (r->repairs.value("components") == 1 ? QObject::tr(" 1 closed shell.") : QObject::tr(" %1 closed shells.").arg(r->repairs.value("components")));
        const bool open = r->repairs.value("open_edges") > 0 || r->repairs.value("nonmanifold_edges") > 0;
        v->addWidget(label(text, open ? "warning" : "muted"));
        if (open) v->addWidget(label(QObject::tr("Open or non-manifold edges mean the STL is not closed; kutta mesh does not fill holes. Close the surfaces in CAD."), "warning"));
    }
    if (!r->missingSurfaceZones) v->addWidget(label(QObject::tr("Missing solids: not reported."), "muted"));
    else if (!r->missingSurfaceZones->isEmpty())
        v->addWidget(label(QObject::tr("These solids have no faces in the mesh, too small for its sizes: %1. Lower the minimum size or add a refinement region.").arg(r->missingSurfaceZones->join(", ")), "warning"));
    if (!r->warnings.isEmpty()) {
        v->addWidget(section(QObject::tr("Mesher notes")));
        v->addWidget(label(r->warnings.join("\n"), "muted"));
    }

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
    v->addWidget(nextButton(c, "mesh/sensitivity", QObject::tr("Mesh Sensitivity")), 0, Qt::AlignLeft);

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

QWidget* meshSensitivity(const InspectorContext& c)
{
    Session* session = c.session;
    auto* mw = session->meshing();
    if (!mw->hasProject()) return noProject(c, QObject::tr("Mesh Sensitivity"));
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Mesh Sensitivity"), &v);
    auto* study = mw->sensitivity();
    const auto result = study->latest(mw->project().dir);
    const auto grids = result.value("grids").toArray();
    const QString state = result.value("state").toString();
    if (!result.isEmpty()) {
        v->addWidget(label(result.value("message").toString(), state == "failed" ? "error" : "muted"));
        auto* table = new QGridLayout;
        table->setHorizontalSpacing(8);
        table->setVerticalSpacing(5);
        const QStringList headings{QObject::tr("Coarse"), QObject::tr("Medium"), QObject::tr("Fine")};
        const QList<std::pair<QString, QString>> metrics{
            {QObject::tr("Cells"), "cells"}, {QObject::tr("Drag (N)"), "drag_n"},
            {QObject::tr("Downforce (N)"), "downforce_n"}, {QObject::tr("Front (%)"), "front_percent"},
            {QObject::tr("Gaps <3 cells"), "under_resolved_gaps"},
            {QObject::tr("Outcome"), "verdict"}};
        for (int i = 0; i < headings.size(); ++i) {
            table->addWidget(label(headings[i], "muted"), 0, i + 1);
            table->setColumnStretch(i + 1, 1);
        }
        const auto verdictText = [](const QString& key) {
            if (key == "converged") return QObject::tr("Converged");
            if (key == "statistically_steady") return QObject::tr("Cycling");
            if (key == "forces_settled_residuals_waived") return QObject::tr("Residuals waived");
            if (key == "iteration_limit") return QObject::tr("Iteration limit");
            return QStringLiteral("—");
        };
        for (int row = 0; row < metrics.size(); ++row) {
            const auto& [name, key] = metrics[row];
            table->addWidget(label(name, "muted"), row + 1, 0);
            for (int i = 0; i < headings.size(); ++i) {
                const auto entry = i < grids.size() ? grids[i].toObject() : QJsonObject{};
                const auto value = entry.value(key);
                QString text = QStringLiteral("—");
                if (key == "verdict") text = verdictText(value.toString());
                else if (value.isDouble())
                    text = key == "cells" ? QLocale().toString(value.toInteger()) : QLocale().toString(value.toDouble(), 'g', 5);
                table->addWidget(label(text), row + 1, i + 1);
            }
        }
        v->addLayout(table);
        if (state == "complete" && grids.size() == 3) {
            const auto medium = grids[1].toObject(), fine = grids[2].toObject();
            const auto change = [&](const QString& key) {
                const auto a = medium.value(key), b = fine.value(key);
                return a.isDouble() && b.isDouble() ? QLocale().toString(b.toDouble() - a.toDouble(), 'g', 4) : QStringLiteral("—");
            };
            v->addWidget(label(QObject::tr("Final loads, fine − medium: %1 N drag, %2 N downforce, %3 points front balance.")
                                   .arg(change("drag_n"), change("downforce_n"), change("front_percent")), "muted"));
            QStringList limitations;
            for (const auto& item : grids) {
                const auto grid = item.toObject();
                if (grid.value("verdict") != QStringLiteral("converged")) limitations << QObject::tr("some grids stopped without converging");
                if (grid.value("under_resolved_gaps").toInt() > 0) limitations << QObject::tr("some sampled gaps have fewer than three cells");
                if (grid.value("approximate").toBool()) limitations << QObject::tr("the solver used approved approximations");
            }
            limitations.removeDuplicates();
            if (!limitations.isEmpty()) v->addWidget(label(QObject::tr("Diagnostic comparison: %1.").arg(limitations.join(" · ")), "warning"));
            v->addWidget(label(QObject::tr("A numerical uncertainty interval needs verified asymptotic refinement and a bound on iterative load error. This study reports the observed grid changes."), "muted"));
            if (fine.value("front_percent").isNull())
                v->addWidget(label(QObject::tr("Front balance needs a vehicle reference and positive net downforce."), "muted"));
        }
        const QString directory = result.value("directory").toString();
        if (!directory.isEmpty()) {
            auto* open = button(QObject::tr("Open study folder"), "folder-open");
            QObject::connect(open, &QPushButton::clicked, page, [directory] { QDesktopServices::openUrl(QUrl::fromLocalFile(directory)); });
            v->addWidget(open, 0, Qt::AlignLeft);
        }
    }

    v->addWidget(section(QObject::tr("Grid setup")));
    v->addWidget(label(QObject::tr("The project sizing defines the coarse grid. Medium and fine halve the cell sizes and first-layer height in turn. "
                                  "All three share the same geometry, fixed tunnel, case setup and stopping criteria."), "muted"));
    v->addWidget(label(QObject::tr("The fine grid can require much more memory and disk space. Each grid keeps the project’s cell limit."), "muted"));
    if (mw->hasMesh() && mw->meshStale())
        v->addWidget(label(QObject::tr("The study uses your current geometry and sizing, including edits made since the project mesh was generated."), "muted"));

    auto* error = label({}, "error");
    error->hide();
    v->addWidget(error);
    if (study->running()) {
        auto* stop = button(QObject::tr("Stop study"), "circle-stop", "danger");
        QObject::connect(stop, &QPushButton::clicked, study, &MeshSensitivity::stop);
        v->addWidget(stop, 0, Qt::AlignLeft);
    } else {
        auto* run = button(QObject::tr("Generate and solve three grids"), "play", "primary");
        QObject::connect(run, &QPushButton::clicked, page, [session, error] {
            QString why;
            if (!session->startMeshSensitivity(&why)) {
                error->setText(why);
                error->show();
            }
        });
        v->addWidget(run, 0, Qt::AlignLeft);
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
