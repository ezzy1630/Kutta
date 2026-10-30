// Inspector pages for Results, Runs and All Case Settings.

#include "core/ChartScale.h"
#include "ui/ChartsPanel.h"
#include "ui/GraphicsArea.h"
#include "ui/Inspector.h"
#include "ui/PageKit.h"
#include "ui/SettingWidgets.h"

#include <QAbstractTableModel>
#include <QButtonGroup>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QPointer>
#include <QSettings>
#include <QFormLayout>
#include <QHeaderView>
#include <QInputDialog>
#include <QJsonArray>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QRadioButton>
#include <QSortFilterProxyModel>
#include <QSpinBox>
#include <QTableView>
#include <QTableWidget>
#include <QTimer>
#include <QUuid>

#include <cmath>
#include <numeric>

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
    v->addWidget(content);
    v->addStretch(1);
    return page;
}

// The run a results page acts on: the active viewport's run, else the latest.
QString currentRunId(const InspectorContext& c)
{
    if (auto* p = c.graphics->active()) {
        const QString k = p->datasetKey();
        if (k.startsWith("run:")) return k.mid(4);
    }
    if (const RunRecord* r = c.session->history().latestFinished()) return r->id;
    return {};
}

QLabel* targetLine(const InspectorContext& c)
{
    auto* p = c.graphics->active();
    const QString text = p && !p->datasetKey().isEmpty()
                             ? QObject::tr("Viewport %1 · %2").arg(p->index() + 1).arg(p->datasetLabel())
                             : QObject::tr("The active viewport has no dataset. Run the calculation or open a run from Run History.");
    return label(text, "muted");
}

// Applies an edit to the active viewport's display.
void editDisplay(const InspectorContext& c, const std::function<void(ViewDisplay&)>& fn)
{
    auto* p = c.graphics->active();
    if (!p) return;
    ViewDisplay d = p->display();
    fn(d);
    p->setDisplay(d);
}

QLineEdit* numberEdit(double v, const QString& name)
{
    auto* e = new QLineEdit(QString::number(v, 'g', 6));
    e->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    e->setAccessibleName(name);
    return e;
}

std::optional<double> numberOf(QLineEdit* e)
{
    bool ok = false;
    double v = QLocale::c().toDouble(e->text().trimmed(), &ok);
    if (!ok) v = QLocale().toDouble(e->text().trimmed(), &ok);
    return ok && std::isfinite(v) ? std::optional<double>(v) : std::nullopt;
}

struct Vec3Edit { QWidget* widget; std::array<QLineEdit*, 3> e; std::optional<QVector3D> value() const; };
std::optional<QVector3D> Vec3Edit::value() const
{
    QVector3D v;
    for (int i = 0; i < 3; ++i) {
        auto n = numberOf(e[std::size_t(i)]);
        if (!n) return std::nullopt;
        v[i] = float(*n);
    }
    return v;
}
Vec3Edit vec3Edit(QVector3D v, const QString& name)
{
    Vec3Edit out{new QWidget, {}};
    auto* h = new QHBoxLayout(out.widget);
    h->setContentsMargins(0, 0, 0, 0);
    h->setSpacing(4);
    for (int i = 0; i < 3; ++i) {
        out.e[std::size_t(i)] = numberEdit(v[i], QStringLiteral("%1 %2").arg(name, QString("XYZ"[i])));
        out.e[std::size_t(i)]->setPlaceholderText(QString("xyz"[i]));
        h->addWidget(out.e[std::size_t(i)]);
    }
    return out;
}

QFormLayout* formLayout()
{
    auto* f = new QFormLayout;
    f->setLabelAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    f->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    f->setHorizontalSpacing(10);
    f->setVerticalSpacing(6);
    return f;
}

// Checklist of the dataset's surfaces; edits the active viewport's visibility.
QListWidget* surfaceList(const InspectorContext& c)
{
    auto* list = new QListWidget;
    list->setAccessibleName(QObject::tr("Surfaces shown"));
    auto* p = c.graphics->active();
    const auto s = p ? p->scene() : nullptr;
    if (!s) {
        list->setEnabled(false);
        return list;
    }
    const QStringList vis = p->display().visible;
    auto addGroup = [&](const QString& title, auto pred) {
        bool header = false;
        for (const auto& it : s->items) {
            if (!pred(it)) continue;
            if (!header) {
                auto* h = new QListWidgetItem(title.toUpper());
                h->setFlags(Qt::NoItemFlags);
                QFont f = list->font();
                f.setPointSizeF(f.pointSizeF() * 0.85);
                f.setBold(true);
                h->setFont(f);
                h->setForeground(tokens().muted);
                list->addItem(h);
                header = true;
            }
            auto* li = new QListWidgetItem(it.label + (it.zoneType.isEmpty() || it.kind != ItemKind::Zone ? QString() : "  ·  " + it.zoneType));
            li->setData(Qt::UserRole, it.id);
            li->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable | Qt::ItemIsSelectable);
            li->setCheckState(vis.contains(it.id) ? Qt::Checked : Qt::Unchecked);
            list->addItem(li);
        }
    };
    addGroup(QObject::tr("Walls"), [](const SceneItem& i) { return i.kind == ItemKind::Zone && i.zoneType == "wall"; });
    addGroup(QObject::tr("Boundaries"), [](const SceneItem& i) { return i.kind == ItemKind::Zone && i.zoneType != "wall"; });
    addGroup(QObject::tr("Planes"), [](const SceneItem& i) { return i.kind == ItemKind::Slice; });
    addGroup(QObject::tr("Iso-surfaces and pathlines"), [](const SceneItem& i) { return i.kind == ItemKind::IsoSurface || i.kind == ItemKind::Streamlines; });
    list->setMinimumHeight(std::min(320, 12 + list->count() * 24));
    QObject::connect(list, &QListWidget::itemChanged, list, [c, list] {
        QStringList ids;
        for (int i = 0; i < list->count(); ++i)
            if (list->item(i)->checkState() == Qt::Checked) ids << list->item(i)->data(Qt::UserRole).toString();
        editDisplay(c, [&](ViewDisplay& d) { d.visible = ids; d.visibleInitialized = true; });
    });
    return list;
}

QIcon colormapIcon(ColormapId id)
{
    QImage img = colormap(id).image().scaled(48, 12);
    return QIcon(QPixmap::fromImage(img));
}

QComboBox* fieldCombo(const SceneData* s, const FieldRef& current, bool scalarOnly = false)
{
    auto* cb = new QComboBox;
    cb->setAccessibleName(QObject::tr("Field"));
    cb->setMaxVisibleItems(24);
    if (!s) return cb;
    std::optional<FieldRole> last;
    for (const auto& f : s->availableFields()) {
        if (scalarOnly && f.component >= 0) continue;
        if (last && *last != f.role) cb->insertSeparator(cb->count());
        last = f.role;
        cb->addItem(displayLabel(f), f.key());
    }
    cb->setCurrentIndex(std::max(0, cb->findData(current.key())));
    return cb;
}

void reportError(QLabel* l, const QString& message)
{
    l->setText(message);
    l->setVisible(!message.isEmpty());
}

QString newPostId(const QString& kind) { return kind + "-" + QUuid::createUuid().toString(QUuid::Id128).left(8); }

} // namespace

// ---- Graphics ------------------------------------------------------------------------

QWidget* contours(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Contours"), QObject::tr("Results ▸ Graphics"), &v);
    v->addWidget(targetLine(c));
    auto* p = c.graphics->active();
    const auto s = p ? p->scene() : nullptr;
    if (!s) return page;
    const ViewDisplay d = p->display();

    auto* f = formLayout();
    auto* field = fieldCombo(s.get(), d.field);
    f->addRow(QObject::tr("Field"), field);
    QObject::connect(field, &QComboBox::activated, field, [c, field](int i) {
        editDisplay(c, [&](ViewDisplay& d) { d.field = FieldRef::fromKey(field->itemData(i).toString()); d.contours = true; });
    });
    auto* cmap = new QComboBox;
    cmap->setAccessibleName(QObject::tr("Colormap"));
    cmap->setIconSize(QSize(48, 12));
    for (ColormapId id : allColormaps()) cmap->addItem(colormapIcon(id), colormap(id).name, int(id));
    cmap->setCurrentIndex(cmap->findData(int(d.colormap)));
    QObject::connect(cmap, &QComboBox::activated, cmap, [c, cmap](int i) {
        editDisplay(c, [&](ViewDisplay& d) { d.colormap = ColormapId(cmap->itemData(i).toInt()); d.colormapFollowsRole = false; });
    });
    f->addRow(QObject::tr("Colormap"), cmap);
    if (d.symmetricRange())
        f->addRow({}, label(QObject::tr("Differences use a diverging map centred on zero unless you pick one."), "faint"));

    auto* style = new QWidget;
    auto* sh = new QHBoxLayout(style);
    sh->setContentsMargins(0, 0, 0, 0);
    auto* smooth = new QRadioButton(QObject::tr("Smooth"));
    auto* banded = new QRadioButton(QObject::tr("Banded"));
    auto* bands = new QSpinBox;
    bands->setRange(4, 40);
    bands->setValue(d.bands > 0 ? d.bands : 20);
    bands->setAccessibleName(QObject::tr("Number of bands"));
    bands->setEnabled(d.bands > 0);
    (d.bands > 0 ? banded : smooth)->setChecked(true);
    sh->addWidget(smooth);
    sh->addWidget(banded);
    sh->addWidget(bands);
    sh->addStretch(1);
    f->addRow(QObject::tr("Style"), style);
    auto applyBands = [c, banded, bands] {
        bands->setEnabled(banded->isChecked());
        editDisplay(c, [&](ViewDisplay& d) { d.bands = banded->isChecked() ? bands->value() : 0; });
    };
    QObject::connect(banded, &QRadioButton::toggled, banded, applyBands);
    QObject::connect(bands, &QSpinBox::valueChanged, bands, applyBands);

    auto* rangeW = new QWidget;
    auto* rh = new QHBoxLayout(rangeW);
    rh->setContentsMargins(0, 0, 0, 0);
    auto* autoR = new QCheckBox(d.symmetricRange() ? QObject::tr("Auto (symmetric)") : QObject::tr("Auto"));
    autoR->setChecked(d.rangeMode == RangeMode::Auto);
    const ColorRange shown = p->view()->resolvedRange();
    auto* lo = numberEdit(d.rangeMode == RangeMode::Fixed ? d.fixedMin : shown.min, QObject::tr("Minimum"));
    auto* hi = numberEdit(d.rangeMode == RangeMode::Fixed ? d.fixedMax : shown.max, QObject::tr("Maximum"));
    lo->setEnabled(!autoR->isChecked());
    hi->setEnabled(!autoR->isChecked());
    rh->addWidget(autoR);
    rh->addWidget(lo, 1);
    rh->addWidget(label("–"));
    rh->addWidget(hi, 1);
    f->addRow(QObject::tr("Range"), rangeW);
    auto* rangeErr = label({}, "error");
    rangeErr->hide();
    f->addRow({}, rangeErr);
    auto applyRange = [c, autoR, lo, hi, rangeErr] {
        lo->setEnabled(!autoR->isChecked());
        hi->setEnabled(!autoR->isChecked());
        if (autoR->isChecked()) {
            reportError(rangeErr, {});
            editDisplay(c, [](ViewDisplay& d) { d.rangeMode = RangeMode::Auto; });
            return;
        }
        auto a = numberOf(lo), b = numberOf(hi);
        if (!a || !b || *b <= *a) return reportError(rangeErr, QObject::tr("Enter a minimum below the maximum."));
        reportError(rangeErr, {});
        editDisplay(c, [&](ViewDisplay& d) { d.rangeMode = RangeMode::Fixed; d.fixedMin = *a; d.fixedMax = *b; });
    };
    QObject::connect(autoR, &QCheckBox::toggled, autoR, applyRange);
    QObject::connect(lo, &QLineEdit::editingFinished, lo, applyRange);
    QObject::connect(hi, &QLineEdit::editingFinished, hi, applyRange);
    v->addLayout(f);

    auto* opts = new QHBoxLayout;
    auto* edges = new QCheckBox(QObject::tr("Mesh edges"));
    edges->setChecked(d.edges);
    QObject::connect(edges, &QCheckBox::toggled, edges, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.edges = on; }); });
    auto* light = new QCheckBox(QObject::tr("Lighting"));
    light->setChecked(d.lighting);
    QObject::connect(light, &QCheckBox::toggled, light, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.lighting = on; }); });
    opts->addWidget(edges);
    opts->addWidget(light);
    opts->addStretch(1);
    v->addLayout(opts);

    v->addWidget(section(QObject::tr("Surfaces")));
    v->addWidget(surfaceList(c));
    return page;
}

QWidget* vectors(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Vectors"), QObject::tr("Results ▸ Graphics"), &v);
    v->addWidget(targetLine(c));
    auto* p = c.graphics->active();
    const auto s = p ? p->scene() : nullptr;
    if (!s) return page;
    const ViewDisplay d = p->display();
    auto* on = new QCheckBox(QObject::tr("Show velocity vectors"));
    on->setChecked(d.vectors);
    v->addWidget(on);
    auto* f = formLayout();
    auto* item = new QComboBox;
    item->setAccessibleName(QObject::tr("Vectors on"));
    for (const auto& it : s->items)
        if (it.mesh && (it.mesh->hasField({"velocity", FieldRole::Kutta, -1}) || it.mesh->hasField({"wall_shear", FieldRole::Kutta, -1})))
            item->addItem(it.label, it.id);
    item->setCurrentIndex(std::max(0, item->findData(d.vectorItem)));
    f->addRow(QObject::tr("Surface"), item);
    auto* len = new QLineEdit(d.vectorLength > 0 ? QString::number(d.vectorLength) : QString());
    len->setPlaceholderText(QObject::tr("automatic"));
    len->setAccessibleName(QObject::tr("Arrow length"));
    f->addRow(QObject::tr("Arrow length (m)"), len);
    auto* stride = new QSpinBox;
    stride->setRange(1, 64);
    stride->setValue(d.vectorStride);
    stride->setAccessibleName(QObject::tr("Skip"));
    f->addRow(QObject::tr("Every Nth vertex"), stride);
    v->addLayout(f);
    if (item->count() == 0) v->addWidget(label(QObject::tr("No surface in this dataset carries a velocity field. Add a plane in Results ▸ Graphics ▸ Planes."), "muted"));
    auto apply = [c, on, item, len, stride] {
        editDisplay(c, [&](ViewDisplay& d) {
            d.vectors = on->isChecked() && item->count() > 0;
            d.vectorItem = item->currentData().toString();
            d.vectorLength = float(numberOf(len).value_or(0));
            d.vectorStride = stride->value();
            if (d.vectors && !d.visible.contains(d.vectorItem)) d.visible << d.vectorItem;
        });
    };
    QObject::connect(on, &QCheckBox::toggled, on, apply);
    QObject::connect(item, &QComboBox::activated, item, apply);
    QObject::connect(len, &QLineEdit::editingFinished, len, apply);
    QObject::connect(stride, &QSpinBox::valueChanged, stride, apply);
    return page;
}

namespace {

// Shared shape of pages that request post-processing from the CLI.
struct PostForm {
    QWidget* page;
    QVBoxLayout* body;
    QPushButton* go;
    QLabel* error;
    QLabel* status;
};

PostForm postPage(const InspectorContext& c, const QString& title, const QString& crumb, const QString& action)
{
    PostForm pf{};
    pf.page = makePage(title, crumb, &pf.body);
    const QString run = currentRunId(c);
    const RunRecord* r = c.session->history().find(run);
    pf.body->addWidget(label(r ? QObject::tr("From %1").arg(r->title()) : QObject::tr("Needs a finished run."), "muted"));
    pf.go = button(action, "play", "primary");
    pf.go->setEnabled(r != nullptr);
    pf.error = label({}, "error");
    pf.error->hide();
    pf.status = label({}, "muted");
    pf.status->hide();
    return pf;
}

void finishPostPage(PostForm& pf)
{
    auto* h = new QHBoxLayout;
    h->addWidget(pf.go);
    h->addStretch(1);
    pf.body->addLayout(h);
    pf.body->addWidget(pf.status);
    pf.body->addWidget(pf.error);
}

void runPost(const InspectorContext& c, PostForm pf, QJsonObject item, std::function<void(const QJsonObject&)> onDone)
{
    const QString run = currentRunId(c);
    if (run.isEmpty()) return reportError(pf.error, QObject::tr("Run the calculation first."));
    reportError(pf.error, {});
    pf.go->setEnabled(false);
    pf.status->setText(QObject::tr("Computing in kutta…"));
    pf.status->show();
    QPointer<QPushButton> go = pf.go;
    QPointer<QLabel> err = pf.error, st = pf.status;
    c.session->postProcess(run, item, [go, err, st, onDone](Result<QJsonObject> r) {
        if (go) go->setEnabled(true);
        if (st) st->hide();
        if (!r) {
            if (err) reportError(err, r.error());
            return;
        }
        onDone(r.value());
    });
}

} // namespace

QWidget* pathlines(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("Pathlines"), QObject::tr("Results ▸ Graphics"), QObject::tr("Compute pathlines"));
    pf.body->addWidget(label(QObject::tr("Streamlines of the converged velocity field, seeded along a straight rake."), "muted"));
    const CaseSummary* cs = c.session->hasCase() ? &c.session->caseSummary() : nullptr;
    QVector3D from(-0.6f, -0.5f, 0.12f), to(-0.6f, 0.5f, 0.12f);
    if (cs && cs->bounds.valid && !std::any_of(cs->zones.begin(), cs->zones.end(), [](const Zone& z) { return z.name.contains("ground"); })) {
        const QVector3D a = cs->bounds.min, b = cs->bounds.max, m = cs->bounds.center();
        from = QVector3D(a.x(), m.y(), m.z());
        to = QVector3D(a.x(), b.y(), m.z());
    }
    auto* f = formLayout();
    auto rakeFrom = vec3Edit(from, QObject::tr("Rake start"));
    auto rakeTo = vec3Edit(to, QObject::tr("Rake end"));
    auto* count = new QSpinBox;
    count->setRange(1, 400);
    count->setValue(40);
    count->setAccessibleName(QObject::tr("Number of pathlines"));
    f->addRow(QObject::tr("Rake start (m)"), rakeFrom.widget);
    f->addRow(QObject::tr("Rake end (m)"), rakeTo.widget);
    f->addRow(QObject::tr("Pathlines"), count);
    pf.body->addLayout(f);
    finishPostPage(pf);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, rakeFrom, rakeTo, count] {
        auto a = rakeFrom.value(), b = rakeTo.value();
        if (!a || !b) return reportError(pf.error, QObject::tr("Enter numbers for both rake points."));
        QJsonObject item{{"id", newPostId("pathlines")}, {"kind", "streamlines"}, {"name", QObject::tr("Pathlines (%1)").arg(count->value())},
                         {"seed", QJsonObject{{"from", QJsonArray{a->x(), a->y(), a->z()}}, {"to", QJsonArray{b->x(), b->y(), b->z()}}, {"count", count->value()}}}};
        runPost(c, pf, item, [c](const QJsonObject&) {
            editDisplay(c, [](ViewDisplay& d) { d.field = {"velocity_magnitude", FieldRole::Kutta, -1}; d.contours = true; });
        });
    });
    return pf.page;
}

QWidget* isoSurfaces(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("Iso-Surfaces"), QObject::tr("Results ▸ Graphics"), QObject::tr("Compute iso-surface"));
    auto* f = formLayout();
    auto* field = new QComboBox;
    field->setAccessibleName(QObject::tr("Field"));
    for (const char* q : {"total_pressure_coefficient", "velocity_magnitude", "pressure_coefficient", "pressure", "viscosity_ratio"})
        field->addItem(fluentLabel(q), QString::fromLatin1(q));
    auto* value = numberEdit(0, QObject::tr("Iso-value"));
    f->addRow(QObject::tr("Field"), field);
    f->addRow(QObject::tr("Value"), value);
    pf.body->addLayout(f);
    pf.body->addWidget(label(QObject::tr("Total pressure coefficient = 0 outlines the wake: where the flow has lost its free-stream total pressure."), "faint"));
    finishPostPage(pf);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, field, value] {
        auto val = numberOf(value);
        if (!val) return reportError(pf.error, QObject::tr("Enter a number."));
        QJsonObject item{{"id", newPostId("iso")}, {"kind", "iso-surface"}, {"field", field->currentData().toString()}, {"value", *val},
                         {"name", QStringLiteral("%1 = %2").arg(field->currentText(), formatNumber(*val))}};
        runPost(c, pf, item, [](const QJsonObject&) {});
    });
    return pf.page;
}

QWidget* planes(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Planes"), QObject::tr("Results ▸ Graphics"), &v);
    const QString run = currentRunId(c);
    const RunRecord* r = c.session->history().find(run);
    if (!r) {
        v->addWidget(label(QObject::tr("Planes are cut from a finished run. Run the calculation first."), "muted"));
        return page;
    }
    v->addWidget(label(QObject::tr("Cutting planes exported from %1. Kutta samples the cell solution on each plane.").arg(r->title()), "muted"));
    auto* list = new QListWidget;
    list->setAccessibleName(QObject::tr("Planes"));
    auto planesNow = std::make_shared<std::vector<Plane>>(c.session->slices(run));
    auto refill = [list, planesNow] {
        list->clear();
        for (const Plane& p : *planesNow) {
            ViewPart part;
            part.origin = p.origin;
            part.normal = p.normal;
            list->addItem(itemLabelForSlice(part));
        }
    };
    refill();
    list->setMinimumHeight(110);
    v->addWidget(list);
    auto* f = formLayout();
    auto* axis = new QComboBox;
    axis->addItem(QObject::tr("x = constant"), 0);
    axis->addItem(QObject::tr("y = constant"), 1);
    axis->addItem(QObject::tr("z = constant"), 2);
    axis->setAccessibleName(QObject::tr("Plane orientation"));
    auto* pos = numberEdit(0, QObject::tr("Plane position"));
    f->addRow(QObject::tr("Orientation"), axis);
    f->addRow(QObject::tr("Position (m)"), pos);
    v->addLayout(f);
    auto* row = new QHBoxLayout;
    auto* add = button(QObject::tr("Add plane"), "plus");
    auto* remove = button(QObject::tr("Remove"), "minus");
    auto* apply = button(QObject::tr("Export planes"), "refresh-cw", "primary");
    row->addWidget(add);
    row->addWidget(remove);
    row->addStretch(1);
    row->addWidget(apply);
    v->addLayout(row);
    auto* err = label({}, "error");
    err->hide();
    v->addWidget(err);
    QObject::connect(add, &QPushButton::clicked, page, [planesNow, refill, axis, pos, err, c] {
        auto p = numberOf(pos);
        if (!p) return reportError(err, QObject::tr("Enter a position."));
        reportError(err, {});
        const int a = axis->currentData().toInt();
        QVector3D o = c.session->hasCase() && c.session->caseSummary().bounds.valid ? c.session->caseSummary().bounds.center() : QVector3D();
        o[a] = float(*p);
        QVector3D n;
        n[a] = 1;
        planesNow->push_back({o, n});
        refill();
    });
    QObject::connect(remove, &QPushButton::clicked, page, [planesNow, refill, list] {
        const int r = list->currentRow();
        if (r >= 0 && r < int(planesNow->size())) planesNow->erase(planesNow->begin() + r);
        refill();
    });
    QObject::connect(apply, &QPushButton::clicked, page, [c, run, planesNow, apply, err] {
        apply->setEnabled(false);
        QPointer<QPushButton> b = apply;
        QPointer<QLabel> e = err;
        c.session->setSlices(run, *planesNow, [b, e](const QString&, Result<std::shared_ptr<SceneData>> s) {
            if (b) b->setEnabled(true);
            if (e) reportError(e, s ? QString() : s.error());
        });
    });
    return page;
}

QWidget* mesh(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Mesh"), QObject::tr("Results ▸ Graphics"), &v);
    v->addWidget(targetLine(c));
    auto* p = c.graphics->active();
    if (!p || !p->scene()) return page;
    const ViewDisplay d = p->display();
    auto* edges = new QCheckBox(QObject::tr("Show mesh edges"));
    edges->setChecked(d.edges);
    QObject::connect(edges, &QCheckBox::toggled, edges, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.edges = on; }); });
    auto* solid = new QCheckBox(QObject::tr("Colour by field"));
    solid->setChecked(d.contours);
    QObject::connect(solid, &QCheckBox::toggled, solid, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.contours = on; }); });
    v->addWidget(edges);
    v->addWidget(solid);
    std::size_t tris = 0;
    for (const auto& it : p->scene()->items)
        if (it.mesh) tris += it.mesh->triangleCount();
    v->addWidget(label(QObject::tr("%1 surface triangles in this dataset (boundary faces split into triangles).").arg(QLocale().toString(qulonglong(tris))), "faint"));
    v->addWidget(section(QObject::tr("Surfaces")));
    v->addWidget(surfaceList(c));
    return page;
}

// ---- Plots ---------------------------------------------------------------------------

QWidget* xyPlot(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("XY Plot"), QObject::tr("Results ▸ Plots"), QObject::tr("Plot along line"));
    pf.body->addWidget(section(QObject::tr("Line probe")));
    auto* f = formLayout();
    QVector3D from(-1, 0, 0.1f), to(4, 0, 0.1f);
    if (c.session->hasCase() && c.session->caseSummary().bounds.valid) {
        const Bounds b = c.session->caseSummary().bounds;
        const bool ext = std::any_of(c.session->caseSummary().zones.begin(), c.session->caseSummary().zones.end(), [](const Zone& z) { return z.name.contains("ground"); });
        if (!ext) { from = QVector3D(b.min.x(), b.center().y(), b.center().z()); to = QVector3D(b.max.x(), b.center().y(), b.center().z()); }
    }
    auto a = vec3Edit(from, QObject::tr("Line start"));
    auto b = vec3Edit(to, QObject::tr("Line end"));
    auto* samples = new QSpinBox;
    samples->setRange(10, 5000);
    samples->setValue(200);
    samples->setAccessibleName(QObject::tr("Samples"));
    auto* field = new QComboBox;
    field->setAccessibleName(QObject::tr("Field"));
    for (const char* q : {"velocity_magnitude", "pressure", "pressure_coefficient", "total_pressure_coefficient", "k", "omega"})
        field->addItem(fluentLabel(q), QString::fromLatin1(q));
    f->addRow(QObject::tr("Start (m)"), a.widget);
    f->addRow(QObject::tr("End (m)"), b.widget);
    f->addRow(QObject::tr("Samples"), samples);
    f->addRow(QObject::tr("Field"), field);
    pf.body->addLayout(f);
    finishPostPage(pf);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, a, b, samples, field] {
        auto pa = a.value(), pb = b.value();
        if (!pa || !pb) return reportError(pf.error, QObject::tr("Enter numbers for both points."));
        const QString q = field->currentData().toString();
        QJsonObject item{{"id", newPostId("line")}, {"kind", "line-probe"}, {"samples", samples->value()},
                         {"from", QJsonArray{pa->x(), pa->y(), pa->z()}}, {"to", QJsonArray{pb->x(), pb->y(), pb->z()}},
                         {"fields", QJsonArray{q}}};
        const QString runId = currentRunId(c);
        runPost(c, pf, item, [c, q, runId](const QJsonObject& m) {
            std::vector<PlotSeries> series;
            int k = 0;
            QString unit;
            for (const QJsonValue& sv : m.value("series").toArray()) {
                const QJsonObject so = sv.toObject();
                PlotSeries ps;
                ps.name = displayLabel(FieldRef::parse(so.value("name").toString()));
                unit = so.value("unit").toString();
                for (const QJsonValue& x : so.value("x").toArray()) ps.x.push_back(x.toDouble());
                for (const QJsonValue& y : so.value("y").toArray()) ps.y.push_back(y.isDouble() ? y.toDouble() : NAN);
                ps.color = tokens().series[std::size_t(k++) % tokens().series.size()];
                series.push_back(std::move(ps));
            }
            const RunRecord* r = c.session->history().find(runId);
            QStringList badges;
            if (m.value("synthetic").toBool()) badges << "SYNTHETIC";
            const QString prov = m.value("provenance").toString();
            if (prov.contains("APPROXIMATE")) badges << "APPROXIMATE";
            if (prov.contains("MODIFIED")) badges << "MODIFIED";
            c.charts->setXY(QObject::tr("%1 along line · %2").arg(fluentLabel(q), r ? r->title() : QString()),
                            QObject::tr("Distance [m]"), unit.isEmpty() ? fluentLabel(q) : QStringLiteral("%1 [%2]").arg(fluentLabel(q), prettyUnit(unit)),
                            std::move(series), badges, prov);
        });
    });

    // Section through surfaces: Cp along a wing chord, the elbow's wall pressure.
    pf.body->addSpacing(8);
    pf.body->addWidget(section(QObject::tr("Section through surfaces")));
    auto* sf = formLayout();
    auto originEdit = vec3Edit(QVector3D(0, 0.3f, 0), QObject::tr("Section origin"));
    auto normalEdit = vec3Edit(QVector3D(0, 1, 0), QObject::tr("Section normal"));
    auto* sfield = new QComboBox;
    for (const char* q : {"pressure_coefficient", "pressure", "wall_shear_magnitude", "y_plus"}) sfield->addItem(fluentLabel(q), QString::fromLatin1(q));
    sfield->setAccessibleName(QObject::tr("Section field"));
    auto* zones = new QListWidget;
    zones->setAccessibleName(QObject::tr("Section surfaces"));
    if (c.session->hasCase())
        for (const auto& z : c.session->caseSummary().zones)
            if (z.type == "wall" && !z.name.contains("ground")) {
                auto* it = new QListWidgetItem(z.name);
                it->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
                it->setCheckState(Qt::Unchecked);
                zones->addItem(it);
            }
    if (zones->count()) zones->item(0)->setCheckState(Qt::Checked);
    zones->setMaximumHeight(130);
    auto* xAxis = new QComboBox;
    xAxis->addItems({"x", "y", "z"});
    xAxis->setAccessibleName(QObject::tr("Horizontal axis"));
    sf->addRow(QObject::tr("Origin (m)"), originEdit.widget);
    sf->addRow(QObject::tr("Normal"), normalEdit.widget);
    sf->addRow(QObject::tr("Field"), sfield);
    sf->addRow(QObject::tr("Plot against"), xAxis);
    sf->addRow(QObject::tr("Surfaces"), zones);
    pf.body->addLayout(sf);
    auto* secGo = button(QObject::tr("Plot section"), "spline");
    auto* secErr = label({}, "error");
    secErr->hide();
    pf.body->addWidget(secGo, 0, Qt::AlignLeft);
    pf.body->addWidget(secErr);
    secGo->setEnabled(pf.go->isEnabled());
    QObject::connect(secGo, &QPushButton::clicked, pf.page, [c, pf, originEdit, normalEdit, sfield, zones, xAxis, secGo, secErr] {
        auto o = originEdit.value(), n = normalEdit.value();
        if (!o || !n || n->length() == 0) return reportError(secErr, QObject::tr("Enter an origin and a non-zero normal."));
        QJsonArray zl;
        for (int i = 0; i < zones->count(); ++i)
            if (zones->item(i)->checkState() == Qt::Checked) zl.append(zones->item(i)->text());
        if (zl.isEmpty()) return reportError(secErr, QObject::tr("Choose at least one surface."));
        reportError(secErr, {});
        const QString q = sfield->currentData().toString();
        QJsonObject item{{"id", newPostId("section")}, {"kind", "section"}, {"field", q}, {"zones", zl}, {"x_axis", xAxis->currentIndex()},
                         {"origin", QJsonArray{o->x(), o->y(), o->z()}}, {"normal", QJsonArray{n->x(), n->y(), n->z()}}};
        PostForm sp = pf;
        sp.go = secGo;
        sp.error = secErr;
        runPost(c, sp, item, [c, q](const QJsonObject& m) {
            std::vector<PlotSeries> series;
            int k = 0;
            QString unit;
            for (const QJsonValue& sv : m.value("series").toArray()) {
                const QJsonObject so = sv.toObject();
                PlotSeries ps;
                ps.name = so.value("name").toString();
                unit = so.value("unit").toString();
                for (const QJsonValue& x : so.value("x").toArray()) ps.x.push_back(x.toDouble());
                for (const QJsonValue& y : so.value("y").toArray()) ps.y.push_back(y.isDouble() ? y.toDouble() : NAN);
                // Section samples are unordered segment end points: plot as markers.
                std::vector<std::size_t> idx(ps.x.size());
                std::iota(idx.begin(), idx.end(), 0);
                std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return ps.x[a] < ps.x[b]; });
                PlotSeries sorted = ps;
                for (std::size_t i = 0; i < idx.size(); ++i) { sorted.x[i] = ps.x[idx[i]]; sorted.y[i] = ps.y[idx[i]]; }
                sorted.markers = true;
                sorted.dashed = false;
                sorted.color = tokens().series[std::size_t(k++) % tokens().series.size()];
                series.push_back(std::move(sorted));
            }
            QStringList badges;
            if (m.value("synthetic").toBool()) badges << "SYNTHETIC";
            const QString prov = m.value("provenance").toString();
            if (prov.contains("APPROXIMATE")) badges << "APPROXIMATE";
            c.charts->setXY(QObject::tr("%1 on section").arg(fluentLabel(q)), m.value("x_label").toString() + " [m]",
                            unit.isEmpty() || unit == "1" ? fluentLabel(q) : QStringLiteral("%1 [%2]").arg(fluentLabel(q), prettyUnit(unit)),
                            std::move(series), badges, prov);
        });
    });
    return pf.page;
}

QWidget* residuals(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Residuals"), QObject::tr("Results ▸ Plots"), &v);
    c.charts->showTab(ChartsPanel::Tab::Residuals);
    v->addWidget(label(QObject::tr("Scaled residuals per iteration, in Fluent's order. The chart follows the live run; select runs in Run History to overlay them."), "muted"));
    auto* log = new QCheckBox(QObject::tr("Logarithmic axis"));
    log->setChecked(c.charts->residuals()->yScale() == Axis::Scale::Log);
    QObject::connect(log, &QCheckBox::toggled, log, [c](bool on) { c.charts->residuals()->setYScale(on ? Axis::Scale::Log : Axis::Scale::Linear); });
    v->addWidget(log);
    auto* row = new QHBoxLayout;
    auto* csv = button(QObject::tr("Export CSV…"), "file-down");
    auto* png = button(QObject::tr("Export PNG…"), "image-down");
    row->addWidget(csv);
    row->addWidget(png);
    row->addStretch(1);
    v->addLayout(row);
    QObject::connect(csv, &QPushButton::clicked, page, [c, csv] {
        const QString f = QFileDialog::getSaveFileName(csv->window(), QObject::tr("Export residuals"), "residuals.csv", "CSV (*.csv)");
        if (!f.isEmpty()) c.charts->residuals()->exportCsv(f);
    });
    QObject::connect(png, &QPushButton::clicked, page, [c, png] {
        const QString f = QFileDialog::getSaveFileName(png->window(), QObject::tr("Export residuals"), "residuals.png", "PNG (*.png)");
        if (!f.isEmpty()) c.charts->residuals()->exportPng(f);
    });
    v->addWidget(label(QObject::tr("Drag to pan, scroll to zoom the iteration axis (Shift: values), Shift-drag to zoom a box, double-click to fit. Hover to read every residual at an iteration; the report plot follows the same cursor."), "faint"));
    return page;
}

QWidget* reportPlots(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Report Plots"), QObject::tr("Results ▸ Plots"), &v);
    c.charts->showTab(ChartsPanel::Tab::Reports);
    v->addWidget(label(QObject::tr("Force reports defined in the case, per iteration. Choose force or coefficient above the chart."), "muted"));
    if (c.session->hasCase() && c.session->caseSummary().reports.empty())
        v->addWidget(label(QObject::tr("This case defines no reports, so there is nothing to plot. Wall forces are in Results ▸ Reports ▸ Forces."), "faint"));
    return page;
}

// ---- Reports -----------------------------------------------------------------------

QWidget* forces(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Forces"), QObject::tr("Results ▸ Reports"), &v);
    const QString run = currentRunId(c);
    const RunRecord* r = c.session->history().find(run);
    if (!r) {
        v->addWidget(label(QObject::tr("Forces are written when a run finishes."), "muted"));
        return page;
    }
    auto ft = c.session->history().forces(run);
    auto* top = new QHBoxLayout;
    top->addWidget(label(QObject::tr("From %1 (forces.csv)").arg(r->title()), "muted"), 1);
    if (r->approximate()) top->addWidget(badge("APPROXIMATE", "approximate"));
    if (r->modified()) top->addWidget(badge("MODIFIED", "modified"));
    if (r->synthetic) top->addWidget(badge("SYNTHETIC", "modified"));
    else if (r->fixture) top->addWidget(badge("FIXTURE REPLAY", "modified"));
    v->addLayout(top);
    if (!ft) {
        v->addWidget(label(QObject::tr("This run wrote no forces.csv."), "error"));
        return page;
    }
    if (!ft->notice.isEmpty()) v->addWidget(label(ft->notice, "faint"));
    auto* axis = new QComboBox;
    axis->addItems({QObject::tr("Wall forces along X"), QObject::tr("Wall forces along Y"), QObject::tr("Wall forces along Z")});
    axis->setAccessibleName(QObject::tr("Force direction for walls"));
    axis->setCurrentIndex(0);
    v->addWidget(axis, 0, Qt::AlignLeft);
    auto* table = new QTableWidget;
    table->setColumnCount(5);
    table->setHorizontalHeaderLabels({QObject::tr("Zone / report"), QObject::tr("Pressure [N]"), QObject::tr("Viscous [N]"), QObject::tr("Total [N]"), QObject::tr("Coefficient")});
    table->verticalHeader()->hide();
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setShowGrid(false);
    table->setAccessibleName(QObject::tr("Forces"));
    QFont tab = table->font();
    tab.setFeature(QFont::Tag("tnum"), 1);
    table->setFont(tab);
    auto fill = [table, rows = ft->rows](int ax) {
        table->setRowCount(0);
        double sp = 0, sv = 0;
        int walls = 0;
        auto num = [](double x) {
            auto* it = new QTableWidgetItem(QLocale().toString(x, 'f', 3));
            it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
            return it;
        };
        for (const auto& row : rows) {
            const int r = table->rowCount();
            table->insertRow(r);
            const bool report = row.kind == "report";
            const double p = report && row.pressureProjection ? *row.pressureProjection : row.pressure[ax];
            const double vv = report && row.viscousProjection ? *row.viscousProjection : row.viscous[ax];
            auto* name = new QTableWidgetItem(report ? row.name + QObject::tr("  (report)") : row.name);
            if (report) { QFont f = name->font(); f.setBold(true); name->setFont(f); }
            table->setItem(r, 0, name);
            table->setItem(r, 1, num(p));
            table->setItem(r, 2, num(vv));
            table->setItem(r, 3, num(p + vv));
            table->setItem(r, 4, row.coefficient ? num(*row.coefficient) : new QTableWidgetItem(QStringLiteral("—")));
            if (!report) { sp += p; sv += vv; ++walls; }
        }
        if (walls > 1) {
            const int r = table->rowCount();
            table->insertRow(r);
            auto* name = new QTableWidgetItem(QObject::tr("Net (all walls)"));
            QFont f = name->font();
            f.setBold(true);
            name->setFont(f);
            table->setItem(r, 0, name);
            table->setItem(r, 1, num(sp));
            table->setItem(r, 2, num(sv));
            table->setItem(r, 3, num(sp + sv));
            table->setItem(r, 4, new QTableWidgetItem(QStringLiteral("—")));
        }
        table->resizeColumnsToContents();
        table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
        table->setMinimumHeight(std::min(460, 34 + table->rowCount() * 26));
    };
    fill(0);
    QObject::connect(axis, &QComboBox::currentIndexChanged, table, fill);
    v->addWidget(table);
    auto* exportBtn = button(QObject::tr("Save forces.csv as…"), "file-down");
    const QString src = QDir(r->dir).filePath("forces.csv");
    QObject::connect(exportBtn, &QPushButton::clicked, exportBtn, [src, exportBtn] {
        const QString f = QFileDialog::getSaveFileName(exportBtn->window(), QObject::tr("Save forces"), "forces.csv", "CSV (*.csv)");
        if (!f.isEmpty()) { QFile::remove(f); QFile::copy(src, f); }
    });
    v->addWidget(exportBtn, 0, Qt::AlignLeft);
    v->addWidget(label(QObject::tr("Report rows are projected on the report's direction; wall rows on the axis chosen above."), "faint"));
    return page;
}

QWidget* integrals(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("Surface Integrals"), QObject::tr("Results ▸ Reports"), QObject::tr("Compute"));
    auto* f = formLayout();
    auto* report = new QComboBox;
    report->addItem(QObject::tr("Area-weighted average"), "area-weighted-average");
    report->addItem(QObject::tr("Integral"), "integral");
    report->addItem(QObject::tr("Area"), "area");
    report->addItem(QObject::tr("Minimum"), "minimum");
    report->addItem(QObject::tr("Maximum"), "maximum");
    report->setAccessibleName(QObject::tr("Report type"));
    auto* field = new QComboBox;
    field->setAccessibleName(QObject::tr("Field"));
    for (const char* q : {"pressure", "pressure_coefficient", "velocity_magnitude", "wall_shear_magnitude", "y_plus"})
        field->addItem(fluentLabel(q), QString::fromLatin1(q));
    auto* withFluent = new QCheckBox(QObject::tr("Also compute on Fluent's field (needs a data file)"));
    withFluent->setEnabled(c.session->hasCase() && !c.session->dataPath().isEmpty());
    auto* zones = new QListWidget;
    zones->setAccessibleName(QObject::tr("Surfaces"));
    if (c.session->hasCase())
        for (const auto& z : c.session->caseSummary().zones) {
            if (z.type == "interior") continue;
            auto* it = new QListWidgetItem(QStringLiteral("%1  ·  %2").arg(z.name, z.type));
            it->setData(Qt::UserRole, z.name);
            it->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
            it->setCheckState(z.type.contains("inlet") || z.type.contains("outlet") ? Qt::Checked : Qt::Unchecked);
            zones->addItem(it);
        }
    zones->setMinimumHeight(std::min(220, 10 + zones->count() * 24));
    f->addRow(QObject::tr("Report"), report);
    f->addRow(QObject::tr("Field"), field);
    f->addRow(QObject::tr("Surfaces"), zones);
    pf.body->addLayout(f);
    pf.body->addWidget(withFluent);
    finishPostPage(pf);
    auto* table = new QTableWidget(0, 3);
    table->setHorizontalHeaderLabels({QObject::tr("Surface"), QObject::tr("Kutta"), QObject::tr("Fluent")});
    table->verticalHeader()->hide();
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setShowGrid(false);
    table->setAccessibleName(QObject::tr("Surface integral results"));
    table->hide();
    pf.body->addWidget(table);
    auto* note = label({}, "faint");
    pf.body->addWidget(note);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, report, field, zones, withFluent, table, note] {
        QJsonArray zl;
        for (int i = 0; i < zones->count(); ++i)
            if (zones->item(i)->checkState() == Qt::Checked) zl.append(zones->item(i)->data(Qt::UserRole).toString());
        if (zl.isEmpty()) return reportError(pf.error, QObject::tr("Choose at least one surface."));
        const QString q = field->currentData().toString();
        auto showRows = [table, note](const QJsonObject& m, int col) {
            const QJsonArray rows = m.value("rows").toArray();
            table->show();
            if (table->rowCount() < rows.size() + 1) table->setRowCount(int(rows.size()) + 1);
            const QString unit = prettyUnit(m.value("unit").toString());
            auto num = [](const QJsonValue& x) {
                auto* it = new QTableWidgetItem(x.isDouble() ? formatNumber(x.toDouble(), 6) : QStringLiteral("—"));
                it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
                return it;
            };
            for (int i = 0; i < rows.size(); ++i) {
                const QJsonObject r = rows[i].toObject();
                table->setItem(i, 0, new QTableWidgetItem(r.value("zone").toString()));
                table->setItem(i, col, num(r.value("value")));
            }
            auto* tot = new QTableWidgetItem(QObject::tr("Net"));
            QFont bf = tot->font();
            bf.setBold(true);
            tot->setFont(bf);
            table->setItem(int(rows.size()), 0, tot);
            table->setItem(int(rows.size()), col, num(m.value("total")));
            table->setHorizontalHeaderItem(1, new QTableWidgetItem(unit.isEmpty() ? QObject::tr("Kutta") : QObject::tr("Kutta [%1]").arg(unit)));
            table->resizeColumnsToContents();
            table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
            table->setMinimumHeight(34 + table->rowCount() * 26);
            const QString prov = m.value("provenance").toString();
            note->setText(QObject::tr("Integrated over the exported surface triangles by kutta post.") + (prov.isEmpty() ? QString() : "\n" + prov));
        };
        table->setRowCount(0);
        table->setColumnHidden(2, !withFluent->isChecked());
        QJsonObject item{{"id", newPostId("integral")}, {"kind", "surface-integral"}, {"field", q}, {"report", report->currentData().toString()}, {"zones", zl}};
        const bool fluent = withFluent->isChecked();
        runPost(c, pf, item, [c, pf, showRows, item, q, fluent](const QJsonObject& m) {
            showRows(m, 1);
            if (!fluent) return;
            QJsonObject ref = item;
            ref["id"] = newPostId("integral");
            ref["field"] = q + "@fluent";
            runPost(c, pf, ref, [showRows](const QJsonObject& m2) { showRows(m2, 2); });
        });
    });
    return pf.page;
}

// ---- Runs -------------------------------------------------------------------------

namespace {

QString duration(double s)
{
    if (s <= 0) return QStringLiteral("—");
    if (s < 90) return QStringLiteral("%1 s").arg(s, 0, 'f', 1);
    return QStringLiteral("%1 min %2 s").arg(int(s) / 60).arg(int(s) % 60);
}

QString flags(const RunRecord& r)
{
    QStringList f;
    if (r.approximate()) f << "APPROXIMATE";
    if (r.modified()) f << "MODIFIED";
    if (r.synthetic) f << "SYNTHETIC";
    else if (r.fixture) f << "FIXTURE";
    return f.join(" · ");
}

} // namespace

QWidget* runHistory(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Run History"), QObject::tr("Runs"), &v);
    if (!c.session->hasCase()) {
        v->addWidget(label(QObject::tr("Open a case to see its runs."), "muted"));
        return page;
    }
    const auto& runs = c.session->history().runs();
    if (runs.empty()) {
        v->addWidget(label(QObject::tr("No runs yet. Every calculation is kept here with its settings, residuals and results."), "muted"));
        return page;
    }
    auto* table = new QTableWidget(int(runs.size()), 6);
    table->setHorizontalHeaderLabels({QObject::tr("Run"), QObject::tr("Status"), QObject::tr("Iter."), QObject::tr("Time"), QObject::tr("Marks"), QObject::tr("Started")});
    table->verticalHeader()->hide();
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setSelectionMode(QAbstractItemView::ExtendedSelection);
    table->setShowGrid(false);
    table->setAccessibleName(QObject::tr("Runs"));
    QFont tab = table->font();
    tab.setFeature(QFont::Tag("tnum"), 1);
    table->setFont(tab);
    const auto& t = tokens();
    for (int i = 0; i < int(runs.size()); ++i) {
        const RunRecord& r = runs[runs.size() - 1 - std::size_t(i)];
        auto* name = new QTableWidgetItem(r.title());
        name->setData(Qt::UserRole, r.id);
        name->setToolTip(QStringLiteral("%1\n%2").arg(QLocale().toString(r.started, QLocale::ShortFormat), QDir::toNativeSeparators(r.dir)));
        table->setItem(i, 0, name);
        auto* st = new QTableWidgetItem(runStatusLabel(r.status) + (r.status == RunStatus::Stopped && r.stopReason == "stopped" ? QObject::tr(" by you") : QString()));
        if (r.status == RunStatus::Failed || r.status == RunStatus::Interrupted) st->setForeground(t.red);
        if (r.status == RunStatus::Running) st->setForeground(t.dark ? t.accent : t.accent.darker(115));
        if (!r.error.isEmpty()) st->setToolTip(r.error);
        table->setItem(i, 1, st);
        auto* it = new QTableWidgetItem(QString::number(r.iterationsDone));
        it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        table->setItem(i, 2, it);
        auto* el = new QTableWidgetItem(duration(r.elapsed));
        el->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        table->setItem(i, 3, el);
        auto* fl = new QTableWidgetItem(flags(r));
        if (r.approximate() || r.modified() || r.synthetic || r.fixture) fl->setForeground(t.dark ? t.amber : t.amberText);
        table->setItem(i, 4, fl);
        table->setItem(i, 5, new QTableWidgetItem(QLocale().toString(r.started, QLocale::ShortFormat)));
    }
    table->resizeColumnsToContents();
    table->horizontalHeader()->setSectionResizeMode(4, QHeaderView::Stretch);
    table->setColumnHidden(5, true);  // start time is in the row tooltip
    table->setFixedHeight(std::min(380, 34 + table->rowCount() * 26));
    table->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    table->selectRow(0);
    v->addWidget(table);

    auto selected = [table] {
        QStringList ids;
        for (const auto& idx : table->selectionModel()->selectedRows()) ids << table->item(idx.row(), 0)->data(Qt::UserRole).toString();
        return ids;
    };
    auto* detail = new QWidget;
    auto* dv = new QVBoxLayout(detail);
    dv->setContentsMargins(0, 4, 0, 0);
    dv->setSpacing(6);
    v->addWidget(detail);
    auto* row1 = new QHBoxLayout;
    auto* show = button(QObject::tr("Show"), "eye");
    auto* plot = button(QObject::tr("Residuals"), "chart-line");
    auto* comp = button(QObject::tr("Compare"), "git-compare-arrows");
    row1->addWidget(show);
    row1->addWidget(plot);
    row1->addWidget(comp);
    row1->addStretch(1);
    v->addLayout(row1);
    auto* row2 = new QHBoxLayout;
    auto* restore = button(QObject::tr("Restore settings"), "undo-2");
    auto* rename = button(QObject::tr("Label…"), "pencil");
    auto* del = button(QObject::tr("Delete…"), "trash-2");
    del->setProperty("flat", true);
    row2->addWidget(restore);
    row2->addWidget(rename);
    row2->addStretch(1);
    row2->addWidget(del);
    v->addLayout(row2);

    auto update = [c, selected, show, plot, comp, restore, rename, del, dv, detail] {
        const QStringList ids = selected();
        const RunRecord* r = ids.size() == 1 ? c.session->history().find(ids[0]) : nullptr;
        const bool done = r && (r->status == RunStatus::Completed || r->status == RunStatus::Stopped);
        show->setEnabled(done || (ids.size() == 2));
        show->setText(ids.size() == 2 ? QObject::tr("Side by side") : QObject::tr("Show"));
        plot->setEnabled(!ids.isEmpty());
        comp->setEnabled(ids.size() == 2);
        restore->setEnabled(r != nullptr);
        rename->setEnabled(r != nullptr);
        del->setEnabled(!ids.isEmpty() && !(r && r->status == RunStatus::Running));
        while (QLayoutItem* it = dv->takeAt(0)) {
            delete it->widget();
            delete it;
        }
        if (!r) return;
        auto* grid = new QFormLayout;
        grid->setHorizontalSpacing(12);
        grid->setVerticalSpacing(3);
        auto add = [&](const QString& k, const QString& v2) {
            auto* l = label(v2);
            l->setTextInteractionFlags(Qt::TextSelectableByMouse);
            grid->addRow(label(k, "muted"), l);
        };
        add(QObject::tr("Device"), r->device.isEmpty() ? r->backend : r->device);
        add(QObject::tr("Settings"), r->overrides.isEmpty() ? QObject::tr("as in the case") : QObject::tr("%n changed", nullptr, int(r->overrides.size())));
        for (auto it = r->overrides.begin(); it != r->overrides.end(); ++it)
            if (const Setting* s = c.session->settings()->find(it.key())) add("  " + settingLabel(*s), displayValue(*s, it.value()));
        if (r->approximate()) {
            QStringList names;
            for (const auto& i : c.session->caseSummary().issues)
                if (r->approximations.contains(i.approximation) && !names.contains(i.label)) names << i.label;
            add(QObject::tr("Approximations"), names.isEmpty() ? r->approximations.join(", ") : names.join(", "));
        }
        if (!r->dataPath.isEmpty()) add(QObject::tr("Initial data"), QFileInfo(r->dataPath).fileName());
        if (!r->error.isEmpty()) add(QObject::tr("Error"), r->error);
        add(QObject::tr("Folder"), QDir::toNativeSeparators(r->dir));
        auto* w = new QWidget;
        w->setLayout(grid);
        dv->addWidget(w);
        detail->setVisible(true);
    };
    QObject::connect(table->selectionModel(), &QItemSelectionModel::selectionChanged, page, update);
    update();
    QObject::connect(show, &QPushButton::clicked, page, [c, selected] {
        const QStringList ids = selected();
        if (ids.size() == 2) {
            c.graphics->setViewLayout(GraphicsArea::Layout::SideBySide);
            c.graphics->setLinkedCameras(true);
            c.graphics->show("run:" + ids[1], 0);
            c.graphics->show("run:" + ids[0], 1);
        } else if (ids.size() == 1) {
            c.graphics->show("run:" + ids[0]);
        }
    });
    QObject::connect(plot, &QPushButton::clicked, page, [c, selected] {
        QStringList ids = selected();
        std::reverse(ids.begin(), ids.end());
        c.charts->showRuns(ids);
        c.charts->showTab(ChartsPanel::Tab::Residuals);
    });
    QObject::connect(comp, &QPushButton::clicked, page, [c, selected] {
        const QStringList ids = selected();
        if (ids.size() != 2) return;
        QSettings().setValue("compare/a", ids[1]);
        QSettings().setValue("compare/b", ids[0]);
        c.navigate("runs/compare");
    });
    QObject::connect(restore, &QPushButton::clicked, page, [c, selected] {
        const QStringList ids = selected();
        if (ids.size() == 1) c.session->restoreRunSettings(ids[0]);
    });
    QObject::connect(rename, &QPushButton::clicked, page, [c, selected, rename] {
        const QStringList ids = selected();
        const RunRecord* r = ids.size() == 1 ? c.session->history().find(ids[0]) : nullptr;
        if (!r) return;
        bool ok = false;
        const QString note = QInputDialog::getText(rename->window(), QObject::tr("Label run"), QObject::tr("Label for %1:").arg(QObject::tr("Run %1").arg(r->number)),
                                                   QLineEdit::Normal, r->note, &ok);
        if (!ok) return;
        c.session->history().setNote(ids[0], note.trimmed());
        emit c.session->historyChanged();
    });
    QObject::connect(del, &QPushButton::clicked, page, [c, selected, del] {
        const QStringList ids = selected();
        if (ids.isEmpty()) return;
        const auto answer = QMessageBox::question(del->window(), QObject::tr("Delete runs"),
                                                  QObject::tr("Delete %n run(s) and their result files? This cannot be undone.", nullptr, int(ids.size())),
                                                  QMessageBox::Cancel | QMessageBox::Yes, QMessageBox::Cancel);
        if (answer != QMessageBox::Yes) return;
        for (const QString& id : ids) c.session->deleteRun(id);
    });
    return page;
}

QWidget* compare(const InspectorContext& c)
{
    QVBoxLayout* v;
    QWidget* page = makePage(QObject::tr("Compare Runs"), QObject::tr("Runs"), &v);
    if (!c.session->hasCase()) {
        v->addWidget(label(QObject::tr("Open a case to compare its runs."), "muted"));
        return page;
    }
    std::vector<const RunRecord*> done;
    for (const auto& r : c.session->history().runs())
        if (r.status == RunStatus::Completed || r.status == RunStatus::Stopped) done.push_back(&r);
    if (done.size() < 2) {
        v->addWidget(label(QObject::tr("Compare needs two finished runs. Change a setting and run again: the differences in settings, forces and fields show here."), "muted"));
        return page;
    }
    auto* pick = new QHBoxLayout;
    auto* a = new QComboBox;
    auto* b = new QComboBox;
    a->setAccessibleName(QObject::tr("Run A"));
    b->setAccessibleName(QObject::tr("Run B"));
    for (const auto* r : done) {
        a->addItem(r->title(), r->id);
        b->addItem(r->title(), r->id);
    }
    const QString sa = QSettings().value("compare/a").toString(), sb = QSettings().value("compare/b").toString();
    a->setCurrentIndex(a->findData(sa) >= 0 ? a->findData(sa) : a->count() - 1);
    b->setCurrentIndex(b->findData(sb) >= 0 ? b->findData(sb) : std::max(0, b->count() - 2));
    pick->addWidget(label("A", "muted"));
    pick->addWidget(a, 1);
    pick->addWidget(label("B", "muted"));
    pick->addWidget(b, 1);
    v->addLayout(pick);

    auto* content = new QWidget;
    auto* cv = new QVBoxLayout(content);
    cv->setContentsMargins(0, 0, 0, 0);
    cv->setSpacing(8);
    v->addWidget(content);

    auto build = [c, a, b, cv] {
        while (QLayoutItem* it = cv->takeAt(0)) {
            delete it->widget();
            delete it;
        }
        const RunRecord* ra = c.session->history().find(a->currentData().toString());
        const RunRecord* rb = c.session->history().find(b->currentData().toString());
        if (!ra || !rb) return;
        QSettings().setValue("compare/a", ra->id);
        QSettings().setValue("compare/b", rb->id);
        if (ra->id == rb->id) {
            cv->addWidget(label(QObject::tr("Choose two different runs."), "muted"));
            return;
        }
        cv->addWidget(section(QObject::tr("Settings that differ")));
        const auto diffs = diffRuns(*ra, *rb, c.session->settings()->settings());
        if (diffs.empty()) {
            cv->addWidget(label(QObject::tr("None. Both runs used the same settings."), "muted"));
        } else {
            auto* t = new QTableWidget(int(diffs.size()), 3);
            t->setHorizontalHeaderLabels({QObject::tr("Setting"), "A", "B"});
            t->verticalHeader()->hide();
            t->setEditTriggers(QAbstractItemView::NoEditTriggers);
            t->setShowGrid(false);
            t->setAccessibleName(QObject::tr("Settings that differ"));
            for (int i = 0; i < int(diffs.size()); ++i) {
                const auto& d = diffs[std::size_t(i)];
                const Setting* s = c.session->settings()->find(d.path);
                auto fmt = [&](const QJsonValue& x) { return s ? displayValue(*s, x) : (x.isString() ? x.toString() : QString::number(x.toDouble())); };
                t->setItem(i, 0, new QTableWidgetItem(s ? settingLabel(*s) + "  ·  " + s->owner : d.label));
                t->setItem(i, 1, new QTableWidgetItem(fmt(d.a)));
                t->setItem(i, 2, new QTableWidgetItem(fmt(d.b)));
            }
            t->resizeColumnsToContents();
            t->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
            t->setMinimumHeight(34 + t->rowCount() * 26);
            cv->addWidget(t);
        }

        cv->addWidget(section(QObject::tr("Forces")));
        auto fa = c.session->history().forces(ra->id), fb = c.session->history().forces(rb->id);
        if (!fa || !fb) {
            cv->addWidget(label(QObject::tr("One of the runs wrote no forces."), "muted"));
        } else {
            auto* t = new QTableWidget(0, 4);
            t->setHorizontalHeaderLabels({QObject::tr("Zone / report"), QObject::tr("A [N]"), QObject::tr("B [N]"), QObject::tr("A − B")});
            t->verticalHeader()->hide();
            t->setEditTriggers(QAbstractItemView::NoEditTriggers);
            t->setShowGrid(false);
            t->setAccessibleName(QObject::tr("Force differences"));
            QFont tab = t->font();
            tab.setFeature(QFont::Tag("tnum"), 1);
            t->setFont(tab);
            auto val = [](const ForceRow& r) { return r.kind == "report" && r.totalProjection ? *r.totalProjection : r.total.length(); };
            for (const auto& rowA : fa->rows) {
                auto it = std::find_if(fb->rows.begin(), fb->rows.end(), [&](const ForceRow& x) { return x.name == rowA.name && x.kind == rowA.kind; });
                if (it == fb->rows.end()) continue;
                const int r = t->rowCount();
                t->insertRow(r);
                const double va = val(rowA), vb = val(*it);
                t->setItem(r, 0, new QTableWidgetItem(rowA.kind == "report" ? rowA.name : rowA.name + QObject::tr(" (|F|)")));
                auto num = [](const QString& s) {
                    auto* x = new QTableWidgetItem(s);
                    x->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
                    return x;
                };
                t->setItem(r, 1, num(QLocale().toString(va, 'f', 3)));
                t->setItem(r, 2, num(QLocale().toString(vb, 'f', 3)));
                const double pct = vb != 0 ? (va - vb) / std::abs(vb) * 100 : 0;
                t->setItem(r, 3, num(QStringLiteral("%1 (%2%3 %)").arg(QLocale().toString(va - vb, 'f', 3), pct >= 0 ? "+" : "", QLocale().toString(pct, 'f', 1))));
            }
            t->resizeColumnsToContents();
            t->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
            t->setMinimumHeight(std::min(360, 34 + t->rowCount() * 26));
            cv->addWidget(t);
        }

        cv->addWidget(section(QObject::tr("Fields")));
        auto* row = new QHBoxLayout;
        auto* diff = button(QObject::tr("Show A − B"), "blend", "primary");
        auto* side = button(QObject::tr("Side by side"), "columns-2");
        auto* res = button(QObject::tr("Overlay residuals"), "chart-line");
        row->addWidget(diff);
        row->addWidget(side);
        row->addWidget(res);
        row->addStretch(1);
        cv->addLayout(row);
        cv->addWidget(label(QObject::tr("A − B is computed on surfaces with identical geometry; others are listed in the console and left out. Differences use a diverging map: white means agreement."), "faint"));
        const QString ia = ra->id, ib = rb->id;
        const QString titleA = ra->title(), titleB = rb->title();
        QObject::connect(diff, &QPushButton::clicked, diff, [c, ia, ib, titleA, titleB] {
            const QString key = QStringLiteral("diff:%1:%2").arg(ia, ib);
            c.graphics->addDataset(key, QStringLiteral("%1 − %2").arg(titleA, titleB));
            c.graphics->setViewLayout(GraphicsArea::Layout::SideBySide);
            c.graphics->setLinkedCameras(true);
            c.graphics->show("run:" + ia, 0);
            c.graphics->show(key, 1);
        });
        QObject::connect(side, &QPushButton::clicked, side, [c, ia, ib] {
            c.graphics->setViewLayout(GraphicsArea::Layout::SideBySide);
            c.graphics->setLinkedCameras(true);
            c.graphics->show("run:" + ia, 0);
            c.graphics->show("run:" + ib, 1);
        });
        QObject::connect(res, &QPushButton::clicked, res, [c, ia, ib] {
            c.charts->showRuns({ia, ib});
            c.charts->showTab(ChartsPanel::Tab::Residuals);
        });
    };
    QObject::connect(a, &QComboBox::activated, page, build);
    QObject::connect(b, &QComboBox::activated, page, build);
    build();
    return page;
}

// ---- All settings -----------------------------------------------------------------

namespace {

class SettingsTableModel : public QAbstractTableModel {
public:
    explicit SettingsTableModel(SettingsModel* m, QObject* parent) : QAbstractTableModel(parent), m_model(m)
    {
        QObject::connect(m, &SettingsModel::valueChanged, this, [this](const QString& path) {
            const auto& all = m_model->settings();
            for (int i = 0; i < int(all.size()); ++i)
                if (all[std::size_t(i)].path == path) emit dataChanged(index(i, 0), index(i, 3));
        });
    }
    int rowCount(const QModelIndex& p = {}) const override { return p.isValid() ? 0 : int(m_model->settings().size()); }
    int columnCount(const QModelIndex& = {}) const override { return 4; }
    QVariant headerData(int s, Qt::Orientation o, int role) const override
    {
        if (o != Qt::Horizontal || role != Qt::DisplayRole) return {};
        static const char* h[] = {"Setting", "Value", "Unit", "Status"};
        return QObject::tr(h[s]);
    }
    QVariant data(const QModelIndex& idx, int role) const override
    {
        const Setting& s = m_model->settings()[std::size_t(idx.row())];
        const bool mod = m_model->isModified(s.path);
        if (role == Qt::DisplayRole) {
            switch (idx.column()) {
            case 0: return s.path;
            case 1: {
                QString t = displayValue(s, m_model->value(s.path));
                const DisplayUnit du = displayUnit(s);
                if (!du.unit.isEmpty() && t.endsWith(" " + du.unit)) t.chop(du.unit.size() + 1);
                return t.left(160);
            }
            case 2: return displayUnit(s).unit;
            case 3: return mod ? QObject::tr("Modified") : s.editable ? QObject::tr("Editable") : s.reason;
            }
        }
        if (role == Qt::ToolTipRole) return QStringLiteral("%1\n%2").arg(s.path, s.editable ? QObject::tr("Editable") : s.reason);
        if (role == Qt::ForegroundRole) {
            const auto& t = tokens();
            if (idx.column() == 3) return mod ? (t.dark ? t.amber : t.amberText) : s.editable ? t.text : t.faint;
            if (idx.column() == 0 && !s.editable) return t.muted;
        }
        if (role == Qt::UserRole) return s.path;
        if (role == Qt::UserRole + 1) return QStringLiteral("%1 %2 %3").arg(s.path, s.label, s.group);
        return {};
    }

private:
    SettingsModel* m_model;
};

class SettingsFilter : public QSortFilterProxyModel {
public:
    using QSortFilterProxyModel::QSortFilterProxyModel;
    QString text;
    bool editableOnly = false, modifiedOnly = false;
    SettingsModel* model = nullptr;
    void apply() { beginFilterChange(); endFilterChange(); }

protected:
    bool filterAcceptsRow(int row, const QModelIndex&) const override
    {
        const Setting& s = model->settings()[std::size_t(row)];
        if (editableOnly && !s.editable) return false;
        if (modifiedOnly && !model->isModified(s.path)) return false;
        if (text.isEmpty()) return true;
        const QString hay = sourceModel()->index(row, 0).data(Qt::UserRole + 1).toString();
        for (const QString& w : text.split(' ', Qt::SkipEmptyParts))
            if (!hay.contains(w, Qt::CaseInsensitive)) return false;
        return true;
    }
};

// The inspector page that holds an editable setting.
QString pageFor(const Setting& s)
{
    if (s.section == "boundary") return "setup/boundary/" + s.owner;
    if (s.section == "materials") return "setup/material/" + s.owner;
    static const QHash<QString, QString> m = {{"reference", "setup/reference"}, {"methods", "solution/methods"}, {"controls", "solution/controls"},
                                              {"run", "solution/run"}, {"models", "setup/models"}, {"cell-zone", "setup/cell-zones"},
                                              {"monitors", "solution/monitors"}, {"initialization", "solution/initialization"}};
    return m.value(s.section, "setup/general");
}

} // namespace

QWidget* allSettings(const InspectorContext& c, const QString& filter)
{
    auto* page = new QWidget;
    auto* v = new QVBoxLayout(page);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    v->addWidget(header(QObject::tr("All Case Settings"), {}, page));
    auto* content = new QWidget(page);
    auto* cv = body(content);
    v->addWidget(content, 1);
    SettingsModel* sm = c.session->settings();
    const int total = int(sm->settings().size());
    int editable = 0;
    for (const auto& s : sm->settings()) editable += s.editable;
    cv->addWidget(label(QObject::tr("Every setting the case contains: %1 in total, %2 editable in Kutta. Read-only settings show why; Kutta reads them to decide whether it can reproduce the case.")
                            .arg(QLocale().toString(total), QLocale().toString(editable)), "muted"));
    auto* search = new QLineEdit(filter);
    search->setPlaceholderText(QObject::tr("Filter by path or name, e.g. inlet turb"));
    search->setClearButtonEnabled(true);
    search->setAccessibleName(QObject::tr("Filter settings"));
    cv->addWidget(search);
    auto* opts = new QHBoxLayout;
    auto* ed = new QCheckBox(QObject::tr("Editable only"));
    auto* mo = new QCheckBox(QObject::tr("Modified only"));
    auto* countLabel = label({}, "faint");
    opts->addWidget(ed);
    opts->addWidget(mo);
    opts->addStretch(1);
    opts->addWidget(countLabel);
    cv->addLayout(opts);
    auto* view = new QTableView;
    auto* model = new SettingsTableModel(sm, view);
    auto* proxy = new SettingsFilter(view);
    proxy->model = sm;
    proxy->setSourceModel(model);
    proxy->text = filter;
    view->setModel(proxy);
    view->verticalHeader()->hide();
    view->verticalHeader()->setDefaultSectionSize(24);
    view->setSelectionBehavior(QAbstractItemView::SelectRows);
    view->setShowGrid(false);
    view->setWordWrap(false);
    view->setAccessibleName(QObject::tr("Case settings"));
    view->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    view->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Interactive);
    view->setColumnWidth(1, 150);
    view->setColumnWidth(2, 60);
    view->setColumnWidth(3, 170);
    view->setMinimumHeight(420);
    view->setTextElideMode(Qt::ElideMiddle);
    cv->addWidget(view, 1);
    auto refresh = [proxy, search, ed, mo, countLabel] {
        proxy->text = search->text().trimmed();
        proxy->editableOnly = ed->isChecked();
        proxy->modifiedOnly = mo->isChecked();
        proxy->apply();
        countLabel->setText(QObject::tr("%1 shown").arg(QLocale().toString(proxy->rowCount())));
    };
    auto* debounce = new QTimer(page);
    debounce->setSingleShot(true);
    debounce->setInterval(150);
    QObject::connect(debounce, &QTimer::timeout, page, refresh);
    QObject::connect(search, &QLineEdit::textChanged, debounce, qOverload<>(&QTimer::start));
    QObject::connect(ed, &QCheckBox::toggled, page, refresh);
    QObject::connect(mo, &QCheckBox::toggled, page, refresh);
    QObject::connect(view, &QTableView::activated, page, [c, proxy, sm](const QModelIndex& idx) {
        const QString path = proxy->data(proxy->index(idx.row(), 0), Qt::UserRole).toString();
        if (const Setting* s = sm->find(path); s && s->editable) c.navigate(pageFor(*s));
    });
    refresh();
    return page;
}

} // namespace kapp::pages
