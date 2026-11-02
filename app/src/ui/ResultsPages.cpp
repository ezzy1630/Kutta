// Inspector pages for Results, Runs and All Case Settings.

#include "core/ChartScale.h"
#include "ui/ChartsPanel.h"
#include "ui/GraphicsArea.h"
#include "ui/Inspector.h"
#include "ui/PageKit.h"
#include "ui/SettingWidgets.h"

#include <QButtonGroup>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHeaderView>
#include <QInputDialog>
#include <QJsonArray>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPointer>
#include <QRadioButton>
#include <QSettings>
#include <QSortFilterProxyModel>
#include <QSpinBox>
#include <QStandardItemModel>
#include <QTableWidget>
#include <QTreeView>
#include <QTimer>

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

// A table exactly as tall as its rows (up to maxRows, then it scrolls).
void fitRows(QTableWidget* t, int maxRows)
{
    const int rows = std::min(t->rowCount(), maxRows);
    const int rowH = t->verticalHeader()->defaultSectionSize();
    t->setFixedHeight(t->horizontalHeader()->sizeHint().height() + rows * rowH + 2 * t->frameWidth() + 2);
}

std::optional<double> numberOf(QLineEdit* e)
{
    bool ok = false;
    double v = QLocale::c().toDouble(e->text().trimmed(), &ok);
    if (!ok) v = QLocale().toDouble(e->text().trimmed(), &ok);
    return ok && std::isfinite(v) ? std::optional<double>(v) : std::nullopt;
}

struct Vec3Edit {
    QWidget* widget;
    std::array<QLineEdit*, 3> e;
    std::optional<QVector3D> value() const;
    void set(const QVector3D& v) const
    {
        for (int i = 0; i < 3; ++i) e[std::size_t(i)]->setText(QString::number(v[i], 'g', 6));
    }
};
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
    addGroup(QObject::tr("Probes"), [](const SceneItem& i) { return i.kind == ItemKind::Probe; });
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

void runPost(const InspectorContext& c, PostForm pf, PostKind kind, const QJsonObject& op, std::function<void(const PostOutput&)> onDone)
{
    const QString run = currentRunId(c);
    if (run.isEmpty()) return reportError(pf.error, QObject::tr("Run the calculation first."));
    reportError(pf.error, {});
    pf.go->setEnabled(false);
    pf.status->setText(QObject::tr("Computing in kutta…"));
    pf.status->show();
    QPointer<QPushButton> go = pf.go;
    QPointer<QLabel> err = pf.error, st = pf.status;
    c.session->postProcess(run, kind, op, [go, err, st, onDone](Result<PostOutput> r) {
        if (go) go->setEnabled(true);
        if (st) st->hide();
        if (!r) {
            if (err) reportError(err, r.error());
            return;
        }
        onDone(r.value());
    });
}

// The preferred quantities kutta can sample for this run. Before the run's
// view is first published the list is unknown and the preferences stand.
QStringList sampleable(const InspectorContext& c, const QStringList& preferred, bool volumeOnly)
{
    const auto avail = c.session->availableFields(currentRunId(c));
    if (avail.empty()) return preferred;
    QStringList out;
    for (const QString& q : preferred)
        if (std::any_of(avail.begin(), avail.end(), [&](const AvailableField& a) { return a.name == q && (a.volume || !volumeOnly); }))
            out << q;
    return out;
}

bool hasFluentField(const InspectorContext& c, const QString& q)
{
    const auto avail = c.session->availableFields(currentRunId(c));
    return std::any_of(avail.begin(), avail.end(), [&](const AvailableField& a) { return a.name == q + "@fluent"; });
}

QJsonArray json(const QVector3D& v) { return {v.x(), v.y(), v.z()}; }

std::optional<QVector3D> vec(const QJsonValue& v)
{
    const QJsonArray a = v.toArray();
    if (a.size() != 3) return std::nullopt;
    return QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble()));
}

// The current run's previous request from this form, so the form reopens as left.
QJsonObject lastOperation(const InspectorContext& c, const QString& name)
{
    return c.session->postOperation(currentRunId(c), name);
}

void selectData(QComboBox* box, const QJsonValue& v)
{
    if (const int i = box->findData(v.toString()); i >= 0) box->setCurrentIndex(i);
}

} // namespace

QWidget* pathlines(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("Pathlines"), QObject::tr("Results ▸ Graphics"), QObject::tr("Compute pathlines"));
    pf.body->addWidget(label(QObject::tr("Streamlines of the solved velocity field, released along a line or from a boundary."), "muted"));
    const CaseSummary* cs = c.session->hasCase() ? &c.session->caseSummary() : nullptr;
    QVector3D from(-0.6f, -0.5f, 0.12f), to(-0.6f, 0.5f, 0.12f);
    if (cs && cs->bounds.valid && !std::any_of(cs->zones.begin(), cs->zones.end(), [](const Zone& z) { return z.name.contains("ground"); })) {
        const QVector3D a = cs->bounds.min, b = cs->bounds.max, m = cs->bounds.center();
        from = QVector3D(a.x(), m.y(), m.z());
        to = QVector3D(a.x(), b.y(), m.z());
    }
    auto* f = formLayout();
    auto* release = new QComboBox;
    release->setAccessibleName(QObject::tr("Release from"));
    release->addItem(QObject::tr("Line (rake)"), QString());
    if (cs)
        for (const auto& z : cs->zones)
            if (z.type.contains("inlet")) release->addItem(z.name, z.name);
    auto rakeFrom = vec3Edit(from, QObject::tr("Rake start"));
    auto rakeTo = vec3Edit(to, QObject::tr("Rake end"));
    auto* count = new QSpinBox;
    count->setRange(1, 2000);
    count->setValue(40);
    count->setAccessibleName(QObject::tr("Number of pathlines"));
    f->addRow(QObject::tr("Release from"), release);
    f->addRow(QObject::tr("Rake start (m)"), rakeFrom.widget);
    f->addRow(QObject::tr("Rake end (m)"), rakeTo.widget);
    f->addRow(QObject::tr("Pathlines"), count);
    QObject::connect(release, &QComboBox::currentIndexChanged, release, [release, rakeFrom, rakeTo] {
        const bool rake = release->currentData().toString().isEmpty();
        rakeFrom.widget->setEnabled(rake);
        rakeTo.widget->setEnabled(rake);
    });
    if (const QJsonObject last = lastOperation(c, QObject::tr("Pathlines")); !last.isEmpty()) {
        const QJsonObject seeds = last.value("seeds").toObject();
        count->setValue(seeds.value("count").toInt(count->value()));
        if (seeds.value("type").toString() == "surface") selectData(release, seeds.value("zones").toArray().first());
        if (auto a = vec(seeds.value("start"))) rakeFrom.set(*a);
        if (auto b = vec(seeds.value("end"))) rakeTo.set(*b);
    }
    pf.body->addLayout(f);

    // How the lines look in the active viewport (Fluent's pathline display
    // options): the colouring field, width and direction arrows.
    pf.body->addWidget(section(QObject::tr("Display")));
    auto* df = formLayout();
    auto* colour = new QComboBox;
    colour->setAccessibleName(QObject::tr("Colour pathlines by"));
    const QStringList lineFields = sampleable(c, {"velocity_magnitude", "pressure", "pressure_coefficient", "total_pressure_coefficient"}, true);
    for (const QString& q : lineFields) colour->addItem(fluentLabel(q), q);
    colour->addItem(QObject::tr("Solid colour"), QString());
    auto* width = new QSpinBox;
    width->setRange(1, 10);
    width->setSuffix(QObject::tr(" px"));
    width->setAccessibleName(QObject::tr("Line width"));
    auto* arrows = new QCheckBox(QObject::tr("Arrows along the flow"));
    if (auto* vp = c.graphics->active()) {
        const ViewDisplay& d = vp->display();
        width->setValue(int(std::lround(d.ribbonWidth)));
        arrows->setChecked(d.lineArrows);
        colour->setCurrentIndex(d.lineSolid ? colour->count() - 1 : std::max(0, colour->findData(d.field.quantity)));
    }
    df->addRow(QObject::tr("Colour by"), colour);
    df->addRow(QObject::tr("Line width"), width);
    df->addRow({}, arrows);
    pf.body->addLayout(df);
    // Also runs when a computation finishes, possibly after the page is gone.
    auto applyStyle = [c, colour = QPointer<QComboBox>(colour), width = QPointer<QSpinBox>(width), arrows = QPointer<QCheckBox>(arrows)] {
        if (!colour || !width || !arrows) return;
        editDisplay(c, [&](ViewDisplay& d) {
            const QString q = colour->currentData().toString();
            d.lineSolid = q.isEmpty();
            if (!q.isEmpty()) {
                d.field = {q, FieldRole::Kutta, -1};
                d.contours = true;
            }
            d.ribbonWidth = float(width->value());
            d.lineArrows = arrows->isChecked();
        });
    };
    QObject::connect(colour, &QComboBox::activated, colour, applyStyle);
    QObject::connect(width, &QSpinBox::valueChanged, width, applyStyle);
    QObject::connect(arrows, &QCheckBox::toggled, arrows, applyStyle);

    finishPostPage(pf);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, release, rakeFrom, rakeTo, count, applyStyle] {
        const QString zone = release->currentData().toString();
        QJsonObject seeds;
        if (zone.isEmpty()) {
            auto a = rakeFrom.value(), b = rakeTo.value();
            if (!a || !b) return reportError(pf.error, QObject::tr("Enter numbers for both rake points."));
            seeds = {{"type", "rake"}, {"start", json(*a)}, {"end", json(*b)}, {"count", count->value()}};
        } else {
            seeds = {{"type", "surface"}, {"zones", QJsonArray{zone}}, {"count", count->value()}, {"distribution", "uniform"}, {"offset", 0}};
        }
        const QStringList fields = sampleable(c, {"velocity_magnitude", "pressure", "pressure_coefficient", "total_pressure_coefficient"}, true);
        const QJsonObject op{{"name", QObject::tr("Pathlines")}, {"mode", "streamline"}, {"seeds", seeds}, {"direction", "forward"},
                             {"fields", QJsonArray::fromStringList(fields)}};
        runPost(c, pf, PostKind::Streamlines, op, [applyStyle](const PostOutput&) { applyStyle(); });
    });
    return pf.page;
}

QWidget* isoSurfaces(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("Iso-Surfaces"), QObject::tr("Results ▸ Graphics"), QObject::tr("Compute iso-surface"));
    auto* f = formLayout();
    auto* field = new QComboBox;
    field->setAccessibleName(QObject::tr("Field"));
    for (const QString& q : sampleable(c, {"total_pressure_coefficient", "q_criterion", "velocity_magnitude", "velocity_x", "pressure_coefficient",
                                           "pressure", "vorticity_magnitude", "viscosity_ratio"}, true))
        field->addItem(fluentLabel(q), q);
    auto* value = numberEdit(0, QObject::tr("Iso-value"));
    f->addRow(QObject::tr("Field"), field);
    f->addRow(QObject::tr("Value"), value);
    pf.body->addLayout(f);
    pf.body->addWidget(label(QObject::tr("Total pressure coefficient = 0 outlines the wake: where the flow has lost its free-stream total pressure. "
                                         "A positive Q criterion shows vortex cores."), "faint"));
    finishPostPage(pf);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, field, value] {
        auto val = numberOf(value);
        if (!val) return reportError(pf.error, QObject::tr("Enter a number."));
        const QString q = field->currentData().toString();
        QStringList extra = sampleable(c, {"velocity_magnitude", "pressure"}, true);
        extra.removeAll(q);
        const QJsonObject op{{"name", QStringLiteral("%1 = %2").arg(field->currentText(), formatNumber(*val))}, {"field", q}, {"value", *val},
                             {"fields", QJsonArray::fromStringList(extra)}};
        runPost(c, pf, PostKind::IsoSurface, op, [](const PostOutput&) {});
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
    if (c.session->hasCase() && (!p || !p->scene() || !p->scene()->geometryOnly)) {
        // The case's mesh is always one click away, run or not.
        auto* show = button(QObject::tr("Show the case mesh"), "box");
        QObject::connect(show, &QPushButton::clicked, show, [c] { c.graphics->show(c.session->meshSceneKey()); });
        v->addWidget(show, 0, Qt::AlignLeft);
    }
    if (!p || !p->scene()) return page;
    const ViewDisplay d = p->display();
    const bool meshOnly = p->scene()->geometryOnly;
    auto* edges = new QCheckBox(QObject::tr("Edges"));
    edges->setChecked(d.edges);
    QObject::connect(edges, &QCheckBox::toggled, edges, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.edges = on; }); });
    v->addWidget(edges);
    if (!meshOnly) {
        auto* field = new QCheckBox(QObject::tr("Colour by field"));
        field->setChecked(d.contours);
        QObject::connect(field, &QCheckBox::toggled, field, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.contours = on; }); });
        v->addWidget(field);
    }
    auto* zones = new QCheckBox(meshOnly ? QObject::tr("Colour by zone") : QObject::tr("Colour by zone when not colouring by field"));
    zones->setChecked(d.zoneColors);
    QObject::connect(zones, &QCheckBox::toggled, zones, [c](bool on) { editDisplay(c, [&](ViewDisplay& d) { d.zoneColors = on; }); });
    v->addWidget(zones);
    std::size_t tris = 0;
    bool faceEdges = true;
    for (const auto& it : p->scene()->items)
        if (it.mesh && it.kind == ItemKind::Zone) {
            tris += it.mesh->triangleCount();
            faceEdges = faceEdges && it.mesh->hasFaceEdges();
        }
    v->addWidget(label(faceEdges ? QObject::tr("%1 surface triangles. Edges are the mesh's own face edges.").arg(QLocale().toString(qulonglong(tris)))
                                 : QObject::tr("%1 surface triangles. Edges include the diagonals that split polygon faces into triangles (this view predates face edges).")
                                       .arg(QLocale().toString(qulonglong(tris))), "faint"));
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
    samples->setRange(2, 100000);
    samples->setValue(200);
    samples->setAccessibleName(QObject::tr("Samples"));
    auto* field = new QComboBox;
    field->setAccessibleName(QObject::tr("Field"));
    for (const QString& q : sampleable(c, {"velocity_magnitude", "velocity_x", "pressure", "pressure_coefficient", "total_pressure_coefficient",
                                           "vorticity_magnitude", "k", "omega"}, true))
        field->addItem(fluentLabel(q), q);
    if (const QJsonObject last = lastOperation(c, QObject::tr("XY line")); !last.isEmpty()) {
        if (auto p = vec(last.value("start"))) a.set(*p);
        if (auto p = vec(last.value("end"))) b.set(*p);
        samples->setValue(last.value("count").toInt(samples->value()));
        selectData(field, last.value("fields").toArray().first());
    }
    f->addRow(QObject::tr("Start (m)"), a.widget);
    f->addRow(QObject::tr("End (m)"), b.widget);
    f->addRow(QObject::tr("Samples"), samples);
    f->addRow(QObject::tr("Field"), field);
    pf.body->addLayout(f);
    pf.body->addWidget(label(QObject::tr("Fluent's values are plotted dashed when the case was opened with its data file."), "faint"));
    finishPostPage(pf);
    QObject::connect(pf.go, &QPushButton::clicked, pf.page, [c, pf, a, b, samples, field] {
        auto pa = a.value(), pb = b.value();
        if (!pa || !pb) return reportError(pf.error, QObject::tr("Enter numbers for both points."));
        const QString q = field->currentData().toString();
        QJsonArray fields{q};
        if (hasFluentField(c, q)) fields.append(q + "@fluent");
        const QJsonObject op{{"name", QObject::tr("XY line")}, {"type", "line"}, {"start", json(*pa)}, {"end", json(*pb)},
                             {"count", samples->value()}, {"fields", fields}};
        const QString runId = currentRunId(c);
        runPost(c, pf, PostKind::Probe, op, [c, q, runId](const PostOutput& out) {
            const RunRecord* r = c.session->history().find(runId);
            c.charts->plotProbe(out, QObject::tr("%1 along line · %2").arg(fluentLabel(q), r ? r->title() : QString()), -1);
        });
    });

    // Section through surfaces: Cp around a wing section, the elbow's wall pressure.
    pf.body->addSpacing(8);
    pf.body->addWidget(section(QObject::tr("Section through surfaces")));
    auto* sf = formLayout();
    // Default: the plane the viewports slice by (the car's centre plane, the elbow's mid-depth).
    const std::vector<Plane> planes = c.session->hasCase() ? c.session->defaultSlices() : std::vector<Plane>{};
    auto originEdit = vec3Edit(planes.empty() ? QVector3D() : planes[0].origin, QObject::tr("Section origin"));
    auto normalEdit = vec3Edit(planes.empty() ? QVector3D(0, 1, 0) : planes[0].normal, QObject::tr("Section normal"));
    auto* sfield = new QComboBox;
    for (const QString& q : sampleable(c, {"pressure_coefficient", "pressure", "skin_friction", "wall_shear_magnitude", "y_plus"}, false))
        sfield->addItem(fluentLabel(q), q);
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
    const QJsonObject lastSection = lastOperation(c, QObject::tr("Section"));
    const QJsonArray lastZones = lastSection.value("zones").toArray();
    for (int i = 0; i < zones->count(); ++i)
        if (lastZones.contains(zones->item(i)->text())) zones->item(i)->setCheckState(Qt::Checked);
    if (lastZones.isEmpty() && c.session->hasCase()) {
        // The largest wall: usually the body of interest.
        const Zone* biggest = nullptr;
        for (const auto& z : c.session->caseSummary().zones)
            if (z.type == "wall" && !z.name.contains("ground") && (!biggest || z.faces > biggest->faces)) biggest = &z;
        for (int i = 0; i < zones->count(); ++i)
            zones->item(i)->setCheckState(biggest && zones->item(i)->text() == biggest->name ? Qt::Checked : Qt::Unchecked);
    }
    if (!lastSection.isEmpty()) {
        if (auto p = vec(lastSection.value("origin"))) originEdit.set(*p);
        if (auto p = vec(lastSection.value("normal"))) normalEdit.set(*p);
        selectData(sfield, lastSection.value("fields").toArray().first());
    }
    zones->setMaximumHeight(std::min(130, 8 + zones->count() * 28));
    auto* xAxis = new QComboBox;
    xAxis->addItem(QObject::tr("x"), 0);
    xAxis->addItem(QObject::tr("y"), 1);
    xAxis->addItem(QObject::tr("z"), 2);
    xAxis->addItem(QObject::tr("Distance along curve"), -1);
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
        QJsonArray fields{q};
        if (hasFluentField(c, q)) fields.append(q + "@fluent");
        const QJsonObject op{{"name", QObject::tr("Section")}, {"type", "section"}, {"zones", zl}, {"origin", json(*o)}, {"normal", json(*n)},
                             {"fields", fields}};
        PostForm sp = pf;
        sp.go = secGo;
        sp.error = secErr;
        const int axis = xAxis->currentData().toInt();
        runPost(c, sp, PostKind::Probe, op, [c, q, axis](const PostOutput& out) {
            c.charts->plotProbe(out, QObject::tr("%1 on section").arg(fluentLabel(q)), axis);
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
    if (r->approximate()) top->addWidget(badge(r->approximateBadge(), "approximate"));
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
    table->setHorizontalHeaderLabels({QObject::tr("Zone / report"), QObject::tr("Pressure [N]"), QObject::tr("Viscous [N]"), QObject::tr("Total [N]"), QObject::tr("Coeff.")});
    table->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    table->setWordWrap(false);
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
            auto* it = new QTableWidgetItem(QLocale().toString(x, 'f', 2));
            it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
            return it;
        };
        for (const auto& row : rows) {
            const int r = table->rowCount();
            table->insertRow(r);
            const bool report = row.kind == "report";
            const double p = report && row.pressureProjection ? *row.pressureProjection : row.pressure[ax];
            const double vv = report && row.viscousProjection ? *row.viscousProjection : row.viscous[ax];
            auto* name = new QTableWidgetItem(report ? row.name + QObject::tr(" (report)") : row.name);
            name->setToolTip(name->text());
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
        // Names get the room the numbers leave; numbers never elide.
        table->resizeColumnsToContents();
        for (int col = 1; col < table->columnCount(); ++col) table->horizontalHeader()->setSectionResizeMode(col, QHeaderView::ResizeToContents);
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
    v->addWidget(label(QObject::tr("Report rows are projected on the report's direction; wall rows on the axis chosen above. Coefficients, groups and balance are under Aerodynamics."), "faint"));
    return page;
}

namespace {

// One surface-report value from a reports.json zone or total entry.
std::optional<double> reportValue(const QJsonObject& z, const QString& report, const QString& field)
{
    auto num = [](const QJsonValue& v) { return v.isDouble() ? std::optional<double>(v.toDouble()) : std::nullopt; };
    if (report == "area") return num(z.value("area"));
    if (report == "mass_flow") return num(z.value("mass_flow"));
    const QJsonObject a = z.value("averages").toObject().value(field).toObject();
    return num(a.value(report));
}

} // namespace

QWidget* integrals(const InspectorContext& c)
{
    PostForm pf = postPage(c, QObject::tr("Surface Integrals"), QObject::tr("Results ▸ Reports"), QObject::tr("Compute"));
    auto* f = formLayout();
    auto* report = new QComboBox;
    report->addItem(QObject::tr("Area-weighted average"), "area_weighted_average");
    report->addItem(QObject::tr("Mass-weighted average"), "mass_flow_weighted_average");
    report->addItem(QObject::tr("Integral"), "area_integral");
    report->addItem(QObject::tr("Area"), "area");
    report->addItem(QObject::tr("Mass flow rate"), "mass_flow");
    report->setAccessibleName(QObject::tr("Report type"));
    auto* field = new QComboBox;
    field->setAccessibleName(QObject::tr("Field"));
    for (const QString& q : sampleable(c, {"pressure", "total_pressure", "pressure_coefficient", "total_pressure_coefficient", "velocity_magnitude",
                                           "wall_shear_magnitude", "skin_friction", "y_plus"}, false))
        field->addItem(fluentLabel(q), q);
    QObject::connect(report, &QComboBox::currentIndexChanged, field, [report, field] {
        const QString r = report->currentData().toString();
        field->setEnabled(r != "area" && r != "mass_flow");
    });
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
    zones->setFixedHeight(std::min(220, 6 + zones->count() * 28));
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
        const QString kind = report->currentData().toString();
        const bool fluent = withFluent->isChecked() && field->isEnabled();
        QJsonArray fields{q};
        if (fluent) fields.append(q + "@fluent");
        const QJsonObject op{{"name", QObject::tr("Surface integrals")}, {"zones", zl}, {"fields", fields}};
        runPost(c, pf, PostKind::Integral, op, [table, note, q, kind, fluent](const PostOutput& out) {
            const QJsonObject rep = out.report;
            const QJsonArray rows = rep.value("zones").toArray();
            const QJsonObject units = rep.value("units").toObject();
            QString unit = kind == "area" ? units.value("area").toString() : kind == "mass_flow" ? units.value("mass_flow").toString()
                                                                                                 : rep.value("total").toObject().value("averages").toObject().value(q).toObject().value("unit").toString();
            if (kind == "area_integral") unit = unit == "1" ? QStringLiteral("m^2") : unit + " m^2";
            unit = prettyUnit(unit);
            auto cell = [](std::optional<double> v) {
                auto* it = new QTableWidgetItem(v ? formatNumber(*v, 6) : QStringLiteral("—"));
                it->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
                return it;
            };
            table->setRowCount(int(rows.size()) + 1);
            table->setColumnHidden(2, !fluent);
            for (int i = 0; i < rows.size(); ++i) {
                const QJsonObject z = rows[i].toObject();
                table->setItem(i, 0, new QTableWidgetItem(z.value("zone").toString()));
                table->setItem(i, 1, cell(reportValue(z, kind, q)));
                table->setItem(i, 2, cell(reportValue(z, kind, q + "@fluent")));
            }
            const QJsonObject total = rep.value("total").toObject();
            auto* tot = new QTableWidgetItem(QObject::tr("Net"));
            QFont bf = tot->font();
            bf.setBold(true);
            tot->setFont(bf);
            table->setItem(int(rows.size()), 0, tot);
            table->setItem(int(rows.size()), 1, cell(reportValue(total, kind, q)));
            table->setItem(int(rows.size()), 2, cell(reportValue(total, kind, q + "@fluent")));
            table->setHorizontalHeaderItem(1, new QTableWidgetItem(unit.isEmpty() ? QObject::tr("Kutta") : QObject::tr("Kutta [%1]").arg(unit)));
            table->setHorizontalHeaderItem(2, new QTableWidgetItem(unit.isEmpty() ? QObject::tr("Fluent") : QObject::tr("Fluent [%1]").arg(unit)));
            table->resizeColumnsToContents();
            table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
            table->setMinimumHeight(34 + table->rowCount() * 26);
            table->show();
            QString text = QObject::tr("Discrete face integrals of the saved solution, computed by kutta export-view.");
            if (kind == "mass_flow") text += " " + QObject::tr("Mass flow is positive out of the domain.");
            QStringList marks;
            if (out.provenance.contains("APPROXIMATE")) marks << "APPROXIMATE";
            if (out.provenance.contains("MODIFIED")) marks << "MODIFIED";
            if (!marks.isEmpty()) text += " " + QObject::tr("These values are %1; the console and this note's tooltip give the details.").arg(marks.join(" and "));
            note->setText(text);
            note->setToolTip(out.provenance);
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
        auto* st = new QTableWidgetItem(runStatusText(r) + (r.status == RunStatus::Stopped && r.stopReason == "stopped" ? QObject::tr(" by you") : QString()));
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
        if (r->status == RunStatus::Completed || r->status == RunStatus::Stopped) {
            const auto summary = c.session->history().summary(r->id);
            const QString why = stopReasonText(r->stopReason, r->convergedBy, summary ? summary->convergence : ConvergenceSetup{});
            if (!why.isEmpty()) add(QObject::tr("Stopped"), why);
            if (const QString poor = summary ? poorMeshNumericsText(*summary) : QString(); !poor.isEmpty()) add(QObject::tr("Mesh numerics"), poor);
        }
        if (!r->coupling.isEmpty()) add(QObject::tr("Coupling"), couplingLabel(r->coupling));
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
        auto* folder = pathLabel(r->dir);
        folder->setTextInteractionFlags(Qt::TextSelectableByMouse);
        grid->addRow(label(QObject::tr("Folder"), "muted"), folder);
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
            fitRows(t, 8);
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
            fitRows(t, 12);
            cv->addWidget(t);
        }

        cv->addWidget(section(QObject::tr("Fields")));
        auto* row = new QHBoxLayout;
        auto* diff = button(QObject::tr("Show A − B"), "blend", "primary");
        auto* side = button(QObject::tr("Side by side"), "columns-2");
        auto* res = button(QObject::tr("Residuals"), "chart-line");
        res->setToolTip(QObject::tr("Overlay both runs' residuals"));
        row->addWidget(diff);
        row->addWidget(side);
        row->addWidget(res);
        row->addStretch(1);
        cv->addLayout(row);
        cv->addWidget(label(QObject::tr("A − B is computed on surfaces with identical geometry; others are listed in the console and left out. Differences use a diverging map: white means agreement."), "faint"));
        cv->addStretch(1);
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

// Leaf rows carry the setting path; group rows carry none.
constexpr int kPathRole = Qt::UserRole;
constexpr int kSearchRole = Qt::UserRole + 1;
constexpr int kRetainedRole = Qt::UserRole + 2;  // the raw-settings folder

// Settings grouped as Fluent's outline: group, then zone or material, then
// Kutta's controls before the retained raw Fluent settings, which sit in their
// own collapsed folder so 7,000 rp entries never bury the 26 that matter.
QString unusedReason(const SettingsModel& m, const Setting& s);

class SettingsTree : public QStandardItemModel {
public:
    SettingsTree(SettingsModel* m, QObject* parent) : QStandardItemModel(parent), m_model(m)
    {
        setHorizontalHeaderLabels({QObject::tr("Setting"), QObject::tr("Value"), QObject::tr("Unit")});
        static const QStringList order{"General", "Models", "Materials", "Cell Zone Conditions", "Boundary Conditions", "Reference Values",
                                       "Methods", "Controls", "Report Definitions", "Monitors", "Initialization", "Run Calculation"};
        std::map<int, std::vector<const Setting*>> byGroup;
        QStringList extra;
        for (const auto& st : m->settings()) {
            qsizetype g = order.indexOf(st.group);
            if (g < 0) {
                if (!extra.contains(st.group)) extra << st.group;
                g = order.size() + extra.indexOf(st.group);
            }
            byGroup[int(g)].push_back(&st);
        }
        for (const auto& [g, list] : byGroup) {
            const QString name = g < order.size() ? order[g] : extra[g - order.size()];
            int editable = 0;
            for (const Setting* st : list) editable += st->editable;
            QStandardItem* group = folder(QStringLiteral("%1   %2").arg(name, counts(int(list.size()), editable)));
            QHash<QString, QStandardItem*> owners;
            QStandardItem* raw = nullptr;
            int rawCount = 0;
            for (const Setting* st : list) rawCount += st->raw;
            for (const Setting* st : list) {
                QStandardItem* parent = group;
                if (st->raw) {
                    if (!raw) {
                        raw = folder(QObject::tr("Retained Fluent settings   %1").arg(QLocale().toString(rawCount)));
                        raw->setData(true, kRetainedRole);
                    }
                    parent = raw;
                } else if (!st->owner.isEmpty()) {
                    QStandardItem*& o = owners[st->owner];
                    if (!o) group->appendRow({o = folder(st->owner), blank(), blank()});
                    parent = o;
                }
                parent->appendRow(row(*st));
            }
            if (raw) group->appendRow({raw, blank(), blank()});
            appendRow({group, blank(), blank()});
        }
        QObject::connect(m, &SettingsModel::valueChanged, this, [this](const QString& path) { refresh(path); });
    }

    const Setting* setting(const QModelIndex& i) const { return m_model->find(i.siblingAtColumn(0).data(kPathRole).toString()); }

private:
    static QString counts(int n, int editable)
    {
        return editable ? QObject::tr("%1 · %2 editable").arg(QLocale().toString(n)).arg(editable) : QLocale().toString(n);
    }
    static QStandardItem* blank()
    {
        auto* it = new QStandardItem;
        it->setEditable(false);
        return it;
    }
    static QStandardItem* folder(const QString& text)
    {
        auto* it = new QStandardItem(text);
        it->setEditable(false);
        QFont f = it->font();
        f.setWeight(QFont::DemiBold);
        it->setFont(f);
        return it;
    }
    QList<QStandardItem*> row(const Setting& st)
    {
        auto* name = new QStandardItem(settingLabel(st));
        name->setData(st.path, kPathRole);
        name->setData(QStringLiteral("%1 %2 %3 %4").arg(st.path, st.label, settingLabel(st), st.owner), kSearchRole);
        name->setEditable(false);
        auto* value = blank();
        auto* unit = new QStandardItem(displayUnit(st).unit);
        unit->setEditable(false);
        m_rows.insert(st.path, {name, value});
        QList<QStandardItem*> r{name, value, unit};
        paint(st, r);
        return r;
    }
    void paint(const Setting& st, const QList<QStandardItem*>& r)
    {
        const auto& t = tokens();
        const bool mod = m_model->isModified(st.path);
        QString v = displayValue(st, m_model->value(st.path));
        const DisplayUnit du = displayUnit(st);
        if (!du.unit.isEmpty() && v.endsWith(" " + du.unit)) v.chop(du.unit.size() + 1);
        r[1]->setText(v.left(160));
        // Fields the selected method does not use are listed, dimmed, as Fluent keeps them too.
        const bool active = m_model->isActive(st.path);
        const QColor fg = !active ? t.faint : mod ? (t.dark ? t.amber : t.amberText) : st.editable ? t.text : t.muted;
        for (auto* it : r) it->setForeground(fg);
        QString tip = QStringLiteral("%1\n%2").arg(st.path, mod ? QObject::tr("Modified") : st.editable ? QObject::tr("Editable") : st.reason);
        if (!active) tip += "\n" + unusedReason(*m_model, st);
        for (auto* it : r) it->setToolTip(tip);
    }
    void refresh(const QString& path)
    {
        auto it = m_rows.find(path);
        const Setting* st = m_model->find(path);
        if (it == m_rows.end() || !st) return;
        QStandardItem* name = it->first;
        QStandardItem* parent = name->parent() ? name->parent() : invisibleRootItem();
        paint(*st, {name, it->second, parent->child(name->row(), 2)});
    }

    SettingsModel* m_model;
    QHash<QString, std::pair<QStandardItem*, QStandardItem*>> m_rows;
};

class SettingsFilter : public QSortFilterProxyModel {
public:
    using QSortFilterProxyModel::QSortFilterProxyModel;
    QString text;
    bool editableOnly = false, modifiedOnly = false;
    SettingsModel* model = nullptr;
    void apply() { beginFilterChange(); endFilterChange(); }
    bool filtering() const { return !text.isEmpty() || editableOnly || modifiedOnly; }

protected:
    // Recursive filtering keeps the folders of every matching setting.
    bool filterAcceptsRow(int row, const QModelIndex& parent) const override
    {
        const QModelIndex idx = sourceModel()->index(row, 0, parent);
        const QString path = idx.data(kPathRole).toString();
        if (path.isEmpty()) return !filtering();
        const Setting* s = model->find(path);
        if (!s) return false;
        if (editableOnly && !s->editable) return false;
        if (modifiedOnly && !model->isModified(s->path)) return false;
        const QString hay = idx.data(kSearchRole).toString();
        for (const QString& w : text.split(' ', Qt::SkipEmptyParts))
            if (!hay.contains(w, Qt::CaseInsensitive)) return false;
        return true;
    }
};

// Why a listed setting has no effect, in the words of its controlling choice.
QString unusedReason(const SettingsModel& m, const Setting& s)
{
    if (const Setting* c = s.dependsOn.isEmpty() ? nullptr : m.find(s.dependsOn))
        return QObject::tr("Not used: %1 is %2.").arg(settingLabel(*c), displayValue(*c, m.value(c->path)));
    return QObject::tr("Not used by this zone or model.");
}

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
    cv->addWidget(label(QObject::tr("Every setting saved in the case: %1, of which %2 are editable in Kutta. Kutta reads the rest to decide whether it can reproduce the case; select one to see why it is read-only.")
                            .arg(QLocale().toString(total), QLocale().toString(editable)), "muted"));
    auto* search = new QLineEdit(filter);
    search->setPlaceholderText(QObject::tr("Search names and paths, e.g. inlet turb"));
    search->setClearButtonEnabled(true);
    search->setAccessibleName(QObject::tr("Search settings"));
    cv->addWidget(search);
    auto* opts = new QHBoxLayout;
    auto* ed = new QCheckBox(QObject::tr("Editable only"));
    auto* mo = new QCheckBox(QObject::tr("Changed only"));
    auto* countLabel = label({}, "faint");
    countLabel->setWordWrap(false);
    opts->addWidget(ed);
    opts->addWidget(mo);
    opts->addStretch(1);
    opts->addWidget(countLabel);
    cv->addLayout(opts);
    auto* view = new QTreeView;
    auto* model = new SettingsTree(sm, view);
    auto* proxy = new SettingsFilter(view);
    proxy->model = sm;
    proxy->setRecursiveFilteringEnabled(true);
    proxy->setSourceModel(model);
    proxy->text = filter;
    view->setModel(proxy);
    view->setUniformRowHeights(true);  // keeps scrolling through 10,000 rows cheap
    view->setIndentation(14);
    view->setSelectionBehavior(QAbstractItemView::SelectRows);
    view->setWordWrap(false);
    view->setAccessibleName(QObject::tr("Case settings"));
    view->header()->setStretchLastSection(false);
    view->header()->setSectionResizeMode(0, QHeaderView::Stretch);
    view->header()->setSectionResizeMode(1, QHeaderView::Interactive);
    view->header()->setSectionResizeMode(2, QHeaderView::Fixed);
    view->setColumnWidth(1, 96);
    view->setColumnWidth(2, 50);
    view->setMinimumHeight(380);
    view->setTextElideMode(Qt::ElideRight);
    cv->addWidget(view, 1);

    // Details of the selected setting: its identity, value and why it is (not) editable.
    auto* detail = new QFrame;
    detail->setProperty("callout", "neutral");
    auto* dv = new QVBoxLayout(detail);
    dv->setContentsMargins(10, 8, 10, 10);
    dv->setSpacing(3);
    auto* dTitle = label({});
    QFont tf = dTitle->font();
    tf.setWeight(QFont::DemiBold);
    dTitle->setFont(tf);
    auto* dPath = label({}, "faint");
    dPath->setTextInteractionFlags(Qt::TextSelectableByMouse);
    auto* dValue = label({});
    auto* dWhy = label({}, "muted");
    auto* dGo = button(QObject::tr("Edit on its page"), "chevron-right");
    dv->addWidget(dTitle);
    dv->addWidget(dPath);
    dv->addWidget(dValue);
    dv->addWidget(dWhy);
    dv->addWidget(dGo, 0, Qt::AlignLeft);
    detail->hide();
    cv->addWidget(detail);
    auto showDetail = [=](const QModelIndex& idx) {
        const Setting* s = model->setting(proxy->mapToSource(idx));
        detail->setVisible(s != nullptr);
        if (!s) return;
        dTitle->setText(settingLabel(*s));
        dPath->setText(s->path);
        dValue->setText(QObject::tr("Value: %1").arg(displayValue(*s, sm->value(s->path))));
        QString why = sm->isModified(s->path) ? QObject::tr("Changed from the case value %1.").arg(displayValue(*s, s->value))
                      : s->editable           ? QObject::tr("Editable in Kutta.")
                                              : QObject::tr("Read-only: %1").arg(s->reason);
        if (!sm->isActive(s->path)) why = unusedReason(*sm, *s) + " " + why;
        dWhy->setText(why);
        dGo->setVisible(s->editable);
        dGo->setProperty("target", pageFor(*s));
    };
    QObject::connect(dGo, &QPushButton::clicked, page, [c, dGo] { c.navigate(dGo->property("target").toString()); });
    QObject::connect(view->selectionModel(), &QItemSelectionModel::currentRowChanged, page, [showDetail](const QModelIndex& i) { showDetail(i); });

    auto refresh = [proxy, view, search, ed, mo, countLabel] {
        proxy->text = search->text().trimmed();
        proxy->editableOnly = ed->isChecked();
        proxy->modifiedOnly = mo->isChecked();
        proxy->apply();
        int shown = 0;
        std::function<void(const QModelIndex&)> countLeaves = [&](const QModelIndex& parent) {
            for (int r = 0; r < proxy->rowCount(parent); ++r) {
                const QModelIndex i = proxy->index(r, 0, parent);
                if (i.data(kPathRole).toString().isEmpty()) countLeaves(i);
                else ++shown;
            }
        };
        countLeaves({});
        countLabel->setText(proxy->filtering() ? QObject::tr("%1 match").arg(QLocale().toString(shown)) : QString());
        // Folder rows span the columns so their names and counts stay readable.
        std::function<void(const QModelIndex&)> span = [&](const QModelIndex& parent) {
            for (int r = 0; r < proxy->rowCount(parent); ++r) {
                const QModelIndex i = proxy->index(r, 0, parent);
                if (!i.data(kPathRole).toString().isEmpty()) continue;
                view->setFirstColumnSpanned(r, parent, true);
                span(i);
            }
        };
        span({});
        // Open the folders of a focused result, except retained Fluent settings
        // when Kutta's own controls also match.
        if (!proxy->filtering() || shown > 400) return view->collapseAll();
        view->expandAll();
        if (shown == 0) return;
        std::function<int(const QModelIndex&)> controls = [&](const QModelIndex& parent) {
            int n = 0;
            for (int r = 0; r < proxy->rowCount(parent); ++r) {
                const QModelIndex i = proxy->index(r, 0, parent);
                if (i.data(kRetainedRole).toBool()) continue;
                n += i.data(kPathRole).toString().isEmpty() ? controls(i) : 1;
            }
            return n;
        };
        if (controls({}) == 0) return;
        std::function<void(const QModelIndex&)> fold = [&](const QModelIndex& parent) {
            for (int r = 0; r < proxy->rowCount(parent); ++r) {
                const QModelIndex i = proxy->index(r, 0, parent);
                if (i.data(kRetainedRole).toBool()) view->collapse(i);
                else if (i.data(kPathRole).toString().isEmpty()) fold(i);
            }
        };
        fold({});
    };
    auto* debounce = new QTimer(page);
    debounce->setSingleShot(true);
    debounce->setInterval(150);
    QObject::connect(debounce, &QTimer::timeout, page, refresh);
    QObject::connect(search, &QLineEdit::textChanged, debounce, qOverload<>(&QTimer::start));
    QObject::connect(ed, &QCheckBox::toggled, page, refresh);
    QObject::connect(mo, &QCheckBox::toggled, page, refresh);
    QObject::connect(view, &QTreeView::activated, page, [c, model, proxy](const QModelIndex& idx) {
        if (const Setting* s = model->setting(proxy->mapToSource(idx)); s && s->editable) c.navigate(pageFor(*s));
    });
    refresh();
    return page;
}

} // namespace kapp::pages
