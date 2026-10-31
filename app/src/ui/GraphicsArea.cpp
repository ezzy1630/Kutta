#include "ui/GraphicsArea.h"
#include "core/ChartScale.h"

#include "ui/Theme.h"

#include <QComboBox>
#include <QFileDialog>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLocale>
#include <QMenu>
#include <QStandardItemModel>
#include <QTimer>
#include <QToolButton>
#include <QVBoxLayout>

namespace kapp {

namespace {

QToolButton* toolButton(const QString& icon, const QString& name, QWidget* parent)
{
    auto* b = new QToolButton(parent);
    b->setIcon(Theme::instance().icon(icon));
    b->setIconSize(QSize(16, 16));
    b->setToolTip(name);
    b->setAccessibleName(name);
    b->setAutoRaise(true);
    return b;
}

// First view of a dataset: along a clearly thin axis (a planar or symmetric
// internal flow such as the elbow shows its shape), else isometric (a car).
Camera::Preset initialView(const Bounds& b)
{
    if (!b.valid) return Camera::Preset::Isometric;
    const QVector3D e = b.max - b.min;
    int thin = 0;
    for (int a = 1; a < 3; ++a)
        if (e[a] < e[thin]) thin = a;
    float other = 0;
    for (int a = 0; a < 3; ++a)
        if (a != thin) other = std::max(other, e[a]);
    if (e[thin] > 0.15f * other) return Camera::Preset::Isometric;
    return thin == 2 ? Camera::Preset::Top : thin == 1 ? Camera::Preset::Left : Camera::Preset::Front;
}

QStringList badgesFor(const SceneData& s, const QString& key, const Session* session)
{
    QStringList b;
    const QString p = s.provenance.toUpper();
    if (p.contains("APPROXIMATE")) {
        // Name the approximations when the view is a run's (or the live run's).
        const auto& runs = session->history().runs();
        const RunRecord* r = key.startsWith("run:") ? session->history().find(key.mid(4))
                             : key == "live" && !runs.empty() ? &runs.back() : nullptr;
        b << (r && r->approximate() ? r->approximateBadge() : QStringLiteral("APPROXIMATE"));
    }
    if (p.contains("MODIFIED")) b << "MODIFIED";
    if (s.synthetic) b << "SYNTHETIC";
    else if (p.contains("FIXTURE")) b << "FIXTURE REPLAY";
    if (key == "live" && session->running()) b << QStringLiteral("LIVE · ITER %1").arg(session->live().lastSnapshotIter);
    else if (key == "live") b << QStringLiteral("SNAPSHOT · ITER %1").arg(session->live().lastSnapshotIter);
    return b;
}

} // namespace

// ---- ViewportPanel -----------------------------------------------------------------

ViewportPanel::ViewportPanel(Session* session, int index, QWidget* parent)
    : QWidget(parent), m_session(session), m_index(index)
{
    auto* v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    auto* header = new QWidget(this);
    header->setObjectName("ViewportHeader");
    header->setAttribute(Qt::WA_StyledBackground);
    auto* h = new QHBoxLayout(header);
    h->setContentsMargins(6, 3, 4, 3);
    h->setSpacing(4);
    m_dataset = new QComboBox(header);
    m_dataset->setAccessibleName(tr("Viewport %1 dataset").arg(index + 1));
    m_dataset->setToolTip(tr("Dataset shown in this viewport"));
    m_dataset->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    m_field = new QComboBox(header);
    m_field->setAccessibleName(tr("Viewport %1 field").arg(index + 1));
    m_field->setToolTip(tr("Field coloured on the surfaces"));
    m_field->setSizeAdjustPolicy(QComboBox::AdjustToContents);
    m_field->setMaxVisibleItems(24);
    h->addWidget(m_dataset);
    h->addWidget(m_field);
    h->addStretch(1);
    m_camera = toolButton("camera", tr("Camera views"), header);
    m_camera->setPopupMode(QToolButton::InstantPopup);
    auto* camMenu = new QMenu(m_camera);
    for (auto [preset, key] : std::initializer_list<std::pair<Camera::Preset, const char*>>{
             {Camera::Preset::Front, "1"}, {Camera::Preset::Back, "2"}, {Camera::Preset::Left, "3"}, {Camera::Preset::Right, "4"},
             {Camera::Preset::Top, "5"}, {Camera::Preset::Bottom, "6"}, {Camera::Preset::Isometric, "7"}}) {
        auto* a = camMenu->addAction(Camera::presetName(preset), this, [this, p = preset] { m_view->applyPreset(p); });
        a->setShortcut(QKeySequence(QString::fromLatin1(key)));
        a->setShortcutContext(Qt::WidgetShortcut);
    }
    camMenu->addSeparator();
    camMenu->addAction(tr("Fit to visible"), this, [this] { m_view->fitView(); })->setShortcut(QKeySequence("F"));
    auto* proj = camMenu->addAction(tr("Orthographic"));
    proj->setCheckable(true);
    connect(camMenu, &QMenu::aboutToShow, this, [this, proj] { proj->setChecked(m_view->display().camera.fovY <= 0); });
    connect(proj, &QAction::toggled, this, [this](bool ortho) {
        Camera c = m_view->display().camera;
        c.fovY = ortho ? 0 : 30;
        m_view->setCamera(c);
    });
    m_camera->setMenu(camMenu);
    auto* shot = toolButton("image-down", tr("Save image of this viewport"), header);
    connect(shot, &QToolButton::clicked, this, [this] { emit exportRequested(this); });
    h->addWidget(m_camera);
    h->addWidget(shot);
    v->addWidget(header);

    m_view = new SceneView(this);
    m_view->setAccessibleName(tr("Viewport %1").arg(index + 1));
    m_view->setFocusPolicy(Qt::StrongFocus);
    v->addWidget(m_view, 1);
    const auto& t = tokens();
    m_view->setTheme(t.dark, t.viewport, t.accent);
    connect(&Theme::instance(), &Theme::changed, this, [this] {
        const auto& tk = tokens();
        m_view->setTheme(tk.dark, tk.viewport, tk.accent);
        updateOverlay();
    });

    connect(m_view, &SceneView::activated, this, [this] { emit activated(this); });
    connect(m_view, &SceneView::picked, this, [this](const PickResult& r, QPointF) {
        m_probe = r.hit ? std::optional<PickResult>(r) : std::nullopt;
        updateOverlay();
        emit probed(this, r);
    });
    connect(m_view, &SceneView::cameraChanged, this, [this] { emit displayChanged(this); });
    connect(m_dataset, &QComboBox::activated, this, [this](int i) {
        setDataset(m_dataset->itemData(i).toString());
        emit activated(this);
    });
    connect(m_field, &QComboBox::activated, this, [this](int i) {
        ViewDisplay d = m_view->display();
        d.field = FieldRef::fromKey(m_field->itemData(i).toString());
        d.contours = true;
        setDisplay(d);
        emit displayChanged(this);
    });
    connect(session, &Session::sceneChanged, this, [this](const QString& key) {
        if (key == m_key) sceneArrived(m_session->sceneIfLoaded(key));
    });
    connect(session, &Session::liveSceneChanged, this, [this] {
        if (m_key == "live") sceneArrived(m_session->sceneIfLoaded("live"));
    });
    connect(session, &Session::runStateChanged, this, [this] { if (m_key == "live") updateOverlay(); });
    updateOverlay();
}

QString ViewportPanel::datasetLabel() const
{
    const int i = m_dataset->findData(m_key);
    return i >= 0 ? m_dataset->itemText(i) : m_key;
}

void ViewportPanel::refreshDatasets(const QList<QPair<QString, QString>>& datasets)
{
    const QSignalBlocker b(m_dataset);
    m_dataset->clear();
    for (const auto& [key, label] : datasets) m_dataset->addItem(label, key);
    if (!m_key.isEmpty() && m_dataset->findData(m_key) < 0) m_dataset->addItem(m_key, m_key);
    m_dataset->setCurrentIndex(m_dataset->findData(m_key));
    m_dataset->setEnabled(m_dataset->count() > 0);
}

void ViewportPanel::setDataset(const QString& key)
{
    m_pendingKey = key;
    m_probe.reset();
    if (key.isEmpty()) {
        m_key.clear();
        m_view->setScene(nullptr);
        updateOverlay();
        return;
    }
    if (m_key != key) {
        m_key = key;
        m_dataset->setCurrentIndex(m_dataset->findData(key));
    }
    if (auto s = m_session->sceneIfLoaded(key)) {
        sceneArrived(s);
        return;
    }
    OverlayState o = m_view->overlayState();
    o.emptyText = tr("Loading %1…").arg(datasetLabel());
    m_view->setScene(nullptr);
    m_view->setOverlay(o);
    QPointer<ViewportPanel> self(this);
    m_session->loadScene(key, [self, key](const QString&, Result<std::shared_ptr<SceneData>> r) {
        if (!self || self->m_pendingKey != key) return;
        if (!r) {
            OverlayState o = self->m_view->overlayState();
            // An older CLI without export-geometry is not a failure of the case.
            o.emptyText = key.startsWith("mesh:") && r.error().contains("export-geometry")
                              ? tr("No mesh preview: this kutta CLI cannot draw a case before a run.\nRun the calculation to see the case and its results.")
                              : tr("Could not load %1\n%2").arg(self->datasetLabel(), r.error());
            o.badges.clear();
            o.legend = false;
            self->m_view->setOverlay(o);
            return;
        }
        self->sceneArrived(r.value());
    });
    emit datasetChanged(this);
}

void ViewportPanel::sceneArrived(std::shared_ptr<const SceneData> scene)
{
    const bool first = !m_view->scene() || m_view->scene()->label.section(" · ", 0, 0) != (scene ? scene->label.section(" · ", 0, 0) : QString());
    ViewDisplay d = m_view->display();
    // The mesh and a solution are displayed differently (zones coloured and
    // edged versus a field), so moving between them starts from each one's defaults.
    const bool wasMesh = m_showedMesh;
    const bool isMesh = scene && scene->geometryOnly;
    if (scene) m_showedMesh = isMesh;
    if (scene && (wasMesh != isMesh || (isMesh && first))) {
        d.visibleInitialized = false;
        d.zoneColors = isMesh;
        d.contours = !isMesh;
        std::size_t tris = 0;
        for (const auto& it : scene->items) tris += it.mesh ? it.mesh->triangleCount() : 0;
        // Edges on a multi-million-face body are a solid black blur; Fluent
        // users switch them on when they want them (Results ▸ Graphics ▸ Mesh).
        d.edges = isMesh && tris < 2'000'000;
    }
    if (scene) {
        QStringList ids;
        for (const auto& it : scene->items) ids << it.id;
        if (!d.visibleInitialized || !std::all_of(d.visible.begin(), d.visible.end(), [&](const QString& v) { return ids.contains(v); })) {
            d.visible = isMesh ? defaultMeshItems(*scene) : defaultVisibleItems(*scene);
            d.visibleInitialized = true;
        } else {
            // New post-processing parts become visible as they arrive.
            for (const auto& it : scene->items)
                if (it.kind != ItemKind::Zone && it.kind != ItemKind::Slice && !m_knownItems.contains(it.id)) d.visible << it.id;
        }
        m_knownItems = ids;
        const auto fields = scene->availableFields();
        // A viewport's first dataset gets the field that suits the case; later
        // datasets keep what the user chose when they have it.
        if (!isMesh && (!m_fieldChosen || std::find(fields.begin(), fields.end(), d.field) == fields.end())) d.field = defaultField(*scene);
        if (scene->runDelta && d.field.role != FieldRole::Run) {
            FieldRef f = d.field;
            f.role = FieldRole::Run;
            if (std::find(fields.begin(), fields.end(), f) != fields.end()) d.field = f;
        }
    }
    // The first solution after the mesh preview shows its contour plane face on,
    // not from the mesh view's angle (a nozzle's axial slice seen obliquely).
    if (scene && wasMesh && !isMesh && !scene->planar && commonPlaneAxis(*scene, d.visible)) d.cameraInitialized = false;
    if (scene && scene->planar) {
        // A 2D case is always seen face on; a camera from a 3D view is replaced.
        if (d.camera.up != QVector3D(0, 1, 0) || d.camera.pitch != 0 || d.camera.yaw != 90) d.cameraInitialized = false;
        if (!d.cameraInitialized) d.camera.planView();
    } else if (scene && !d.cameraInitialized) {
        if (d.camera.up != QVector3D(0, 0, 1)) d.camera = Camera{};  // leaving a 2D case's view
        // Contours on one plane (a slice or symmetry plane) are seen face on;
        // otherwise the dataset's shape decides.
        const auto axis = isMesh ? std::nullopt : commonPlaneAxis(*scene, d.visible);
        d.camera.preset(axis ? (*axis == 2 ? Camera::Preset::Top : *axis == 1 ? Camera::Preset::Left : Camera::Preset::Front)
                             : initialView(scene->bounds));
    }
    m_view->setScene(scene);
    m_view->setDisplay(d);
    if (scene && (!d.cameraInitialized || wasMesh != isMesh)) m_view->fitView();
    else if (first && scene && !d.cameraInitialized) m_view->fitView();
    rebuildFieldCombo();
    updateOverlay();
    emit displayChanged(this);
}

void ViewportPanel::setDisplay(const ViewDisplay& d)
{
    m_fieldChosen = true;
    m_view->setDisplay(d);
    rebuildFieldCombo();
    updateOverlay();
}

void ViewportPanel::rebuildFieldCombo()
{
    const QSignalBlocker b(m_field);
    m_field->clear();
    const auto s = m_view->scene();
    if (!s) {
        m_field->setEnabled(false);
        return;
    }
    m_field->setEnabled(true);
    std::optional<FieldRole> lastRole;
    for (const FieldRef& f : s->availableFields()) {
        if (lastRole && *lastRole != f.role) m_field->insertSeparator(m_field->count());
        lastRole = f.role;
        m_field->addItem(displayLabel(f), f.key());
    }
    m_field->setCurrentIndex(m_field->findData(m_view->display().field.key()));
}

QString ViewportPanel::probeText(const PickResult& r) const
{
    const auto s = m_view->scene();
    const ViewDisplay& d = m_view->display();
    QStringList lines;
    if (s && d.contours) {
        const QString unit = prettyUnit(s->unitOf(d.field.quantity));
        lines << QStringLiteral("%1  %2%3").arg(displayLabel(d.field),
                                                 r.value ? formatNumber(*r.value, 5) : tr("no data here"),
                                                 r.value && !unit.isEmpty() ? " " + unit : QString());
    }
    lines << r.itemLabel;
    lines << QStringLiteral("x %1  y %2  z %3 m").arg(formatNumber(r.world.x(), 4), formatNumber(r.world.y(), 4), formatNumber(r.world.z(), 4));
    return lines.join('\n');
}

void ViewportPanel::showProbe(const std::optional<PickResult>& r, bool keepValue)
{
    m_probe = r;
    if (m_probe && m_probe->hit && !keepValue) {
        // Re-sample this viewport's own field at the shared world point.
        const auto s = m_view->scene();
        m_probe->value.reset();
        if (s)
            if (const SceneItem* it = s->item(m_probe->itemId); it && it->mesh)
                m_probe->value = it->mesh->sampleNear(m_probe->world, m_view->display().field);
    }
    updateOverlay();
}

void ViewportPanel::setActive(bool active)
{
    m_active = active;
    updateOverlay();
}

void ViewportPanel::updateOverlay()
{
    OverlayState o;
    o.dark = tokens().dark;
    o.active = m_active;
    const auto s = m_view->scene();
    const ViewDisplay& d = m_view->display();
    o.view = d.camera.view();
    m_dataset->setVisible(m_dataset->count() > 0);
    m_field->setVisible(s != nullptr);
    if (!s) {
        o.triad = false;
        o.emptyText = m_key.isEmpty() ? tr("No results yet.\nRun the calculation, or open a run from Run History.") : m_view->overlayState().emptyText;
        m_view->setOverlay(o);
        return;
    }
    if (s->geometryOnly) {
        // Fluent's mesh display: zones by colour, with a key, before any run.
        m_field->setVisible(false);
        o.title = tr("Mesh");
        std::size_t faces = 0;
        for (const auto& z : m_session->hasCase() ? m_session->caseSummary().zones : std::vector<Zone>{})
            if (z.type != "interior") faces += std::size_t(z.faces);
        o.subtitle = !m_session->hasCase() ? QString()
                     : s->planar ? tr("%1 · 2D · %2 cells").arg(m_session->caseName(), QLocale().toString(m_session->caseSummary().cells))
                                 : tr("%1 · %2 boundary faces").arg(m_session->caseName(), QLocale().toString(qulonglong(faces)));
        for (const auto& it : s->items)
            if (it.kind == ItemKind::Zone && d.visible.contains(it.id)) o.zoneKey.push_back({it.label, d.zoneColors ? it.color : QColor()});
        if (m_probe && m_probe->hit) {
            o.probe = m_probe;
            o.probePixel = m_view->projectToLogical(m_probe->world);
            o.probeText = probeText(*m_probe);
        }
        if (o.zoneKey.empty()) o.emptyText = tr("No zones are shown. Choose them in Results ▸ Graphics ▸ Mesh.");
        m_view->setOverlay(o);
        return;
    }
    const QString unit = prettyUnit(s->unitOf(d.field.quantity));
    o.unit = unit;
    const bool fieldShown = d.contours && !d.field.quantity.isEmpty();
    // Fluent names the display by what is drawn: pathlines alone are
    // "Pathlines Coloured by", surfaces are "Contours of".
    bool surfaces = false, pathlines = false;
    for (const auto& it : s->items)
        if (d.visible.contains(it.id)) {
            if (it.kind == ItemKind::Streamlines) pathlines = true;
            else if (it.kind != ItemKind::Probe) surfaces = true;
        }
    const QString drawn = pathlines && !surfaces ? tr("Pathlines Coloured by %1") : tr("Contours of %1");
    o.title = fieldShown ? drawn.arg(displayLabel(d.field)) : tr("Surfaces");
    if (fieldShown && !unit.isEmpty()) o.title += QStringLiteral(" [%1]").arg(unit);
    if (fieldShown && pathlines && surfaces) o.title += tr(" · pathlines");
    if (d.vectors) o.title += fieldShown ? tr(" · velocity vectors") : tr(" with velocity vectors");
    o.subtitle = datasetLabel() + (m_session->hasCase() ? " · " + m_session->caseName() : QString());
    o.badges = badgesFor(*s, m_key, m_session);
    o.provenance = s->provenance;
    if (s->synthetic && !o.provenance.contains("SYNTHETIC", Qt::CaseInsensitive))
        o.provenance = (o.provenance.isEmpty() ? QString() : o.provenance + " | ") + "SYNTHETIC fixture data, not a Kutta solution";
    o.colormap = d.effectiveColormap();
    o.range = m_view->resolvedRange();
    const auto dataRange = s->range(d.field, d.visible);
    if (!o.range.valid()) {
        if (dataRange) o.range = autoRange(dataRange->first, dataRange->second, d.symmetricRange());
        if (d.rangeMode == RangeMode::Fixed) o.range = {d.fixedMin, d.fixedMax};
    }
    // Walls carry no velocity, inlets no wall shear: say so rather than show a legend for grey.
    o.legend = fieldShown && dataRange;
    if (fieldShown && !dataRange) o.subtitle += tr(" · %1 is not defined on the surfaces shown").arg(displayLabel(d.field));
    o.bands = d.bands;
    if (m_probe && m_probe->hit) {
        o.probe = m_probe;
        o.probePixel = m_view->projectToLogical(m_probe->world);
        o.probeText = probeText(*m_probe);
    }
    bool anyVisible = false;
    for (const auto& it : s->items) anyVisible |= d.visible.contains(it.id);
    if (!anyVisible) o.emptyText = tr("Nothing is shown. Choose surfaces in Results ▸ Graphics.");
    m_view->setOverlay(o);
}

// ---- GraphicsArea ------------------------------------------------------------------

GraphicsArea::GraphicsArea(Session* session, QWidget* parent) : QWidget(parent), m_session(session)
{
    m_grid = new QGridLayout(this);
    m_grid->setContentsMargins(0, 0, 0, 0);
    m_grid->setSpacing(1);
    setAutoFillBackground(true);
    for (int i = 0; i < 4; ++i) {
        auto* p = new ViewportPanel(session, i, this);
        m_panels[std::size_t(i)] = p;
        connect(p, &ViewportPanel::activated, this, &GraphicsArea::setActive);
        connect(p, &ViewportPanel::displayChanged, this, [this](ViewportPanel* src) {
            if (m_linkCameras && !m_syncing) {
                m_syncing = true;
                for (auto* o : m_panels)
                    if (o != src && o->isVisible() && o->display().camera != src->display().camera) o->view()->setCamera(src->display().camera);
                m_syncing = false;
            }
            emit displayChanged(src);
        });
        connect(p, &ViewportPanel::probed, this, [this](ViewportPanel* src, const PickResult& r) {
            for (auto* o : m_panels) {
                if (o == src) continue;
                // Same dataset and field: show the picked value itself, not a re-sample.
                const bool same = o->scene() == src->scene() && o->display().field == src->display().field;
                o->showProbe(r.hit ? std::optional<PickResult>(r) : std::nullopt, same);
            }
        });
        connect(p, &ViewportPanel::exportRequested, this, [this](ViewportPanel* src) {
            const QString f = QFileDialog::getSaveFileName(this, tr("Save viewport image"),
                                                           QStringLiteral("%1.png").arg(src->datasetLabel()), tr("PNG image (*.png)"));
            if (!f.isEmpty()) src->exportImage().save(f, "PNG");
        });
    }
    connect(session, &Session::historyChanged, this, &GraphicsArea::refreshDatasets);
    connect(session, &Session::liveSceneChanged, this, &GraphicsArea::refreshDatasets);
    connect(session, &Session::caseOpened, this, [this] {
        m_extra.clear();
        for (auto* p : m_panels) p->setDataset({});
        refreshDatasets();
    });
    auto bg = [this] {
        QPalette pal = palette();
        pal.setColor(QPalette::Window, tokens().border);
        setPalette(pal);
    };
    bg();
    connect(&Theme::instance(), &Theme::changed, this, bg);
    relayout();
    setActive(m_panels[0]);
}

int GraphicsArea::visibleCount() const
{
    switch (m_layout) {
    case Layout::Single: return 1;
    case Layout::SideBySide:
    case Layout::Stacked: return 2;
    case Layout::Quad: return 4;
    }
    return 1;
}

void GraphicsArea::setViewLayout(Layout l)
{
    m_layout = l;
    relayout();
    // New viewports start on the active viewport's dataset and display.
    for (int i = 1; i < visibleCount(); ++i) {
        ViewportPanel* p = m_panels[std::size_t(i)];
        if (p->datasetKey().isEmpty() && m_active && !m_active->datasetKey().isEmpty()) {
            p->setDataset(m_active->datasetKey());
        }
    }
    if (m_active && m_active->index() >= visibleCount()) setActive(m_panels[0]);
    // Views changed shape: frame their contents again once the layout settles.
    QTimer::singleShot(0, this, [this] {
        for (int i = 0; i < visibleCount(); ++i)
            if (m_panels[std::size_t(i)]->scene()) m_panels[std::size_t(i)]->view()->fitView();
    });
}

void GraphicsArea::relayout()
{
    for (auto* p : m_panels) {
        m_grid->removeWidget(p);
        p->hide();
    }
    const int n = visibleCount();
    for (int i = 0; i < n; ++i) {
        const int row = m_layout == Layout::Stacked ? i : m_layout == Layout::Quad ? i / 2 : 0;
        const int col = m_layout == Layout::SideBySide ? i : m_layout == Layout::Quad ? i % 2 : 0;
        m_grid->addWidget(m_panels[std::size_t(i)], row, col);
        m_panels[std::size_t(i)]->show();
    }
    for (int i = 0; i < 2; ++i) {
        m_grid->setRowStretch(i, (m_layout == Layout::Stacked || m_layout == Layout::Quad || i == 0) ? 1 : 0);
        m_grid->setColumnStretch(i, (m_layout == Layout::SideBySide || m_layout == Layout::Quad || i == 0) ? 1 : 0);
    }
    for (auto* p : m_panels) p->setActive(p == m_active && n > 1);
}

void GraphicsArea::setActive(ViewportPanel* p)
{
    if (!p) return;
    const bool changed = p != m_active;
    m_active = p;
    for (auto* o : m_panels) o->setActive(o == p && visibleCount() > 1);
    if (changed) emit activeChanged(p);
}

void GraphicsArea::addDataset(const QString& key, const QString& label)
{
    for (const auto& e : m_extra)
        if (e.first == key) return;
    m_extra.append({key, label});
    refreshDatasets();
}

void GraphicsArea::refreshDatasets()
{
    QList<QPair<QString, QString>> list;
    if (m_session->sceneIfLoaded("live")) list.append({"live", m_session->running() ? tr("Live solution") : tr("Last live snapshot")});
    const auto& runs = m_session->history().runs();
    for (auto it = runs.rbegin(); it != runs.rend(); ++it)
        if (it->status == RunStatus::Completed || it->status == RunStatus::Stopped) list.append({"run:" + it->id, it->title()});
    list.append(m_extra);
    if (m_session->hasCase()) list.append({m_session->meshSceneKey(), tr("Mesh")});
    for (auto* p : m_panels) p->refreshDatasets(list);
}

void GraphicsArea::show(const QString& key, int panelIndex)
{
    ViewportPanel* p = panelIndex >= 0 ? m_panels[std::size_t(panelIndex)] : m_active;
    if (!p) return;
    p->setDataset(key);
}

} // namespace kapp
