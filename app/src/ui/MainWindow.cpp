#include "ui/MainWindow.h"

#include "ui/ChartsPanel.h"
#include "ui/CommandPalette.h"
#include "ui/ConsoleView.h"
#include "ui/GraphicsArea.h"
#include "ui/Inspector.h"
#include "ui/Outline.h"
#include "ui/PageKit.h"
#include "ui/SettingWidgets.h"

#include <DockAreaWidget.h>
#include <DockManager.h>
#include <DockWidget.h>

#include <QApplication>
#include <QCloseEvent>
#include <QComboBox>
#include <QDragEnterEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QMenuBar>
#include <QMessageBox>
#include <QMimeData>
#include <QPushButton>
#include <QSettings>
#include <QStatusBar>
#include <QToolButton>

namespace kapp {

using namespace kit;

namespace {

bool isCaseFile(const QString& p)
{
    const QString n = QFileInfo(p).fileName().toLower();
    return n.endsWith(".cas.h5") || n.endsWith(".cas") || n.endsWith(".cas.gz");
}

QToolButton* iconButton(const QString& icon, const QString& name, QWidget* parent)
{
    auto* b = new QToolButton(parent);
    b->setIcon(Theme::instance().icon(icon));
    b->setIconSize(QSize(16, 16));
    b->setToolTip(name);
    b->setAccessibleName(name);
    b->setAutoRaise(true);
    return b;
}

} // namespace

MainWindow::MainWindow(Session* session, const QStringList& samples, QWidget* parent)
    : QMainWindow(parent), m_session(session), m_samples(samples)
{
    setWindowTitle("Kutta");
    setAcceptDrops(true);
    setUnifiedTitleAndToolBarOnMac(true);
    auto* central = new QWidget(this);
    auto* v = new QVBoxLayout(central);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);
    v->addWidget(buildToolbar());
    buildDocks();
    v->addWidget(m_docks, 1);
    setCentralWidget(central);
    buildMenus();

    m_cliLabel = new QLabel(this);
    const CliInfo& cli = session->cliInfo();
    m_cliLabel->setText(cli.fixture ? tr("Solver: fixture CLI (recorded and synthetic data)") : tr("Solver: %1").arg(QDir::toNativeSeparators(cli.program)));
    m_cliLabel->setToolTip(cli.describe());
    statusBar()->addPermanentWidget(m_cliLabel);
    statusBar()->setSizeGripEnabled(false);

    m_palette = new CommandPalette(this);
    m_palette->setSource([this] { return paletteEntries(); });

    connect(session, &Session::consoleLine, m_console, &ConsoleView::append);
    connect(session, &Session::caseOpened, this, [this] {
        refreshTitle();
        refreshToolbar();
        statusBar()->showMessage(tr("Opened %1").arg(m_session->caseName()), 5000);
        // Show the newest finished run's results, as Fluent shows the last data;
        // without one, the mesh, as Fluent does when it reads a case.
        const QString key = m_session->latestSceneKey();
        m_graphics->show(key.isEmpty() ? m_session->meshSceneKey() : key, 0);
        m_charts->showRuns({});
        navigate("setup/general");
    });
    // Shown in place on the welcome or General page (no modal dialog); the
    // case that was open stays open.
    connect(session, &Session::caseFailed, this, [this](const QString& m) {
        statusBar()->showMessage(tr("Could not open the case: %1").arg(m.section('\n', 0, 0)), 10000);
        navigate("setup/general");
    });
    connect(session, &Session::runStateChanged, this, &MainWindow::refreshToolbar);
    connect(session, &Session::iteration, this, &MainWindow::refreshToolbar);
    connect(session, &Session::devicesChanged, this, &MainWindow::refreshToolbar);
    connect(session, &Session::approvalsChanged, this, &MainWindow::refreshToolbar);
    connect(session->settings(), &SettingsModel::overridesChanged, this, &MainWindow::refreshToolbar);
    connect(session, &Session::runRefused, this, [this](const QString& why) {
        statusBar()->showMessage(why, 8000);
        if (why.contains("approval")) navigate("setup/general");
    });
    connect(session, &Session::runStarted, this, [this] {
        m_charts->showTab(ChartsPanel::Tab::Residuals);
        if (m_session->liveViewInterval() > 0) m_graphics->show("live", 0);
    });
    connect(session, &Session::liveSceneChanged, this, [this] {
        // First snapshot of a run: put it in the active viewport if it shows nothing.
        if (m_graphics->active() && m_graphics->active()->datasetKey().isEmpty()) m_graphics->show("live");
    });
    connect(session, &Session::runFinished, this, [this](const QString& id) {
        const RunRecord* r = m_session->history().find(id);
        if (!r) return;
        statusBar()->showMessage(tr("%1 %2").arg(r->title(), runStatusLabel(r->status).toLower()), 8000);
        if (r->status == RunStatus::Completed || r->status == RunStatus::Stopped) {
            for (int i = 0; i < m_graphics->visibleCount(); ++i) {
                auto* p = m_graphics->panel(i);
                if (p->datasetKey() == "live" || p->datasetKey().isEmpty() || i == 0) {
                    m_graphics->show("run:" + id, i);
                    break;
                }
            }
        }
    });
    connect(&Theme::instance(), &Theme::changed, this, [this] {
        m_themeButton->setIcon(Theme::instance().icon(tokens().dark ? "sun" : "moon"));
        refreshToolbar();
    });
    refreshTitle();
    refreshToolbar();
    QSettings s;
    if (!restoreGeometry(s.value("window/geometry").toByteArray())) resize(1560, 960);
}

MainWindow::~MainWindow() = default;

QWidget* MainWindow::buildToolbar()
{
    auto* bar = new QWidget(this);
    bar->setObjectName("Toolbar");
    bar->setAttribute(Qt::WA_StyledBackground);
    auto* h = new QHBoxLayout(bar);
    h->setContentsMargins(12, 6, 10, 6);
    h->setSpacing(8);
    auto* mark = new QLabel("Kutta", bar);
    QFont mf = mark->font();
    mf.setWeight(QFont::Bold);
    mf.setPointSizeF(mf.pointSizeF() * 1.08);
    mf.setLetterSpacing(QFont::AbsoluteSpacing, 0.3);
    mark->setFont(mf);
    h->addWidget(mark);
    m_caseLabel = new QLabel(bar);
    m_caseLabel->setProperty("role", "muted");
    h->addWidget(m_caseLabel);
    h->addSpacing(14);

    m_runButton = button(tr("Run"), "play", "primary", bar);
    m_runButton->setMinimumWidth(92);
    connect(m_runButton, &QPushButton::clicked, this, [this] {
        if (m_session->running()) m_session->stop();
        else m_session->run();
    });
    h->addWidget(m_runButton);
    auto* itLabel = new QLabel(tr("Iterations"), bar);
    itLabel->setProperty("role", "muted");
    h->addWidget(itLabel);
    m_iterations = new QLineEdit(bar);
    m_iterations->setValidator(new QIntValidator(1, 10'000'000, m_iterations));
    m_iterations->setFixedWidth(72);
    m_iterations->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    m_iterations->setAccessibleName(tr("Number of iterations"));
    itLabel->setBuddy(m_iterations);
    connect(m_iterations, &QLineEdit::editingFinished, this, [this] {
        if (!m_session->hasCase()) return;
        const QString err = m_session->settings()->setValue("run/iterations", m_iterations->text().toInt());
        if (!err.isEmpty()) statusBar()->showMessage(err, 6000);
        refreshToolbar();
    });
    h->addWidget(m_iterations);
    m_device = new QComboBox(bar);
    m_device->setAccessibleName(tr("Device"));
    m_device->setToolTip(tr("Device kutta solves on"));
    connect(m_device, &QComboBox::activated, this, [this](int i) { m_session->setBackend(m_device->itemData(i).toString()); });
    h->addWidget(m_device);
    m_stateLabel = new QLabel(bar);
    m_stateLabel->setProperty("role", "muted");
    m_stateLabel->setMinimumWidth(240);
    m_stateLabel->setTextFormat(Qt::RichText);
    connect(m_stateLabel, &QLabel::linkActivated, this, [this](const QString& href) { navigate(href); });
    QFont sf = m_stateLabel->font();
    sf.setFeature(QFont::Tag("tnum"), 1);
    m_stateLabel->setFont(sf);
    h->addWidget(m_stateLabel);
    h->addStretch(1);

    const char* icons[3] = {"square", "columns-2", "layout-grid"};
    const QString names[3] = {tr("One viewport"), tr("Two viewports"), tr("Four viewports")};
    for (int i = 0; i < 3; ++i) {
        m_layoutButtons[i] = iconButton(icons[i], names[i], bar);
        m_layoutButtons[i]->setCheckable(true);
        connect(m_layoutButtons[i], &QToolButton::clicked, this, [this, i] {
            m_graphics->setViewLayout(i == 0 ? GraphicsArea::Layout::Single : i == 1 ? GraphicsArea::Layout::SideBySide : GraphicsArea::Layout::Quad);
            refreshToolbar();
        });
        h->addWidget(m_layoutButtons[i]);
    }
    h->addSpacing(6);
    auto* search = new QPushButton(bar);
    search->setIcon(Theme::instance().icon("search"));
    search->setText(tr("Search   %1").arg(QKeySequence(QKeySequence::Find).toString(QKeySequence::NativeText).replace("F", "K")));
    search->setProperty("flat", true);
    search->setAccessibleName(tr("Command palette"));
    connect(search, &QPushButton::clicked, this, [this] { showPalette(); });
    h->addWidget(search);
    m_themeButton = iconButton(tokens().dark ? "sun" : "moon", tr("Switch light or dark theme"), bar);
    connect(m_themeButton, &QToolButton::clicked, this, [this] { setTheme(!tokens().dark); });
    h->addWidget(m_themeButton);
    return bar;
}

void MainWindow::buildDocks()
{
    using ads::CDockManager;
    CDockManager::setConfigFlags(CDockManager::DefaultOpaqueConfig);
    CDockManager::setConfigFlag(CDockManager::DisableStylesheet, true);
    CDockManager::setConfigFlag(CDockManager::FocusHighlighting, false);
    CDockManager::setConfigFlag(CDockManager::DockAreaHasUndockButton, false);
    CDockManager::setConfigFlag(CDockManager::DockAreaHasTabsMenuButton, false);
    CDockManager::setConfigFlag(CDockManager::AllTabsHaveCloseButton, false);
    CDockManager::setConfigFlag(CDockManager::ActiveTabHasCloseButton, false);
    CDockManager::setConfigFlag(CDockManager::DockAreaHasCloseButton, false);
    CDockManager::setConfigFlag(CDockManager::EqualSplitOnInsertion, false);
    m_docks = new CDockManager(nullptr);

    m_graphics = new GraphicsArea(m_session);
    m_charts = new ChartsPanel(m_session);
    m_outline = new Outline(m_session);
    m_inspector = new Inspector(m_session, m_graphics, m_charts);
    m_console = new ConsoleView;
    InspectorContext& ctx = m_inspector->context();
    ctx.openDialog = [this] { openCaseDialog(); };
    ctx.openPath = [this](const QString& p) { openCase(p); };
    ctx.samples = m_samples;
    m_inspector->refresh();

    auto dock = [](const QString& title, QWidget* w, const QString& name) {
        auto* d = new ads::CDockWidget(title);
        d->setObjectName(name);
        d->setWidget(w, ads::CDockWidget::ForceNoScrollArea);
        d->setFeature(ads::CDockWidget::DockWidgetClosable, false);
        d->setFeature(ads::CDockWidget::DockWidgetFloatable, true);
        d->setMinimumSizeHintMode(ads::CDockWidget::MinimumSizeHintFromDockWidgetMinimumSize);
        return d;
    };
    auto* graphicsDock = dock(tr("Graphics"), m_graphics, "graphics");
    graphicsDock->setFeature(ads::CDockWidget::NoTab, false);
    auto* centralArea = m_docks->setCentralWidget(graphicsDock);
    auto* outlineDock = dock(tr("Outline"), m_outline, "outline");
    outlineDock->setMinimumSize(180, 200);
    auto* outlineArea = m_docks->addDockWidget(ads::LeftDockWidgetArea, outlineDock);
    auto* inspectorDock = dock(tr("Inspector"), m_inspector, "inspector");
    inspectorDock->setMinimumSize(320, 200);
    m_docks->addDockWidget(ads::RightDockWidgetArea, inspectorDock, outlineArea);
    auto* consoleDock = dock(tr("Console"), m_console, "console");
    auto* bottomArea = m_docks->addDockWidget(ads::BottomDockWidgetArea, consoleDock, centralArea);
    auto* chartsDock = dock(tr("Charts"), m_charts, "charts");
    chartsDock->setMinimumSize(320, 180);
    m_docks->addDockWidget(ads::RightDockWidgetArea, chartsDock, bottomArea);
    // One size per splitter child: outline | inspector | graphics column.
    m_docks->setSplitterSizes(outlineArea, {205, 395, 1000});
    m_docks->setSplitterSizes(centralArea, {640, 260});
    m_docks->setSplitterSizes(bottomArea, {420, 560});

    connect(m_outline, &Outline::nodeSelected, m_inspector, &Inspector::show);
    connect(m_inspector, &Inspector::navigate, this, &MainWindow::navigate);
}

void MainWindow::buildMenus()
{
    QMenuBar* mb = menuBar();
    auto add = [this](QMenu* m, const QString& text, const QString& icon, QKeySequence key, auto fn) {
        QAction* a = m->addAction(text);
        if (!icon.isEmpty()) a->setIcon(Theme::instance().icon(icon));
        a->setIconVisibleInMenu(false);
        if (!key.isEmpty()) a->setShortcut(key);
        connect(a, &QAction::triggered, this, fn);
        addAction(a);
        m_commands << a;
        return a;
    };

    QMenu* file = mb->addMenu(tr("&File"));
    add(file, tr("Open Case…"), "folder-open", QKeySequence::Open, [this] { openCaseDialog(); });
    add(file, tr("Open Case with Data…"), "file-input", QKeySequence(), [this] {
        const QString cas = QFileDialog::getOpenFileName(this, tr("Fluent case"), {}, tr("Fluent case (*.cas.h5 *.cas)"));
        if (cas.isEmpty()) return;
        const QString dat = QFileDialog::getOpenFileName(this, tr("Fluent data"), QFileInfo(cas).absolutePath(), tr("Fluent data (*.dat.h5 *.dat)"));
        openCase(cas, dat);
    });
    m_recentMenu = file->addMenu(tr("Open Recent"));
    connect(m_recentMenu, &QMenu::aboutToShow, this, [this] {
        m_recentMenu->clear();
        for (const QString& p : m_session->recentCases()) {
            QAction* a = m_recentMenu->addAction(QFileInfo(p).fileName(), this, [this, p] { openCase(p); });
            a->setToolTip(p);
            a->setEnabled(QFileInfo::exists(p));
        }
        if (m_recentMenu->isEmpty()) m_recentMenu->addAction(tr("No recent cases"))->setEnabled(false);
    });
    file->addSeparator();
    add(file, tr("Save Viewport Image…"), "image-down", QKeySequence("Ctrl+Shift+S"), [this] { exportActiveViewport(); });
    file->addSeparator();
    QAction* quit = add(file, tr("Quit"), {}, QKeySequence::Quit, [this] { close(); });
    quit->setMenuRole(QAction::QuitRole);

    QMenu* edit = mb->addMenu(tr("&Edit"));
    m_undoAction = m_session->settings()->undoStack()->createUndoAction(this, tr("Undo"));
    m_undoAction->setShortcut(QKeySequence::Undo);
    m_redoAction = m_session->settings()->undoStack()->createRedoAction(this, tr("Redo"));
    m_redoAction->setShortcut(QKeySequence::Redo);
    edit->addAction(m_undoAction);
    edit->addAction(m_redoAction);
    m_commands << m_undoAction << m_redoAction;
    edit->addSeparator();
    add(edit, tr("Reset All Settings to the Case"), "rotate-ccw", QKeySequence(), [this] { m_session->settings()->resetAll(); });
    add(edit, tr("Show Modified Settings"), "list-filter", QKeySequence(), [this] {
        navigate("settings");
        m_inspector->showSettingsFiltered({});
    });

    QMenu* view = mb->addMenu(tr("&View"));
    add(view, tr("Command Palette"), "command", QKeySequence("Ctrl+K"), [this] { showPalette(); });
    view->addSeparator();
    add(view, tr("One Viewport"), "square", QKeySequence("Ctrl+1"), [this] { m_graphics->setViewLayout(GraphicsArea::Layout::Single); refreshToolbar(); });
    add(view, tr("Two Viewports Side by Side"), "columns-2", QKeySequence("Ctrl+2"), [this] { m_graphics->setViewLayout(GraphicsArea::Layout::SideBySide); refreshToolbar(); });
    add(view, tr("Two Viewports Stacked"), "rows-2", QKeySequence("Ctrl+3"), [this] { m_graphics->setViewLayout(GraphicsArea::Layout::Stacked); refreshToolbar(); });
    add(view, tr("Four Viewports"), "layout-grid", QKeySequence("Ctrl+4"), [this] { m_graphics->setViewLayout(GraphicsArea::Layout::Quad); refreshToolbar(); });
    QAction* link = add(view, tr("Link Viewport Cameras"), "focus", QKeySequence(), [this](bool on) { m_graphics->setLinkedCameras(on); });
    link->setCheckable(true);
    connect(view, &QMenu::aboutToShow, this, [this, link] { link->setChecked(m_graphics->linkedCameras()); });
    view->addSeparator();
    add(view, tr("Fit View"), "scan", QKeySequence(), [this] { if (auto* p = m_graphics->active()) p->view()->fitView(); });
    view->addSeparator();
    add(view, tr("Show Getting Started"), "info", QKeySequence(), [this] {
        QSettings().remove("ui/hideGettingStarted");
        navigate("setup/general");
        m_inspector->refresh();
    });
    add(view, tr("Dark Theme"), "moon", QKeySequence(), [this] { setTheme(true); });
    add(view, tr("Light Theme"), "sun", QKeySequence(), [this] { setTheme(false); });

    QMenu* run = mb->addMenu(tr("&Run"));
    m_runAction = add(run, tr("Run Calculation"), "play", QKeySequence("Ctrl+R"), [this] { m_session->run(); });
    m_stopAction = add(run, tr("Stop Calculation"), "circle-stop", QKeySequence("Ctrl+."), [this] { m_session->stop(); });
    run->addSeparator();
    add(run, tr("Run Calculation Settings"), "settings-2", QKeySequence(), [this] { navigate("solution/run"); });
    add(run, tr("Refresh Devices"), "refresh-cw", QKeySequence(), [this] { m_session->refreshDevices(); });

    QMenu* go = mb->addMenu(tr("&Go"));
    add(go, tr("General"), {}, QKeySequence("Ctrl+Shift+G"), [this] { navigate("setup/general"); });
    add(go, tr("Boundary Conditions"), {}, QKeySequence("Ctrl+Shift+B"), [this] { navigate("setup/boundaries"); });
    add(go, tr("Contours"), {}, QKeySequence("Ctrl+Shift+C"), [this] { navigate("results/contours"); });
    add(go, tr("Run History"), {}, QKeySequence("Ctrl+Shift+H"), [this] { navigate("runs/history"); });
    add(go, tr("All Case Settings"), {}, QKeySequence("Ctrl+Shift+A"), [this] { navigate("settings"); });

    QMenu* help = mb->addMenu(tr("&Help"));
    QAction* about = add(help, tr("About Kutta"), "info", QKeySequence(), [this] {
        QMessageBox::about(this, tr("About Kutta"),
                           tr("<b>Kutta %1</b><p>GPU solver workbench for Fluent cases.</p><p>Solver: %2</p><p>Third-party notices are in the licences folder next to the app.</p>")
                               .arg(QApplication::applicationVersion(), m_session->cliInfo().describe().toHtmlEscaped()));
    });
    about->setMenuRole(QAction::AboutRole);
    add(help, tr("Keyboard Shortcuts"), "command", QKeySequence::HelpContents, [this] {
        QMessageBox::information(this, tr("Keyboard shortcuts"),
                                 tr("Ctrl/Cmd-K  command palette\nCtrl/Cmd-R  run   Ctrl/Cmd-.  stop\nCtrl/Cmd-1…4  viewport layouts\n\nIn a viewport:\n"
                                    "drag  orbit   Shift-drag or middle-drag  pan   right-drag or wheel  zoom to cursor\n"
                                    "click  probe (all viewports show the same point)\nF  fit   1–7  front/back/left/right/top/bottom/iso\n\n"
                                    "In charts:\ndrag  pan   wheel  zoom iterations (Shift: values)   Shift-drag  box zoom\ndouble-click  fit   L  log/linear"));
    });
}

void MainWindow::refreshTitle()
{
    if (m_session->hasCase()) {
        const auto& c = m_session->caseSummary();
        setWindowTitle(QStringLiteral("%1 — Kutta").arg(m_session->caseName()));
        setWindowFilePath(m_session->casePath());
        m_caseLabel->setText(QStringLiteral("%1 · %2 cells").arg(m_session->caseName(), QLocale().toString(c.cells)));
    } else {
        setWindowTitle("Kutta");
        m_caseLabel->setText(tr("No case open"));
    }
}

void MainWindow::refreshToolbar()
{
    Session* s = m_session;
    const bool running = s->running();
    QString why;
    const bool can = s->canRun(&why);
    m_runButton->setText(running ? tr("Stop") : tr("Run"));
    m_runButton->setIcon(Theme::instance().icon(running ? "circle-stop" : "play", running ? tokens().red : tokens().accentText));
    m_runButton->setProperty("primary", !running);
    m_runButton->setProperty("danger", running);
    repolish(m_runButton);
    m_runButton->setEnabled(running ? s->state() != RunState::Stopping : can);
    m_runButton->setToolTip(running ? tr("Stop after the current iteration (Ctrl+.)") : can ? tr("Run the calculation (Ctrl+R)") : why);
    m_runAction->setEnabled(can);
    m_stopAction->setEnabled(running && s->state() != RunState::Stopping);
    m_iterations->setEnabled(s->hasCase() && !running);
    if (s->hasCase() && !m_iterations->hasFocus()) m_iterations->setText(QString::number(s->requestedIterations()));
    if (!s->hasCase()) m_iterations->clear();
    m_iterations->setProperty("modified", s->hasCase() && s->settings()->isModified("run/iterations"));
    repolish(m_iterations);

    const QSignalBlocker b(m_device);
    m_device->clear();
    m_device->addItem(Theme::instance().icon("cpu"), tr("Automatic device"), "auto");
    for (const auto& d : s->devices()) m_device->addItem(Theme::instance().icon(d.backend == "cpu" ? "cpu" : "microchip"), d.label(), d.backend);
    m_device->setCurrentIndex(std::max(0, m_device->findData(s->backend())));
    m_device->setEnabled(!running);

    QString state;
    const auto& l = s->live();
    switch (s->state()) {
    case RunState::Starting: state = tr("Starting kutta…"); break;
    case RunState::Stopping: state = tr("Stopping after this iteration…"); break;
    case RunState::Running:
        state = l.last ? tr("Iteration %1 / %2").arg(QLocale().toString(l.last->iter), QLocale().toString(l.requestedIterations))
                       : tr("Setting up on %1").arg(l.setup ? l.setup->device : QString());
        if (l.last && l.last->secPerIter) state += tr(" · %1 ms/iter").arg(*l.last->secPerIter * 1000, 0, 'f', 1);
        break;
    case RunState::Idle:
        if (!s->hasCase()) state.clear();
        else if (!can) state = why;
        else {
            QStringList marks;
            if (!s->approvedApproximations().isEmpty()) marks << tr("approximate");
            if (const int n = s->settings()->caseModifiedCount()) marks << (n == 1 ? tr("1 setting modified") : tr("%1 settings modified").arg(n));
            state = marks.isEmpty() ? tr("Ready") : tr("Ready · ") + marks.join(", ");
        }
        break;
    }
    // A run that only needs approval is a decision, not an error: amber, with a way to it.
    const bool blocked = s->state() == RunState::Idle && s->hasCase() && !can && !s->running();
    const bool approval = blocked && why.startsWith(tr("Needs your approval"));
    // The newest run failing stays visible until the next run: a status-bar
    // message that fades is easy to miss on a long solve.
    const auto& runs = s->history().runs();
    const RunRecord* failed = s->state() == RunState::Idle && s->hasCase() && !blocked && !runs.empty()
                                      && (runs.back().status == RunStatus::Failed || runs.back().status == RunStatus::Interrupted)
                                  ? &runs.back() : nullptr;
    auto withLink = [](const QString& text, const QString& href, const QString& link) {
        return QStringLiteral("%1 &nbsp;<a href=\"%2\" style=\"color:%3\">%4</a>").arg(text.toHtmlEscaped(), href, tokens().text.name(), link);
    };
    if (blocked) {
        const QString text = approval ? tr("Approval needed: %1").arg(why.section(':', 1).trimmed().chopped(1)) : state;
        m_stateLabel->setText(withLink(text, "setup/general", tr("Review")));
        m_stateLabel->setToolTip(state);
    } else if (failed) {
        m_stateLabel->setText(withLink(tr("%1 %2 · Ready to run again").arg(failed->title(), runStatusLabel(failed->status).toLower()), "solution/run", tr("Details")));
        m_stateLabel->setToolTip(failed->error);
    } else {
        m_stateLabel->setText(state.toHtmlEscaped());
        m_stateLabel->setToolTip(state);
    }
    m_stateLabel->setProperty("role", blocked ? (approval ? "warning" : "error") : failed ? "error" : "muted");
    repolish(m_stateLabel);

    const int idx = m_graphics->viewLayout() == GraphicsArea::Layout::Single ? 0 : m_graphics->viewLayout() == GraphicsArea::Layout::Quad ? 2 : 1;
    for (int i = 0; i < 3; ++i) m_layoutButtons[i]->setChecked(i == idx);
}

void MainWindow::navigate(const QString& id)
{
    m_outline->select(id);
    if (m_inspector->current() != id) m_inspector->show(id);
}

void MainWindow::openCaseDialog()
{
    QString start = m_session->hasCase() ? QFileInfo(m_session->casePath()).absolutePath() : QString();
    const QString f = QFileDialog::getOpenFileName(this, tr("Open Fluent case"), start, tr("Fluent case (*.cas.h5 *.cas *.cas.gz);;All files (*)"));
    if (!f.isEmpty()) openCase(f);
}

void MainWindow::openCase(const QString& path, const QString& data)
{
    QString dataPath = data;
    if (dataPath.isEmpty()) {
        // Pair case.cas.h5 with case.dat.h5 when it sits beside it, as Fluent does.
        QString guess = path;
        guess.replace(".cas.h5", ".dat.h5").replace(".cas", ".dat");
        if (guess != path && QFileInfo::exists(guess)) dataPath = guess;
    }
    statusBar()->showMessage(tr("Reading %1…").arg(QFileInfo(path).fileName()));
    m_session->openCase(path, dataPath);
}

void MainWindow::showPalette(const QString& query) { m_palette->open(query); }

void MainWindow::setTheme(bool dark)
{
    Theme::instance().apply(*qApp, dark);
    QSettings().setValue("theme/dark", dark);
}

void MainWindow::exportActiveViewport()
{
    auto* p = m_graphics->active();
    if (!p || !p->scene()) return;
    const QString f = QFileDialog::getSaveFileName(this, tr("Save viewport image"), p->datasetLabel() + ".png", tr("PNG image (*.png)"));
    if (!f.isEmpty()) p->exportImage().save(f, "PNG");
}

std::vector<PaletteEntry> MainWindow::paletteEntries()
{
    std::vector<PaletteEntry> out;
    for (QAction* a : m_commands) {
        if (!a->isEnabled() || a->text().isEmpty()) continue;
        QString text = a->text();
        text.remove('&');
        out.push_back({tr("Command"), text, a->shortcut().toString(QKeySequence::NativeText), [a] { a->trigger(); }});
    }
    for (const auto& n : m_outline->nodes()) {
        const QString id = n.id;
        out.push_back({tr("Go to"), n.path.section(" ▸ ", -1), n.path, [this, id] { navigate(id); }});
    }
    if (m_session->hasCase()) {
        for (const auto& s : m_session->settings()->settings()) {
            if (!s.editable) continue;
            const QString text = settingLabel(s) + (s.owner.isEmpty() ? QString() : " · " + s.owner);
            const QString path = s.path;
            QString page = s.section == "boundary" ? "setup/boundary/" + s.owner
                         : s.section == "materials" ? "setup/material/" + s.owner
                         : s.section == "reference" ? "setup/reference"
                         : s.section == "methods" ? "solution/methods"
                         : s.section == "controls" ? "solution/controls"
                         : s.section == "run" ? "solution/run" : "settings";
            out.push_back({tr("Setting"), text, displayValue(s, m_session->settings()->value(path)), [this, page] { navigate(page); }});
        }
        const auto& runs = m_session->history().runs();
        for (auto it = runs.rbegin(); it != runs.rend(); ++it) {
            if (it->status != RunStatus::Completed && it->status != RunStatus::Stopped) continue;
            const QString id = it->id;
            out.push_back({tr("Run"), tr("Show %1").arg(it->title()), runStatusLabel(it->status), [this, id] { m_graphics->show("run:" + id); }});
        }
    }
    for (const QString& p : m_samples) {
        out.push_back({tr("Sample"), tr("Open %1").arg(QFileInfo(p).fileName()), {}, [this, p] { openCase(p); }});
    }
    return out;
}

void MainWindow::closeEvent(QCloseEvent* e)
{
    if (m_session->running()) {
        const auto r = QMessageBox::question(this, tr("Calculation running"), tr("Stop the calculation and quit? Kutta finishes the current iteration and saves the run."),
                                             QMessageBox::Cancel | QMessageBox::Yes, QMessageBox::Cancel);
        if (r != QMessageBox::Yes) {
            e->ignore();
            return;
        }
        m_session->stop();
    }
    QSettings().setValue("window/geometry", saveGeometry());
    QMainWindow::closeEvent(e);
}

void MainWindow::dragEnterEvent(QDragEnterEvent* e)
{
    for (const QUrl& u : e->mimeData()->urls())
        if (u.isLocalFile() && isCaseFile(u.toLocalFile())) {
            e->acceptProposedAction();
            return;
        }
}

void MainWindow::dropEvent(QDropEvent* e)
{
    QString cas, dat;
    for (const QUrl& u : e->mimeData()->urls()) {
        const QString f = u.toLocalFile();
        if (isCaseFile(f)) cas = f;
        else if (f.contains(".dat")) dat = f;
    }
    if (!cas.isEmpty()) openCase(cas, dat);
}

} // namespace kapp
