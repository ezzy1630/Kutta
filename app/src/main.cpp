// Kutta desktop app. Normal start opens the workbench; --scenario drives it to a
// named state and --screenshot saves the window (including GPU viewports) for
// visual review and the CI smoke test. --benchmark measures frame times and
// idle cost.

#include "core/CliClient.h"
#include "core/Session.h"
#include "ui/ChartsPanel.h"
#include "ui/CommandPalette.h"
#include "ui/GraphicsArea.h"
#include "ui/Inspector.h"
#include "ui/MainWindow.h"
#include "ui/Outline.h"
#include "ui/ScreenChecks.h"
#include "ui/Theme.h"

#include <QAbstractButton>
#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QScrollBar>
#include <QSet>
#include <QSettings>
#include <QStandardPaths>
#include <QStyleHints>
#include <QTableWidget>
#include <QTemporaryDir>
#include <QTimer>

#include <cstdio>
#include <deque>
#include <functional>

#if defined(Q_OS_UNIX)
#include <sys/resource.h>
#endif

using namespace kapp;

namespace {

// Steps run one after another; each calls next() when its condition holds.
class Scenario : public QObject {
public:
    using Step = std::function<void(std::function<void()> next)>;
    explicit Scenario(QObject* parent) : QObject(parent) {}
    void add(Step s) { m_steps.push_back(std::move(s)); }
    void start() { QTimer::singleShot(0, this, [this] { advance(); }); }
    void fail(const QString& why)
    {
        std::fprintf(stderr, "scenario failed: %s\n", qPrintable(why));
        QCoreApplication::exit(2);
    }

    // Poll a condition (scenario mode only) with a timeout.
    void waitFor(std::function<bool()> cond, std::function<void()> next, const QString& what, int timeoutMs = 120000)
    {
        auto* t = new QTimer(this);
        auto* clock = new QElapsedTimer;
        clock->start();
        t->setInterval(40);
        connect(t, &QTimer::timeout, this, [=, this] {
            if (cond()) {
                t->deleteLater();
                delete clock;
                next();
            } else if (clock->elapsed() > timeoutMs) {
                t->deleteLater();
                delete clock;
                fail("timed out waiting for " + what);
            }
        });
        t->start();
    }

private:
    void advance()
    {
        if (m_steps.empty()) return;
        Step s = std::move(m_steps.front());
        m_steps.pop_front();
        s([this] { QTimer::singleShot(0, this, [this] { advance(); }); });
    }
    std::deque<Step> m_steps;
};

QSize parseSize(const QString& s, QSize fallback)
{
    const QStringList p = s.split('x');
    if (p.size() != 2) return fallback;
    return {p[0].toInt(), p[1].toInt()};
}

// Placeholder case files for the fixture CLI's samples: it identifies cases by
// file name and never reads them. A real CLI needs real files.
QStringList makeSamples(const QString& root, bool fixture)
{
    if (!fixture) return {};
    QDir().mkpath(root);
    QStringList out;
    // Scenarios address the first two by position (elbow, car).
    for (const char* name : {"mixing_elbow.cas.h5", "mixing_elbow.dat.h5", "fsae_car.cas.h5", "fsae_car.dat.h5",
                             "naca4412_ground.cas", "elbow_param.cas.h5", "cylinder_tutorial.cas", "naca4412_ground.msh", "fsae_car.msh"}) {
        const QString p = QDir(root).filePath(name);
        if (!QFileInfo::exists(p)) {
            QFile f(p);
            if (f.open(QIODevice::WriteOnly)) f.write("Kutta fixture sample placeholder: the fixture CLI serves this case by name.\n");
        }
        if (QString(name).contains(".cas")) out << p;
    }
    return out;
}

double processCpuSeconds()
{
#if defined(Q_OS_UNIX)
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    return double(u.ru_utime.tv_sec + u.ru_stime.tv_sec) + double(u.ru_utime.tv_usec + u.ru_stime.tv_usec) / 1e6;
#else
    return 0;
#endif
}

} // namespace

int main(int argc, char** argv)
{
    QApplication::setHighDpiScaleFactorRoundingPolicy(Qt::HighDpiScaleFactorRoundingPolicy::PassThrough);
    QApplication app(argc, argv);
    QApplication::setApplicationName("Kutta");
    QApplication::setOrganizationName("Kutta");
    QApplication::setOrganizationDomain("kutta.app");
    QApplication::setApplicationVersion(KUTTA_APP_VERSION);

    QCommandLineParser cli;
    cli.setApplicationDescription("Kutta: GPU solver workbench for Fluent cases");
    cli.addHelpOption();
    cli.addVersionOption();
    cli.addPositionalArgument("case", "Fluent case to open (.cas.h5)");
    QCommandLineOption cliOpt("cli", "Path to the kutta CLI (default: next to the app, $KUTTA_CLI, or the dev fixture).", "path");
    QCommandLineOption dataOpt("data-dir", "Where runs are stored (default: the app's data folder).", "dir");
    QCommandLineOption scenarioOpt("scenario", "Drive the app to a named state (see --list-scenarios).", "name");
    QCommandLineOption shotOpt("screenshot", "Save a PNG of the window after the scenario, then quit.", "file");
    QCommandLineOption themeOpt("theme", "dark or light.", "theme");
    QCommandLineOption sizeOpt("size", "Window size in logical pixels, e.g. 1600x1000.", "WxH");
    QCommandLineOption listOpt("list-scenarios", "List scenario names.");
    QCommandLineOption benchOpt("benchmark", "Measure viewport frame times over N frames and idle cost, print JSON, quit.", "frames");
    QCommandLineOption caseOpt("open-data", "Data file to open with the case.", "dat");
    cli.addOptions({cliOpt, dataOpt, scenarioOpt, shotOpt, themeOpt, sizeOpt, listOpt, benchOpt, caseOpt});
    cli.process(app);

    // case-* scenarios take any case (the positional argument; no data file
    // unless --open-data) and are meant for the real CLI: the mesh, the
    // approvals, the Methods page, a run with convergence checking on
    // (KUTTA_SCENARIO_ITERATIONS, default 300) and a run stopped by the user;
    // KUTTA_SCENARIO_PAGE shows another inspector page at the end.
    // elbow-* scenarios are meant for the real CLI and the Ansys mixing elbow
    // (pass the case path; --open-data defaults to the .dat.h5 beside it).
    // Post-processing ones need the real CLI: the fixture has no volume fields.
    // car-* and bench use the fixture CLI's synthetic FSAE car.
    // new-case* start from a bare mesh (the positional argument, else the
    // fixture's NACA wing mesh): the New Case page, the created case after a
    // run, and that case saved with a changed setting and reopened.
    // new-case-car is the New Case page for the fixture's synthetic car mesh
    // (3D, with wheels). KUTTA_SCENARIO_CASE_NAME names the created case
    // (default: the mesh's name), e.g. to try a long one.
    const QStringList scenarios = {"start", "open-failed", "elbow-open", "elbow-setup", "palette", "settings", "elbow-live", "elbow-run", "elbow-fluent",
                                   "elbow-post", "elbow-xy", "elbow-integrals", "elbow-history", "elbow-compare", "elbow-monitors",
                                   "car-results", "car-xy", "car-compare", "car-failed", "car-aero", "car-monitors", "car-study", "param-run", "cylinder-run",
                                   "wing-models", "wing-aero", "wing-monitors", "bench",
                                   "case-mesh", "case-approvals", "case-methods", "case-run", "case-stop",
                                   "new-case", "new-case-car", "new-case-run", "new-case-save"};
    if (cli.isSet(listOpt)) {
        for (const QString& s : scenarios) std::printf("%s\n", qPrintable(s));
        return 0;
    }
    const QString scenario = cli.value(scenarioOpt);
    const bool scripted = !scenario.isEmpty() || cli.isSet(benchOpt);
    if (!scenario.isEmpty() && !scenarios.contains(scenario)) {
        std::fprintf(stderr, "unknown scenario %s (try --list-scenarios)\n", qPrintable(scenario));
        return 2;
    }

    // Scripted runs use throwaway settings and run folders so every capture
    // starts from the same state.
    std::unique_ptr<QTemporaryDir> temp;
    QString dataRoot = cli.value(dataOpt);
    if (scripted) {
        QApplication::setOrganizationName("KuttaScenario");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        if (dataRoot.isEmpty()) {
            temp = std::make_unique<QTemporaryDir>();
            dataRoot = temp->path();
        }
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, QDir(dataRoot).filePath("settings"));
        QSettings().clear();
        // Fixture replays run quickly in scripted mode.
        if (qEnvironmentVariableIsEmpty("KUTTA_FIXTURE_SPEED")) qputenv("KUTTA_FIXTURE_SPEED", "120");
    }
    if (dataRoot.isEmpty()) dataRoot = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation);

    bool dark = app.styleHints()->colorScheme() != Qt::ColorScheme::Light;
    if (!scripted && QSettings().contains("theme/dark")) dark = QSettings().value("theme/dark").toBool();
    if (cli.isSet(themeOpt)) dark = cli.value(themeOpt) != "light";
    Theme::instance().apply(app, dark);

#ifdef KUTTA_DEV_FIXTURE_CLI
    const QString devFixture = QStringLiteral(KUTTA_DEV_FIXTURE_CLI);
#else
    const QString devFixture;
#endif
#ifdef KUTTA_DEV_REAL_CLI
    const QString devReal = QStringLiteral(KUTTA_DEV_REAL_CLI);
#else
    const QString devReal;
#endif
    // Scenarios are reproducible: they use the fixture unless --cli or KUTTA_CLI says otherwise.
    const bool fixtureScenario = scripted && !cli.isSet(cliOpt) && qEnvironmentVariableIsEmpty("KUTTA_CLI") && QFileInfo(devFixture).isExecutable();
    const CliInfo info = fixtureScenario ? locateCli(devFixture) : locateCli(cli.value(cliOpt), devReal, devFixture);
    if (info.program.isEmpty()) {
        std::fprintf(stderr, "kutta CLI not found; pass --cli or set KUTTA_CLI\n");
    }
    auto* session = new Session(info, dataRoot, &app);
    const QStringList samples = makeSamples(QDir(dataRoot).filePath("samples"), info.fixture);
    MainWindow w(session, samples);
    w.resize(parseSize(cli.value(sizeOpt), QSize(1600, 1000)));
    session->refreshDevices();

    if (!scripted) {
        w.show();
        if (!cli.positionalArguments().isEmpty()) w.openCase(cli.positionalArguments().first(), cli.value(caseOpt));
        return app.exec();
    }

    // ---- Scenarios -------------------------------------------------------------------
    auto* sc = new Scenario(&app);
    const QString elbow = samples.value(0), car = samples.value(1);
    auto sample = [&samples](const QString& stem) {
        for (const QString& p : samples)
            if (QFileInfo(p).fileName().startsWith(stem + ".")) return p;
        return QString();
    };
    QString caseArg = cli.positionalArguments().value(0);
    auto open = [session, sc](const QString& path, const QString& data = {}) {
        return [session, sc, path, data](std::function<void()> next) {
            if (path.isEmpty()) return sc->fail("no case for this scenario (fixture CLI samples need the fixture CLI; pass a case path)");
            bool* done = new bool(false);
            QObject::connect(session, &Session::caseOpened, sc, [done] { *done = true; }, Qt::SingleShotConnection);
            QObject::connect(session, &Session::caseFailed, sc, [sc](const QString& m) { sc->fail(m); }, Qt::SingleShotConnection);
            // Scenarios say exactly which data file to use (none means none).
            session->openCase(path, data);
            sc->waitFor([done] { return *done; }, [done, next] { delete done; next(); }, "case to open");
        };
    };
    // Runs in scenarios opt into whatever approximations the case needs; the
    // results then carry the APPROXIMATE marks that this exercises.
    auto approveAll = [session](std::function<void()> next) {
        for (const QString& a : session->caseSummary().approximations()) session->setApproximationApproved(a, true);
        next();
    };
    auto runToEnd = [session, sc](std::function<void()> next) {
        const QString id = session->run();
        if (id.isEmpty()) return sc->fail("run refused");
        sc->waitFor([session] { return !session->running(); }, next, "run to finish", 300000);
    };
    auto showScene = [&w, sc](int panel, const QString& key) {
        return [&w, sc, panel, key](std::function<void()> next) {
            w.graphics()->show(key, panel);
            sc->waitFor([&w, panel] { return w.graphics()->panel(panel)->scene() != nullptr; }, next, "scene " + key);
        };
    };
    auto display = [&w](int panel, std::function<void(ViewDisplay&)> fn) {
        return [&w, panel, fn](std::function<void()> next) {
            ViewDisplay d = w.graphics()->panel(panel)->display();
            fn(d);
            w.graphics()->panel(panel)->setDisplay(d);
            next();
        };
    };
    auto go = [&w](const QString& id) { return [&w, id](std::function<void()> next) { w.navigate(id); next(); }; };
    auto layout = [&w](GraphicsArea::Layout l) { return [&w, l](std::function<void()> next) { w.graphics()->setViewLayout(l); next(); }; };
    auto preset = [&w](int panel, Camera::Preset p) { return [&w, panel, p](std::function<void()> next) { w.graphics()->panel(panel)->view()->applyPreset(p); next(); }; };
    auto settle = [](int ms) { return [ms](std::function<void()> next) { QTimer::singleShot(ms, next); }; };
    auto field = [](const QString& q, FieldRole role = FieldRole::Kutta) { return FieldRef{q, role, -1}; };
    auto post = [session, sc](PostKind kind, QJsonObject op, std::function<void(const PostOutput&)> then = {}) {
        return [session, sc, kind, op, then](std::function<void()> next) {
            const QString run = session->history().latestFinished() ? session->history().latestFinished()->id : QString();
            session->postProcess(run, kind, op, [sc, next, then](Result<PostOutput> r) {
                if (!r) return sc->fail("post-processing: " + r.error());
                if (then) then(r.value());
                next();
            });
        };
    };

    const QString elbowCase = caseArg.isEmpty() ? elbow : caseArg;
    QString elbowData = cli.value(caseOpt);
    if (elbowData.isEmpty() && QFileInfo::exists(QString(elbowCase).replace(".cas.h5", ".dat.h5")))
        elbowData = QString(elbowCase).replace(".cas.h5", ".dat.h5");
    // Find a widget on the current inspector page by its accessible name or text.
    auto inspectorButton = [&w](const QString& text) -> QAbstractButton* {
        for (auto* b : w.inspector()->findChildren<QAbstractButton*>())
            if (b->isVisible() && (b->text() == text || b->accessibleName() == text)) return b;
        return nullptr;
    };
    auto runWith = [session, sc](int iterations, double coldInlet = 0) {
        return [session, sc, iterations, coldInlet](std::function<void()> next) {
            session->settings()->setValue("run/iterations", iterations);
            if (coldInlet > 0) session->settings()->setValue("boundary/cold-inlet/velocity-magnitude", coldInlet);
            const QString id = session->run();
            if (id.isEmpty()) return sc->fail("run refused");
            sc->waitFor([session] { return !session->running(); }, next, "run to finish", 300000);
        };
    };
    if (scenario == "start") {
        sc->add(go("setup/general"));
    } else if (scenario == "open-failed") {
        // A file that is not a readable case: the reason shows in place.
        sc->add([session, sc, dataRoot](auto next) {
            session->openCase(QDir(dataRoot).filePath("front_wing.cas.h5"));
            sc->waitFor([session] { return !session->openError().message.isEmpty(); }, next, "open to fail");
        });
        sc->add(go("setup/general"));
    } else if (scenario == "elbow-open") {
        // What a first-time user sees right after opening: the case's mesh
        // and the getting-started steps.
        sc->add(open(elbowCase, elbowData));
        sc->add([&w, sc](auto next) {
            sc->waitFor([&w] { return w.graphics()->panel(0)->scene() != nullptr; }, next, "mesh to show");
        });
        sc->add(go("setup/general"));
    } else if (scenario == "elbow-setup") {
        sc->add(open(elbowCase));
        sc->add([session](auto next) {
            session->settings()->setValue("boundary/cold-inlet/velocity-magnitude", 1.4);
            session->settings()->setValue("boundary/cold-inlet/turb-intensity", 0.08);
            next();
        });
        sc->add(go("setup/boundary/cold-inlet"));
    } else if (scenario == "palette") {
        sc->add(open(elbowCase));
        sc->add([&w](auto next) { w.showPalette("inlet"); next(); });
    } else if (scenario == "elbow-run" || scenario == "elbow-fluent") {
        const bool withData = scenario == "elbow-fluent";
        sc->add(open(elbowCase, withData ? elbowData : QString()));
        sc->add(approveAll);
        sc->add(runWith(300));
        sc->add(showScene(0, "run:run-0001"));
        if (withData) {
            sc->add(layout(GraphicsArea::Layout::SideBySide));
            sc->add(showScene(1, "run:run-0001"));
            sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
            sc->add(display(1, [&](ViewDisplay& d) { d.field = field("velocity_magnitude", FieldRole::Delta); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
            sc->add(preset(0, Camera::Preset::Top));
            sc->add(preset(1, Camera::Preset::Top));
            sc->add(go("results/contours"));
        } else {
            sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.bands = 16; d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
            sc->add(preset(0, Camera::Preset::Top));
            sc->add(go("solution/run"));
        }
    } else if (scenario == "elbow-monitors") {
        // The seeded elbow stopping on its own residual criteria.
        sc->add(open(elbowCase, elbowData));
        sc->add(approveAll);
        sc->add([session](auto next) {
            session->settings()->setValue("monitors/check-convergence", true);
            next();
        });
        sc->add(runWith(1000));
        sc->add(showScene(0, "run:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.bands = 16; d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
        sc->add(preset(0, Camera::Preset::Top));
        sc->add(go("solution/monitors"));
    } else if (scenario == "elbow-post") {
        // Volume post-processing needs the real CLI: pass --cli, the case and --open-data.
        sc->add(open(elbowCase, elbowData));
        sc->add(approveAll);
        sc->add(runWith(150));
        sc->add(showScene(0, "run:run-0001"));
        sc->add(post(PostKind::Streamlines,
                     {{"name", "Pathlines"}, {"mode", "streamline"}, {"direction", "forward"},
                      {"seeds", QJsonObject{{"type", "surface"}, {"zones", QJsonArray{"cold-inlet", "hot-inlet"}}, {"count", 36}, {"distribution", "uniform"}, {"offset", 0}}},
                      {"fields", QJsonArray{"velocity_magnitude", "pressure"}}}));
        sc->add(post(PostKind::IsoSurface,
                     {{"name", "Velocity Magnitude = 0.9"}, {"field", "velocity_magnitude"}, {"value", 0.9}, {"fields", QJsonArray{"pressure"}}}));
        sc->add(post(PostKind::Probe,
                     {{"name", "Section"}, {"type", "section"}, {"zones", QJsonArray{"wall-elbow"}}, {"origin", QJsonArray{0, 0, 0.025}},
                      {"normal", QJsonArray{0, 0, 1}}, {"fields", QJsonArray{"pressure", "pressure@fluent"}}},
                     [&w](const PostOutput& out) { w.charts()->plotProbe(out, "Static Pressure on the elbow wall, z = 0.025 m · Run 1", -1); }));
        sc->add(layout(GraphicsArea::Layout::Quad));
        for (int i = 1; i < 4; ++i) sc->add(showScene(i, "run:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.lineArrows = true; d.visible = {"post:Pathlines", "zone:cold-inlet", "zone:hot-inlet", "zone:outlet"}; d.visibleInitialized = true; }));
        sc->add(display(1, [&](ViewDisplay& d) { d.field = field("pressure"); d.visible = {"post:Velocity Magnitude = 0.9"}; d.visibleInitialized = true; }));
        sc->add(display(2, [&](ViewDisplay& d) { d.field = field("pressure", FieldRole::Delta); d.visible = {"zone:symmetry-xyplane", "post:Section"}; d.visibleInitialized = true; }));
        sc->add(display(3, [&](ViewDisplay& d) { d.field = field("total_pressure_coefficient"); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
        for (int i = 0; i < 4; ++i) sc->add(preset(i, i == 1 ? Camera::Preset::Isometric : Camera::Preset::Top));
        sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::XY); next(); });
        sc->add(go("results/pathlines"));
    } else if (scenario == "elbow-live") {
        qputenv("KUTTA_FIXTURE_HOLD_AT", "150");
        qputenv("KUTTA_FIXTURE_SPEED", "20");
        sc->add(open(elbowCase));
        sc->add(approveAll);
        sc->add([session](auto next) {
            session->setLiveViewInterval(20);
            session->settings()->setValue("run/iterations", 2000);  // still running at the screenshot
            next();
        });
        sc->add([session, sc](auto next) {
            session->run();
            sc->waitFor([session] { return session->live().lastSnapshotIter >= 120 && session->live().last && session->live().last->iter >= 150; }, next, "live snapshot");
        });
        sc->add(showScene(0, "live"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
        sc->add(preset(0, Camera::Preset::Top));
        sc->add(go("solution/run"));
    } else if (scenario == "elbow-xy") {
        sc->add(open(elbowCase, elbowData));
        sc->add(approveAll);
        sc->add(runWith(150));
        sc->add(showScene(0, "run:run-0001"));
        sc->add(go("results/xy"));
        // Across the outlet leg, where the hot and cold streams mix.
        sc->add(post(PostKind::Probe,
                     {{"name", "XY line"}, {"type", "line"}, {"count", 200}, {"start", QJsonArray{0.1, 0.1, 0.012}}, {"end", QJsonArray{0.2, 0.1, 0.012}},
                      {"fields", QJsonArray{"velocity_magnitude", "velocity_magnitude@fluent"}}},
                     [&w](const PostOutput& out) { w.charts()->plotProbe(out, "Velocity Magnitude across the outlet leg, y = 0.1 m · Run 1", 0); }));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.visible = {"zone:symmetry-xyplane", "post:XY line"}; d.visibleInitialized = true; }));
        sc->add(preset(0, Camera::Preset::Top));
        sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::XY); next(); });
    } else if (scenario == "elbow-integrals") {
        sc->add(open(elbowCase, elbowData));
        sc->add(approveAll);
        sc->add(runWith(150));
        sc->add(showScene(0, "run:run-0001"));
        sc->add(go("results/integrals"));
        // Through the page itself: Fluent column on, Compute.
        sc->add([&w, sc, inspectorButton](auto next) {
            auto* fluent = inspectorButton(QObject::tr("Also compute on Fluent's field"));
            auto* go = inspectorButton(QObject::tr("Compute"));
            if (!fluent || !go) return sc->fail("surface integrals page controls not found");
            fluent->setChecked(true);
            go->click();
            sc->waitFor([&w] {
                for (auto* t : w.inspector()->findChildren<QTableWidget*>())
                    if (t->accessibleName() == QObject::tr("Surface integral results") && t->isVisible()) return true;
                return false;
            }, next, "surface integrals");
        });
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure"); d.visible = {"zone:cold-inlet", "zone:hot-inlet", "zone:outlet", "zone:wall-elbow"}; d.visibleInitialized = true; }));
        sc->add(preset(0, Camera::Preset::Isometric));
    } else if (scenario == "elbow-history" || scenario == "elbow-compare") {
        sc->add(open(elbowCase));
        sc->add(approveAll);
        sc->add(runWith(300));
        sc->add(runWith(300, 1.4));  // a faster cold stream
        if (scenario == "elbow-compare") {
            sc->add([&w](auto next) {
                QSettings().setValue("compare/a", "run-0002");
                QSettings().setValue("compare/b", "run-0001");
                w.graphics()->addDataset("diff:run-0002:run-0001", "Run 2 − Run 1");
                next();
            });
            sc->add(layout(GraphicsArea::Layout::SideBySide));
            sc->add(showScene(0, "run:run-0002"));
            sc->add(showScene(1, "diff:run-0002:run-0001"));
            sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
            sc->add(display(1, [&](ViewDisplay& d) { d.field = field("velocity_magnitude", FieldRole::Run); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
            for (int i = 0; i < 2; ++i) sc->add(preset(i, Camera::Preset::Top));
            sc->add(go("runs/compare"));
            sc->add([&w](auto next) { w.charts()->showRuns({"run-0002", "run-0001"}); next(); });
        } else {
            sc->add(showScene(0, "run:run-0002"));
            sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
            sc->add(preset(0, Camera::Preset::Top));
            sc->add(go("runs/history"));
            sc->add([&w](auto next) { w.charts()->showRuns({"run-0002", "run-0001"}); next(); });
        }
    } else if (scenario == "car-results" || scenario == "car-xy") {
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add(runToEnd);
        sc->add(showScene(0, "run:run-0001"));
        if (scenario == "car-results") {
            sc->add(post(PostKind::Streamlines,
                         {{"name", "Pathlines"}, {"mode", "streamline"}, {"direction", "forward"},
                          {"seeds", QJsonObject{{"type", "rake"}, {"start", QJsonArray{-0.6, -0.7, 0.06}}, {"end", QJsonArray{-0.6, 0.7, 0.5}}, {"count", 48}}},
                          {"fields", QJsonArray{"velocity_magnitude", "pressure_coefficient"}}}));
            sc->add(layout(GraphicsArea::Layout::Quad));
            for (int i = 1; i < 4; ++i) sc->add(showScene(i, "run:run-0001"));
            sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); d.colormap = ColormapId::FluentRainbow; d.colormapFollowsRole = false; d.bands = 0; }));
            sc->add(display(1, [&](ViewDisplay& d) {
                d.field = field("velocity_magnitude");
                d.visible.clear();
                for (const auto& it : w.graphics()->panel(1)->scene()->items)
                    if (it.kind == ItemKind::Slice || (it.zoneType == "wall" && !it.id.contains("ground"))) d.visible << it.id;
                d.visibleInitialized = true;
            }));
            sc->add(preset(1, Camera::Preset::Left));
            sc->add(display(2, [&](ViewDisplay& d) {
                d.field = field("velocity_magnitude");
                d.contours = true;
                d.visible.clear();
                for (const auto& it : w.graphics()->panel(2)->scene()->items)
                    if (it.kind == ItemKind::Streamlines || (it.zoneType == "wall" && !it.id.contains("ground"))) d.visible << it.id;
                d.visibleInitialized = true;
            }));
            sc->add(display(3, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); d.bands = 14; d.edges = true; }));
            sc->add(preset(3, Camera::Preset::Bottom));
            sc->add(go("results/contours"));
        } else {
            sc->add(go("results/xy"));
            // Plot through the same path the XY page uses.
            sc->add(post(PostKind::Probe,
                         {{"name", "XY line"}, {"type", "line"}, {"count", 300}, {"start", QJsonArray{-1, 0, 0.1}}, {"end", QJsonArray{4, 0, 0.1}},
                          {"fields", QJsonArray{"velocity_magnitude", "total_pressure_coefficient"}}},
                         [&w](const PostOutput& out) { w.charts()->plotProbe(out, "Along y = 0, z = 0.1 m · Run 1", -1); }));
        }
    } else if (scenario == "car-failed") {
        // A run that fails part-way, as when the disk fills.
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add([session, sc](auto next) {
            qputenv("KUTTA_FIXTURE_FAIL_AT", "60");
            qputenv("KUTTA_FIXTURE_FAIL_MESSAGE", "cannot write fields.bin in the run folder: No space left on device");
            if (session->run().isEmpty()) return sc->fail("run refused");
            sc->waitFor([session] { return !session->running(); }, next, "run to fail", 300000);
        });
        sc->add(go("solution/run"));
    } else if (scenario == "car-compare") {
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add(runToEnd);
        sc->add([session](auto next) {
            session->settings()->setValue("boundary/inlet/velocity-magnitude", 20.0);
            next();
        });
        sc->add(runToEnd);
        sc->add([&w, session](auto next) {
            QSettings().setValue("compare/a", "run-0002");
            QSettings().setValue("compare/b", "run-0001");
            w.graphics()->addDataset("diff:run-0002:run-0001", "Run 2 − Run 1");
            next();
        });
        sc->add(layout(GraphicsArea::Layout::SideBySide));
        sc->add(showScene(0, "run:run-0002"));
        sc->add(showScene(1, "diff:run-0002:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); }));
        sc->add(display(1, [&](ViewDisplay& d) { d.field = field("pressure_coefficient", FieldRole::Run); }));
        sc->add(go("runs/compare"));
        sc->add([&w](auto next) { w.charts()->showRuns({"run-0002", "run-0001"}); next(); });
    } else if (scenario == "car-aero") {
        // Coefficients, breakdown by group and balance of one finished run.
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add([session](auto next) {
            AeroSpec spec = session->aeroSpec();
            spec.frontAxle = QVector3D(-0.78f, 0, 0);  // the synthetic car's front tyre contact
            spec.wheelbase = 1.58;
            spec.momentCentre = *spec.frontAxle;
            session->setAeroSpec(spec);
            next();
        });
        sc->add(runToEnd);
        sc->add(showScene(0, "run:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); }));
        sc->add(preset(0, Camera::Preset::Isometric));
        sc->add(go("results/aero"));
        sc->add([sc, session](auto next) {
            sc->waitFor([session] { return session->aeroIfComputed("run-0001").has_value(); }, next, "aero coefficients");
        });
    } else if (scenario == "car-monitors") {
        // A run stopped by its convergence criteria: residuals and two reports.
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add([session](auto next) {
            for (const char* p : {"monitors/check-convergence", "monitors/report/downforce/active", "monitors/report/drag/active"})
                session->settings()->setValue(p, true);
            next();
        });
        sc->add(runToEnd);
        sc->add(showScene(0, "run:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); }));
        sc->add(go("solution/monitors"));
        sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::Reports); next(); });
    } else if (scenario == "car-study") {
        // An inlet speed sweep, run back to back, with the reference velocity following.
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add([session, sc](auto next) {
            AeroSpec spec = session->aeroSpec();
            spec.frontAxle = QVector3D(-0.78f, 0, 0);
            spec.wheelbase = 1.58;
            session->setAeroSpec(spec);
            session->settings()->setValue("run/iterations", 600);
            const QString p = "boundary/inlet/velocity-magnitude";
            auto study = makeSweep(*session->settings(), p, {12, 16, 20, 24, 28}, referenceFollowers(*session->settings(), p));
            if (!study) return sc->fail(study.error());
            study->name = "Speed sweep";
            QString why;
            const QString id = session->addStudy(study.value(), &why);
            if (id.isEmpty() || !session->startStudy(id, &why)) return sc->fail(why);
            QSettings().setValue("studies/plot", "cd");
            sc->waitFor([session] { return session->activeStudy().isEmpty() && !session->running(); }, next, "study to finish", 300000);
        });
        sc->add(showScene(0, "run:run-0005"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); }));
        sc->add(go("runs/studies"));
        sc->add([session, sc](auto next) {
            sc->waitFor([session] {
                for (int i = 1; i <= 5; ++i)
                    if (!session->aeroIfComputed(QStringLiteral("run-%1").arg(i, 4, 10, QChar('0')))) return false;
                return true;
            }, next, "study results");
        });
        sc->add(go("runs/studies"));
    } else if (scenario == "param-run" || scenario == "cylinder-run") {
        // Standard k-epsilon (3D, SIMPLE) and laminar 2D recordings: their own residual sets.
        sc->add(open(caseArg.isEmpty() ? sample(scenario == "param-run" ? "elbow_param" : "cylinder_tutorial") : caseArg));
        sc->add(approveAll);
        sc->add(runWith(300));
        sc->add(showScene(0, "run:run-0001"));
        sc->add(go("solution/run"));
        sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::Residuals); next(); });
    } else if (scenario == "wing-models") {
        sc->add(open(caseArg.isEmpty() ? sample("naca4412_ground") : caseArg));
        sc->add([&w, sc](auto next) { sc->waitFor([&w] { return w.graphics()->panel(0)->scene() != nullptr; }, next, "mesh to show"); });
        sc->add(go("setup/models"));
    } else if (scenario == "wing-aero" || scenario == "wing-monitors") {
        // The ground-effect wing stopping when downforce and drag settle; its
        // continuity residual levels off above the case's 1e-4 criteria.
        sc->add(open(caseArg.isEmpty() ? sample("naca4412_ground") : caseArg));
        sc->add(approveAll);
        sc->add([session](auto next) {
            for (const char* p : {"monitors/check-convergence", "monitors/report/downforce/active", "monitors/report/drag/active"})
                session->settings()->setValue(p, true);
            for (const Setting& st : session->caseSummary().settings)
                if (st.path.startsWith("monitors/residual/") && st.path.endsWith("/check")) session->settings()->setValue(st.path, false);
            next();
        });
        sc->add(runWith(600));
        sc->add(showScene(0, "run:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); }));
        if (scenario == "wing-aero") {
            sc->add(go("results/aero"));
            sc->add([sc, session](auto next) {
                sc->waitFor([session] { return session->aeroIfComputed("run-0001").has_value(); }, next, "aero coefficients");
            });
        } else {
            sc->add(go("solution/monitors"));
        }
        sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::Reports); next(); });
    } else if (scenario.startsWith("new-case")) {
        const QString mesh = !caseArg.isEmpty()               ? caseArg
                             : scenario == "new-case-car" ? QDir(dataRoot).filePath("samples/fsae_car.msh")
                                                          : QDir(dataRoot).filePath("samples/naca4412_ground.msh");
        sc->add([session, sc, mesh](auto next) {
            session->startNewCase(mesh);
            sc->waitFor([session] { return !session->newCaseBusy(); }, [session, sc, next] {
                if (!session->newCaseDraft()) return sc->fail("new case: " + session->newCaseError());
                next();
            }, "mesh to read");
        });
        if (scenario == "new-case-car") {
            sc->add(go("new-case"));
        } else if (scenario == "new-case") {
            // A role and the speed changed from the template, as a user checking the suggestions would.
            sc->add([session](auto next) {
                session->setNewCaseSpeed(20);
                for (const NewCaseZone* z : session->newCaseDraft()->boundaries())
                    if (z->role == "symmetry") {
                        session->setNewCaseRole(z->name, "wall");
                        break;
                    }
                next();
            });
            sc->add(go("new-case"));
        } else {
            const QString name = qEnvironmentVariable("KUTTA_SCENARIO_CASE_NAME", QFileInfo(mesh).completeBaseName());
            const QString out = QDir(dataRoot).filePath("cases/" + name + ".kutta.json");
            QDir().mkpath(QFileInfo(out).absolutePath());
            sc->add([session, sc, out](auto next) {
                bool* done = new bool(false);
                QObject::connect(session, &Session::caseOpened, sc, [done] { *done = true; }, Qt::SingleShotConnection);
                session->createNewCase(out);
                sc->waitFor([done] { return *done; }, [done, next] { delete done; next(); }, "new case to open");
            });
            sc->add(approveAll);
            const int iterations = qEnvironmentVariableIntValue("KUTTA_SCENARIO_ITERATIONS");
            sc->add(runWith(iterations > 0 ? iterations : 300));
            sc->add(showScene(0, "run:run-0001"));
            sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); }));
            if (scenario == "new-case-save") {
                // Save with a changed reference area; the reopened case carries it.
                sc->add([session, sc, out](auto next) {
                    session->settings()->setValue("reference/area", 0.5);
                    bool* done = new bool(false);
                    QObject::connect(session, &Session::caseSaved, sc, [session, sc, done] {
                        QObject::connect(session, &Session::caseOpened, sc, [done] { *done = true; }, Qt::SingleShotConnection);
                    }, Qt::SingleShotConnection);
                    QObject::connect(session, &Session::saveFailed, sc, [sc](const QString& m) { sc->fail("save: " + m); }, Qt::SingleShotConnection);
                    session->saveCase(out);
                    sc->waitFor([done] { return *done; }, [done, next] { delete done; next(); }, "saved case to reopen");
                });
                sc->add(go("setup/reference"));
            } else {
                sc->add(go("setup/general"));
            }
            sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::Reports); next(); });
        }
    } else if (scenario == "bench") {
        // Frame timing on the car's surfaces; KUTTA_FIXTURE_CAR_DETAIL scales
        // the tessellation (5 gives about 10M triangles).
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add(runToEnd);
        sc->add(showScene(0, "run:run-0001"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("pressure_coefficient"); }));
        sc->add([&w](auto next) { w.graphics()->panel(0)->view()->fitView(); next(); });
    } else if (scenario.startsWith("case-")) {
        const QString path = caseArg.isEmpty() ? elbow : caseArg;
        sc->add(open(path, cli.value(caseOpt)));
        const int iterations = qEnvironmentVariableIsSet("KUTTA_SCENARIO_ITERATIONS") ? qEnvironmentVariableIntValue("KUTTA_SCENARIO_ITERATIONS") : 300;
        if (scenario == "case-mesh") {
            sc->add([&w, sc](auto next) { sc->waitFor([&w] { return w.graphics()->panel(0)->scene() != nullptr; }, next, "mesh to show"); });
            sc->add(go("setup/general"));
        } else if (scenario == "case-approvals") {
            sc->add([&w, sc](auto next) { sc->waitFor([&w] { return w.graphics()->panel(0)->scene() != nullptr; }, next, "mesh to show"); });
            sc->add([session](auto next) {
                // The first approximation approved, the others still waiting: both states show.
                QSettings().setValue("ui/hideGettingStarted", true);
                const QStringList a = session->caseSummary().approximations();
                if (!a.isEmpty()) session->setApproximationApproved(a.first(), true);
                next();
            });
            sc->add(go("setup/general"));
        } else if (scenario == "case-methods") {
            sc->add([&w, sc](auto next) { sc->waitFor([&w] { return w.graphics()->panel(0)->scene() != nullptr; }, next, "mesh to show"); });
            sc->add(go("solution/methods"));
        } else if (scenario == "case-run") {
            sc->add(approveAll);
            sc->add([session](auto next) {
                session->settings()->setValue("monitors/check-convergence", true);
                next();
            });
            sc->add(runWith(iterations));
            sc->add(showScene(0, "run:run-0001"));
            sc->add(go("solution/monitors"));
            sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::Residuals); next(); });
        } else {  // case-stop
            sc->add(approveAll);
            sc->add([session, sc, iterations](auto next) {
                session->settings()->setValue("monitors/check-convergence", true);
                session->settings()->setValue("run/iterations", 100000);
                if (session->run().isEmpty()) return sc->fail("run refused");
                // Stop as the user would, once the run is under way.
                sc->waitFor([session, iterations] { return session->live().last && session->live().last->iter >= iterations; }, [session, sc, next] {
                    session->stop();
                    sc->waitFor([session] { return !session->running(); }, next, "run to stop", 120000);
                }, "run to reach the stop iteration", 300000);
            });
            sc->add(showScene(0, "run:run-0001"));
            sc->add(go("solution/run"));
            sc->add([&w](auto next) { w.charts()->showTab(ChartsPanel::Tab::Residuals); next(); });
        }
        // Another page of the same state, e.g. KUTTA_SCENARIO_PAGE=runs/history.
        if (const QString page = qEnvironmentVariable("KUTTA_SCENARIO_PAGE"); !page.isEmpty()) sc->add(go(page));
    } else if (scenario == "settings") {
        sc->add(open(elbowCase));
        sc->add([&w](auto next) {
            // The page is built from the whole inventory (9,844 settings on the elbow).
            QElapsedTimer t;
            t.start();
            w.inspector()->showSettingsFiltered("cold-inlet");
            std::fprintf(stderr, "settings page: %lld ms for %zu settings\n", t.elapsed(), w.session()->settings()->settings().size());
            next();
        });
    }

    const QString shot = cli.value(shotOpt);
    const int benchFrames = cli.value(benchOpt).toInt();
    sc->add(settle(600));
    sc->add([&](auto next) {
        if (benchFrames <= 0) return next();
        auto* p = w.graphics()->active();
        p->view()->measureFrames(benchFrames, [&, next, p](SceneView::FrameTimes t) {
            // Idle cost: no input for 5 s. Nothing should render or use CPU.
            const int framesBefore = p->view()->framesRendered();
            const double cpu0 = processCpuSeconds();
            QElapsedTimer wall;
            wall.start();
            QTimer::singleShot(5000, [&, t, framesBefore, cpu0, wall, next, p] {
                const double cpu = processCpuSeconds() - cpu0;
                QJsonObject o{{"frames", t.frames}, {"mean_ms", t.meanMs}, {"p95_ms", t.p95Ms}, {"max_ms", t.maxMs},
                              {"triangles", double(t.triangles)}, {"fps_mean", t.meanMs > 0 ? 1000.0 / t.meanMs : 0},
                              {"idle_seconds", wall.elapsed() / 1000.0}, {"idle_cpu_percent", cpu / (wall.elapsed() / 1000.0) * 100.0},
                              {"idle_frames_rendered", p->view()->framesRendered() - framesBefore},
                              {"viewport_px", QStringLiteral("%1x%2@%3").arg(p->view()->width()).arg(p->view()->height()).arg(p->view()->devicePixelRatio())}};
                std::printf("%s\n", QJsonDocument(o).toJson(QJsonDocument::Compact).constData());
                std::fflush(stdout);
                next();
            });
        });
    });
    // What a reviewer would catch in the shot: the inspector page laid out
    // badly at any panel width (text overlapping, clipped or cut short), or
    // the home folder anywhere on screen (core/PathDisplay.h).
    sc->add([&](auto next) {
        if (benchFrames > 0) return next();
        QList<int> widths;
        for (int px = 320; px <= 720; px += 20) widths << px;
        QStringList problems = w.inspector()->layoutProblemsAtWidths(widths);
        problems += homePathsOnScreen(&w);
        if (problems.isEmpty()) return next();
        for (const QString& p : problems) std::fprintf(stderr, "screen check: %s\n", qPrintable(p));
        std::fprintf(stderr, "scenario failed: %lld screen check problem(s)\n", qlonglong(problems.size()));
        QCoreApplication::exit(4);
    });
    sc->add(settle(100));
    // KUTTA_SCENARIO_SCROLL=bottom shows the end of a long inspector page.
    sc->add([&](auto next) {
        if (qEnvironmentVariable("KUTTA_SCENARIO_SCROLL") == "bottom") w.inspector()->verticalScrollBar()->setValue(w.inspector()->verticalScrollBar()->maximum());
        next();
    });
    sc->add([&](auto) {
        if (!shot.isEmpty()) {
            const QImage img = w.grab().toImage();
            if (!img.save(shot)) {
                std::fprintf(stderr, "could not write %s\n", qPrintable(shot));
                QCoreApplication::exit(1);
                return;
            }
            std::printf("wrote %s (%dx%d)\n", qPrintable(shot), img.width(), img.height());
            // A GPU viewport showing a dataset must have drawn something: the
            // smoke test relies on this to catch a blank or failed backend.
            for (int i = 0; i < w.graphics()->visibleCount(); ++i) {
                ViewportPanel* p = w.graphics()->panel(i);
                if (!p->scene()) continue;
                const QImage v = p->exportImage();
                QSet<QRgb> colours;
                for (int y = 0; y < v.height(); y += std::max(1, v.height() / 64))
                    for (int x = 0; x < v.width(); x += std::max(1, v.width() / 64)) colours.insert(v.pixel(x, y));
                std::printf("viewport %d: %lld sampled colours\n", i, qlonglong(colours.size()));
                if (colours.size() < 16) {
                    std::fprintf(stderr, "viewport %d rendered blank\n", i);
                    QCoreApplication::exit(3);
                    return;
                }
                // Every scenario colours results by a field; grey results mean
                // the display lost it (once: the mesh's display carried over).
                if (!p->scene()->geometryOnly && !p->view()->fieldShown()) {
                    std::fprintf(stderr, "viewport %d shows results without their field\n", i);
                    QCoreApplication::exit(3);
                    return;
                }
            }
        }
        session->stop();
        QTimer::singleShot(100, [] { QCoreApplication::exit(0); });
    });
    // Scripted windows appear without taking focus: on a machine in use, keys
    // typed elsewhere would otherwise land in the window (the command palette
    // opens on typing) and end up in the screenshot.
    w.setAttribute(Qt::WA_ShowWithoutActivating);
    w.show();
    sc->start();
    return app.exec();
}
