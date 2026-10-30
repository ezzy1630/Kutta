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
#include "ui/Theme.h"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QElapsedTimer>
#include <QFile>
#include <QJsonDocument>
#include <QJsonArray>
#include <QJsonObject>
#include <QSettings>
#include <QStandardPaths>
#include <QStyleHints>
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
    for (const char* name : {"mixing_elbow.cas.h5", "mixing_elbow.dat.h5", "fsae_car.cas.h5", "fsae_car.dat.h5"}) {
        const QString p = QDir(root).filePath(name);
        if (!QFileInfo::exists(p)) {
            QFile f(p);
            if (f.open(QIODevice::WriteOnly)) f.write("Kutta fixture sample placeholder: the fixture CLI serves this case by name.\n");
        }
        if (QString(name).contains(".cas.")) out << p;
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

    const QStringList scenarios = {"start", "elbow-setup", "palette", "elbow-run", "elbow-live", "elbow-fluent", "car-results",
                                   "car-xy", "car-compare", "history", "settings"};
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
    const CliInfo info = locateCli(cli.value(cliOpt), devFixture);
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
    auto post = [&w, session, sc](QJsonObject item) {
        return [&w, session, sc, item](std::function<void()> next) {
            const QString run = session->history().latestFinished() ? session->history().latestFinished()->id : QString();
            session->postProcess(run, item, [sc, next](Result<QJsonObject> r) {
                if (!r) return sc->fail("post-processing: " + r.error());
                next();
            });
        };
    };

    const QString elbowCase = caseArg.isEmpty() ? elbow : caseArg;
    if (scenario == "start") {
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
        QString data = cli.value(caseOpt);
        if (withData && data.isEmpty()) data = QString(elbowCase).replace(".cas.h5", ".dat.h5");
        sc->add(open(elbowCase, withData ? data : QString()));
        sc->add(approveAll);
        sc->add(runToEnd);
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
    } else if (scenario == "elbow-live") {
        qputenv("KUTTA_FIXTURE_HOLD_AT", "150");
        qputenv("KUTTA_FIXTURE_SPEED", "20");
        sc->add(open(elbowCase));
        sc->add(approveAll);
        sc->add([session](auto next) { session->setLiveViewInterval(20); next(); });
        sc->add([session, sc](auto next) {
            session->run();
            sc->waitFor([session] { return session->live().lastSnapshotIter >= 120 && session->live().last && session->live().last->iter >= 150; }, next, "live snapshot");
        });
        sc->add(showScene(0, "live"));
        sc->add(display(0, [&](ViewDisplay& d) { d.field = field("velocity_magnitude"); d.visible = {"zone:symmetry-xyplane"}; d.visibleInitialized = true; }));
        sc->add(preset(0, Camera::Preset::Top));
        sc->add(go("solution/run"));
    } else if (scenario == "car-results" || scenario == "car-xy") {
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add(runToEnd);
        sc->add(showScene(0, "run:run-0001"));
        if (scenario == "car-results") {
            sc->add(post({{"id", "pathlines-a"}, {"kind", "streamlines"}, {"name", "Pathlines"},
                          {"seed", QJsonObject{{"from", QJsonArray{-0.6, -0.7, 0.06}}, {"to", QJsonArray{-0.6, 0.7, 0.5}}, {"count", 48}}}}));
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
            sc->add(post({{"id", "line-a"}, {"kind", "line-probe"}, {"samples", 300}, {"from", QJsonArray{-1, 0, 0.1}}, {"to", QJsonArray{4, 0, 0.1}},
                          {"fields", QJsonArray{"velocity_magnitude", "total_pressure_coefficient"}}}));
            sc->add([&w, session](auto next) {
                // Plot the probe through the same path the XY page uses.
                QFile f(QDir(session->history().latestFinished()->dir).filePath("post/line-a/manifest.json"));
                f.open(QIODevice::ReadOnly);
                const QJsonObject m = QJsonDocument::fromJson(f.readAll()).object();
                std::vector<PlotSeries> series;
                int k = 0;
                for (const QJsonValue& sv : m.value("series").toArray()) {
                    PlotSeries ps;
                    ps.name = displayLabel(FieldRef::parse(sv.toObject().value("name").toString()));
                    for (const QJsonValue& x : sv.toObject().value("x").toArray()) ps.x.push_back(x.toDouble());
                    for (const QJsonValue& y : sv.toObject().value("y").toArray()) ps.y.push_back(y.isDouble() ? y.toDouble() : NAN);
                    ps.color = tokens().series[std::size_t(k++)];
                    series.push_back(std::move(ps));
                }
                w.charts()->setXY("Along y = 0, z = 0.1 m · Run 1", "Distance [m]", "Value", std::move(series), {"SYNTHETIC"},
                                  m.value("provenance").toString());
                next();
            });
        }
    } else if (scenario == "car-compare" || scenario == "history") {
        sc->add(open(caseArg.isEmpty() ? car : caseArg));
        sc->add(approveAll);
        sc->add(runToEnd);
        sc->add([session](auto next) {
            session->settings()->setValue("boundary/inlet/velocity-magnitude", 20.0);
            next();
        });
        sc->add(runToEnd);
        if (scenario == "car-compare") {
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
        } else {
            sc->add(showScene(0, "run:run-0002"));
            sc->add(go("runs/history"));
        }
    } else if (scenario == "settings") {
        sc->add(open(elbowCase));
        sc->add([&w](auto next) { w.inspector()->showSettingsFiltered("cold-inlet turb"); next(); });
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
    sc->add([&](auto) {
        if (!shot.isEmpty()) {
            const QImage img = w.grab().toImage();
            if (!img.save(shot)) {
                std::fprintf(stderr, "could not write %s\n", qPrintable(shot));
                QCoreApplication::exit(1);
                return;
            }
            std::printf("wrote %s (%dx%d)\n", qPrintable(shot), img.width(), img.height());
        }
        session->stop();
        QTimer::singleShot(100, [] { QCoreApplication::exit(0); });
    });
    w.show();
    sc->start();
    return app.exec();
}
