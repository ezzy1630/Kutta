// The session against the fixture CLI, through the real process boundary:
// open, refuse without approval, run, record, stop, and load results.

#include "core/Session.h"

#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

class TestSession : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;
    QString m_elbow;

    std::unique_ptr<Session> open(QString casePath = {})
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        auto s = std::make_unique<Session>(cli, m_dir.filePath("data"));
        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->openCase(casePath.isEmpty() ? m_elbow : casePath);
        if (!opened.wait(20000)) return nullptr;
        return s;
    }

private slots:
    void initTestCase()
    {
        QCoreApplication::setOrganizationName("KuttaTest");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_dir.filePath("settings"));
        qputenv("KUTTA_FIXTURE_SPEED", "400");
        // The fixture CLI serves cases by file name.
        m_elbow = m_dir.filePath("mixing_elbow.cas.h5");
        QFile f(m_elbow);
        QVERIFY(f.open(QIODevice::WriteOnly));
    }

    // The case's mesh is a dataset as soon as the case is read, with no run.
    void meshIsShownWithoutARun()
    {
        auto s = open();
        QVERIFY(s);
        const QString key = s->meshSceneKey();
        QVERIFY(key.startsWith("mesh:"));
        bool done = false;
        std::shared_ptr<SceneData> scene;
        s->loadScene(key, [&](const QString&, Result<std::shared_ptr<SceneData>> r) {
            QVERIFY2(r.ok(), qPrintable(r ? QString() : r.error()));
            scene = r.value();
            done = true;
        });
        QTRY_VERIFY_WITH_TIMEOUT(done, 20000);
        QVERIFY(scene && scene->geometryOnly);
        QStringList zones;
        for (const auto& it : scene->items) {
            QCOMPARE(it.kind, ItemKind::Zone);
            QVERIFY(it.color.isValid());
            zones << it.label;
        }
        QVERIFY(zones.contains("wall-elbow") && zones.contains("cold-inlet"));
        // The mesh view leaves out the symmetry plane, which would hide the pipe.
        const QStringList shown = defaultMeshItems(*scene);
        QVERIFY(shown.contains("zone:wall-elbow"));
        QVERIFY(!shown.contains("zone:symmetry-xyplane"));
        // Reopening reuses the export instead of reading the mesh again.
        QVERIFY(QFileInfo::exists(m_dir.filePath("data/geometry/" + key.mid(5) + "/geometry/manifest.json")));
    }

    // A case that cannot be read leaves the open case open and says why,
    // for the page to show in place (no modal dialog).
    void failedOpenKeepsTheCaseAndTheReason()
    {
        auto s = open();
        QVERIFY(s);
        QSignalSpy failed(s.get(), &Session::caseFailed);
        s->openCase(m_dir.filePath("moved.cas.h5"));  // refused before kutta runs
        QCOMPARE(failed.count(), 1);
        QVERIFY(s->openError().message.contains("does not exist"));
        QFile junk(m_dir.filePath("not-a-case.txt"));
        QVERIFY(junk.open(QIODevice::WriteOnly) && junk.write("hello") > 0);
        junk.close();
        s->openCase(junk.fileName());
        QVERIFY(failed.wait(20000));
        QVERIFY(s->hasCase());
        QCOMPARE(s->caseName(), QString("mixing_elbow"));
        QCOMPARE(QFileInfo(s->openError().casePath).fileName(), QString("not-a-case.txt"));
        QVERIFY(!s->openError().message.isEmpty());
        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->openCase(m_elbow);
        QVERIFY(opened.wait(20000));
        QVERIFY(s->openError().message.isEmpty());
    }

    void refusesApproximateRunWithoutApproval()
    {
        auto s = open();
        QVERIFY(s);
        QString why;
        QVERIFY(!s->canRun(&why));
        QVERIFY(why.contains("Warped-face"));
        QSignalSpy refused(s.get(), &Session::runRefused);
        QVERIFY(s->run().isEmpty());
        QCOMPARE(refused.count(), 1);
    }

    void runIsRecordedWithItsSettings()
    {
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        QCOMPARE(s->settings()->setValue("run/iterations", 200), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QCOMPARE(id, QString("run-0001"));
        QVERIFY(finished.wait(60000));
        const RunRecord* r = s->history().find(id);
        QVERIFY(r);
        QCOMPARE(r->status, RunStatus::Completed);
        QCOMPARE(r->iterationsDone, 200);
        QCOMPARE(r->approximations, QStringList{"wfgc"});
        QCOMPARE(r->iterationsRequested, 200);
        QVERIFY(!r->modified());  // the iteration count is a run control, not a case change
        QVERIFY(r->fixture);
        QCOMPARE(r->coupling, QString("coupled"));  // what the setup event says ran
        QCOMPARE(r->stopReason, QString("iteration_limit"));
        // The events the app received are kept with the run.
        const RunProgress p = s->history().progress(id);
        QVERIFY(p.done);
        QCOMPARE(p.residuals.iterations.back(), 200);
        QVERIFY(s->history().forces(id));

        // Its results load as a dataset.
        bool loaded = false;
        s->loadScene("run:" + id, [&](const QString&, Result<std::shared_ptr<SceneData>> scene) {
            QVERIFY2(scene.ok(), qPrintable(scene ? QString() : scene.error()));
            QVERIFY(scene.value()->item("zone:symmetry-xyplane"));
            QVERIFY(scene.value()->provenance.contains("APPROXIMATE"));
            loaded = true;
        });
        QTRY_VERIFY_WITH_TIMEOUT(loaded, 20000);
    }

    // Each operation joins the run's batch and the whole view is republished;
    // a rejected operation leaves the batch and the published view as they were.
    void postProcessingRepublishesTheRunView()
    {
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        QCOMPARE(s->settings()->setValue("run/iterations", 40), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(finished.wait(60000));

        std::optional<Result<PostOutput>> section;
        s->postProcess(id, PostKind::Probe,
                       {{"name", "Section"}, {"type", "section"}, {"zones", QJsonArray{"wall-elbow"}},
                        {"origin", QJsonArray{0, 0, 0.025}}, {"normal", QJsonArray{0, 0, 1}}, {"fields", QJsonArray{"pressure"}}},
                       [&](Result<PostOutput> r) { section = r; });
        QTRY_VERIFY_WITH_TIMEOUT(section.has_value(), 20000);
        QVERIFY2(section->ok(), qPrintable(*section ? QString() : section->error()));
        QVERIFY(section->value().part.lineCount >= 1);
        QVERIFY(section->value().scene->item("post:Section"));

        std::optional<Result<PostOutput>> report;
        s->postProcess(id, PostKind::Integral, {{"name", "Ports"}, {"zones", QJsonArray{"cold-inlet", "outlet"}}, {"fields", QJsonArray{"pressure"}}},
                       [&](Result<PostOutput> r) { report = r; });
        QTRY_VERIFY_WITH_TIMEOUT(report.has_value(), 20000);
        QVERIFY2(report->ok(), qPrintable(*report ? QString() : report->error()));
        QCOMPARE(report->value().report.value("zones").toArray().size(), 2);
        QVERIFY(report->value().scene->item("post:Section"));  // earlier results persist

        // The recorded elbow has no volume fields: tracing fails and changes nothing.
        std::optional<Result<PostOutput>> traces;
        s->postProcess(id, PostKind::Streamlines,
                       {{"name", "Pathlines"}, {"seeds", QJsonObject{{"type", "rake"}, {"start", QJsonArray{0, 0, 0}}, {"end", QJsonArray{0, 0.1, 0}}, {"count", 4}}}},
                       [&](Result<PostOutput> r) { traces = r; });
        QTRY_VERIFY_WITH_TIMEOUT(traces.has_value(), 20000);
        QVERIFY(!traces->ok());
        QFile batch(QDir(s->history().find(id)->dir).filePath("app-post.json"));
        QVERIFY(batch.open(QIODevice::ReadOnly));
        const QJsonObject b = QJsonDocument::fromJson(batch.readAll()).object();
        QVERIFY(!b.contains("streamlines"));
        QCOMPARE(b.value("probes").toArray().size(), 1);
        auto view = loadViewManifest(QDir(s->history().find(id)->dir).filePath("view/manifest.json"));
        QVERIFY(view && view.value().probes.size() == 1);
    }

    void liveSnapshotsFollowThisRun()
    {
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        s->setLiveViewInterval(20);
        QCOMPARE(s->settings()->setValue("run/iterations", 150), QString());
        QSignalSpy live(s.get(), &Session::liveSceneChanged);
        QSignalSpy finished(s.get(), &Session::runFinished);
        s->run();
        QVERIFY(finished.wait(60000));
        QVERIFY(live.count() >= 2);  // recorded snapshots at 40 and 120
        QCOMPARE(s->live().lastSnapshotSequence, live.count());
        QVERIFY(!s->live().setup->liveRunId.isEmpty());
        QVERIFY(s->sceneIfLoaded("live"));
    }

    void stopEndsTheRunCleanly()
    {
        qputenv("KUTTA_FIXTURE_HOLD_AT", "60");
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QTRY_VERIFY_WITH_TIMEOUT(s->live().last && s->live().last->iter >= 60, 20000);
        s->stop();
        QVERIFY(finished.wait(30000));
        qunsetenv("KUTTA_FIXTURE_HOLD_AT");
        const RunRecord* r = s->history().find(id);
        QCOMPARE(r->status, RunStatus::Stopped);
        QCOMPARE(r->stopReason, QString("stopped"));
        QVERIFY(r->iterationsDone >= 60 && r->iterationsDone < 600);
        QVERIFY(runOutcomeAdvice(*r).contains("Stopped by you"));
    }

    // With convergence checking on, the run stops on the case's own residual
    // criteria and says so; the last sample's condition states are kept.
    void convergedRunSaysWhichConditionsWereMet()
    {
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        QCOMPARE(s->settings()->setValue("monitors/check-convergence", true), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(finished.wait(60000));
        const RunRecord* r = s->history().find(id);
        QVERIFY(r->converged());
        QCOMPARE(r->convergedBy, QStringList{"residuals"});
        QVERIFY(r->iterationsDone < 1000);
        QVERIFY(!r->modified());  // convergence controls are run controls
        const RunProgress p = s->history().progress(id);
        QVERIFY(p.convergence && p.convergence->met);
        QCOMPARE(p.convergenceIter, r->iterationsDone);
        const auto marks = convergenceMarks(p, false);
        QVERIFY(std::any_of(marks.begin(), marks.end(), [](const ConvergenceMark& m) { return m.kind == ConvergenceMark::Kind::Converged; }));
    }

    // Recorded k-epsilon, laminar 2D and transition cases replay with their
    // own equations, approximations and early iterations.
    void recordedCasesKeepTheirOwnEquations_data()
    {
        QTest::addColumn<QString>("file");
        QTest::addColumn<QStringList>("residuals");
        QTest::addColumn<QStringList>("labels");
        QTest::addColumn<bool>("planar");
        QTest::addColumn<QString>("transition");
        QTest::newRow("k-epsilon") << "elbow_param.cas.h5" << QStringList{"continuity", "x-velocity", "y-velocity", "z-velocity", "k", "epsilon"}
                                   << QStringList{"No compressibility correction", "Documented k-epsilon only"} << false << QString();
        QTest::newRow("laminar 2D") << "cylinder_tutorial.cas" << QStringList{"continuity", "x-velocity", "y-velocity"}
                                    << QStringList{"No reverse-flow suppression"} << true << QString();
        QTest::newRow("transition") << "naca4412_ground.cas" << QStringList{"continuity", "x-velocity", "y-velocity", "k", "omega", "intermit"}
                                    << QStringList{} << true << QString("gamma-transport-eqn");
    }
    void recordedCasesKeepTheirOwnEquations()
    {
        QFETCH(QString, file);
        QFETCH(QStringList, residuals);
        QFETCH(QStringList, labels);
        QFETCH(bool, planar);
        QFETCH(QString, transition);
        const QString path = m_dir.filePath(file);
        QFile placeholder(path);
        QVERIFY(placeholder.open(QIODevice::WriteOnly));
        placeholder.close();
        auto s = open(path);
        QVERIFY(s);
        QCOMPARE(s->caseSummary().planar(), planar);
        QCOMPARE(s->caseSummary().transition, transition);
        for (const QString& a : s->caseSummary().approximations()) s->setApproximationApproved(a, true);
        QCOMPARE(s->settings()->setValue("run/iterations", 300), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(!id.isEmpty());
        QVERIFY(finished.wait(60000));
        const RunRecord* r = s->history().find(id);
        QCOMPARE(r->status, RunStatus::Completed);
        QCOMPARE(runStatusText(*r), QString("Iteration limit"));
        QCOMPARE(r->approximationLabels, labels);
        const RunProgress p = s->history().progress(id);
        QCOMPARE(QSet<QString>(p.residuals.names.begin(), p.residuals.names.end()), QSet<QString>(residuals.begin(), residuals.end()));
        // Iterations 1-5 come before the first monitor sample.
        QCOMPARE(std::vector<int>(p.residuals.iterations.begin(), p.residuals.iterations.begin() + 6), (std::vector<int>{1, 2, 3, 4, 5, 20}));
        if (transition.isEmpty()) return;
        // The wing's downforce report is a lift with a downward vector; up stays +y in the plane.
        QCOMPARE(s->aeroSpec().lift, QVector3D(0, 1, 0));
        std::optional<AeroResult> aero;
        s->aero(id, [&](Result<AeroResult> a) {
            QVERIFY2(a.ok(), qPrintable(a ? QString() : a.error()));
            aero = a.value();
        });
        QTRY_VERIFY_WITH_TIMEOUT(aero.has_value(), 20000);
        QVERIFY(aero->total()->downforce > 30);
        QVERIFY(aero->total()->total.drag > 0);
        // The replay stopped at 300 of the recording's 600: aero says 300.
        QCOMPARE(aero->iterations, r->iterationsDone);
        // A zone named wing is suggested into a wing group, not body.
        QCOMPARE(int(s->aeroSpec().groups.size()), 1);
        QCOMPARE(s->aeroSpec().groups[0].name, QString("Wing"));
    }

    // A run that fails says why and what to do: disk full and a crash.
    void failedRunSaysWhatToDo()
    {
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        qputenv("KUTTA_FIXTURE_FAIL_AT", "40");
        qputenv("KUTTA_FIXTURE_FAIL_MESSAGE", "cannot write fields.bin: No space left on device");
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString full = s->run();
        QVERIFY(finished.wait(30000));
        qunsetenv("KUTTA_FIXTURE_FAIL_MESSAGE");
        const QString crash = s->run();
        QVERIFY(finished.wait(30000));
        qunsetenv("KUTTA_FIXTURE_FAIL_AT");

        const RunRecord* r = s->history().find(full);
        QCOMPARE(r->status, RunStatus::Failed);
        QVERIFY(r->error.contains("No space left on device"));
        QVERIFY(runOutcomeAdvice(*r).contains("disk is full"));
        r = s->history().find(crash);
        QCOMPARE(r->status, RunStatus::Failed);
        QVERIFY(runOutcomeAdvice(*r).contains("stopped before finishing"));
    }

    // A run Kutta cannot record would leave no trace (or run without the
    // user's changes): it is refused before kutta starts.
    void unwritableRunFolderRefusesTheRun()
    {
        auto s = open();
        QVERIFY(s);
        s->setApproximationApproved("wfgc", true);
        QCOMPARE(s->settings()->setValue("boundary/cold-inlet/velocity-magnitude", 1.2), QString());
        // A file where the case's run folder belongs blocks the run on every
        // OS; a read-only directory does not on Windows.
        const QString caseDir = s->history().caseDir();
        QVERIFY(QDir(caseDir).removeRecursively());
        QFile blocker(caseDir);
        QVERIFY(blocker.open(QIODevice::WriteOnly));
        blocker.close();
        QSignalSpy refused(s.get(), &Session::runRefused);
        const QString id = s->run();
        QVERIFY(blocker.remove());
        QVERIFY(id.isEmpty());
        QCOMPARE(refused.count(), 1);
        QVERIFY(refused.first().first().toString().contains("Could not start the run"));
        QVERIFY(!s->running());
    }

    // File ▸ New Case from Mesh: suggestions, a changed role and speed, the
    // created Kutta case opened and run, then saved elsewhere and reopened.
    void newCaseFromMeshCreatesRunsAndSaves()
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        Session s(cli, m_dir.filePath("data-new"));
        QDir().mkpath(m_dir.filePath("mesh"));
        const QString mesh = m_dir.filePath("mesh/naca4412_ground.msh");
        QFile f(mesh);
        QVERIFY(f.open(QIODevice::WriteOnly));
        f.close();

        s.startNewCase(mesh);
        QVERIFY(s.newCaseBusy());
        QTRY_VERIFY_WITH_TIMEOUT(!s.newCaseBusy(), 20000);
        QVERIFY2(s.newCaseDraft(), qPrintable(s.newCaseError()));
        const NewCaseDraft& d = *s.newCaseDraft();
        QCOMPARE(d.dimension, 2);
        QMap<QString, QString> suggested;
        for (const NewCaseZone* z : d.boundaries()) suggested[z->name] = z->suggestedRole;
        QCOMPARE(suggested, (QMap<QString, QString>{{"wing", "wall"}, {"inlet", "velocity-inlet"}, {"outlet", "pressure-outlet"},
                                                    {"top", "symmetry"}, {"ground", "moving-ground"}}));
        s.setNewCaseRole("top", "wall");
        s.setNewCaseRole("top", "no-such-role");  // ignored
        s.setNewCaseSpeed(20);
        for (const NewCaseZone* z : s.newCaseDraft()->boundaries())
            if (z->name == "top") QCOMPARE(z->role, QString("wall"));
        QCOMPARE(s.newCaseDraft()->speed, 20.0);

        const QString out = m_dir.filePath("cases/wing.kutta.json");
        QDir().mkpath(m_dir.filePath("cases"));
        QSignalSpy opened(&s, &Session::caseOpened);
        s.createNewCase(out);
        QVERIFY(opened.wait(20000));
        QVERIFY(!s.newCaseDraft());
        QVERIFY(s.caseSummary().kuttaCase());
        QCOMPARE(s.caseName(), QString("wing"));
        QCOMPARE(QFileInfo(s.caseSummary().meshFile).canonicalFilePath(), QFileInfo(mesh).canonicalFilePath());
        QCOMPARE(s.caseSummary().zone("top")->type, QString("wall"));
        const QJsonObject written = QJsonDocument::fromJson([&] { QFile c(out); c.open(QIODevice::ReadOnly); return c.readAll(); }()).object();
        QCOMPARE(written.value("mesh").toString(), QString("../mesh/naca4412_ground.msh"));
        QCOMPARE(written.value("boundaries").toObject().value("inlet").toObject().value("speed").toDouble(), 20.0);
        QCOMPARE(written.value("reports").toObject().value("drag").toObject().value("zones").toArray(), (QJsonArray{"wing", "top"}));

        QCOMPARE(s.settings()->setValue("run/iterations", 20), QString());
        QSignalSpy finished(&s, &Session::runFinished);
        const QString id = s.run();
        QVERIFY(!id.isEmpty());
        QVERIFY(finished.wait(60000));
        QCOMPARE(s.history().find(id)->status, RunStatus::Completed);

        // Save As with no changes: the copy opens from its new place.
        const QString copy = m_dir.filePath("cases/saved/wing-copy.kutta.json");
        QDir().mkpath(QFileInfo(copy).absolutePath());
        QSignalSpy saved(&s, &Session::caseSaved);
        s.saveCase(copy);
        QVERIFY(opened.wait(20000));
        QCOMPARE(saved.count(), 1);
        QCOMPARE(s.casePath(), copy);
        QCOMPARE(QFileInfo(s.caseSummary().meshFile).canonicalFilePath(), QFileInfo(mesh).canonicalFilePath());

        // The starter's values, not the recording's: Fluent's default area, the chosen speed.
        QCOMPARE(s.settings()->value("reference/area").toDouble(), 1.0);
        QCOMPARE(s.settings()->value("reference/velocity").toDouble(), 20.0);
        QCOMPARE(s.settings()->value("boundary/inlet/velocity-magnitude").toDouble(), 20.0);

        // Save with changed settings: written into the case, reopened with them, no longer changes.
        QCOMPARE(s.settings()->setValue("reference/area", 0.5), QString());
        QCOMPARE(s.settings()->setValue("boundary/inlet/velocity-magnitude", 25.0), QString());
        s.saveCase(copy);
        QVERIFY(opened.wait(20000));
        QCOMPARE(saved.count(), 2);
        QCOMPARE(s.settings()->caseOverrides(), QJsonObject());
        QCOMPARE(s.settings()->value("reference/area").toDouble(), 0.5);
        QCOMPARE(s.settings()->value("boundary/inlet/velocity-magnitude").toDouble(), 25.0);
        const QJsonObject resaved = QJsonDocument::fromJson([&] { QFile c(copy); c.open(QIODevice::ReadOnly); return c.readAll(); }()).object();
        QCOMPARE(resaved.value("reference_values").toObject().value("area").toDouble(), 0.5);
        QCOMPARE(resaved.value("boundaries").toObject().value("inlet").toObject().value("speed").toDouble(), 25.0);

        // A change the fixture cannot write: the session says which, and keeps the case and the change.
        QCOMPARE(s.settings()->setValue("boundary/outlet/gauge-pressure", 10.0), QString());
        QSignalSpy failed(&s, &Session::saveFailed);
        s.saveCase(copy);
        QVERIFY(failed.wait(20000));
        QVERIFY(failed.first().first().toString().contains("real kutta CLI"));
        QVERIFY(failed.first().first().toString().contains("boundary/outlet/gauge-pressure"));
        QCOMPARE(s.casePath(), copy);
        QVERIFY(s.settings()->caseOverrides().contains("boundary/outlet/gauge-pressure"));
    }

    // A 3D starter with wheels (the fixture's synthetic car mesh): wheels are
    // grouped, stay walls, and the CLI's note says they will not turn until the
    // case gives them rotating walls.
    void newCaseNotesStationaryWheels()
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        Session s(cli, m_dir.filePath("data-car"));
        const QString mesh = m_dir.filePath("mesh/fsae_car.msh");
        QDir().mkpath(m_dir.filePath("mesh"));
        QFile f(mesh);
        QVERIFY(f.open(QIODevice::WriteOnly));
        f.close();
        s.startNewCase(mesh);
        QTRY_VERIFY_WITH_TIMEOUT(!s.newCaseBusy(), 20000);
        QVERIFY2(s.newCaseDraft(), qPrintable(s.newCaseError()));
        QCOMPARE(s.newCaseDraft()->dimension, 3);
        auto wheelNote = [&s] {
            for (const QString& n : s.newCaseDraft()->notes)
                if (n.startsWith("Wheels are stationary walls")) return true;
            return false;
        };
        QStringList wheels;
        for (const NewCaseZone* z : s.newCaseDraft()->boundaries())
            if (z->group == "Wheels") wheels << z->name + "=" + z->role;
        QCOMPARE(wheels, (QStringList{"tyre-front=wall", "tyre-rear=wall"}));
        QVERIFY(wheelNote());

        // The starter opens as a Kutta case on the synthetic car.
        QSignalSpy opened(&s, &Session::caseOpened);
        s.createNewCase(m_dir.filePath("cases/car.kutta.json"));
        QVERIFY(opened.wait(20000));
        QVERIFY(s.caseSummary().kuttaCase());
        QVERIFY(s.caseSummary().synthetic);
        QCOMPARE(s.caseSummary().zone("tyre-rear")->type, QString("wall"));
    }

    void missingMeshSaysSo()
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        Session s(cli, m_dir.filePath("data-missing"));
        s.startNewCase(m_dir.filePath("nowhere.msh"));
        QVERIFY(!s.newCaseDraft());
        QVERIFY(s.newCaseError().contains("does not exist"));
    }
};

QTEST_GUILESS_MAIN(TestSession)
#include "tst_session.moc"
