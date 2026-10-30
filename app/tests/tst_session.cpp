// The session against the fixture CLI, through the real process boundary:
// open, refuse without approval, run, record, stop, and load results.

#include "core/Session.h"

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

    std::unique_ptr<Session> open()
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        auto s = std::make_unique<Session>(cli, m_dir.filePath("data"));
        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->openCase(m_elbow);
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
    }
};

QTEST_GUILESS_MAIN(TestSession)
#include "tst_session.moc"
