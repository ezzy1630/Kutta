// The session against the fixture CLI, through the real process boundary:
// open, refuse without approval, run, record, stop, and load results.

#include "core/Session.h"

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
        QCOMPARE(r->overrides.value("run/iterations").toInt(), 200);
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
