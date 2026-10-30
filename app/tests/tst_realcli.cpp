// The session against the real kutta CLI on the Ansys mixing elbow. Skipped
// unless KUTTA_REAL_CLI (the binary) and KUTTA_ELBOW_CASE (mixing_elbow.cas.h5,
// with its .dat.h5 beside it) are set: CI has no solver build or case files.

#include "core/Session.h"

#include <QFileInfo>
#include <QJsonArray>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

class TestRealCli : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;
    QString m_cli, m_case, m_data;

    std::unique_ptr<Session> open(bool withData)
    {
        CliInfo cli;
        cli.program = m_cli;
        auto s = std::make_unique<Session>(cli, m_dir.filePath("data"));
        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->openCase(m_case, withData ? m_data : QString());
        if (!opened.wait(60000)) return nullptr;
        for (const QString& a : s->caseSummary().approximations()) s->setApproximationApproved(a, true);
        return s;
    }

    static std::optional<Result<PostOutput>> post(Session& s, const QString& run, PostKind kind, const QJsonObject& op)
    {
        std::optional<Result<PostOutput>> out;
        s.postProcess(run, kind, op, [&](Result<PostOutput> r) { out = r; });
        QTest::qWaitFor([&] { return out.has_value(); }, 120000);
        return out ? out : std::optional<Result<PostOutput>>(Result<PostOutput>::failure("timed out"));
    }

private slots:
    void initTestCase()
    {
        m_cli = qEnvironmentVariable("KUTTA_REAL_CLI");
        m_case = qEnvironmentVariable("KUTTA_ELBOW_CASE");
        m_data = QString(m_case).replace(".cas.h5", ".dat.h5");
        if (m_cli.isEmpty() || m_case.isEmpty()) QSKIP("set KUTTA_REAL_CLI and KUTTA_ELBOW_CASE to run against the real CLI");
        QVERIFY(QFileInfo(m_cli).isExecutable());
        QVERIFY(QFileInfo::exists(m_data));
        QCoreApplication::setOrganizationName("KuttaRealCliTest");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_dir.filePath("settings"));
    }

    // The whole Fluent rp table arrives; only the adapter's controls are editable.
    void inventoryIsComplete()
    {
        auto s = open(true);
        QVERIFY(s);
        const auto& all = s->settings()->settings();
        QVERIFY(all.size() > 5000);
        int editable = 0;
        for (const auto& st : all) editable += st.editable;
        QVERIFY(editable >= 20 && editable < 100);
        QVERIFY(s->settings()->find("boundary/cold-inlet/velocity-magnitude")->editable);
    }

    // Run with live surfaces, then post-process the finished run with Fluent's
    // fields alongside Kutta's.
    void runLiveThenPostProcess()
    {
        auto s = open(true);
        QVERIFY(s);
        s->setLiveViewInterval(20);
        QCOMPARE(s->settings()->setValue("run/iterations", 60), QString());
        QSignalSpy live(s.get(), &Session::liveSceneChanged);
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(finished.wait(120000));
        const RunRecord* r = s->history().find(id);
        QCOMPARE(r->status, RunStatus::Completed);
        QCOMPARE(r->iterationsDone, 60);
        QVERIFY(!r->modified());
        QVERIFY(live.count() >= 1);
        QVERIFY(!s->live().setup->liveRunId.isEmpty());

        auto section = post(*s, id, PostKind::Probe,
                            {{"name", "Section"}, {"type", "section"}, {"zones", QJsonArray{"wall-elbow"}}, {"origin", QJsonArray{0, 0, 0.025}},
                             {"normal", QJsonArray{0, 0, 1}}, {"fields", QJsonArray{"pressure", "pressure@fluent"}}});
        QVERIFY2(section->ok(), qPrintable(*section ? QString() : section->error()));
        QVERIFY(section->value().part.lineCount >= 1);
        auto traces = post(*s, id, PostKind::Streamlines,
                           {{"name", "Pathlines"}, {"mode", "streamline"},
                            {"seeds", QJsonObject{{"type", "surface"}, {"zones", QJsonArray{"cold-inlet"}}, {"count", 10}, {"distribution", "uniform"}, {"offset", 0}}},
                            {"fields", QJsonArray{"velocity_magnitude"}}});
        QVERIFY2(traces->ok(), qPrintable(*traces ? QString() : traces->error()));
        QVERIFY(traces->value().scene->item("post:Section"));  // earlier results persist
        auto report = post(*s, id, PostKind::Integral,
                           {{"name", "Ports"}, {"zones", QJsonArray{"cold-inlet", "hot-inlet", "outlet"}}, {"fields", QJsonArray{"pressure", "pressure@fluent"}}});
        QVERIFY2(report->ok(), qPrintable(*report ? QString() : report->error()));
        // Mass is conserved through the ports (outward-positive flows sum to ~0).
        const double net = report->value().report.value("total").toObject().value("mass_flow").toDouble();
        const double through = report->value().report.value("total").toObject().value("absolute_mass_flow").toDouble();
        QVERIFY(through > 0 && std::abs(net) < 1e-3 * through);
        QVERIFY(!s->availableFields(id).empty());
    }

    void stopFinishesTheCurrentIteration()
    {
        auto s = open(false);
        QVERIFY(s);
        QCOMPARE(s->settings()->setValue("run/iterations", 5000), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QTRY_VERIFY_WITH_TIMEOUT(s->live().last && s->live().last->iter >= 20, 60000);
        s->stop();
        QVERIFY(finished.wait(60000));
        const RunRecord* r = s->history().find(id);
        QCOMPARE(r->status, RunStatus::Stopped);
        QVERIFY(r->iterationsDone >= 20 && r->iterationsDone < 5000);
        QVERIFY(s->history().forces(id));  // outputs were written
    }
};

QTEST_GUILESS_MAIN(TestRealCli)
#include "tst_realcli.moc"
