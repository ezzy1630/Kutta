// The session against the real kutta CLI on the Ansys mixing elbow. Skipped
// unless KUTTA_REAL_CLI (the binary) and KUTTA_ELBOW_CASE (mixing_elbow.cas.h5,
// with its .dat.h5 beside it) are set: CI has no solver build or case files.
// The 2D, laminar and k-epsilon cases also need KUTTA_FLUENT_EXAMPLES.

#include "core/SceneData.h"
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

    std::unique_ptr<Session> open(bool withData, const QString& other = {})
    {
        CliInfo cli;
        cli.program = m_cli;
        auto s = std::make_unique<Session>(cli, m_dir.filePath("data"));
        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->openCase(other.isEmpty() ? m_case : other, withData ? m_data : QString());
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

    // With checking on, the seeded elbow stops on the case's own residual
    // criteria; its wall loads give coefficients without a re-solve.
    void convergesOnTheCaseCriteriaThenAero()
    {
        auto s = open(true);
        QVERIFY(s);
        QCOMPARE(s->settings()->setValue("monitors/check-convergence", true), QString());
        QCOMPARE(s->settings()->setValue("run/iterations", 1000), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(finished.wait(300000));
        const RunRecord* r = s->history().find(id);
        QVERIFY(r->converged());
        QVERIFY(r->iterationsDone < 1000);
        QCOMPARE(r->convergedBy, QStringList{"residuals"});
        QVERIFY(!r->modified());
        const auto summary = s->history().summary(id);
        QVERIFY(summary && summary->convergence.enabled);
        QCOMPARE(summary->convergence.residuals.value("continuity"), 1e-3);

        std::optional<Result<AeroResult>> aero;
        s->aero(id, [&](Result<AeroResult> a) { aero = std::move(a); });
        QTRY_VERIFY_WITH_TIMEOUT(aero.has_value(), 60000);
        QVERIFY2(aero->ok(), qPrintable(*aero ? QString() : aero->error()));
        const AeroResult& a = aero->value();
        QVERIFY(a.approximate);
        QVERIFY(a.total() && a.total()->coefficients);
        QCOMPARE(a.balance.status, QString("needs_vehicle"));
        // The total is the case's walls, as in forces.csv (read as floats).
        const auto forces = s->history().forces(id);
        QVERIFY(forces);
        double drag = 0;
        for (const auto& row : forces->rows)
            if (row.kind == "wall") drag += row.total.x();
        QVERIFY(std::abs(a.total()->total.drag - drag) <= 1e-5 * std::abs(drag));
    }

    // The 2D laminar cylinder: a planar case seen as its plane, saved Coupled
    // without pseudo time (so it runs SIMPLEC), Green-Gauss gradients, no
    // turbulence, opt-in reverse flow with its measured effect.
    void planarLaminarCylinder()
    {
        const QString corpus = qEnvironmentVariable("KUTTA_FLUENT_EXAMPLES");
        const QString cylinder = corpus + "/adjoint_cylinder_tutorial/cylinder_tutorial.cas";
        if (corpus.isEmpty() || !QFileInfo::exists(cylinder)) QSKIP("set KUTTA_FLUENT_EXAMPLES for the 2D cylinder");
        auto s = open(false, cylinder);
        QVERIFY(s);
        const CaseSummary& cs = s->caseSummary();
        QVERIFY(cs.planar());
        QCOMPARE(cs.viscousModel, QString("laminar"));
        auto issue = std::find_if(cs.issues.begin(), cs.issues.end(), [](const SupportIssue& i) { return i.approximation == "reverse-flow"; });
        QVERIFY(issue != cs.issues.end() && !issue->measuredEffect.isEmpty());
        SettingsModel* m = s->settings();
        QCOMPARE(m->find("methods/gradient")->valueLabel, QString("Green-Gauss Cell Based"));
        QVERIFY(!m->isActive("methods/k-scheme") && !m->isActive("boundary/inlet/k"));
        QCOMPARE(m->value("methods/pressure-velocity-coupling").toString(), QString("coupled"));
        QVERIFY(!m->value("methods/pseudo-transient").toBool());
        QCOMPARE(m->setValue("run/iterations", 40), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(finished.wait(300000));
        const RunRecord* r = s->history().find(id);
        QCOMPARE(r->status, RunStatus::Completed);
        QCOMPARE(r->coupling, QString("simplec"));
        const RunProgress p = s->history().progress(id);
        QVERIFY(!p.residuals.names.contains("k"));
        bool loaded = false;
        s->loadScene("run:" + id, [&](const QString&, Result<std::shared_ptr<SceneData>> scene) {
            QVERIFY2(scene.ok(), qPrintable(scene ? QString() : scene.error()));
            QVERIFY(scene.value()->planar);
            QCOMPARE(defaultVisibleItems(*scene.value()), QStringList{planarFace(*scene.value())});
            loaded = true;
        });
        QTRY_VERIFY_WITH_TIMEOUT(loaded, 60000);
    }

    // Standard k-epsilon with SIMPLE: epsilon residuals and scheme, and no
    // Coupled choice without pseudo time.
    void kEpsilonElbow()
    {
        const QString corpus = qEnvironmentVariable("KUTTA_FLUENT_EXAMPLES");
        const QString param = corpus + "/mixing_elbow/elbow_param.cas.h5";
        if (corpus.isEmpty() || !QFileInfo::exists(param)) QSKIP("set KUTTA_FLUENT_EXAMPLES for elbow_param");
        auto s = open(false, param);
        QVERIFY(s);
        SettingsModel* m = s->settings();
        QVERIFY(m->isActive("methods/epsilon-scheme") && !m->isActive("methods/omega-scheme"));
        QVERIFY(!m->find("methods/pressure-velocity-coupling")->options.contains(QJsonValue("coupled")));
        QCOMPARE(m->setValue("run/iterations", 40), QString());
        QSignalSpy finished(s.get(), &Session::runFinished);
        const QString id = s->run();
        QVERIFY(finished.wait(300000));
        const RunRecord* r = s->history().find(id);
        QCOMPARE(r->coupling, QString("simple"));
        QVERIFY(s->history().progress(id).residuals.names.contains("epsilon"));
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
