// Run records persist, number themselves and diff against each other.

#include "TestSupport.h"
#include "core/RunHistory.h"

#include <QJsonDocument>
#include <QTemporaryDir>

#include <algorithm>
#include <QtTest>

using namespace kapp;

class TestRunHistory : public QObject {
    Q_OBJECT
private slots:
    void createSaveReload()
    {
        QTemporaryDir root;
        RunHistory h;
        h.setRoot(root.path());
        h.setCase("/cases/Mixing Elbow.cas.h5", "8145813f50f989a521ed8d1b");
        QVERIFY(h.caseDir().endsWith("Mixing_Elbow-8145813f50f9"));
        RunRecord a = h.create();
        QCOMPARE(a.id, QString("run-0001"));
        a.status = RunStatus::Completed;
        a.overrides = QJsonObject{{"controls/mom-relaxation", 0.3}};
        h.save(a);
        RunRecord b = h.create();
        QCOMPARE(b.number, 2);
        h.save(b);  // still running when the app went away

        RunHistory again;
        again.setRoot(root.path());
        again.setCase("/cases/Mixing Elbow.cas.h5", "8145813f50f989a521ed8d1b");
        QCOMPARE(int(again.runs().size()), 2);
        QCOMPARE(again.find("run-0001")->overrides, a.overrides);
        QCOMPARE(again.find("run-0002")->status, RunStatus::Interrupted);
        QCOMPARE(again.latestFinished()->id, QString("run-0001"));
    }

    void removeDeletesOnlyTheRun()
    {
        QTemporaryDir root;
        RunHistory h;
        h.setRoot(root.path());
        h.setCase("/c/x.cas.h5", "abc");
        RunRecord r = h.create();
        r.status = RunStatus::Completed;
        h.save(r);
        QVERIFY(QFileInfo::exists(r.dir));
        QVERIFY(h.remove(r.id));
        QVERIFY(!QFileInfo::exists(r.dir));
        QVERIFY(QFileInfo::exists(h.caseDir()));
        QVERIFY(!h.remove("run-9999"));
    }

    void progressLogFromRealRun()
    {
        const RunProgress p = parseProgressLog(fixtureBytes("cases/mixing_elbow/progress-cold.jsonl"));
        QVERIFY(p.setup);
        QCOMPARE(p.setup->approximations, QStringList{"wfgc"});
        QVERIFY(p.done);
        // Iterations 1-5 (residuals only), then every monitor interval.
        QCOMPARE(int(p.residuals.iterations.size()), 5 + p.done->iterations / 20);
        QCOMPARE(p.residuals.iterations.front(), 1);
        QVERIFY(p.residuals.names.contains("continuity"));
        for (const auto& series : p.residuals.values) QCOMPARE(series.size(), p.residuals.iterations.size());
        QVERIFY(!p.warnings.isEmpty());
        QVERIFY(p.warnings.first().startsWith("APPROXIMATE RESULT"));
    }

    void diffShowsOnlyDifferences()
    {
        Setting s;
        s.path = "boundary/inlet/velocity-magnitude";
        s.label = "Velocity magnitude";
        s.unit = "m/s";
        s.value = 16.0;
        RunRecord a, b;
        a.overrides = QJsonObject{{s.path, 20.0}};
        a.approximations = QStringList{"wfgc"};
        b.approximations = QStringList{"wfgc"};
        a.device = b.device = "M4";
        const auto d = diffRuns(a, b, {s});
        QCOMPARE(int(d.size()), 1);
        QCOMPARE(d[0].a.toDouble(), 20.0);
        QCOMPARE(d[0].b.toDouble(), 16.0);  // B uses the case value
        b.overrides = a.overrides;
        QVERIFY(diffRuns(a, b, {s}).empty());
    }

    // Why a run stopped, in words, for each of the CLI's stop reasons.
    void stopReasonsInWords()
    {
        ConvergenceSetup conv;
        conv.enabled = true;
        conv.residuals = {{"continuity", 1e-3}};
        conv.reports.push_back({"drag", 1e-3, 5, 100});
        QCOMPARE(stopReasonText("converged", {"residuals"}, conv), QString("Converged: every checked residual reached its criterion"));
        QVERIFY(stopReasonText("converged", {"report:drag"}, conv).contains("drag changed less than 0.1 %"));
        QCOMPARE(stopReasonText("stopped", {}, conv), QString("Stopped by you; not converged"));
        QCOMPARE(stopReasonText("iteration_limit", {}, conv), QString("Iteration limit; convergence criteria not met"));
        conv.enabled = false;
        QCOMPARE(stopReasonText("iteration_limit", {}, conv), QString("Iteration limit; convergence checking was off"));
        QVERIFY(stopReasonText("", {}, conv).isEmpty());  // an older CLI: no reason given
    }

    // Poor mesh numerics is named only when it touched cells, and an older
    // CLI's summary without the count says nothing.
    void poorMeshNumericsOnlyWhenUsed()
    {
        auto parse = [](const char* json) { return parseRunSummary(QJsonDocument::fromJson(json).object()).value(); };
        QVERIFY(poorMeshNumericsText(parse(R"({"v":1,"cells":10})")).isEmpty());
        QVERIFY(poorMeshNumericsText(parse(R"({"v":1,"cells":10,"poor_quality_cells":0,"poor_cell_quality":0.05})")).isEmpty());
        QCOMPARE(poorMeshNumericsText(parse(R"({"v":1,"cells":10,"poor_quality_cells":4,"poor_cell_quality":0.050000000000000003})")),
                 QString("4 poor-quality cells (orthogonal quality ≤ 0.05) used low-order face pressures, as Fluent's poor mesh numerics does"));
        QVERIFY(poorMeshNumericsText(parse(R"({"v":1,"cells":10,"poor_quality_cells":1})")).startsWith("1 poor-quality cell used"));
    }

    // Run History's Status column names how a completed run ended.
    void statusSaysHowARunEnded()
    {
        RunRecord r;
        r.status = RunStatus::Completed;
        r.stopReason = "iteration_limit";
        QCOMPARE(runStatusText(r), QString("Iteration limit"));
        r.stopReason = "converged";
        QCOMPARE(runStatusText(r), QString("Converged"));
        r.stopReason.clear();  // a record from an older CLI
        QCOMPARE(runStatusText(r), QString("Completed"));
        r.status = RunStatus::Stopped;
        r.stopReason = "stopped";
        QCOMPARE(runStatusText(r), QString("Stopped"));
    }

    // The residual chart's marks come from the run's own events: equal
    // criteria share one line, a 2D case has no z-velocity criterion, and the
    // stop is marked with its reason; the last condition states are kept.
    void convergenceMarksFromProgress()
    {
        const QByteArray log =
            R"({"v":1,"event":"setup","cells":10,"backend":"cpu","device":"cpu","approximations":[],"init":"hybrid","pressure_velocity_coupling":"simplec","convergence":{"enabled":true,"condition":"all","residuals":{"continuity":1e-05,"x-velocity":1e-05,"y-velocity":1e-05,"z-velocity":1e-03},"reports":[]}})" "\n"
            R"({"v":1,"event":"iteration","iter":20,"residuals":{"continuity":1e-3,"x-velocity":1e-6,"y-velocity":1e-6,"z-velocity":1e-11},"convergence":{"met":false,"conditions":[{"kind":"residuals","met":false,"equations":{"continuity":{"value":1e-3,"criterion":1e-5,"met":false}}}]}})" "\n"
            R"({"v":1,"event":"done","converged":false,"iterations":20,"elapsed":1,"out":"x","stop_reason":"stopped"})" "\n";
        const RunProgress p = parseProgressLog(log);
        QCOMPARE(p.setup->coupling, QString("simplec"));
        QVERIFY(p.convergence && !p.convergence->met);
        QCOMPARE(p.convergenceIter, 20);
        auto marks = convergenceMarks(p, true);
        QCOMPARE(int(marks.size()), 2);
        QCOMPARE(marks[0].kind, ConvergenceMark::Kind::Criterion);
        QCOMPARE(marks[0].label, QString("criterion 1e-5"));
        QCOMPARE(marks[1].kind, ConvergenceMark::Kind::Stopped);
        QCOMPARE(marks[1].value, 20.0);
        marks = convergenceMarks(p, false);  // 3D: two distinct criteria, each named
        QCOMPARE(int(marks.size()), 3);
        QCOMPARE(marks[0].label, QString("criterion 1e-5 · continuity, x-velocity, y-velocity"));
        QCOMPARE(marks[1].label, QString("criterion 1e-3 · z-velocity"));
    }

    // The CLI stops a run whose forces settled while its residuals stalled
    // (residuals waived) with converged:true. The app must not call that
    // converged: it says the forces settled and the residuals did not.
    void forcesSettledWithWaivedResidualsIsNotConverged()
    {
        const QByteArray log =
            R"({"v":1,"event":"setup","cells":10,"backend":"metal","device":"gpu","approximations":[],"convergence":{"enabled":true,"condition":"all","residuals":{"continuity":0.001},"reports":[],"forces":{"tolerance":0.001,"window":300,"reports":["drag","downforce"]}}})" "\n"
            R"({"v":1,"event":"iteration","iter":1200,"residuals":{"continuity":0.017},"convergence":{"met":true,"conditions":[{"kind":"residuals","met":false,"waived":true,"equations":{"continuity":{"value":0.017,"criterion":0.001,"met":false}}},{"kind":"forces","met":true,"settled":true,"criterion":0.001,"window":300,"reports":{"drag":0.0004,"downforce":null}}]}})" "\n"
            R"({"v":1,"event":"done","converged":true,"iterations":1200,"elapsed":1,"out":"x","stop_reason":"converged","residuals_waived":true,"convergence":{"met_by":["forces"]}})" "\n";
        const RunProgress p = parseProgressLog(log);
        QVERIFY(p.setup->convergence.forces && p.setup->convergence.forces->window == 300);
        QVERIFY(p.done && p.done->residualsWaived);
        const auto& conditions = p.convergence->conditions;
        QVERIFY(conditions[0].waived && !conditions[0].met);
        QCOMPARE(conditions[1].kind, QString("forces"));
        QVERIFY(conditions[1].settled && conditions[1].ranges.value("drag") && !conditions[1].ranges.value("downforce"));
        const auto marks = convergenceMarks(p, false);
        QVERIFY(std::none_of(marks.begin(), marks.end(), [](const ConvergenceMark& m) { return m.kind == ConvergenceMark::Kind::Converged; }));
        RunRecord r;
        r.id = "run-0001";
        r.status = RunStatus::Completed;
        r.stopReason = "converged";
        r.residualsWaived = true;
        QVERIFY(!runStatusText(r).startsWith("Converged"));
        const auto back = RunRecord::fromJson(r.toJson(), "/x");
        QVERIFY(back && back->residualsWaived);
        const QString why = stopReasonText(r.stopReason, p.done->metBy, p.setup->convergence, true);
        QVERIFY2(why.contains("not converged") && why.contains("drag, downforce"), qPrintable(why));
    }

    void couplingPersistsWithTheRun()
    {
        RunRecord r;
        r.id = "run-0001";
        r.coupling = "simplec";
        const auto back = RunRecord::fromJson(r.toJson(), "/x");
        QVERIFY(back);
        QCOMPARE(back->coupling, QString("simplec"));
    }
};

QTEST_GUILESS_MAIN(TestRunHistory)
#include "tst_runhistory.moc"
