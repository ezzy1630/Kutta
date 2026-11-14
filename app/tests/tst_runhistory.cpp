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

    // Why a run stopped, in words, for each verdict.
    void stopReasonsInWords()
    {
        ConvergenceSetup conv;
        conv.enabled = true;
        conv.residuals = {{"continuity", 1e-3}};
        conv.reports.push_back({"drag", 1e-3, 5, 100});
        QCOMPARE(stopReasonText(Verdict::Converged, {"residuals"}, conv), QString("Converged: every checked residual reached its criterion"));
        QVERIFY(stopReasonText(Verdict::Converged, {"report:drag"}, conv).contains("drag changed less than 0.1 %"));
        const QString waived = stopReasonText(Verdict::ForcesSettledResidualsWaived, {"forces"}, conv);
        QVERIFY2(waived.startsWith("Stopped on settled forces") && waived.contains("not converged"), qPrintable(waived));
        QCOMPARE(stopReasonText(Verdict::Stopped, {}, conv), QString("Stopped by you; not converged"));
        QCOMPARE(stopReasonText(Verdict::IterationLimit, {}, conv), QString("Iteration limit; convergence criteria not met"));
        QVERIFY(stopReasonText(Verdict::StatisticallySteady, {}, conv).startsWith("Statistically steady; not converged"));
        QVERIFY(stopReasonText(Verdict::Diverged, {}, conv).startsWith("Diverged"));
        conv.enabled = false;
        QCOMPARE(stopReasonText(Verdict::IterationLimit, {}, conv), QString("Iteration limit; convergence checking was off"));
        QVERIFY(stopReasonText(Verdict::Unknown, {}, conv).isEmpty());  // an older CLI: no reason given
    }

    // Every verdict the CLI sends, read from done, summary.json and a
    // diverged run's error event. The verdict wins over converged and
    // stop_reason: a waived run is converged:true, stop_reason "converged".
    void verdictReadFromEveryEvent_data()
    {
        QTest::addColumn<QString>("key");
        QTest::addColumn<bool>("converged");
        QTest::addColumn<QString>("stopReason");
        QTest::addColumn<int>("verdict");
        QTest::newRow("converged") << "converged" << true << "converged" << int(Verdict::Converged);
        QTest::newRow("waived") << "forces_settled_residuals_waived" << true << "converged" << int(Verdict::ForcesSettledResidualsWaived);
        QTest::newRow("steady") << "statistically_steady" << false << "statistically_steady" << int(Verdict::StatisticallySteady);
        QTest::newRow("limit") << "iteration_limit" << false << "iteration_limit" << int(Verdict::IterationLimit);
        QTest::newRow("stopped") << "stopped" << false << "stopped" << int(Verdict::Stopped);
    }
    void verdictReadFromEveryEvent()
    {
        QFETCH(QString, key);
        QFETCH(bool, converged);
        QFETCH(QString, stopReason);
        QFETCH(int, verdict);
        // residuals_waived left out: the verdict alone must say a waived run is not converged.
        const QString done = QStringLiteral(R"({"v":1,"event":"done","converged":%1,"iterations":40,"elapsed":1,"out":"x","stop_reason":"%2","verdict":"%3"})")
                                 .arg(converged ? "true" : "false", stopReason, key);
        auto ev = parseProgressLine(done.toUtf8());
        QVERIFY(ev);
        QCOMPARE(int(std::get<DoneEvent>(ev.value()).verdict), verdict);
        QCOMPARE(verdictKey(Verdict(verdict)), key);
        const QString summary = QStringLiteral(R"({"v":1,"cells":10,"converged":%1,"stop_reason":"%2","verdict":"%3"})")
                                    .arg(converged ? "true" : "false", stopReason, key);
        QCOMPARE(int(parseRunSummary(QJsonDocument::fromJson(summary.toUtf8()).object()).value().verdict), verdict);
        // The record keeps it, and only "converged" is converged.
        RunRecord r;
        r.id = "run-0001";
        r.status = Verdict(verdict) == Verdict::Stopped ? RunStatus::Stopped : RunStatus::Completed;
        r.stopReason = stopReason;
        r.verdict = Verdict(verdict);
        const auto back = RunRecord::fromJson(r.toJson(), "/x");
        QVERIFY(back);
        QCOMPARE(int(back->verdict), verdict);
        QCOMPARE(back->converged(), Verdict(verdict) == Verdict::Converged);
        QCOMPARE(runStatusText(*back) == QLatin1String("Converged"), Verdict(verdict) == Verdict::Converged);
    }

    // A diverged run's error event carries verdict "diverged"; other errors have none.
    void divergedErrorIsAVerdict()
    {
        auto ev = parseProgressLine(R"({"v":1,"event":"error","message":"the solution diverged at iteration 37: continuity 3.2e+02","verdict":"diverged"})");
        QVERIFY(ev);
        QCOMPARE(std::get<ErrorEvent>(ev.value()).verdict, Verdict::Diverged);
        ev = parseProgressLine(R"({"v":1,"event":"error","message":"out of memory"})");
        QCOMPARE(std::get<ErrorEvent>(ev.value()).verdict, Verdict::Unknown);
        const QByteArray log =
            R"({"v":1,"event":"setup","cells":10,"backend":"cpu","device":"cpu","approximations":[],"convergence":{"enabled":false}})" "\n"
            R"({"v":1,"event":"iteration","iter":36,"residuals":{"continuity":1e2}})" "\n"
            R"({"v":1,"event":"error","message":"the solution diverged at iteration 37","verdict":"diverged"})" "\n";
        const auto marks = convergenceMarks(parseProgressLog(log), false);
        QCOMPARE(int(marks.size()), 1);
        QCOMPARE(marks[0].kind, ConvergenceMark::Kind::Diverged);
        QCOMPARE(marks[0].value, 36.0);
        RunRecord r;
        r.id = "run-0001";
        r.status = RunStatus::Failed;
        r.verdict = Verdict::Diverged;
        QCOMPARE(runStatusText(r), QString("Diverged"));
        QCOMPARE(RunRecord::fromJson(r.toJson(), "/x")->verdict, Verdict::Diverged);
        r.verdict = Verdict::Unknown;
        QCOMPARE(runStatusText(r), QString("Failed"));
    }

    // Older CLIs send no verdict: it comes from stop_reason and residuals_waived
    // as before, in events, summaries and the app's own older records.
    void verdictFallbackWithoutOne()
    {
        auto done = [](const char* json) { return std::get<DoneEvent>(parseProgressLine(json).value()).verdict; };
        QCOMPARE(done(R"({"v":1,"event":"done","converged":true,"iterations":9,"elapsed":1,"out":"x","stop_reason":"converged"})"), Verdict::Converged);
        QCOMPARE(done(R"({"v":1,"event":"done","converged":true,"iterations":9,"elapsed":1,"out":"x","stop_reason":"converged","residuals_waived":true})"),
                 Verdict::ForcesSettledResidualsWaived);
        QCOMPARE(done(R"({"v":1,"event":"done","converged":true,"iterations":9,"elapsed":1,"out":"x","stop_reason":"converged","convergence":{"residuals_waived":true}})"),
                 Verdict::ForcesSettledResidualsWaived);
        QCOMPARE(done(R"({"v":1,"event":"done","converged":false,"iterations":9,"elapsed":1,"out":"x","stop_reason":"statistically_steady"})"), Verdict::StatisticallySteady);
        QCOMPARE(done(R"({"v":1,"event":"done","converged":false,"iterations":9,"elapsed":1,"out":"x","stop_reason":"iteration_limit"})"), Verdict::IterationLimit);
        QCOMPARE(done(R"({"v":1,"event":"done","converged":false,"iterations":9,"elapsed":1,"out":"x","stop_reason":"stopped"})"), Verdict::Stopped);
        QCOMPARE(done(R"({"v":1,"event":"done","converged":true,"iterations":9,"elapsed":1,"out":"x"})"), Verdict::Unknown);  // before stop_reason
        // A verdict this app does not know is never taken as converged, whatever converged says.
        QCOMPARE(done(R"({"v":1,"event":"done","converged":true,"iterations":9,"elapsed":1,"out":"x","stop_reason":"converged","verdict":"something_new"})"), Verdict::Unknown);
        auto summary = [](const char* json) { return parseRunSummary(QJsonDocument::fromJson(json).object()).value().verdict; };
        QCOMPARE(summary(R"({"v":1,"cells":1,"converged":true,"stop_reason":"converged","convergence":{"residuals_waived":true}})"), Verdict::ForcesSettledResidualsWaived);
        QCOMPARE(summary(R"({"v":1,"cells":1,"converged":true,"stop_reason":"converged"})"), Verdict::Converged);
        // An app-run.json from before verdicts.
        auto record = [](const char* json) { return RunRecord::fromJson(QJsonDocument::fromJson(json).object(), "/x")->verdict; };
        QCOMPARE(record(R"({"v":1,"id":"run-0001","status":"completed","stop_reason":"converged","residuals_waived":true})"), Verdict::ForcesSettledResidualsWaived);
        QCOMPARE(record(R"({"v":1,"id":"run-0001","status":"completed","stop_reason":"converged"})"), Verdict::Converged);
        QCOMPARE(record(R"({"v":1,"id":"run-0001","status":"completed","stop_reason":"iteration_limit"})"), Verdict::IterationLimit);
        QCOMPARE(record(R"({"v":1,"id":"run-0001","status":"completed"})"), Verdict::Unknown);
        // The chart marks follow the verdict: a waived run is not marked converged.
        const QByteArray log =
            R"({"v":1,"event":"iteration","iter":9,"residuals":{"continuity":1e-2}})" "\n"
            R"({"v":1,"event":"done","converged":true,"iterations":9,"elapsed":1,"out":"x","stop_reason":"converged","verdict":"forces_settled_residuals_waived"})" "\n";
        const auto marks = convergenceMarks(parseProgressLog(log), false);
        QCOMPARE(int(marks.size()), 1);
        QVERIFY(marks[0].kind != ConvergenceMark::Kind::Converged);
        QVERIFY(marks[0].label.startsWith("Forces settled"));
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

    // Run History's Status column names how a run ended, from its verdict.
    void statusSaysHowARunEnded()
    {
        RunRecord r;
        r.status = RunStatus::Completed;
        r.verdict = Verdict::IterationLimit;
        QCOMPARE(runStatusText(r), QString("Iteration limit"));
        r.verdict = Verdict::Converged;
        QCOMPARE(runStatusText(r), QString("Converged"));
        r.verdict = Verdict::ForcesSettledResidualsWaived;
        QCOMPARE(runStatusText(r), QString("Forces settled; residuals not converged"));
        r.verdict = Verdict::StatisticallySteady;
        QCOMPARE(runStatusText(r), QString("Statistically steady"));
        r.verdict = Verdict::Unknown;  // a record from an older CLI
        QCOMPARE(runStatusText(r), QString("Completed"));
        r.status = RunStatus::Stopped;
        r.verdict = Verdict::Stopped;
        QCOMPARE(runStatusText(r), QString("Stopped"));
        r.status = RunStatus::Failed;
        r.verdict = Verdict::Diverged;
        QCOMPARE(runStatusText(r), QString("Diverged"));
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
        r.verdict = p.done->verdict;  // no verdict in this older log: from residuals_waived
        QCOMPARE(r.verdict, Verdict::ForcesSettledResidualsWaived);
        QVERIFY(!runStatusText(r).startsWith("Converged"));
        QVERIFY(!r.converged());
        const auto back = RunRecord::fromJson(r.toJson(), "/x");
        QVERIFY(back && back->verdict == Verdict::ForcesSettledResidualsWaived);
        const QString why = stopReasonText(r.verdict, p.done->metBy, p.setup->convergence);
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
