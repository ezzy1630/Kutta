// Run records persist, number themselves and diff against each other.

#include "TestSupport.h"
#include "core/RunHistory.h"

#include <QTemporaryDir>
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
        QCOMPARE(int(p.residuals.iterations.size()), p.done->iterations / 20);
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
};

QTEST_GUILESS_MAIN(TestRunHistory)
#include "tst_runhistory.moc"
