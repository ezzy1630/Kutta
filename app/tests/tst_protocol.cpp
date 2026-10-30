// Protocol parsing against recorded output of the real kutta CLI.

#include "TestSupport.h"
#include "core/Protocol.h"

#include <QtTest>

using namespace kapp;

class TestProtocol : public QObject {
    Q_OBJECT
private slots:
    void inventoryFromRealInspect()
    {
        auto c = parseCaseSummary(fixtureJson("cases/mixing_elbow/inspect.json"));
        QVERIFY2(c.ok(), qPrintable(c ? QString() : c.error()));
        const CaseSummary& cs = c.value();
        QCOMPARE(cs.cells, 17822);
        QCOMPARE(int(cs.settings.size()), 9844);
        QCOMPARE(int(std::count_if(cs.settings.begin(), cs.settings.end(), [](const Setting& s) { return s.editable; })), 26);
        QVERIFY(cs.zone("cold-inlet"));
        QCOMPARE(cs.zone("cold-inlet")->type, QString("velocity-inlet"));
    }

    void approximationsBlockUntilApproved()
    {
        const CaseSummary cs = parseCaseSummary(fixtureJson("cases/mixing_elbow/inspect.json")).value();
        QCOMPARE(cs.status, SupportStatus::Approximate);
        QVERIFY(cs.approximations().contains("wfgc"));
        QVERIFY(!cs.blocking({}).empty());
        QVERIFY(cs.blocking({"wfgc"}).empty());
    }

    void settingPathsSplitIntoOwnerAndLeaf()
    {
        const CaseSummary cs = parseCaseSummary(fixtureJson("cases/mixing_elbow/inspect.json")).value();
        auto find = [&](const QString& p) {
            auto it = std::find_if(cs.settings.begin(), cs.settings.end(), [&](const Setting& s) { return s.path == p; });
            return it == cs.settings.end() ? nullptr : &*it;
        };
        const Setting* v = find("boundary/cold-inlet/velocity-magnitude");
        QVERIFY(v);
        QCOMPARE(v->section, QString("boundary"));
        QCOMPARE(v->owner, QString("cold-inlet"));
        QCOMPARE(v->leaf, QString("velocity-magnitude"));
        QVERIFY(!v->raw);
        QVERIFY(v->editable);
        const Setting* raw = find("boundary/cold-inlet/raw/ai");
        QVERIFY(raw && raw->raw && !raw->editable && !raw->reason.isEmpty());
        const Setting* flag = find("general/case-config/rp-3d%3f");
        QVERIFY(flag && flag->raw);
        QCOMPARE(decodePathSegment("%43ortex%20%56ariables"), QString("Cortex Variables"));
    }

    void progressEvents()
    {
        auto setup = parseProgressLine(R"({"v":1,"event":"setup","cells":17822,"backend":"metal","device":"M4","approximations":["wfgc"],"init":"hybrid","overrides":{},"modified":false})");
        QVERIFY(setup.ok());
        QCOMPARE(std::get<SetupEvent>(setup.value()).approximations, QStringList{"wfgc"});
        auto it = parseProgressLine(R"({"v":1,"event":"iteration","iter":20,"residuals":{"continuity":0.1,"k":null},"reports":{"drag":{"force":1.5,"coefficient":null}},"sec_per_iter":0.01})");
        QVERIFY(it.ok());
        const auto& e = std::get<IterationEvent>(it.value());
        QCOMPARE(e.iter, 20);
        QCOMPARE(e.residuals.value("continuity"), 0.1);
        QVERIFY(!e.residuals.contains("k"));  // null means not reported
        QCOMPARE(*e.reports.value("drag").force, 1.5);
        QVERIFY(!e.reports.value("drag").coefficient);
        auto done = parseProgressLine(R"({"v":1,"event":"done","converged":false,"iterations":600,"elapsed":3.5,"out":"/x","stop_reason":"stopped"})");
        QCOMPARE(std::get<DoneEvent>(done.value()).stopReason, QString("stopped"));
        auto snap = parseProgressLine(R"({"v":1,"event":"snapshot","iter":40,"view":"/tmp/v/manifest.json"})");
        QCOMPARE(std::get<SnapshotEvent>(snap.value()).iter, 40);
        QVERIFY(!parseProgressLine("not json"));
        QVERIFY(!parseProgressLine(R"({"v":2,"event":"done"})"));
    }

    void rejectedRunCarriesIssues()
    {
        const QByteArray line = fixtureBytes("cases/mixing_elbow/rejected.jsonl").split('\n').first();
        auto r = parseProgressLine(line);
        QVERIFY(r.ok());
        const auto& e = std::get<ErrorEvent>(r.value());
        QVERIFY(e.message.contains("--approximate wfgc"));
        QVERIFY(std::any_of(e.issues.begin(), e.issues.end(), [](const SupportIssue& i) { return i.cls == IssueClass::Approximable; }));
    }

    void forcesCsvKeepsNotice()
    {
        auto f = parseForcesCsv(fixtureBytes("cases/mixing_elbow/run-cold/forces.csv"));
        QVERIFY(f.ok());
        QVERIFY(f.value().notice.startsWith("APPROXIMATE RESULT"));
        QVERIFY(!f.value().rows.empty());
        QCOMPARE(f.value().rows.front().name, QString("wall-elbow"));
        QCOMPARE(f.value().rows.front().zoneId, 34);
    }

    void residualsInFluentOrder()
    {
        QCOMPARE(residualOrder({"omega", "k", "y-velocity", "continuity", "x-velocity"}),
                 (QStringList{"continuity", "x-velocity", "y-velocity", "k", "omega"}));
    }
};

QTEST_GUILESS_MAIN(TestProtocol)
#include "tst_protocol.moc"
