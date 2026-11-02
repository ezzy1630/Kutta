// Protocol parsing against recorded output of the real kutta CLI.

#include "TestSupport.h"
#include "core/Protocol.h"

#include <QJsonArray>
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
        QCOMPARE(int(cs.settings.size()), 9857);  // with controls/courant-number
        // Includes each inlet's and outlet's turbulence method and the values
        // of the methods it does not use, so switching method is possible,
        // the editable pressure scheme, the Coupled flow Courant number and
        // the limiter type.
        QCOMPARE(int(std::count_if(cs.settings.begin(), cs.settings.end(), [](const Setting& s) { return s.editable; })), 60);
        for (const char* path : {"controls/courant-number", "methods/limiter"})
            QVERIFY2(std::any_of(cs.settings.begin(), cs.settings.end(),
                                 [&](const Setting& s) { return s.path == path && s.editable; }), path);
        QVERIFY(cs.zone("cold-inlet"));
        QCOMPARE(cs.zone("cold-inlet")->type, QString("velocity-inlet"));

        // Enum codes carry Fluent's words, and fields follow their method.
        auto find = [&](const QString& p) -> const Setting& {
            static const Setting none;
            auto it = std::find_if(cs.settings.begin(), cs.settings.end(), [&](const Setting& s) { return s.path == p; });
            return it == cs.settings.end() ? none : *it;
        };
        const Setting& spec = find("boundary/cold-inlet/turbulence-specification");
        QCOMPARE(spec.value.toInt(), 3);
        QCOMPARE(spec.valueLabel, QString("Intensity and Hydraulic Diameter"));
        QCOMPARE(spec.labelFor(1), QString("Intensity and Length Scale"));
        QVERIFY(spec.editable);
        const Setting& k = find("boundary/cold-inlet/k");
        QCOMPARE(k.dependsOn, spec.path);
        QVERIFY(!k.active);
        QVERIFY(find("boundary/cold-inlet/turb-hydraulic-diam").active);
        QCOMPARE(find("boundary/cold-inlet/velocity-specification").valueLabel, QString("Magnitude, Normal to Boundary"));

        // Fluent's Methods panel: the three couplings, the gradient, and the
        // schemes of the equations SST solves (epsilon's is listed, unused).
        const Setting& coupling = find("methods/pressure-velocity-coupling");
        QCOMPARE(coupling.value.toString(), QString("coupled"));
        QCOMPARE(coupling.optionLabels, (QStringList{"SIMPLE", "SIMPLEC", "Coupled"}));
        QVERIFY(coupling.editable);
        QVERIFY(find("methods/pseudo-transient").value.toBool());
        QCOMPARE(find("methods/gradient").valueLabel, QString("Least Squares Cell Based"));
        QVERIFY(find("methods/omega-scheme").active);
        QVERIFY(!find("methods/epsilon-scheme").active);
        QCOMPARE(cs.dimension, std::optional<int>(3));
        QVERIFY(!cs.planar());
    }

    // Approximable issues carry their measured effect separately (newer CLIs);
    // the explanation then drops the copy the CLI appends to it.
    void approximationMeasuredEffect()
    {
        const CaseSummary cs = parseCaseSummary(fixtureJson("cases/mixing_elbow/inspect.json")).value();
        auto it = std::find_if(cs.issues.begin(), cs.issues.end(), [](const SupportIssue& i) { return i.approximation == "wfgc"; });
        QVERIFY(it != cs.issues.end());
        QVERIFY(it->measuredEffect.contains("3.366748%"));
        QVERIFY(it->explanation.endsWith(it->measuredEffect + '.'));
        QVERIFY(!it->explanationOnly().contains("3.366748%"));
        QVERIFY(it->explanationOnly().endsWith("--approximate wfgc."));
        SupportIssue older;  // no measured_effect: the explanation stands as sent
        older.explanation = "Text. effect.";
        QCOMPARE(older.explanationOnly(), older.explanation);
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

    // Source locations from kutta's internal checks move to the end, without
    // the build machine's directories; other text is untouched.
    void cliMessagesDropBuildPaths()
    {
        QCOMPARE(cliMessage("/build/kutta/src/kutta/post/request.cpp:30: unknown post-processing property: x"),
                 QString("unknown post-processing property: x (request.cpp:30)"));
        QCOMPARE(cliMessage("kutta export-view: C:\\b\\src\\view.cpp:7: bad view\nsecond line"),
                 QString("kutta export-view: bad view (view.cpp:7)\nsecond line"));
        QCOMPARE(cliMessage("cannot open HDF5 file: /tmp/a.cas.h5"), QString("cannot open HDF5 file: /tmp/a.cas.h5"));
    }

    void progressEvents()
    {
        auto setup = parseProgressLine(R"({"v":1,"event":"setup","cells":17822,"backend":"metal","device":"M4","approximations":["wfgc"],"init":"hybrid","overrides":{},"modified":false,"pressure_velocity_coupling":"simplec"})");
        QVERIFY(setup.ok());
        QCOMPARE(std::get<SetupEvent>(setup.value()).approximations, QStringList{"wfgc"});
        QCOMPARE(std::get<SetupEvent>(setup.value()).coupling, QString("simplec"));
        QCOMPARE(couplingLabel("simplec"), QString("SIMPLEC"));
        // The recorded run states the coupling that ran.
        auto recorded = parseProgressLine(fixtureBytes("cases/mixing_elbow/progress-cold.jsonl").split('\n').first());
        QCOMPARE(std::get<SetupEvent>(recorded.value()).coupling, QString("coupled"));
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
        QCOMPARE(residualOrder({"epsilon", "k", "z-velocity", "continuity"}), (QStringList{"continuity", "z-velocity", "k", "epsilon"}));
        // A 2D case has no z-velocity, whatever the extrusion reports.
        QVERIFY(!residualShown("z-velocity", true));
        QVERIFY(residualShown("z-velocity", false));
        QVERIFY(residualShown("y-velocity", true));
    }

    // kutta new-case --json (recorded from the real CLI on the NACA wing mesh).
    void wheelGeometryValidation()
    {
        WheelGeometry w{0.2, {0.4, 0.25, 0}, {0, 0, -2}};
        QVERIFY(w.error(2).isEmpty());
        QCOMPARE(w.rollingAxis(2), (std::array<double, 3>{0, 0, 1}));
        w.axis = {0, 2, 0};
        QVERIFY(w.error(3).isEmpty());
        QCOMPARE(w.rollingAxis(3), (std::array<double, 3>{0, -1, 0}));
        w.axis = {0, -2, 0.1};  // camber does not change rolling sense
        QVERIFY(w.error(3).isEmpty());
        QVERIFY(w.rollingAxis(3)[1] < 0);
        QVERIFY(w.rollingAxis(3)[2] > 0);
        w.axis = {1, 0, 0};
        QVERIFY(!w.error(3).isEmpty());
        w.axis = {0, 0, 0};
        QVERIFY(w.error(3).contains("cannot be zero"));
        w.radius = -0.2;
        QVERIFY(w.error(3).contains("radius greater than zero"));
    }

    void newCaseReply()
    {
        auto d = parseNewCase(fixtureJson("cases/naca4412_ground/new-case.json"));
        QVERIFY2(d.ok(), qPrintable(d.error()));
        QCOMPARE(d->dimension, 2);
        QCOMPARE(d->templateName, QString("fsae-external"));
        QCOMPARE(d->speed, 15.0);
        QCOMPARE(int(d->boundaries().size()), 5);
        QVERIFY(d->roles.contains("moving-ground"));
        const NewCaseZone* wing = d->zone("wing");
        QVERIFY(wing);
        QCOMPARE(wing->group, QString("Wing"));
        QCOMPARE(d->zone("ground")->suggestedRole, QString("moving-ground"));
        QVERIFY(d->kuttaCase.value("boundaries").isObject());
        QVERIFY(!parseNewCase(QJsonObject{{"v", 2}}));
        QVERIFY(!parseNewCase(QJsonObject{{"v", 1}}));
    }

    void kuttaCaseSummaryAndConvertReply()
    {
        QJsonObject o = fixtureJson("cases/naca4412_ground/inspect.json");
        auto fluent = parseCaseSummary(o);
        QVERIFY(fluent && !fluent->kuttaCase());  // recordings from before the field
        QJsonObject file = o.value("file").toObject();
        file["format"] = "kutta";
        file["mesh_file"] = QJsonObject{{"path", "/m/wing.msh"}, {"sha256", "ab"}};
        o["file"] = file;
        auto kutta = parseCaseSummary(o);
        QVERIFY(kutta && kutta->kuttaCase());
        QCOMPARE(kutta->meshFile, QString("/m/wing.msh"));

        auto r = parseConvert({{"v", 1}, {"out", "/c/wing.kutta.json"}, {"replaced_input", true},
                               {"omitted", QJsonArray{QJsonObject{{"path", "case-config/rf-energy?"}, {"reason", "energy"}, {"class", "informational"}}}},
                               {"run_controls_not_saved", QJsonArray{"run/iterations"}}});
        QVERIFY(r);
        QCOMPARE(r->out, QString("/c/wing.kutta.json"));
        QCOMPARE(int(r->omitted.size()), 1);
        QVERIFY(!r->omitted[0].changesSolution);
        QCOMPARE(r->runControlsNotSaved, QStringList{"run/iterations"});
        QVERIFY(!parseConvert({{"v", 1}}));
    }
};

QTEST_GUILESS_MAIN(TestProtocol)
#include "tst_protocol.moc"
