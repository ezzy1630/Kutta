// Aero coefficients, balance and convergence criteria: the fixture CLI's
// `aero` against the hand-computed cases tests/post/test_aero.cpp checks the
// real CLI with, the app's parsing of both, and a session run that stops on
// its criteria.

#include "core/Protocol.h"
#include "core/Session.h"
#include "core/SettingsModel.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QProcess>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

namespace {

QJsonArray v3(double x, double y, double z) { return {x, y, z}; }
QVector3D cross(QVector3D a, QVector3D b) { return QVector3D::crossProduct(a, b); }
QJsonArray v3(QVector3D v) { return {v.x(), v.y(), v.z()}; }

// A point load: moments about the origin are r x F.
QJsonObject zone(const QString& name, QVector3D pressure, QVector3D viscous, QVector3D at)
{
    return {{"name", name}, {"zone_id", 7}, {"area", 1.0}, {"pressure", v3(pressure)}, {"viscous", v3(viscous)},
            {"pressure_moment", v3(cross(at, pressure))}, {"viscous_moment", v3(cross(at, viscous))}};
}

QJsonObject loads(const QJsonArray& zones, QJsonObject reference = {})
{
    if (reference.isEmpty()) reference = {{"area", 1.0}, {"density", 1.2}, {"velocity", 20.0}, {"length", 2.5}, {"pressure", 0.0}};
    return {{"v", 1}, {"origin", v3(0, 0, 0)}, {"units", QJsonObject{{"force", "N"}, {"moment", "N m"}, {"area", "m^2"}}},
            {"reference", reference}, {"report_zones", QJsonArray{}}, {"zones", zones},
            {"provenance", "APPROXIMATE RESULT: test"}, {"approximate", true}, {"modified", false}, {"iterations", 600}};
}

const AeroRow* row(const AeroResult& a, const QString& name)
{
    for (const AeroRow& r : a.rows)
        if (r.name == name) return &r;
    return nullptr;
}

} // namespace

class TestAero : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;
    int m_runs = 0;

    // Runs the fixture's `kutta aero DIR --aero spec --json` on a loads.json.
    Result<AeroResult> aero(const QJsonObject& loadsJson, const QJsonObject& spec, int* exitCode = nullptr)
    {
        const QString dir = m_dir.filePath(QStringLiteral("run-%1").arg(++m_runs));
        QDir().mkpath(dir);
        auto write = [](const QString& path, const QJsonObject& o) {
            QFile f(path);
            if (f.open(QIODevice::WriteOnly)) f.write(QJsonDocument(o).toJson());
        };
        write(dir + "/loads.json", loadsJson);
        write(dir + "/spec.json", spec);
        QProcess p;
        p.start(QStringLiteral(KUTTA_FIXTURE_CLI), {"aero", dir, "--aero", dir + "/spec.json", "--json"});
        p.waitForFinished(20000);
        if (exitCode) *exitCode = p.exitCode();
        const QJsonObject event = QJsonDocument::fromJson(p.readAllStandardOutput().trimmed()).object();
        if (event.value("event").toString() != "aero")
            return Result<AeroResult>::failure(event.value("message").toString("no aero event"));
        return parseAero(event);
    }

private slots:
    void initTestCase()
    {
        QCoreApplication::setOrganizationName("KuttaTest");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_dir.filePath("settings"));
        qputenv("KUTTA_FIXTURE_SPEED", "400");
    }

    // F=(100,0,-1000) N at (1,0,0.3): M=(0,1030,0); q=240 Pa, A=1, L=2.5,
    // front axle at the origin, wheelbase 2.5.
    void handCaseCoefficientsAndBalance()
    {
        const auto r = aero(loads({zone("car", {100, 0, -1000}, {}, {1, 0, 0.3f})}),
                            {{"v", 1}, {"vehicle", QJsonObject{{"front_axle", v3(0, 0, 0)}, {"wheelbase", 2.5}}}});
        QVERIFY2(r.ok(), qPrintable(r ? QString() : r.error()));
        const AeroResult& a = r.value();
        QCOMPARE(*a.dynamicPressure, 240.0);
        const AeroRow* t = a.total();
        QVERIFY(t && t->coefficients);
        QCOMPARE(t->total.drag, 100.0);
        QCOMPARE(t->downforce, 1000.0);
        QVERIFY(qAbs(t->coefficients->cd - 0.416667) < 1e-6);
        QVERIFY(qAbs(t->coefficients->cl + 4.166667) < 1e-6);
        QVERIFY(qAbs(*t->cla + 4.166667) < 1e-6);
        QVERIFY(qAbs(*t->cda - 0.416667) < 1e-6);
        QVERIFY(qAbs(t->pitch - 1030) < 1e-3);
        QVERIFY(qAbs(*t->coefficients->cmPitch - 1.716667) < 1e-6);
        QCOMPARE(*t->downforceShare, 1.0);
        QVERIFY(qAbs(*a.liftToDrag - 10) < 1e-9);
        QCOMPARE(a.balance.status, QString("defined"));
        QVERIFY(qAbs(*a.balance.rearLoad - 412) < 1e-3);
        QVERIFY(qAbs(*a.balance.frontLoad - 588) < 1e-3);
        QVERIFY(qAbs(*a.balance.frontPercent - 58.8) < 1e-4);
        QVERIFY(qAbs(*a.balance.copFromFrontAxle - 1.03) < 1e-6);
        QVERIFY(a.approximate);
        QCOMPARE(a.provenance, QString("APPROXIMATE RESULT: test"));
    }

    // The CLI's group case: pressure/viscous split, shares, moments moved to
    // the moment centre, and a balance that does not depend on it.
    void groupsSplitAndMoveMoments()
    {
        const QJsonArray zones{zone("front-wing", {2, 0, -300}, {1, 0, 0}, {-0.5f, 0, 0.2f}),
                               zone("rear-wing", {5, 0, -400}, {0.5f, 0, 0}, {2, 0, 0.8f}),
                               zone("floor", {0, 0, -200}, {3, 0, 0}, {1, 0, 0})};
        const QJsonObject spec{{"moment_centre", v3(1, 0, 0)},
                               {"groups", QJsonArray{QJsonObject{{"name", "Wings"}, {"zones", QJsonArray{"front-wing", "rear-wing"}}}}},
                               {"vehicle", QJsonObject{{"front_axle", v3(0, 0, 0)}, {"wheelbase", 2.5}}}};
        const auto r = aero(loads(zones), spec);
        QVERIFY2(r.ok(), qPrintable(r ? QString() : r.error()));
        const AeroResult& a = r.value();
        QStringList names;
        for (const AeroRow& x : a.rows) names << x.name;
        QCOMPARE(names, (QStringList{"front-wing", "rear-wing", "floor", "Wings", "Total"}));
        QCOMPARE(a.ungrouped, QStringList{"floor"});
        const AeroRow* w = row(a, "Wings");
        QVERIFY(w && w->kind == "group" && w->coefficients);
        QVERIFY(qAbs(w->coefficients->cdPressure - 7. / 240) < 1e-9);
        QVERIFY(qAbs(w->coefficients->cdViscous - 1.5 / 240) < 1e-9);
        QVERIFY(qAbs(w->viscous.drag - 1.5) < 1e-6);
        QVERIFY(qAbs(*w->downforceShare - 7. / 9) < 1e-6);
        QVERIFY(qAbs(row(a, "front-wing")->pitch + 449.4) < 1e-3);
        QVERIFY(qAbs(w->pitch + 45) < 1e-3);
        QVERIFY(qAbs(*a.total()->coefficients->cmPitch + 45. / 600) < 1e-6);
        QVERIFY(qAbs(*a.liftToDrag - 900 / 11.5) < 1e-6);
        QVERIFY(qAbs(*a.balance.rearLoad - 342) < 1e-3);
        QVERIFY(qAbs(*a.balance.frontPercent - 62) < 1e-4);
        QVERIFY(qAbs(*a.balance.copFromFrontAxle - 0.95) < 1e-6);
    }

    // Missing references give nulls the app shows as "—", not zeros.
    void missingReferencesAreNotZero()
    {
        // No area: CL/CD are unknown, ClA/CdA are not (they need only q).
        auto r = aero(loads({zone("car", {100, 0, -1000}, {}, {1, 0, 0})},
                            {{"area", QJsonValue()}, {"density", 1.2}, {"velocity", 20.0}, {"length", 2.5}, {"pressure", 0.0}}),
                      {});
        QVERIFY(r.ok());
        QVERIFY(!r.value().total()->coefficients);
        QVERIFY(qAbs(*r.value().total()->cla + 1000. / 240) < 1e-9);
        QCOMPARE(r.value().balance.status, QString("needs_vehicle"));
        QVERIFY(!r.value().balance.frontPercent);
        // No density: nothing is a coefficient.
        r = aero(loads({zone("car", {100, 0, -1000}, {}, {1, 0, 0})},
                       {{"area", 1.0}, {"density", QJsonValue()}, {"velocity", 20.0}, {"length", 2.5}, {"pressure", 0.0}}),
                 {});
        QVERIFY(r.ok());
        QVERIFY(!r.value().total()->coefficients && !r.value().total()->cla && !r.value().dynamicPressure);
        // Lift instead of downforce: no split.
        r = aero(loads({zone("car", {10, 0, 100}, {}, {1, 0, 0})}),
                 {{"vehicle", QJsonObject{{"front_axle", v3(0, 0, 0)}, {"wheelbase", 2.0}}}});
        QVERIFY(r.ok());
        QCOMPARE(r.value().balance.status, QString("no_downforce"));
        QVERIFY(!r.value().balance.frontPercent && !r.value().balance.copFromFrontAxle);
    }

    void invalidSpecIsAUsageError()
    {
        int code = 0;
        const auto r = aero(loads({zone("a", {0, 0, -1}, {}, {})}),
                            {{"groups", QJsonArray{QJsonObject{{"name", "x"}, {"zones", QJsonArray{"missing"}}}}}}, &code);
        QVERIFY(!r.ok());
        QCOMPARE(code, 2);
    }

    // What the app writes is a spec the CLI accepts: empty groups are dropped
    // and a half-entered vehicle is left out.
    void specRoundTripsAndOmitsIncompleteParts()
    {
        AeroSpec s;
        s.carZones = {"a", "b"};
        s.groups = {{"Wings", {"a"}}, {"Empty", {}}};
        s.frontAxle = QVector3D(-0.78f, 0, 0);
        QJsonObject o = s.toJson();
        QVERIFY(!o.contains("vehicle"));
        QCOMPARE(o.value("groups").toArray().size(), 1);
        s.wheelbase = 1.58;
        o = s.toJson();
        const AeroSpec back = AeroSpec::fromJson(o);
        QVERIFY(back.hasVehicle());
        QCOMPARE(*back.wheelbase, 1.58);
        QCOMPARE(back.carZones, s.carZones);
        for (const QString& key : o.keys())
            QVERIFY2(QStringList({"v", "axes", "moment_centre", "car_zones", "groups", "vehicle"}).contains(key), qPrintable(key));
    }

    void suggestedGroupsFollowZoneNames()
    {
        const auto groups = suggestAeroGroups({"fw-mainplane", "front_wing_flap", "rw-main", "RearWing_EP", "undertray", "diffuser",
                                               "tyre-fl", "wheel_rr", "monocoque", "driver"});
        QMap<QString, QStringList> by;
        for (const auto& g : groups) by[g.name] = g.zones;
        QCOMPARE(by.value("Front wing"), (QStringList{"fw-mainplane", "front_wing_flap"}));
        QCOMPARE(by.value("Rear wing"), (QStringList{"rw-main", "RearWing_EP"}));
        QCOMPARE(by.value("Undertray / diffuser"), (QStringList{"undertray", "diffuser"}));
        QCOMPARE(by.value("Wheels"), (QStringList{"tyre-fl", "wheel_rr"}));
        QCOMPARE(by.value("Body"), (QStringList{"monocoque", "driver"}));
    }

    void convergenceEventsParse()
    {
        const auto setup = parseConvergenceSetup(QJsonDocument::fromJson(R"({"enabled":true,"condition":"any",
            "residuals":{"continuity":0.001},
            "reports":[{"report":"drag","stop_criterion":0.001,"previous_values":5,"ignore_before":100}]})").object());
        QVERIFY(setup.enabled && setup.condition == "any" && setup.hasConditions());
        QCOMPARE(setup.residuals.value("continuity"), 0.001);
        QCOMPARE(int(setup.reports.size()), 1);
        QCOMPARE(setup.reports[0].previousValues, 5);
        QCOMPARE(describeCondition("report:drag", setup), QString("drag changed less than 0.1 % between each of its last 5 values"));
        const auto state = parseConvergenceState(QJsonDocument::fromJson(R"({"met":false,"conditions":[
            {"kind":"residuals","met":false,"equations":{"continuity":{"value":0.0021,"criterion":0.001,"met":false}}},
            {"kind":"report","report":"drag","value":null,"criterion":0.001,"met":false}]})").object());
        QCOMPARE(int(state.conditions.size()), 2);
        QCOMPARE(state.conditions[0].id(), QString("residuals"));
        QCOMPARE(state.conditions[0].equations.value("continuity").value, 0.0021);
        QCOMPARE(state.conditions[1].id(), QString("report:drag"));
        QVERIFY(!state.conditions[1].value);  // too few values yet
        // An older CLI sends no convergence object: checking is off.
        QVERIFY(!parseConvergenceSetup({}).enabled);
    }

    // Convergence controls say when to stop, so they travel as run controls
    // and never make a run MODIFIED.
    void monitorsAreRunControls()
    {
        QVERIFY(SettingsModel::isRunControl("monitors/check-convergence"));
        QVERIFY(SettingsModel::isRunControl("monitors/report/drag/active"));
        QVERIFY(SettingsModel::isRunControl("run/iterations"));
        QVERIFY(!SettingsModel::isRunControl("reference/velocity"));
        const QJsonObject setup = SettingsModel::setupJsonFor({{"reference/velocity", 20.0}},
                                                              {{"run/iterations", 600}, {"monitors/check-convergence", true}});
        const QJsonObject overrides = setup.value("overrides").toObject();
        QVERIFY(overrides.contains("reference/velocity") && overrides.contains("monitors/check-convergence"));
        QVERIFY(!overrides.contains("run/iterations"));  // passed as --iterations
    }

    // The synthetic car with its criteria on stops before the iteration
    // limit; the run says which conditions stopped it and stays unmodified.
    void runStopsOnItsCriteria()
    {
        QFile car(m_dir.filePath("fsae_car.cas.h5"));
        QVERIFY(car.open(QIODevice::WriteOnly));
        car.close();
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        Session s(cli, m_dir.filePath("data"));
        QSignalSpy opened(&s, &Session::caseOpened);
        s.openCase(car.fileName());
        QVERIFY(opened.wait(20000));
        for (const QString& a : s.caseSummary().approximations()) s.setApproximationApproved(a, true);
        for (const char* p : {"monitors/check-convergence", "monitors/report/downforce/active", "monitors/report/drag/active"})
            QCOMPARE(s.settings()->setValue(p, true), QString());
        QCOMPARE(s.settings()->setValue("run/iterations", 2000), QString());
        QSignalSpy finished(&s, &Session::runFinished);
        const QString id = s.run();
        QVERIFY(!id.isEmpty());
        QVERIFY(finished.wait(120000));
        const RunRecord* r = s.history().find(id);
        QVERIFY(r);
        QVERIFY(r->converged());
        QVERIFY(r->iterationsDone < 2000);
        QVERIFY(!r->modified());
        QVERIFY(r->convergedBy.contains("report:downforce") && r->convergedBy.contains("report:drag"));
        QCOMPARE(runStatusText(*r), QString("Converged"));

        // Its aero numbers come from loads.json, without a re-solve.
        std::optional<Result<AeroResult>> aero;
        s.aero(id, [&](Result<AeroResult> a) { aero = std::move(a); });
        QTRY_VERIFY_WITH_TIMEOUT(aero.has_value(), 20000);
        QVERIFY2(aero->ok(), qPrintable(*aero ? QString() : aero->error()));
        QVERIFY(aero->value().total()->downforce > 0);
        QVERIFY(aero->value().provenance.contains("SYNTHETIC"));
    }
};

QTEST_MAIN(TestAero)
#include "tst_aero.moc"
