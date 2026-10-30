// Overrides, validation and undo for the settings model.

#include "TestSupport.h"
#include "core/SettingsModel.h"

#include <QSignalSpy>
#include <QtTest>

using namespace kapp;

class TestSettings : public QObject {
    Q_OBJECT
    SettingsModel m;
private slots:
    void init()
    {
        m.reset(parseCaseSummary(fixtureJson("cases/mixing_elbow/inspect.json")).value().settings);
        m.undoStack()->clear();
    }

    void editMarksModifiedAndWritesSetup()
    {
        QSignalSpy changed(&m, &SettingsModel::valueChanged);
        QCOMPARE(m.setValue("boundary/cold-inlet/velocity-magnitude", 1.4), QString());
        QVERIFY(m.isModified("boundary/cold-inlet/velocity-magnitude"));
        QCOMPARE(m.value("boundary/cold-inlet/velocity-magnitude").toDouble(), 1.4);
        QCOMPARE(changed.count(), 1);
        const QJsonObject setup = m.setupJson();
        QCOMPARE(setup.value("v").toInt(), 1);
        QCOMPARE(setup.value("overrides").toObject().value("boundary/cold-inlet/velocity-magnitude").toDouble(), 1.4);
    }

    void settingTheCaseValueClearsTheOverride()
    {
        m.setValue("boundary/cold-inlet/velocity-magnitude", 1.4);
        m.setValue("boundary/cold-inlet/velocity-magnitude", 1.0);
        QVERIFY(!m.isModified("boundary/cold-inlet/velocity-magnitude"));
        QCOMPARE(m.modifiedCount(), 0);
    }

    void rejectsInvalidValues()
    {
        QVERIFY(!m.setValue("boundary/cold-inlet/velocity-magnitude", "fast").isEmpty());
        QVERIFY(!m.setValue("methods/momentum-scheme", "third-order").isEmpty());
        QVERIFY(!m.setValue("run/iterations", 10.5).isEmpty());
        // Read-only settings explain themselves.
        QCOMPARE(m.setValue("reference/length", 2.0), QString("not yet configurable in Kutta"));
        QVERIFY(!m.setValue("no/such/setting", 1).isEmpty());
        QCOMPARE(m.modifiedCount(), 0);
        QCOMPARE(m.undoStack()->count(), 0);
    }

    void undoRedoAndResetAll()
    {
        m.setValue("controls/mom-relaxation", 0.3);
        m.setValue("methods/momentum-scheme", "first-order-upwind");
        QCOMPARE(m.modifiedCount(), 2);
        m.undoStack()->undo();
        QCOMPARE(m.modifiedCount(), 1);
        m.undoStack()->redo();
        QCOMPARE(m.modifiedCount(), 2);
        m.resetAll();
        QCOMPARE(m.modifiedCount(), 0);
        m.undoStack()->undo();  // reset all is one step
        QCOMPARE(m.modifiedCount(), 2);
    }

    void consecutiveEditsOfOneFieldMerge()
    {
        m.setValue("controls/mom-relaxation", 0.3);
        m.setValue("controls/mom-relaxation", 0.35);
        QCOMPARE(m.undoStack()->count(), 1);
        m.undoStack()->undo();
        QVERIFY(!m.isModified("controls/mom-relaxation"));
    }

    void applyOverridesRestoresARun()
    {
        const QJsonObject run{{"controls/mom-relaxation", 0.4}, {"run/iterations", 200}};
        QCOMPARE(m.applyOverrides(run, "restore"), QString());
        QCOMPARE(m.overrides(), run);
        QVERIFY(!m.applyOverrides(QJsonObject{{"reference/length", 2.0}}, "bad").isEmpty());
    }

    // A method choice decides which of its fields the case uses, live, with
    // Fluent's words for the codes; unused fields never reach the solver.
    void methodChoiceShowsOnlyItsFields()
    {
        const QString spec = "boundary/in/turbulence-specification";
        auto setting = [](QJsonObject o) {
            o.insert("group", "Boundary Conditions");
            o.insert("unit", "");
            if (!o.contains("editable")) o.insert("editable", true);
            return o;
        };
        const QJsonArray settings{
            setting({{"path", spec}, {"label", "Turbulence specification"}, {"kind", "enum"}, {"value", 3},
                     {"options", QJsonArray{1, 3}}, {"option_labels", QJsonArray{"Intensity and Length Scale", "Intensity and Hydraulic Diameter"}},
                     {"value_label", "Intensity and Hydraulic Diameter"}}),
            setting({{"path", "boundary/in/turb-hydraulic-diam"}, {"label", "Hydraulic diameter"}, {"kind", "number"}, {"value", 0.1},
                     {"depends_on", spec}, {"active_when", QJsonArray{3}}, {"active", true}}),
            setting({{"path", "boundary/in/turb-length-scale"}, {"label", "Turbulent length scale"}, {"kind", "number"}, {"value", 0.01},
                     {"depends_on", spec}, {"active_when", QJsonArray{1}}, {"active", false}}),
            setting({{"path", "boundary/in/k"}, {"label", "Turbulent kinetic energy"}, {"kind", "number"}, {"value", 1},
                     {"depends_on", spec}, {"active_when", QJsonArray{0}}, {"active", false}, {"editable", false}, {"reason", "not yet configurable in Kutta"}}),
            setting({{"path", "run/monitor-interval"}, {"label", "Monitor interval"}, {"kind", "integer"}, {"value", 20}, {"min", 1}}),
            setting({{"path", "run/iterations"}, {"label", "Number of iterations"}, {"kind", "integer"}, {"value", 1000}, {"min", 1}}),
        };
        auto c = parseCaseSummary({{"v", 1}, {"support", QJsonObject{{"status", "supported"}}}, {"settings", settings}});
        QVERIFY2(c.ok(), qPrintable(c ? QString() : c.error()));
        SettingsModel s;
        s.reset(c.value().settings);
        const Setting* sp = s.find(spec);
        QVERIFY(sp);
        QCOMPARE(sp->labelFor(3), QString("Intensity and Hydraulic Diameter"));
        QCOMPARE(sp->labelFor(1), QString("Intensity and Length Scale"));
        QVERIFY(s.isActive("boundary/in/turb-hydraulic-diam"));
        QVERIFY(!s.isActive("boundary/in/turb-length-scale"));
        QVERIFY(!s.isActive("boundary/in/k"));
        QCOMPARE(s.dependents(spec).size(), 3);

        // Switching the method re-evaluates its fields and announces them.
        QSignalSpy changed(&s, &SettingsModel::valueChanged);
        QCOMPARE(s.setValue(spec, 1), QString());
        QVERIFY(changed.count() >= 4);  // the choice and its three fields
        QVERIFY(!s.isActive("boundary/in/turb-hydraulic-diam"));
        QVERIFY(s.isActive("boundary/in/turb-length-scale"));
        QVERIFY(!s.setValue(spec, 2).isEmpty());  // not one of the case's options

        // A value edited under the other method stays in the model (as in
        // Fluent) but is not sent while that method is not selected.
        QCOMPARE(s.setValue("boundary/in/turb-length-scale", 0.02), QString());
        QVERIFY(s.caseOverrides().contains("boundary/in/turb-length-scale"));
        s.undoStack()->undo();
        s.undoStack()->undo();  // back to Intensity and Hydraulic Diameter
        QCOMPARE(s.value(spec).toInt(), 3);
        s.setValue("boundary/in/turb-length-scale", 0.02);
        QVERIFY(s.isModified("boundary/in/turb-length-scale"));
        QVERIFY(!s.caseOverrides().contains("boundary/in/turb-length-scale"));

        // Run controls never make a result MODIFIED; the monitor interval
        // still reaches kutta through the setup file.
        s.resetAll();
        QCOMPARE(s.setValue("run/monitor-interval", 1), QString());
        QCOMPARE(s.setValue("run/iterations", 300), QString());
        QCOMPARE(s.caseModifiedCount(), 0);
        QVERIFY(s.needsSetupFile());
        const QJsonObject o = s.setupJson().value("overrides").toObject();
        QCOMPARE(o.value("run/monitor-interval").toInt(), 1);
        QVERIFY(!o.contains("run/iterations"));  // passed as --iterations
    }
};

QTEST_GUILESS_MAIN(TestSettings)
#include "tst_settings.moc"
