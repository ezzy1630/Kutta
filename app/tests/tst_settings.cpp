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
};

QTEST_GUILESS_MAIN(TestSettings)
#include "tst_settings.moc"
