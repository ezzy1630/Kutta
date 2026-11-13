// Geometry to case through the fixture CLI's `mesh` (the real process
// boundary): import STLs, generate, cancel, refuse with advice, set up the
// case from the mesh, and duplicate it as a design variant.

#include "core/MeshWorkflow.h"
#include "core/Session.h"

#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QProcess>
#include <QSignalSpy>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

class TestMeshFlow : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;
    QString m_stl;

    std::unique_ptr<Session> session()
    {
        CliInfo cli;
        cli.program = QStringLiteral(KUTTA_FIXTURE_CLI);
        cli.fixture = true;
        return std::make_unique<Session>(cli, m_dir.filePath("data"));
    }
    QString stl(const QString& name) const { return QDir(m_stl).filePath(name); }

    // Imports into a fresh project and waits for the surfaces.
    bool import(MeshWorkflow* mw, const QStringList& files, const QString& name)
    {
        mw->importGeometry(files, m_dir.filePath("projects/" + name));
        return QTest::qWaitFor([mw] { return mw->hasProject() && !mw->loading(); }, 30000) && mw->error().isEmpty();
    }
    bool generated(MeshWorkflow* mw)
    {
        if (!mw->generate()) return false;
        return QTest::qWaitFor([mw] { return !mw->meshing(); }, 60000);
    }

private slots:
    void initTestCase()
    {
        QCoreApplication::setOrganizationName("KuttaTest");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_dir.filePath("settings"));
        qputenv("KUTTA_FIXTURE_SPEED", "400");
        m_stl = m_dir.filePath("stl");
        QProcess p;
        p.start(QStringLiteral(KUTTA_FIXTURE_CLI), {"sample-stl", m_stl});
        QVERIFY(p.waitForFinished(60000));
        QCOMPARE(p.exitCode(), 0);
        QVERIFY(QFileInfo::exists(stl("rear-wing-v2.stl")));
    }

    // Surfaces, tunnel and refinement regions are in the scene before any CLI runs;
    // Generate writes the mesh and its report, and a settings change makes it stale.
    void importsAndGenerates()
    {
        auto s = session();
        MeshWorkflow* mw = s->meshing();
        QVERIFY(import(mw, {stl("fsae-car.stl"), stl("rear-wing.stl")}, "car"));
        QCOMPARE(mw->project().setup.units, LengthUnit::Millimetre);  // a 3,960-"metre" car
        QCOMPARE(int(mw->project().solids.size()), 6);
        int wheels = 0;
        for (const auto& solid : mw->project().solids) wheels += solid.wheel;
        QCOMPARE(wheels, 2);  // tyre-front and tyre-rear, by name
        QStringList ids;
        for (const auto& it : mw->scene()->items) ids << it.id;
        QVERIFY(ids.contains("zone:rear-wing") && ids.contains("zone:ground") && ids.contains("outline:tunnel"));
        QVERIFY2(ids.contains("box:body") && ids.contains("box:wake"), qPrintable(ids.join(' ')));
        QVERIFY(mw->problems().isEmpty());

        QVERIFY(generated(mw));
        QCOMPARE(mw->outcome(), MeshWorkflow::Outcome::Succeeded);
        QVERIFY(mw->hasMesh() && !mw->meshStale());
        QVERIFY(QFileInfo::exists(mw->project().dir + "/mesh.json"));
        QVERIFY(mw->report() && mw->report()->synthetic && mw->report()->cells > 0);
        QVERIFY(mw->report()->metric("orthogonal_quality"));
        QCOMPARE(mw->progress().fraction, 1.0);

        MeshSetup changed = mw->project().setup;
        changed.layers = 8;
        mw->setSetup(changed);
        QVERIFY(mw->meshStale());

        // The project reopens with its names, roles and mesh.
        const QString dir = mw->project().dir;
        mw->closeProject();
        mw->openProject(dir);
        QVERIFY(QTest::qWaitFor([mw] { return mw->hasProject() && !mw->loading(); }, 30000));
        QCOMPARE(mw->project().setup.layers, 8);
        QVERIFY(mw->hasMesh() && mw->meshStale());
    }

    // Cancel stops the mesher and keeps nothing half-written.
    void cancelsAMeshInProgress()
    {
        auto s = session();
        MeshWorkflow* mw = s->meshing();
        QVERIFY(import(mw, {stl("fsae-car.stl"), stl("rear-wing.stl")}, "cancel"));
        qputenv("KUTTA_FIXTURE_MESH_HOLD", "octree");
        QVERIFY(mw->generate());
        QVERIFY(QTest::qWaitFor([mw] { return mw->progress().stage == "octree"; }, 30000));
        mw->cancel();
        QVERIFY(QTest::qWaitFor([mw] { return !mw->meshing(); }, 15000));
        qunsetenv("KUTTA_FIXTURE_MESH_HOLD");
        QCOMPARE(mw->outcome(), MeshWorkflow::Outcome::Cancelled);
        QVERIFY(!mw->hasMesh());
    }

    // A mesher error arrives with its advice, in place, not as a dialog.
    void showsMesherErrorsWithAdvice()
    {
        auto s = session();
        MeshWorkflow* mw = s->meshing();
        QVERIFY(import(mw, {stl("rear-wing.stl")}, "wing-only"));
        QVERIFY(generated(mw));
        QCOMPARE(mw->outcome(), MeshWorkflow::Outcome::Failed);
        QVERIFY(!mw->failure().isEmpty());
        QVERIFY(!mw->advice().isEmpty());
        QVERIFY(!mw->hasMesh());
    }

    // The mesh becomes a case (wheels flagged from the geometry), and a design
    // variant with a swapped rear wing reuses that case's settings.
    void setsUpTheCaseAndAVariant()
    {
        auto s = session();
        MeshWorkflow* mw = s->meshing();
        QVERIFY(import(mw, {stl("fsae-car.stl"), stl("rear-wing.stl")}, "base"));
        QVERIFY(generated(mw));
        QCOMPARE(mw->outcome(), MeshWorkflow::Outcome::Succeeded);

        s->setUpCaseFromMesh();
        QVERIFY(QTest::qWaitFor([&] { return s->newCaseDraft() && !s->newCaseBusy(); }, 30000));
        const QString target = s->newCaseTarget();
        QCOMPARE(target, mw->project().defaultCaseFile());
        for (const NewCaseZone* z : s->newCaseDraft()->boundaries())
            if (z->wheel && !z->wheelConfirmed) s->confirmNewCaseWheel(z->name, true);
        QVERIFY2(s->newCaseDraft()->wheelError().isEmpty(), qPrintable(s->newCaseDraft()->wheelError()));
        s->createNewCase(target);
        QVERIFY(QTest::qWaitFor([&] { return s->hasCase() && s->casePath() == target; }, 30000));
        QCOMPARE(mw->project().caseFile, target);
        // The case's mesh keeps the geometry view's colours.
        std::shared_ptr<SceneData> meshScene;
        s->loadScene(s->meshSceneKey(), [&](const QString&, Result<std::shared_ptr<SceneData>> r) { meshScene = r ? r.value() : nullptr; });
        QVERIFY(QTest::qWaitFor([&] { return meshScene != nullptr; }, 30000));
        QCOMPARE(meshScene->item("zone:rear-wing")->color, mw->scene()->item("zone:rear-wing")->color);
        QCOMPARE(meshScene->item("zone:tyre-rear")->color, mw->scene()->item("zone:tyre-rear")->color);
        QCOMPARE(MeshWorkflow::projectForCase(target), mw->project().projectFile());

        const QString variant = mw->duplicateAsVariant();
        QVERIFY2(!variant.isEmpty(), qPrintable(mw->error()));
        QCOMPARE(mw->project().baseCase, target);
        QCOMPARE(int(mw->family().size()), 2);
        int wingFile = -1;
        for (int i = 0; i < mw->project().files.size(); ++i)
            if (mw->project().files[i].endsWith("rear-wing.stl")) wingFile = i;
        QVERIFY(wingFile >= 0);
        mw->replaceGeometry(wingFile, stl("rear-wing-v2.stl"));
        QVERIFY(QTest::qWaitFor([&] { return !mw->loading() && mw->project().files.value(wingFile).endsWith("rear-wing-v2.stl"); }, 30000));
        QVERIFY(std::any_of(mw->project().solids.begin(), mw->project().solids.end(), [](const GeometrySolid& g) { return g.name == "rear-wing"; }));
        QVERIFY(generated(mw));
        QCOMPARE(mw->outcome(), MeshWorkflow::Outcome::Succeeded);

        QSignalSpy opened(s.get(), &Session::caseOpened);
        s->setUpCaseFromMesh();
        QVERIFY(opened.wait(30000));
        const QString variantCase = mw->project().defaultCaseFile();
        QCOMPARE(s->casePath(), variantCase);
        QFile in(variantCase), base(target);
        QVERIFY(in.open(QIODevice::ReadOnly) && base.open(QIODevice::ReadOnly));
        const QJsonObject v = QJsonDocument::fromJson(in.readAll()).object();
        const QJsonObject b = QJsonDocument::fromJson(base.readAll()).object();
        // Same settings, its own mesh, and every wall of the new mesh has a boundary.
        QCOMPARE(v.value("setup").toObject().value("models"), b.value("setup").toObject().value("models"));
        QVERIFY(QFileInfo(QDir(QFileInfo(variantCase).absolutePath()).filePath(v.value("mesh").toString())).canonicalFilePath()
                == QFileInfo(mw->meshFile()).canonicalFilePath());
        QVERIFY(!QDir::isAbsolutePath(v.value("mesh").toString()));
        const QString boundaries = QString::fromUtf8(QJsonDocument(v).toJson());
        for (const auto& z : mw->report()->zones)
            QVERIFY2(boundaries.contains('"' + z.name + '"'), qPrintable(z.name));
    }
};

QTEST_MAIN(TestMeshFlow)
#include "tst_meshflow.moc"
