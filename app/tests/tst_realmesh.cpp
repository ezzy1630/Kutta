// New Case from Geometry against the real `kutta mesh`: a small closed body
// in millimetres, meshed coarse through MeshWorkflow, its report read
// strictly, and the case set up from the mesh with kutta new-case. Runs with
// KUTTA_REAL_CLI, or the in-tree CLI when the app is built with the solver;
// skipped otherwise (CI's app jobs have no solver).

#include "core/MeshWorkflow.h"
#include "core/Session.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

class TestRealMesh : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;
    QString m_cli;

    // An ASCII STL box (mm) standing `lift` above the ground plane z = 0.
    QString writeBox(const QString& name, QVector3D lo, QVector3D hi)
    {
        const QString path = m_dir.filePath(name + ".stl");
        QFile f(path);
        if (!f.open(QIODevice::WriteOnly)) return {};
        auto c = [&](int i) { return QVector3D(i & 1 ? hi.x() : lo.x(), i & 2 ? hi.y() : lo.y(), i & 4 ? hi.z() : lo.z()); };
        // Outward-facing quads as corner indices.
        const int quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
        QByteArray s = "solid " + name.toUtf8() + "\n";
        for (const auto& q : quads)
            for (const auto& tri : {std::array<int, 3>{q[0], q[1], q[2]}, std::array<int, 3>{q[0], q[2], q[3]}}) {
                s += " facet normal 0 0 0\n  outer loop\n";
                for (int i : tri) s += QStringLiteral("   vertex %1 %2 %3\n").arg(c(i).x()).arg(c(i).y()).arg(c(i).z()).toUtf8();
                s += "  endloop\n endfacet\n";
            }
        s += "endsolid " + name.toUtf8() + "\n";
        f.write(s);
        return path;
    }

private slots:
    void initTestCase()
    {
        m_cli = qEnvironmentVariable("KUTTA_REAL_CLI");
#ifdef KUTTA_DEV_REAL_CLI
        if (m_cli.isEmpty()) m_cli = QStringLiteral(KUTTA_DEV_REAL_CLI);
#endif
        if (m_cli.isEmpty()) QSKIP("set KUTTA_REAL_CLI, or build the app with the solver, to mesh with the real kutta mesh");
        QVERIFY(QFileInfo(m_cli).isExecutable());
        QCoreApplication::setOrganizationName("KuttaRealMeshTest");
        QSettings::setDefaultFormat(QSettings::IniFormat);
        QSettings::setPath(QSettings::IniFormat, QSettings::UserScope, m_dir.filePath("settings"));
    }

    void meshesAndSetsUpTheCase()
    {
        CliInfo cli;
        cli.program = m_cli;
        Session s(cli, m_dir.filePath("data"));
        MeshWorkflow* mw = s.meshing();
        // 1 m long, 50 mm above the ground: millimetres, read as such from the size.
        const QString stl = writeBox("body", QVector3D(0, -200, 50), QVector3D(1000, 200, 350));
        QVERIFY(!stl.isEmpty());
        mw->importGeometry({stl}, m_dir.filePath("projects/box"));
        QVERIFY(QTest::qWaitFor([mw] { return mw->hasProject() && !mw->loading(); }, 30000));
        QCOMPARE(mw->project().setup.units, LengthUnit::Millimetre);
        MeshSetup m = mw->project().setup;
        m.baseSize = 0.4;
        m.maxSize = 0.1;
        m.minSize = 0.05;
        m.proximityCells = 2;
        m.boxes.clear();
        m.layers = 2;
        m.firstHeight = 0.005;
        m.tunnel = {1, 2, 1, 1, true, {}};
        mw->setSetup(m);
        mw->setSpeed(30);
        QVERIFY2(mw->problems().isEmpty(), qPrintable(mw->problems().join('\n')));
        QVERIFY(mw->generate());
        QVERIFY(QTest::qWaitFor([mw] { return !mw->meshing(); }, 120000));
        QVERIFY2(mw->outcome() == MeshWorkflow::Outcome::Succeeded, qPrintable(mw->failure() + "\n" + mw->advice() + "\n" + mw->reportError()));
        QVERIFY(!mw->meshSynthetic());
        const MeshReport& r = *mw->report();
        QVERIFY(r.cells > 1000 && r.faces > r.cells);
        QStringList zones;
        for (const auto& z : r.zones) zones << z.name;
        for (const auto& f : tunnelFaces(m)) QVERIFY2(zones.contains(f.name), qPrintable(f.name + " not in " + zones.join(',')));
        QVERIFY(zones.contains("body"));
        // The figures the Mesh Report shows, from the real report.
        QVERIFY(r.quality && r.quality->orthogonalMin && *r.quality->orthogonalMin >= 0.15 - 1e-9);
        QVERIFY(!r.quality->orthogonalHistogram.empty());
        QVERIFY(r.layers && r.layers->coverageFull && r.layers->coverageNone);
        QCOMPARE(*r.layers->requested, qint64(2));
        QVERIFY(r.snap && r.snap->distanceMax);
        QVERIFY(r.yPlus && r.yPlus->p50 && r.yPlus->from30to300);
        QCOMPARE(*r.yPlus->speed, 30.0);  // the project's speed reached the mesher
        QVERIFY(r.gaps);
        QVERIFY(r.missingSurfaceZones && r.missingSurfaceZones->isEmpty());
        // Lengths in metres: the tunnel the app drew is the one the mesher built.
        const Bounds t = tunnelBounds(m, mw->geometryBounds());
        const QJsonArray b = r.raw.value("bounds").toArray();
        for (int a = 0; a < 3; ++a) {
            QVERIFY2(qAbs(b[0].toArray()[a].toDouble() - t.min[a]) < 1e-3, qPrintable(QStringLiteral("min %1").arg(a)));
            QVERIFY2(qAbs(b[1].toArray()[a].toDouble() - t.max[a]) < 1e-3, qPrintable(QStringLiteral("max %1").arg(a)));
        }

        // The case from the mesh: the tunnel's roles from its names, the speed carried over.
        s.setUpCaseFromMesh();
        QVERIFY(QTest::qWaitFor([&] { return s.newCaseDraft() && !s.newCaseBusy(); }, 60000));
        const NewCaseDraft& d = *s.newCaseDraft();
        QCOMPARE(d.speed, 30.0);
        QMap<QString, QString> roles;
        for (const NewCaseZone* z : d.boundaries()) roles[z->name] = z->role;
        QCOMPARE(roles.value("inlet"), QString("velocity-inlet"));
        QCOMPARE(roles.value("outlet"), QString("pressure-outlet"));
        QCOMPARE(roles.value("ground"), QString("moving-ground"));
        QCOMPARE(roles.value("side-low"), QString("symmetry"));
        QCOMPARE(roles.value("side-high"), QString("symmetry"));
        QCOMPARE(roles.value("top"), QString("symmetry"));
        QCOMPARE(roles.value("body"), QString("wall"));
        const QString target = s.newCaseTarget();
        s.createNewCase(target);
        QVERIFY(QTest::qWaitFor([&] { return s.hasCase() && s.casePath() == target; }, 60000));
        QFile kase(target);
        QVERIFY(kase.open(QIODevice::ReadOnly));
        const QJsonObject k = QJsonDocument::fromJson(kase.readAll()).object();
        QCOMPARE(k.value("methods").toObject().value("warped_face_gradient_correction").toBool(true), false);
        QCOMPARE(k.value("boundaries").toObject().value("inlet").toObject().value("speed").toDouble(), 30.0);
    }
};

QTEST_GUILESS_MAIN(TestRealMesh)
#include "tst_realmesh.moc"
