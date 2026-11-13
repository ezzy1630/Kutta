// New Case from Geometry without a CLI: STL reading and writing, mesh.json,
// the planning estimate, and the `kutta mesh` adapter's parsers.

#include "core/MeshCli.h"
#include "core/MeshSetup.h"
#include "core/Stl.h"

#include <QDataStream>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

namespace {

void writeText(const QString& path, const QByteArray& text)
{
    QFile f(path);
    QVERIFY(f.open(QIODevice::WriteOnly));
    f.write(text);
}

// Little-endian binary STL: 80-byte header, count, 50 bytes per triangle.
void writeBinary(const QString& path, const std::vector<std::array<float, 9>>& tris, bool truncate = false)
{
    QFile f(path);
    QVERIFY(f.open(QIODevice::WriteOnly));
    QDataStream s(&f);
    s.setByteOrder(QDataStream::LittleEndian);
    s.setFloatingPointPrecision(QDataStream::SinglePrecision);
    QByteArray header(80, ' ');
    header.replace(0, 5, "solid");  // some exporters start binary files with "solid"
    s.writeRawData(header.constData(), 80);
    s << quint32(tris.size());
    for (std::size_t t = 0; t < tris.size() - (truncate ? 1 : 0); ++t) {
        for (int k = 0; k < 3; ++k) s << 0.f;
        for (float x : tris[t]) s << x;
        s << quint16(0);
    }
}

} // namespace

class TestGeometry : public QObject {
    Q_OBJECT
    QTemporaryDir m_dir;

private slots:
    // Named solids become zones; shared corners are welded; mm become m.
    void readsAsciiSolidsInMillimetres()
    {
        const QString p = m_dir.filePath("two.stl");
        writeText(p, "solid wing\n"
                     " facet normal 0 0 1\n  outer loop\n   vertex 0 0 0\n   vertex 1000 0 0\n   vertex 1000 1000 0\n  endloop\n endfacet\n"
                     " facet normal 0 0 1\n  outer loop\n   vertex 0 0 0\n   vertex 1000 1000 0\n   vertex 0 1000 0\n  endloop\n endfacet\n"
                     " facet normal 0 0 1\n  outer loop\n   vertex 5 5 5\n   vertex 5 5 5\n   vertex 6 6 6\n  endloop\n endfacet\n"
                     "endsolid wing\n"
                     "solid nose cone\n"
                     " facet normal 0 0 1\n  outer loop\n   vertex 0 0 1000\n   vertex 1000 0 1000\n   vertex 0 1000 1000\n  endloop\n endfacet\n"
                     "endsolid nose cone\n");
        auto f = readStl(p, LengthUnit::Millimetre);
        QVERIFY2(f.ok(), qPrintable(f ? QString() : f.error()));
        QVERIFY(!f->binary);
        QCOMPARE(int(f->solids.size()), 2);
        const StlSolid& wing = f->solids[0];
        QCOMPARE(wing.name, QString("wing"));
        QCOMPARE(int(wing.triangles), 2);  // the degenerate facet is dropped
        QCOMPARE(int(wing.mesh.positions.size() / 3), 4);  // the diagonal's ends are shared
        QVERIFY(qAbs(wing.area - 1.0) < 1e-9);
        QCOMPARE(f->solids[1].name, QString("nose cone"));
        const Bounds b = stlBounds(f.value());
        QVERIFY(qAbs(b.max.z() - 1.0f) < 1e-6f);
    }

    // Binary files are one solid named after the file, even with a "solid" header.
    void readsBinaryAsOneSolid()
    {
        const QString p = m_dir.filePath("rear-wing.stl");
        writeBinary(p, {{0, 0, 0, 1, 0, 0, 0, 1, 0}, {1, 0, 0, 1, 1, 0, 0, 1, 0}});
        auto f = readStl(p, LengthUnit::Metre);
        QVERIFY2(f.ok(), qPrintable(f ? QString() : f.error()));
        QVERIFY(f->binary);
        QCOMPARE(int(f->solids.size()), 1);
        QCOMPARE(f->solids[0].name, QString("rear-wing"));
        QCOMPARE(int(f->solids[0].mesh.positions.size() / 3), 4);
    }

    // Sharp edges keep their own vertices (flat shading across them); smooth
    // and flat neighbours stay welded.
    void splitsSharpEdges()
    {
        std::vector<std::array<float, 9>> cube;
        const float p[8][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
        const int quads[6][4] = {{0, 2, 3, 1}, {4, 5, 7, 6}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 4, 6, 2}, {1, 3, 7, 5}};
        for (const auto& q : quads)
            for (const auto& t : {std::array<int, 3>{q[0], q[1], q[2]}, std::array<int, 3>{q[0], q[2], q[3]}})
                cube.push_back({p[t[0]][0], p[t[0]][1], p[t[0]][2], p[t[1]][0], p[t[1]][1], p[t[1]][2], p[t[2]][0], p[t[2]][1], p[t[2]][2]});
        const QString path = m_dir.filePath("cube.stl");
        writeBinary(path, cube);
        auto f = readStl(path, LengthUnit::Metre);
        QVERIFY(f.ok());
        QCOMPARE(int(f->solids[0].mesh.positions.size() / 3), 24);
        QVERIFY(qAbs(f->solids[0].area - 6.0) < 1e-9);
    }

    // A file that is neither a whole binary STL nor ASCII is refused, not read short.
    void refusesTruncatedFiles()
    {
        const QString p = m_dir.filePath("cut.stl");
        writeBinary(p, {{0, 0, 0, 1, 0, 0, 0, 1, 0}, {1, 0, 0, 1, 1, 0, 0, 1, 0}}, true);
        QVERIFY(!readStl(p, LengthUnit::Metre).ok());
        QVERIFY(!readStl(m_dir.filePath("missing.stl"), LengthUnit::Metre).ok());
    }

    // Renamed solids go to the mesher as ASCII copies in the project's unit.
    void writesRenamedSolidsBack()
    {
        const QString src = m_dir.filePath("body.stl");
        writeBinary(src, {{0, 0, 0, 1, 0, 0, 0, 1, 0}});
        auto f = readStl(src, LengthUnit::Metre);
        QVERIFY(f.ok());
        const QString out = m_dir.filePath("renamed.stl");
        QString error;
        QVERIFY2(writeAsciiStl(out, {{"monocoque", &f->solids[0].mesh}}, LengthUnit::Millimetre, &error), qPrintable(error));
        auto back = readStl(out, LengthUnit::Millimetre);
        QVERIFY(back.ok());
        QCOMPARE(back->solids[0].name, QString("monocoque"));
        QCOMPARE(back->solids[0].mesh.positions.size(), f->solids[0].mesh.positions.size());
        for (std::size_t k = 0; k < f->solids[0].mesh.positions.size(); ++k)
            QVERIFY(qAbs(back->solids[0].mesh.positions[k] - f->solids[0].mesh.positions[k]) < 1e-6f);
    }

    void guessesMillimetresFromSize()
    {
        Bounds car;
        car.valid = true;
        car.min = QVector3D(0, -700, 0);
        car.max = QVector3D(3900, 700, 1200);
        QCOMPARE(guessUnit(car), LengthUnit::Millimetre);
        car.max = QVector3D(3.9f, 0.7f, 1.2f);
        car.min = QVector3D(0, -0.7f, 0);
        QCOMPARE(guessUnit(car), LengthUnit::Metre);
    }

    // mesh.json keeps lengths in its own unit; reading it back gives the same setup.
    void meshJsonRoundTripsInMillimetres()
    {
        Bounds car;
        car.valid = true;
        car.min = QVector3D(-1, -0.7f, 0);
        car.max = QVector3D(3, 0.7f, 1.2f);
        const MeshSetup s = defaultMeshSetup(car, LengthUnit::Millimetre);
        QVERIFY2(s.problems().isEmpty(), qPrintable(s.problems().join(' ')));
        const QJsonObject j = s.toJson();
        QCOMPARE(j.value("v").toInt(), 1);
        QCOMPARE(j.value("units").toString(), QString("mm"));
        QVERIFY(qAbs(j.value("base_size").toDouble() - s.baseSize * 1000) < 1e-6);
        QCOMPARE(j.value("refinement_boxes").toArray().size(), 2);
        QVERIFY(j.value("tunnel").toObject().value("ground").toBool());
        QVERIFY(MeshSetup::fromJson(j) == s);
        // Through text, as the mesher reads it.
        QVERIFY(MeshSetup::fromJson(QJsonDocument::fromJson(QJsonDocument(j).toJson()).object()) == s);
    }

    void reportsSetupProblems()
    {
        MeshSetup s;
        s.minSize = 0.1;
        s.maxSize = 0.05;
        s.up = "+x";
        s.boxes.push_back({"wake", QVector3D(1, 0, 0), QVector3D(0, 1, 1), 0.05});
        const QString all = s.problems().join('\n');
        QVERIFY(all.contains("larger than the maximum"));
        QVERIFY(all.contains("different axes"));
        QVERIFY(all.contains("wake"));
    }

    // The tunnel is sized in car lengths from the geometry; the ground sits under the wheels.
    void sizesTheTunnelInCarLengths()
    {
        Bounds car;
        car.valid = true;
        car.min = QVector3D(0, -0.7f, 0.05f);
        car.max = QVector3D(4, 0.7f, 1.2f);
        MeshSetup s;
        const Bounds t = tunnelBounds(s, car);
        QVERIFY(t.valid);
        QVERIFY(qAbs(t.min.x() - (0 - 3 * 4)) < 1e-4f);
        QVERIFY(qAbs(t.max.x() - (4 + 6 * 4)) < 1e-4f);
        QVERIFY(qAbs(t.min.z() - 0.05f) < 1e-4f);
        QVERIFY(qAbs(t.max.y() - (0.7f + 2 * 4)) < 1e-4f);
        s.tunnel.symmetry = "y";
        QVERIFY(qAbs(tunnelBounds(s, car).min.y()) < 1e-6f);
        QStringList names;
        for (const auto& f : tunnelFaces(s)) names << f.name;
        QVERIFY(names.contains("inlet") && names.contains("outlet") && names.contains("ground") && names.contains("symmetry"));
    }

    // Finer surfaces cost cells; the estimate turns that into device memory.
    void estimatesCellsAndMemory()
    {
        Bounds car;
        car.valid = true;
        car.min = QVector3D(0, -0.7f, 0);
        car.max = QVector3D(4, 0.7f, 1.2f);
        MeshSetup s = defaultMeshSetup(car, LengthUnit::Metre);
        const CellEstimate coarse = estimateCells(s, car, 12.0);
        s.minSize /= 2;
        const CellEstimate fine = estimateCells(s, car, 12.0);
        QVERIFY(coarse.total() > 1e5);
        QVERIFY(fine.total() > coarse.total());
        QVERIFY(fine.solverBytes() > coarse.solverBytes());
        const double h = firstHeightForYPlus(30, 20, 4);
        QVERIFY(h > 1e-5 && h < 1e-2);
        QVERIFY(qAbs(yPlusForFirstHeight(h, 20, 4) - 30) < 1e-6);
    }

    void buildsTheMeshCommand()
    {
        const QStringList a = meshArguments({{"/g/car.stl", "/g/wing.stl"}, "/p/mesh.json", "/p/mesh"});
        QCOMPARE(a.first(), QString("mesh"));
        QVERIFY(a.indexOf("/g/car.stl") == 1 && a.indexOf("/g/wing.stl") == 2);
        QCOMPARE(a.value(a.indexOf("--setup") + 1), QString("/p/mesh.json"));
        QCOMPARE(a.value(a.indexOf("--out") + 1), QString("/p/mesh"));
        QCOMPARE(a.value(a.indexOf("--progress") + 1), QString("jsonl"));
        QCOMPARE(meshStages().first(), QString("read"));
        QCOMPARE(meshStages().last(), QString("write"));
    }

    void parsesMeshEvents()
    {
        auto p = parseMeshEvent(R"({"v":1,"event":"mesh-progress","stage":"layers","fraction":1.4,"message":"5 layers"})");
        QVERIFY(p.ok());
        const auto* prog = std::get_if<MeshProgressEvent>(&p.value());
        QVERIFY(prog && prog->stage == "layers" && prog->fraction == 1.0);
        auto e = parseMeshEvent(R"({"v":1,"event":"error","message":"open shell","advice":"close the gap"})");
        QVERIFY(e.ok() && std::get_if<MeshErrorEvent>(&e.value())->advice == "close the gap");
        QVERIFY(std::holds_alternative<MeshDoneEvent>(parseMeshEvent(R"({"v":1,"event":"mesh-done","cells":10})").value()));
        QVERIFY(!parseMeshEvent(R"({"v":2,"event":"mesh-done"})").ok());
        QVERIFY(!parseMeshEvent("Reading STL...").ok());
    }

    // The report is read tolerantly: either layout of counts and quality, missing parts empty.
    void readsMeshReportsTolerantly()
    {
        auto a = parseMeshReport(QJsonDocument::fromJson(R"({"v":1,
            "mesh":{"cells":{"total":1000,"hex":900,"prism":100},"faces":3000,"nodes":1100},
            "zones":[{"name":"rear-wing","type":"wall","faces":120}],
            "quality":{"orthogonal_quality":{"min":0.2,"mean":0.9,"histogram":{"edges":[0,0.5,1],"counts":[3,997]}},
                       "skewness":{"max":0.7,"mean":0.1,"histogram":{"edges":[0,1],"counts":[1,2]}}},
            "worst_cells":[{"metric":"skewness","value":0.7,"centroid":[1,2,3],"zone":"rear-wing"}],
            "repairs":["closed 2 holes",{"kind":"merged duplicate vertices","count":14,"solid":"monocoque"}],
            "y_plus":{"mean":28,"speed":20}})").object());
        QVERIFY2(a.ok(), qPrintable(a ? QString() : a.error()));
        QCOMPARE(a->cells, qint64(1000));
        QCOMPARE(a->cellTypes.value("prism"), qint64(100));
        QCOMPARE(int(a->zones.size()), 1);
        const MeshQuality* oq = a->metric("orthogonal_quality");
        QVERIFY(oq && oq->higherIsBetter && *oq->min == 0.2 && oq->counts.size() == 2);
        const MeshQuality* sk = a->metric("skewness");
        QVERIFY(sk && !sk->higherIsBetter && sk->counts.empty());  // edges and counts disagree: dropped
        QCOMPARE(int(a->worst.size()), 1);
        QVERIFY(a->worst[0].position && a->worst[0].position->z() == 3.f);
        QCOMPARE(a->repairs.size(), 2);
        QVERIFY(a->repairs[1].startsWith("monocoque: ") && a->repairs[1].contains("14"));

        auto b = parseMeshReport(QJsonDocument::fromJson(R"({"cells":500,"cell_types":{"hex":400,"poly":100},
            "quality":[{"name":"volume_ratio","max":12,"mean":1.4}]})").object());
        QVERIFY(b.ok());
        QCOMPARE(b->cells, qint64(500));
        QVERIFY(b->metric("volume_ratio") && *b->metric("volume_ratio")->max == 12);
        QVERIFY(b->worst.empty() && b->repairs.isEmpty() && !b->yPlus);
        QVERIFY(!parseMeshReport({}).ok());
        QVERIFY(!parseMeshReport(QJsonObject{{"v", 2}}).ok());
    }
};

QTEST_GUILESS_MAIN(TestGeometry)
#include "tst_geometry.moc"
