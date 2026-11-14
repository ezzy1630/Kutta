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

    // mesh.json lengths are metres whatever the STL unit; `units` scales STL
    // coordinates only. Reading it back gives the same setup.
    void meshJsonIsInMetresForMillimetreGeometry()
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
        QCOMPARE(j.value("base_size").toDouble(), s.baseSize);  // metres, not mm
        QVERIFY(j.value("base_size").toDouble() < 10);
        QCOMPARE(j.value("layers").toObject().value("first_height").toDouble(), s.firstHeight);
        const QJsonArray box = j.value("refinement_boxes").toArray().first().toObject().value("max").toArray();
        QVERIFY(qAbs(box[0].toDouble() - double(s.boxes.front().max.x())) < 1e-6);
        QVERIFY(!j.value("surface").toObject().contains("proximity_min_size"));  // automatic: left to the mesher
        QVERIFY(j.value("tunnel").toObject().value("ground").toBool());
        QVERIFY(MeshSetup::fromJson(j) == s);
        // Through text, as the mesher reads it.
        QVERIFY(MeshSetup::fromJson(QJsonDocument::fromJson(QJsonDocument(j).toJson()).object()) == s);
        // A project saved before the real mesher wrote its lengths in mm.
        QJsonObject old = j;
        old["base_size"] = s.baseSize * 1000;
        QCOMPARE(MeshSetup::fromJson(old, true).baseSize, s.baseSize);
        // Keys a setup omits take kutta mesh's defaults.
        const MeshSetup d = MeshSetup::fromJson(QJsonObject{{"v", 1}});
        QCOMPARE(d.baseSize, 0.25);
        QCOMPARE(d.proximityCells, 3);
        QCOMPARE(d.layers, 4);
        QCOMPARE(d.firstHeight, 0.002);
        QCOMPARE(d.maxCells, qint64(10'000'000));
    }

    // Refused before the mesher would refuse it, in the mesher's ranges.
    void reportsSetupProblems()
    {
        MeshSetup s;
        s.minSize = 0.1;
        s.maxSize = 0.05;
        s.up = "+y";
        s.proximityCells = 0;
        s.growth = 1.0;
        s.boxes.push_back({"wake", QVector3D(1, 0, 0), QVector3D(0, 1, 1), 0.05});
        s.boxes.push_back({"wake", QVector3D(0, 0, 0), QVector3D(1, 1, 1), 0.05});
        const QString all = s.problems().join('\n');
        QVERIFY(all.contains("larger than the maximum"));
        QVERIFY(all.contains("+x and +z up"));
        QVERIFY(all.contains("Cells per gap"));
        QVERIFY(all.contains("Growth rate"));
        QVERIFY(all.contains("unique"));
        QVERIFY(all.contains("each minimum must be below"));
        QVERIFY(MeshSetup().problems().isEmpty());
    }

    // The tunnel as kutta mesh builds it: body x lengths from the bounds, the
    // ground at z = 0, upper limits rounded out to whole base cells.
    void sizesTheTunnelLikeTheMesher()
    {
        Bounds car;
        car.valid = true;
        car.min = QVector3D(0, -0.7f, 0.05f);
        car.max = QVector3D(4, 0.7f, 1.2f);
        MeshSetup s;
        s.baseSize = 0.25;
        const Bounds t = tunnelBounds(s, car);
        QVERIFY(t.valid);
        QVERIFY(qAbs(t.min.x() - (0 - 3 * 4)) < 1e-4f);
        QVERIFY(qAbs(t.max.x() - (4 + 6 * 4)) < 1e-4f);  // 40 m: a whole number of 0.25 m cells
        QCOMPARE(t.min.z(), 0.0f);                        // ground at z = 0, not at the body
        QVERIFY(qAbs(t.min.y() - (-0.7f - 2 * 4)) < 1e-4f);
        // y: 17.4 m from -8.7 rounds out to 17.5 m.
        QVERIFY(qAbs(t.max.y() - (-8.7f + 17.5f)) < 1e-4f);
        s.tunnel.ground = false;
        QVERIFY(qAbs(tunnelBounds(s, car).min.z() - (0.05f - 2 * 4)) < 1e-4f);  // `side` lengths below
        s.tunnel.ground = true;
        QStringList names;
        for (const auto& f : tunnelFaces(s)) names << f.name;
        QCOMPARE(names, (QStringList{"inlet", "outlet", "side-low", "side-high", "ground", "top"}));
        s.tunnel.symmetry = "y";
        QCOMPARE(tunnelBounds(s, car).min.y(), 0.0f);
        names.clear();
        for (const auto& f : tunnelFaces(s)) names << f.name;
        QVERIFY(names.contains("symmetry") && !names.contains("side-low"));
        for (const QString& n : names) QVERIFY2(reservedZoneNames().contains(n), qPrintable(n));
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

    // The estimate against the real coarse meshes it was calibrated on.
    void estimateMatchesRealCoarseMeshes()
    {
        auto check = [](const char* json, Bounds body, double area, double cells) {
            const MeshSetup s = MeshSetup::fromJson(QJsonDocument::fromJson(json).object());
            const double n = estimateCells(s, body, area).total();
            QVERIFY2(qAbs(n / cells - 1) < 0.1, qPrintable(QStringLiteral("estimate %1 for %2 real cells").arg(n).arg(cells)));
        };
        Bounds ahmed;
        ahmed.valid = true;
        ahmed.min = QVector3D(-1.044f, -0.1945f, 0.05f);
        ahmed.max = QVector3D(0, 0.1945f, 0.338f);
        check(R"({"v":1,"base_size":0.2,"surface":{"min_size":0.025,"max_size":0.05,"proximity_min_size":0.00625,"proximity_cells":3},
                  "layers":{"count":3,"first_height":0.002,"growth":1.3}})", ahmed, 1.517, 144877);
        Bounds fsae;
        fsae.valid = true;
        fsae.min = QVector3D(-1.2f, -0.825f, 0);
        fsae.max = QVector3D(1.906f, 0.825f, 1.24f);
        check(R"({"v":1,"base_size":0.5,"surface":{"min_size":0.04,"max_size":0.1,"proximity_min_size":0.003,"proximity_cells":1},
                  "layers":{"count":4,"first_height":0.002,"growth":1.3}})", fsae, 21.38, 1312313);
    }

    void buildsTheMeshCommand()
    {
        const QStringList a = meshArguments({{"/g/car.stl", "/g/wing.stl"}, "/p/mesh.json", "/p/mesh", 25});
        QCOMPARE(a.first(), QString("mesh"));
        QVERIFY(a.indexOf("/g/car.stl") == 1 && a.indexOf("/g/wing.stl") == 2);
        QCOMPARE(a.value(a.indexOf("--setup") + 1), QString("/p/mesh.json"));
        QCOMPARE(a.value(a.indexOf("--out") + 1), QString("/p/mesh"));
        QCOMPARE(a.value(a.indexOf("--speed") + 1), QString("25"));
        QCOMPARE(a.value(a.indexOf("--progress") + 1), QString("jsonl"));
        QVERIFY(!a.contains("--case"));
        QVERIFY(!meshArguments({{"/g/car.stl"}, {}, "/p/mesh"}).contains("--speed"));
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

    // A real report (Ahmed body, coarse setup, kutta mesh at fb44554) reads
    // into every figure the Mesh Report page shows.
    void readsARealMeshReport()
    {
        auto r = loadMeshReport(QStringLiteral(KUTTA_TEST_DATA "/mesh-report-ahmed-coarse.json"));
        QVERIFY2(r.ok(), qPrintable(r ? QString() : r.error()));
        QCOMPARE(r->cells, qint64(144877));
        QCOMPARE(r->faces, qint64(459815));
        QCOMPARE(*r->cellTypes.prism, qint64(51090));
        QVERIFY(qAbs(*r->cellTypes.hexFractionOfCore - 0.99495) < 1e-4);
        QCOMPARE(int(r->zones.size()), 7);
        QCOMPARE(r->zones[2].name, QString("side-low"));
        QVERIFY(r->quality && qAbs(*r->quality->orthogonalMin - 0.15446) < 1e-4);
        QCOMPARE(int(r->quality->orthogonalHistogram.size()), 9);
        QCOMPARE(r->quality->orthogonalHistogram[3].count, qint64(7));
        QCOMPARE(*r->quality->cellsBelow015, qint64(0));
        QVERIFY(r->quality->worstCell && qAbs(r->quality->worstCell->y() - 0.1081f) < 1e-3f);
        QVERIFY(r->quality->volumeRatioCore && qAbs(*r->quality->volumeRatioCore->p999 - 9.556) < 1e-3);
        QCOMPARE(int(r->quality->volumeRatioLayers->histogram.size()), 7);
        QVERIFY(r->quality->volumeRatioLayers->histogram.back().to > 1e299);
        QVERIFY(r->gaps && r->gaps->size() == 1);
        QCOMPARE(r->gaps->front().b, QString("ground"));
        QCOMPARE(r->gaps->front().cellsAcross, qint64(3));
        QVERIFY(r->closedGaps().empty());
        QVERIFY(r->layers && qAbs(*r->layers->coverageNone - 0.03541) < 1e-4);
        QCOMPARE(*r->layers->generated, qint64(3));
        QVERIFY(r->snap && qAbs(*r->snap->distanceMax - 0.007843) < 1e-5);
        QVERIFY(r->yPlus && qAbs(*r->yPlus->p50 - 56.73) < 0.01 && *r->yPlus->speed == 20);
        QVERIFY(qAbs(*r->yPlus->from30to300 - 0.9564) < 1e-4);
        QCOMPARE(r->repairs.value("components"), qint64(1));
        QVERIFY(r->missingSurfaceZones && r->missingSurfaceZones->isEmpty());
        QCOMPARE(r->warnings.size(), 3);
        QVERIFY(r->peakMemoryBytes && *r->peakMemoryBytes > 0);
        QVERIFY(qAbs(solverBytesEstimate(r->cells) - 144877 * 2600.0) < 1);
        QVERIFY(solverBytesEstimate(r->cells) < kTwelveGbUsableBytes);
        QVERIFY(solverBytesEstimate(4'417'972) > kTwelveGbUsableBytes);  // the default FSAE mesh
    }

    // Strict: a key of the wrong type fails with its path; what the report
    // leaves out (or reports null) stays unset, never a default.
    void readsMeshReportsStrictly()
    {
        auto minimal = parseMeshReport(QJsonDocument::fromJson(R"({"v":1,"cells":10,"faces":40,"nodes":20,"zones":[],"peak_memory_bytes":null})").object());
        QVERIFY2(minimal.ok(), qPrintable(minimal ? QString() : minimal.error()));
        QVERIFY(!minimal->quality && !minimal->gaps && !minimal->layers && !minimal->snap && !minimal->yPlus);
        QVERIFY(!minimal->missingSurfaceZones && !minimal->peakMemoryBytes && !minimal->wallSeconds && !minimal->cellTypes.hex);
        QVERIFY(minimal->repairs.isEmpty() && minimal->warnings.isEmpty());

        auto wrong = parseMeshReport(QJsonDocument::fromJson(R"({"v":1,"cells":10,"faces":40,"nodes":20,"zones":[],
            "quality":{"orthogonal_min":"low"},"layers":{"coverage":{"none":"all"}}})").object());
        QVERIFY(!wrong.ok());
        QVERIFY2(wrong.error().contains("quality/orthogonal_min") && wrong.error().contains("layers/coverage/none"), qPrintable(wrong.error()));
        QVERIFY(!parseMeshReport(QJsonDocument::fromJson(R"({"v":1,"cells":10,"zones":[]})").object()).ok());  // faces, nodes required
        QVERIFY(!parseMeshReport(QJsonDocument::fromJson(R"({"v":1,"cells":10.5,"faces":40,"nodes":20,"zones":[]})").object()).ok());
        QVERIFY(!parseMeshReport({}).ok());
        QVERIFY(!parseMeshReport(QJsonObject{{"v", 2}}).ok());
        QVERIFY(!parseMeshReport(QJsonObject{{"cells", 10}}).ok());  // no version

        // A closed slot, against parts that merely touch (a tyre on its plinth).
        auto gaps = parseMeshReport(QJsonDocument::fromJson(R"({"v":1,"cells":10,"faces":40,"nodes":20,"zones":[],
            "castellation":{"gaps":[{"between":["rear-wing-element-1","rear-wing-element-2"],"min_cells_across":0,"width_m":0.0128},
                                    {"between":["wheel-fl","ground"],"min_cells_across":0,"width_m":0},
                                    {"between":["front-wing-element-2","front-wing-element-3"],"min_cells_across":1,"width_m":0.02}]},
            "layers":{"requested":0,"generated":0,"cells":0}})").object());
        QVERIFY(gaps.ok());
        const auto closed = gaps->closedGaps();
        QCOMPARE(int(closed.size()), 2);
        QCOMPARE(int(std::count_if(closed.begin(), closed.end(), [](const MeshGap& g) { return g.width > 0; })), 1);
        QCOMPARE(*gaps->layers->requested, qint64(0));
        QVERIFY(!gaps->layers->coverageNone);
    }
};

QTEST_GUILESS_MAIN(TestGeometry)
#include "tst_geometry.moc"
