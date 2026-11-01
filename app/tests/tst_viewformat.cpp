// The KUTV/KUTL binary readers reject anything malformed.

#include "TestSupport.h"
#include "core/SceneData.h"
#include "core/ViewFormat.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QTemporaryDir>
#include <QtTest>

using namespace kapp;

namespace {
MeshData triangle()
{
    MeshData m;
    m.positions = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    m.indices = {0, 1, 2};
    MeshData::Field f;
    f.desc.name = "pressure";
    f.desc.unit = "Pa";
    f.desc.components = 1;
    f.values = {1.f, 2.f, std::nanf("")};
    m.fields.push_back(f);
    return m;
}
MeshData quad(QVector3D a, QVector3D b, QVector3D c, QVector3D d)
{
    MeshData m;
    for (const QVector3D& v : {a, b, c, d}) m.positions.insert(m.positions.end(), {v.x(), v.y(), v.z()});
    m.indices = {0, 1, 2, 0, 2, 3};
    return m;
}

SceneItem zone(const QString& name, const QString& type, MeshData mesh)
{
    SceneItem it;
    it.id = "zone:" + name;
    it.kind = ItemKind::Zone;
    it.zoneType = type;
    it.mesh = std::make_shared<PreparedMesh>(std::move(mesh));
    return it;
}

std::shared_ptr<SceneData> scene(std::vector<SceneItem> items)
{
    auto s = std::make_shared<SceneData>();
    for (auto& it : items) s->bounds.expand(it.mesh->bounds());
    s->items = std::move(items);
    return s;
}
} // namespace

class TestViewFormat : public QObject {
    Q_OBJECT
private slots:
    // A ground plane is found by its shape, whatever it is called: a flat wall
    // spanning the domain (with float noise) is not the body. Pipe and duct
    // walls are internal flow.
    void enclosingWallsAreNotTheBody()
    {
        const auto car = scene({
            zone("road", "wall", quad({-5, -2, 0}, {10, -2, 0}, {10, 2, 1e-4f}, {-5, 2, 0})),
            zone("body", "wall", quad({0, -0.5f, 0.1f}, {2, -0.5f, 0.1f}, {2, 0.5f, 0.6f}, {0, 0.5f, 0.6f})),
            zone("roof", "wall", quad({-5, -2, 4}, {10, -2, 4}, {10, 2, 4}, {-5, 2, 4})),
            zone("side", "symmetry", quad({-5, 0, 0}, {10, 0, 0}, {10, 0, 4}, {-5, 0, 4})),
        });
        QVERIFY(isExternalFlow(*car));
        QCOMPARE(defaultMeshItems(*car), QStringList{"zone:body"});
        QCOMPARE(defaultVisibleItems(*car), QStringList{"zone:body"});

        // A square duct's walls are all flat and span it: internal, so its symmetry plane shows.
        const auto duct = scene({
            zone("floor", "wall", quad({0, 0, 0}, {10, 0, 0}, {10, 1, 0}, {0, 1, 0})),
            zone("top", "wall", quad({0, 0, 1}, {10, 0, 1}, {10, 1, 1}, {0, 1, 1})),
            zone("mid", "symmetry", quad({0, 0.5f, 0}, {10, 0.5f, 0}, {10, 0.5f, 1}, {0, 0.5f, 1})),
        });
        QVERIFY(!isExternalFlow(*duct));
        QCOMPARE(defaultVisibleItems(*duct), QStringList{"zone:mid"});
    }

    // A 2D case is shown as its plane: the +z extrusion cap, labelled as the
    // domain; the opposite cap (the same cells) is dropped, and the walls and
    // inlets (edge-on strips) are not shown first.
    void planarCaseShowsItsPlane()
    {
        auto s = scene({
            zone("wall", "wall", quad({4, 4, 0}, {5, 4, 0}, {5, 4, 1}, {4, 4, 1})),
            zone("inlet", "velocity-inlet", quad({0, 0, 0}, {0, 10, 0}, {0, 10, 1}, {0, 0, 1})),
            zone("bottom", "symmetry", quad({0, 0, 0}, {10, 0, 0}, {10, 0, 1}, {0, 0, 1})),
            zone("front", "symmetry", quad({0, 0, 0}, {10, 0, 0}, {10, 10, 0}, {0, 10, 0})),
            zone("back", "symmetry", quad({0, 0, 1}, {10, 0, 1}, {10, 10, 1}, {0, 10, 1})),
        });
        QVERIFY(planarFace(*s).isEmpty());  // not marked planar: an ordinary thin domain
        makePlanar(*s);
        QVERIFY(s->planar);
        QCOMPARE(planarFace(*s), QString("zone:back"));  // the cap at z = 1 faces the plan view
        QVERIFY(!s->item("zone:front"));
        QCOMPARE(s->item("zone:back")->label, QString("2D domain"));
        QVERIFY(s->item("zone:bottom"));  // a symmetry boundary of the 2D case stays
        QCOMPARE(defaultVisibleItems(*s), QStringList{"zone:back"});
        QCOMPARE(defaultMeshItems(*s), QStringList{"zone:back"});
    }

    // A 2D body in a far field of inlets, outlets and symmetry (the Ansys
    // cylinder) opens framed on the body, as a wing over a ground does; a
    // channel, whose walls are its edges, opens on the whole domain.
    void planarExternalFlowFramesTheBody()
    {
        auto cylinder = scene({
            zone("wall", "wall", quad({9, -0.5f, 0}, {10, -0.5f, 0}, {10, 0.5f, 1}, {9, 0.5f, 1})),
            zone("inlet", "velocity-inlet", quad({0, -10, 0}, {0, 10, 0}, {0, 10, 1}, {0, -10, 1})),
            zone("outlet", "pressure-outlet", quad({30, -10, 0}, {30, 10, 0}, {30, 10, 1}, {30, -10, 1})),
            zone("top", "symmetry", quad({0, 10, 0}, {30, 10, 0}, {30, 10, 1}, {0, 10, 1})),
            zone("front", "symmetry", quad({0, -10, 0}, {30, -10, 0}, {30, 10, 0}, {0, 10, 0})),
            zone("back", "symmetry", quad({0, -10, 1}, {30, -10, 1}, {30, 10, 1}, {0, 10, 1})),
        });
        QVERIFY(!isExternalFlow(*cylinder));  // no ground or tunnel wall
        makePlanar(*cylinder);
        const auto focus = planarFocus(*cylinder);
        QVERIFY(focus);
        QCOMPARE(focus->min, QVector3D(8, -1.5f, 0));
        QCOMPARE(focus->max, QVector3D(13, 1.5f, 1));  // three diameters of wake: 20 behind against 9 ahead

        // A wing just above its ground: the room above is not a wake, so the
        // frame keeps one chord all round (clipped at the ground).
        auto wing = scene({
            zone("wing", "wall", quad({0, 0.05f, 0}, {1, 0.05f, 0}, {1, 0.15f, 1}, {0, 0.15f, 1})),
            zone("ground", "wall", quad({-4, 0, 0}, {5, 0, 0}, {5, 0, 1}, {-4, 0, 1})),
            zone("top", "symmetry", quad({-4, 4, 0}, {5, 4, 0}, {5, 4, 1}, {-4, 4, 1})),
            zone("front", "symmetry", quad({-4, 0, 0}, {5, 0, 0}, {5, 4, 0}, {-4, 4, 0})),
            zone("back", "symmetry", quad({-4, 0, 1}, {5, 0, 1}, {5, 4, 1}, {-4, 4, 1})),
        });
        makePlanar(*wing);
        const auto wingFocus = planarFocus(*wing);
        QVERIFY(wingFocus);
        QCOMPARE(wingFocus->min, QVector3D(-1, 0, 0));
        QCOMPARE(wingFocus->max, QVector3D(2, 1.15f, 1));

        auto channel = scene({
            zone("floor", "wall", quad({0, 0, 0}, {10, 0, 0}, {10, 0, 1}, {0, 0, 1})),
            zone("roof", "wall", quad({0, 1, 0}, {10, 1, 0}, {10, 1, 1}, {0, 1, 1})),
            zone("inlet", "velocity-inlet", quad({0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1})),
            zone("front", "symmetry", quad({0, 0, 0}, {10, 0, 0}, {10, 1, 0}, {0, 1, 0})),
            zone("back", "symmetry", quad({0, 0, 1}, {10, 0, 1}, {10, 1, 1}, {0, 1, 1})),
        });
        makePlanar(*channel);
        QVERIFY(!planarFocus(*channel));
    }

    void roundTrip()
    {
        ViewPart part;
        const QByteArray bytes = writeMeshBytes(triangle(), &part);
        QCOMPARE(part.vertices, 3u);
        QCOMPARE(*part.fields.front().min, 1.0);  // NaN is a missing sample, not a value
        auto back = readMeshBytes(bytes, part);
        QVERIFY2(back.ok(), qPrintable(back ? QString() : back.error()));
        QCOMPARE(back.value().indices, triangle().indices);
        QVERIFY(std::isnan(back.value().fields.front().values[2]));
    }

    void rejectsCorruptFiles()
    {
        ViewPart part;
        QByteArray bytes = writeMeshBytes(triangle(), &part);
        QVERIFY(!readMeshBytes(bytes.left(bytes.size() - 4), part));     // truncated
        QVERIFY(!readMeshBytes(bytes + QByteArray(4, 0), part));         // trailing bytes
        QByteArray magic = bytes;
        magic[0] = 'X';
        QVERIFY(!readMeshBytes(magic, part));
        QByteArray badIndex = bytes;
        const quint32 big = 7;
        std::memcpy(badIndex.data() + 32 + 9 * 4, &big, 4);               // index beyond vertex count
        QVERIFY(!readMeshBytes(badIndex, part));
    }

    // Zone edge files list the polygon edges (no triangle diagonals) as
    // polylines over the zone's own vertices.
    void readsZoneEdgeLines()
    {
        auto kutl = [](quint32 V, std::vector<quint32> idx, std::vector<quint32> off) {
            QByteArray b("KUTL", 4);
            const quint32 head[] = {1, 32, V, quint32(idx.size()), 0, quint32(off.size() - 1), 0};
            b.append(reinterpret_cast<const char*>(head), sizeof head);
            b.append(QByteArray(int(12 * V), 0));
            b.append(reinterpret_cast<const char*>(idx.data()), int(4 * idx.size()));
            b.append(reinterpret_cast<const char*>(off.data()), int(4 * off.size()));
            return b;
        };
        ViewPart e;
        e.kind = PartKind::Probe;
        e.file = "zone-0.edges.bin";
        e.vertices = 5;
        e.indices = 7;
        e.lineCount = 2;
        // A closed quad 0-1-2-3-0 and a separate edge 1-4.
        const QByteArray bytes = kutl(5, {0, 1, 2, 3, 0, 1, 4}, {0, 5, 7});
        auto lines = readEdgeLines(bytes, e, 5);
        QVERIFY2(lines.ok(), qPrintable(lines ? QString() : lines.error()));
        QCOMPARE(lines.value(), (std::vector<std::uint32_t>{0, 1, 1, 2, 2, 3, 3, 0, 1, 4}));
        QVERIFY(!readEdgeLines(bytes, e, 6));                                    // not the zone's vertices
        QVERIFY(!readEdgeLines(kutl(5, {0, 1, 2, 3, 0, 1, 9}, {0, 5, 7}), e, 5)); // index out of range
        QVERIFY(!readEdgeLines(kutl(5, {0, 1, 2, 3, 0, 1, 4}, {0, 6, 7}).left(bytes.size() - 4), e, 5));

        // The manifest carries them per zone; geometry-only views say so.
        const QJsonObject zone{{"id", 3}, {"name", "wall"}, {"type", "wall"}, {"file", "zone-0.surface.bin"}, {"vertices", 5}, {"indices", 6},
                               {"fields", QJsonArray{}},
                               {"edges", QJsonObject{{"file", "zone-0.edges.bin"}, {"primitive", "polylines"}, {"vertices", 5}, {"indices", 7},
                                                     {"line_count", 2}, {"fields", QJsonArray{}}}}};
        auto m = parseViewManifest({{"v", 1}, {"geometry_only", true}, {"zones", QJsonArray{zone}},
                                    {"unavailable_fields", QJsonArray{QJsonObject{{"name", "skin_friction"}, {"operation", QJsonValue()}}}}}, "/tmp");
        QVERIFY2(m.ok(), qPrintable(m ? QString() : m.error()));
        QVERIFY(m.value().geometryOnly);
        QVERIFY(m.value().zones.front().edges);
        QCOMPARE(m.value().zones.front().edges->lineCount, 2u);
        QCOMPARE(m.value().unavailableFields, QStringList{"skin_friction"});
        QJsonObject unsafe = zone;
        unsafe["edges"] = QJsonObject{{"file", "../x.bin"}, {"vertices", 5}, {"indices", 7}, {"line_count", 2}};
        QVERIFY(!parseViewManifest({{"v", 1}, {"zones", QJsonArray{unsafe}}}, "/tmp"));
    }

    void polylinesRoundTrip()
    {
        MeshData m;
        m.polylines = true;
        m.positions = {0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 1, 1, 0};
        m.indices = {0, 3, 5};  // two lines: three and two points
        m.fields.push_back({{"arc_length", "m", 1, -1, {}, {}}, {0, 1, 2, 0, 1}});
        ViewPart part;
        part.kind = PartKind::Streamlines;
        const QByteArray bytes = writeMeshBytes(m, &part);
        QCOMPARE(part.lineCount, 2u);
        QCOMPARE(part.indices, 5u);
        auto back = readMeshBytes(bytes, part);
        QVERIFY2(back.ok(), qPrintable(back ? QString() : back.error()));
        QCOMPARE(back.value().indices, m.indices);
        QCOMPARE(back.value().positions, m.positions);
        part.lineCount = 3;  // the manifest must agree with the header
        QVERIFY(!readMeshBytes(bytes, part));
    }

    // Probes written by the real CLI (export-view --post on the Ansys mixing
    // elbow): the line probe skips samples outside the volume through its
    // index array; the section has two curves on one wall.
    void readsRealCliProbesAndReports()
    {
        const QString dir = QStringLiteral(KUTTA_TEST_DATA "/real-cli-post");
        auto m = loadViewManifest(dir + "/manifest.json");
        QVERIFY2(m.ok(), qPrintable(m ? QString() : m.error()));
        QCOMPARE(m.value().probes.size(), std::size_t(2));
        QVERIFY(!m.value().availableFields.empty());
        QCOMPARE(m.value().reportsFile, QString("reports.json"));
        QVERIFY(m.value().referenceValues.contains("density"));

        const ViewPart& line = m.value().probes[0];
        QCOMPARE(line.name, QString("centre-line"));
        auto lm = readMeshFile(dir + "/" + line.file, line);
        QVERIFY2(lm.ok(), qPrintable(lm ? QString() : lm.error()));
        QVERIFY(lm.value().polylines);
        QCOMPARE(lm.value().vertexCount(), std::size_t(line.indices));  // gathered valid samples only
        QVERIFY(line.indices < line.vertices);
        const auto* valid = lm.value().field("valid");
        QVERIFY(valid);
        QVERIFY(std::all_of(valid->values.begin(), valid->values.end(), [](float v) { return v == 1; }));
        const auto* arc = lm.value().field("arc_length");
        QVERIFY(arc && std::is_sorted(arc->values.begin(), arc->values.end()));

        const ViewPart& section = m.value().probes[1];
        auto sm = readMeshFile(dir + "/" + section.file, section);
        QVERIFY2(sm.ok(), qPrintable(sm ? QString() : sm.error()));
        QCOMPARE(sm.value().indices.size(), std::size_t(section.lineCount) + 1);
        QCOMPARE(section.meta.value("lines").toArray().size(), qsizetype(section.lineCount));
        QVERIFY(sm.value().field("pressure@fluent"));

        QFile rf(dir + "/reports.json");
        QVERIFY(rf.open(QIODevice::ReadOnly));
        const QJsonObject report = QJsonDocument::fromJson(rf.readAll()).object().value("reports").toArray().first().toObject();
        QCOMPARE(report.value("zones").toArray().size(), 3);
        QVERIFY(report.value("total").toObject().value("averages").toObject().contains("pressure@fluent"));
    }

    void loadsTheRecordedElbowView()
    {
        QTemporaryDir dir;
        const QString src = "cases/mixing_elbow/run-seeded/view/";
        QFile::copy(QStringLiteral(KUTTA_FIXTURE_DIR "/") + src + "manifest.json", dir.filePath("manifest.json"));
        auto m = loadViewManifest(dir.filePath("manifest.json"));
        QVERIFY2(m.ok(), qPrintable(m ? QString() : m.error()));
        QVERIFY(m.value().reference);
        for (const auto& z : m.value().zones) {
            QFile out(dir.filePath(z.file));
            QVERIFY(out.open(QIODevice::WriteOnly));
            out.write(fixtureBytes(src + z.file));
        }
        for (const auto& z : m.value().slices) {
            QFile out(dir.filePath(z.file));
            QVERIFY(out.open(QIODevice::WriteOnly));
            out.write(fixtureBytes(src + z.file));
        }
        auto scene = loadScene(m.value(), "Run 1");
        QVERIFY2(scene.ok(), qPrintable(scene ? QString() : scene.error()));
        const auto fields = scene.value()->availableFields();
        QVERIFY(std::find(fields.begin(), fields.end(), FieldRef{"pressure", FieldRole::Delta, -1}) != fields.end());
        auto r = scene.value()->range({"velocity_magnitude", FieldRole::Kutta, -1}, {"zone:symmetry-xyplane"});
        QVERIFY(r && r->first > 0 && r->second < 2);
    }
};

QTEST_GUILESS_MAIN(TestViewFormat)
#include "tst_viewformat.moc"
