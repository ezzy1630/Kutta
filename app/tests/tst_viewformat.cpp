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
} // namespace

class TestViewFormat : public QObject {
    Q_OBJECT
private slots:
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
