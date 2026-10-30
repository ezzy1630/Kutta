// The KUTV/KUTL binary readers reject anything malformed.

#include "TestSupport.h"
#include "core/SceneData.h"
#include "core/ViewFormat.h"

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
