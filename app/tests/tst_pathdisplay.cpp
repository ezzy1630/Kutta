// The app's rule for showing paths: never the home folder (screenshots and
// screen shares carry what is on screen), and compact in narrow fields.

#include "core/PathDisplay.h"

#include <QDir>
#include <QTest>

using namespace kapp;

class TestPathDisplay : public QObject {
    Q_OBJECT
private slots:
    void homeBecomesTilde()
    {
        const QString home = "/Users/ana";
        QCOMPARE(homeRelativePath("/Users/ana/cases/car.kutta.json", home), QDir::toNativeSeparators("~/cases/car.kutta.json"));
        QCOMPARE(homeRelativePath("/Users/ana", home), QStringLiteral("~"));
        // A sibling whose name starts with home's is not home.
        QCOMPARE(homeRelativePath("/Users/anabel/car.msh", home), QDir::toNativeSeparators("/Users/anabel/car.msh"));
        QCOMPARE(homeRelativePath("/tmp/run/car.msh", home), QDir::toNativeSeparators("/tmp/run/car.msh"));
        // The real home, as the app uses it.
        QVERIFY(!homeRelativePath(QDir::home().filePath("x/y.json")).contains(QDir::toNativeSeparators(QDir::homePath())));
    }

    void compactKeepsTheFileAndItsNearestFolders()
    {
        const QString home = "/Users/ana";
        const auto native = [](const char* s) { return QDir::toNativeSeparators(QString::fromUtf8(s)); };
        QCOMPARE(compactPath("/Users/ana/.t3/worktrees/Kutta/t3-c77b8b28/build/kutta", 2, home), native("~/…/t3-c77b8b28/build/kutta"));
        QCOMPARE(compactPath("/var/folders/p7/mm8q/T/Kutta-MpJino/cases/naca4412_ground.kutta.json", 2, home),
                 native("…/Kutta-MpJino/cases/naca4412_ground.kutta.json"));
        QCOMPARE(compactPath("/Users/ana/cases/car.kutta.json", 2, home), native("~/cases/car.kutta.json"));
        QCOMPARE(compactPath("/tmp/kc77/runs/wing/naca4412_ground.msh", 1, home), native("…/wing/naca4412_ground.msh"));
        QCOMPARE(compactPath("car.msh", 2, home), QStringLiteral("car.msh"));
        QCOMPARE(compactPath("/Users/ana/a/b/c/d/e/f.json", 0, home), native("~/…/f.json"));
    }

    void freeTextLosesTheHomeFolder()
    {
        const QString home = "/Users/ana";
        QCOMPARE(withHomeAsTilde("Reading mesh /Users/ana/cases/car.msh", home), QStringLiteral("Reading mesh ~/cases/car.msh"));
        QCOMPARE(withHomeAsTilde("$ kutta solve /Users/ana/a.cas.h5 --out /Users/ana/runs/r1", home),
                 QStringLiteral("$ kutta solve ~/a.cas.h5 --out ~/runs/r1"));
        QCOMPARE(withHomeAsTilde("'/Users/ana': no such file", home), QStringLiteral("'~': no such file"));
        QCOMPARE(withHomeAsTilde("C:\\Users\\ana\\car.msh and C:/Users/ana/x", "C:/Users/ana"), QStringLiteral("~\\car.msh and ~/x"));
        QCOMPARE(withHomeAsTilde("/Users/anabel/car.msh", home), QStringLiteral("/Users/anabel/car.msh"));
        QCOMPARE(withHomeAsTilde("/Users/ana.bak/x", home), QStringLiteral("/Users/ana.bak/x"));
        QCOMPARE(withHomeAsTilde("/ is root", "/"), QStringLiteral("/ is root"));
    }
};

QTEST_GUILESS_MAIN(TestPathDisplay)
#include "tst_pathdisplay.moc"
