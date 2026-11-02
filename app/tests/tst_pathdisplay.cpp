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
        const QString home = "/Users/user";
        QCOMPARE(homeRelativePath("/Users/user/cases/car.kutta.json", home), QDir::toNativeSeparators("~/cases/car.kutta.json"));
        QCOMPARE(homeRelativePath("/Users/user", home), QStringLiteral("~"));
        // A sibling whose name starts with home's is not home.
        QCOMPARE(homeRelativePath("/Users/username/car.msh", home), QDir::toNativeSeparators("/Users/username/car.msh"));
        QCOMPARE(homeRelativePath("/tmp/run/car.msh", home), QDir::toNativeSeparators("/tmp/run/car.msh"));
        // The real home, as the app uses it.
        QVERIFY(!homeRelativePath(QDir::home().filePath("x/y.json")).contains(QDir::toNativeSeparators(QDir::homePath())));
    }

    void compactKeepsTheFileAndItsNearestFolders()
    {
        const QString home = "/Users/user";
        const auto native = [](const char* s) { return QDir::toNativeSeparators(QString::fromUtf8(s)); };
        QCOMPARE(compactPath("/Users/user/.t3/worktrees/Kutta/t3-c77b8b28/build/kutta", 2, home), native("~/…/t3-c77b8b28/build/kutta"));
        QCOMPARE(compactPath("/var/folders/p7/mm8q/T/Kutta-MpJino/cases/naca4412_ground.kutta.json", 2, home),
                 native("…/Kutta-MpJino/cases/naca4412_ground.kutta.json"));
        QCOMPARE(compactPath("/Users/user/cases/car.kutta.json", 2, home), native("~/cases/car.kutta.json"));
        QCOMPARE(compactPath("/tmp/kc77/runs/wing/naca4412_ground.msh", 1, home), native("…/wing/naca4412_ground.msh"));
        QCOMPARE(compactPath("car.msh", 2, home), QStringLiteral("car.msh"));
        QCOMPARE(compactPath("/Users/user/a/b/c/d/e/f.json", 0, home), native("~/…/f.json"));
    }

    void freeTextLosesTheHomeFolder()
    {
        const QString home = "/Users/user";
        QCOMPARE(withHomeAsTilde("Reading mesh /Users/user/cases/car.msh", home), QStringLiteral("Reading mesh ~/cases/car.msh"));
        QCOMPARE(withHomeAsTilde("$ kutta solve /Users/user/a.cas.h5 --out /Users/user/runs/r1", home),
                 QStringLiteral("$ kutta solve ~/a.cas.h5 --out ~/runs/r1"));
        QCOMPARE(withHomeAsTilde("'/Users/user': no such file", home), QStringLiteral("'~': no such file"));
        // Windows: home in both separator styles, as CLI messages and Qt write it.
        const QString windows = QString("C:/Users/user/car.msh").replace('/', '\\') + " and C:/Users/user/x";
        QCOMPARE(withHomeAsTilde(windows, "C:/Users/user"), QStringLiteral("~\\car.msh and ~/x"));
        QCOMPARE(withHomeAsTilde("/Users/username/car.msh", home), QStringLiteral("/Users/username/car.msh"));
        QCOMPARE(withHomeAsTilde("/ is root", "/"), QStringLiteral("/ is root"));
    }
};

QTEST_GUILESS_MAIN(TestPathDisplay)
#include "tst_pathdisplay.moc"
