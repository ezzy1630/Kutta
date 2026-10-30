// Colormaps, colour ranges, chart axes and palette matching.

#include "core/ChartScale.h"
#include "core/Colormap.h"
#include "core/TextMatch.h"

#include <QtTest>

using namespace kapp;

class TestPlotting : public QObject {
    Q_OBJECT
private slots:
    void colormapEndpointsMatchMatplotlib()
    {
        QCOMPARE(colormap(ColormapId::Viridis).sample(0).rgb(), qRgb(68, 1, 84));
        QCOMPARE(colormap(ColormapId::Viridis).sample(1).rgb(), qRgb(253, 231, 37));
        QCOMPARE(colormap(ColormapId::FluentRainbow).sample(0).rgb(), qRgb(0, 0, 255));
        QCOMPARE(colormap(ColormapId::FluentRainbow).sample(1).rgb(), qRgb(255, 0, 0));
        for (ColormapId id : allColormaps()) QCOMPARE(colormapFromKey(colormap(id).key), id);
    }

    void deltasUseADivergingMapCentredOnZero()
    {
        QVERIFY(colormap(defaultColormap(true)).diverging);
        QVERIFY(!colormap(defaultColormap(false)).diverging);
        const ColorRange r = autoRange(-0.3, 0.2, true);
        QCOMPARE(r.min, -0.3);
        QCOMPARE(r.max, 0.3);
        QCOMPARE(r.normalize(0), 0.5);
        const QColor mid = colormap(ColormapId::CoolWarm).sample(0.5);
        QVERIFY(mid.saturationF() < 0.1);  // agreement reads as neutral
    }

    void bandsQuantizeLikeTheShader()
    {
        QCOMPARE(bandOf(0.0, 10), 0);
        QCOMPARE(bandOf(0.999, 10), 9);
        QCOMPARE(bandOf(1.0, 10), 9);
        QCOMPARE(bandOf(0.55, 10), 5);
    }

    void logAxisTicksOnDecades()
    {
        Axis a;
        a.setScale(Axis::Scale::Log);
        a.setRange(1e-6, 1);
        a.setPixels(300, 0);
        QStringList majors;
        for (const Tick& t : a.ticks(7))
            if (t.major) majors << t.label;
        QVERIFY(majors.contains("1e-6"));
        QVERIFY(majors.contains("1e0"));
        QCOMPARE(a.toPixel(1e-3), 150.0);
        QVERIFY(qAbs(a.fromPixel(a.toPixel(3e-4)) - 3e-4) < 1e-12);
    }

    void linearAxisNiceRangeAndZoom()
    {
        double lo = 0.13, hi = 9.7;
        Axis::niceRange(Axis::Scale::Linear, lo, hi);
        QCOMPARE(lo, 0.0);
        QCOMPARE(hi, 10.0);
        Axis a;
        a.setRange(0, 100);
        a.setPixels(0, 1000);
        a.zoom(50, 0.5);
        QCOMPARE(a.lo(), 25.0);
        QCOMPARE(a.hi(), 75.0);
        QCOMPARE(niceStep(0.37), 0.5);
        QCOMPARE(formatNumber(1500), QString("1500"));
        QCOMPARE(formatNumber(0.0012), QString("0.0012"));
    }

    void paletteRanksPrefixesFirst()
    {
        QVERIFY(fuzzyScore("inl", "inlet") > fuzzyScore("inl", "cold-inlet"));
        QVERIFY(fuzzyScore("bc", "Boundary Conditions") > 0);
        QVERIFY(fuzzyScore("xyz", "Boundary Conditions") < 0);
    }
};

QTEST_GUILESS_MAIN(TestPlotting)
#include "tst_plotting.moc"
