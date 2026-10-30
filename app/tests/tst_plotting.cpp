// Colormaps, colour ranges, chart axes and palette matching.

#include "core/ChartLayout.h"
#include "core/ChartScale.h"
#include "core/Colormap.h"
#include "core/TextMatch.h"

#include <QFontMetricsF>
#include <QtTest>

using namespace kapp;

class TestPlotting : public QObject {
    Q_OBJECT
private slots:
    // Legend, axis titles and tick labels never collide, at any pane size
    // (the XY plot's legend once sat on its rotated y title).
    void chartLayoutNeverOverlaps()
    {
        Axis x, y;
        x.setRange(0.1, 0.2);
        y.setRange(-0.4, 1.4);
        const QStringList shortKeys{"Velocity Magnitude", "Velocity Magnitude · Fluent"};
        const QStringList residuals{"continuity", "x-velocity", "y-velocity", "z-velocity", "k", "omega"};
        QStringList many;
        for (int i = 0; i < 14; ++i) many << QStringLiteral("Static Pressure · Fluent · wall-elbow-section-%1").arg(i);
        QFont font;
        font.setPointSizeF(11);
        const QFontMetricsF fm(font);
        for (const QStringList& keys : {QStringList{}, shortKeys, residuals, many})
            for (int w : {220, 300, 480, 760, 1400})
                for (int h : {140, 200, 320, 640})
                    for (const QString& yl : {QString(), QStringLiteral("Velocity Magnitude [m/s]"), QStringLiteral("Static Pressure Coefficient along the section [1]")}) {
                        ChartLayoutInput in{QRectF(0, 0, w, h), font, 26, keys, "Distance along curve [m]", yl, &x, &y};
                        const ChartLayout l = chartLayout(in);
                        const QString where = QStringLiteral("%1x%2, %3 keys, y '%4'").arg(w).arg(h).arg(keys.size()).arg(yl);
                        QVERIFY2(l.plot.width() > 0 && l.plot.height() > 0, qPrintable(where));
                        QVERIFY2(in.rect.contains(l.plot), qPrintable(where));
                        for (std::size_t i = 0; i < l.legend.size(); ++i) {
                            const QRectF a = l.legend[i].rect;
                            QVERIFY2(a.right() <= in.rect.right() && a.left() >= in.rect.left(), qPrintable(where));
                            QVERIFY2(!a.intersects(l.plot), qPrintable(where));
                            QVERIFY2(l.yTitle.isEmpty() || !a.intersects(l.yTitle), qPrintable(where + " legend/y title"));
                            for (std::size_t j = i + 1; j < l.legend.size(); ++j) QVERIFY2(!a.intersects(l.legend[j].rect), qPrintable(where));
                        }
                        // Every key is either shown or folded into "+N more".
                        int shown = 0;
                        for (const auto& c : l.legend) shown += c.key >= 0;
                        QCOMPARE(shown + l.hiddenLegend.size(), keys.size());
                        if (!yl.isEmpty()) {
                            QVERIFY2(!l.yTitle.isEmpty() && !l.yTitle.intersects(l.plot), qPrintable(where));
                            // Rotated titles fit whole along the plot; flat ones sit above it.
                            if (l.yRotated) QVERIFY2(l.yTitle.height() >= fm.horizontalAdvance(yl), qPrintable(where));
                            else QVERIFY2(l.yTitle.bottom() <= l.plot.top(), qPrintable(where));
                        }
                        QVERIFY2(!l.xTitle.intersects(l.plot) && l.xTitle.bottom() <= in.rect.bottom(), qPrintable(where));
                    }
    }

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

QTEST_MAIN(TestPlotting)
#include "tst_plotting.moc"
