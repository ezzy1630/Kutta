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

    // Long approximation badges wrap below the title instead of cutting it
    // short ("Scaled Resid…"); nothing in the header overlaps the legend or plot.
    void chartHeaderKeepsTitleAndBadges()
    {
        Axis x, y;
        x.setRange(0, 300);
        y.setScale(Axis::Scale::Log);
        y.setRange(1e-9, 1);
        QFont font, titleFont, badgeFont;
        font.setPointSizeF(11);
        titleFont.setPointSizeF(12);
        titleFont.setWeight(QFont::DemiBold);
        badgeFont.setPointSizeF(9.4);
        badgeFont.setWeight(QFont::Bold);
        const QFontMetricsF tfm(titleFont);
        const QString title = "Scaled Residuals · Run 1";
        const QStringList keys{"continuity", "x-velocity", "y-velocity", "z-velocity", "k", "epsilon"};
        for (const QStringList& badges : {QStringList{}, QStringList{"FIXTURE REPLAY"},
                                          QStringList{"APPROXIMATE · No compressibility correction, Documented k-epsilon only", "MODIFIED"}})
            for (int w : {220, 300, 480, 760, 1400})
                for (int h : {140, 200, 320}) {
                    ChartLayoutInput in{QRectF(0, 0, w, h), font, 26, keys, "Iteration", {}, &x, &y, title, titleFont, badges, badgeFont};
                    const ChartLayout l = chartLayout(in);
                    const QString where = QStringLiteral("%1x%2, %3 badges").arg(w).arg(h).arg(badges.size());
                    QVERIFY2(in.rect.contains(l.title), qPrintable(where));
                    if (tfm.horizontalAdvance(title) <= w - 20) QVERIFY2(l.titleText == title, qPrintable(where + ": title elided"));
                    QCOMPARE(int(l.badges.size()), int(badges.size()));
                    for (std::size_t i = 0; i < l.badges.size(); ++i) {
                        const QRectF b = l.badges[i].rect;
                        QVERIFY2(b.left() >= in.rect.left() && b.right() <= in.rect.right(), qPrintable(where));
                        QVERIFY2(!b.intersects(l.title), qPrintable(where + ": badge on title"));
                        QVERIFY2(!b.intersects(l.plot), qPrintable(where + ": badge on plot"));
                        for (const auto& c : l.legend) QVERIFY2(!b.intersects(c.rect), qPrintable(where + ": badge on legend"));
                        for (std::size_t j = i + 1; j < l.badges.size(); ++j) QVERIFY2(!b.intersects(l.badges[j].rect), qPrintable(where));
                    }
                    for (const auto& c : l.legend) QVERIFY2(!c.rect.intersects(l.title), qPrintable(where + ": legend on title"));
                    QVERIFY2(l.plot.height() > 0 && in.rect.contains(l.plot), qPrintable(where));
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

        // SIMPLE's paired state takes residuals to 1e-10: a short chart labels
        // both ends of the span, and every decade keeps a grid line.
        a.setRange(1e-10, 1);
        for (int target : {3, 4, 5}) {
            majors.clear();
            int lines = 0;
            for (const Tick& t : a.ticks(target)) {
                if (t.major) majors << t.label;
                lines += t.value == std::pow(10.0, std::round(std::log10(t.value)));
            }
            QVERIFY2(majors.contains("1e0") && majors.contains("1e-10"), qPrintable(majors.join(' ')));
            QVERIFY(majors.size() <= target + 1);
            QCOMPARE(lines, 11);
        }
        QCOMPARE(formatScientific(7.7148e-6), QString("7.71e-6"));
        QCOMPARE(formatScientific(1e-5), QString("1e-5"));
        QCOMPARE(formatScientific(5.1848e-10), QString("5.18e-10"));
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
