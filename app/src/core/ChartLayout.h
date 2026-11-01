#pragma once

// Where each part of a 2D chart goes. Every label owns its own rectangle:
// the title and its badges (badges wrap below a title they do not fit beside,
// so neither is cut short), legend rows, the y title (rotated beside the axis when the plot is tall
// enough, else on one line above it), tick labels and the x title, so nothing
// overlaps at any pane size. Pure geometry: PlotWidget paints it, tests check it.

#include "core/ChartScale.h"

#include <QFont>
#include <QRectF>
#include <QStringList>

#include <vector>

namespace kapp {

struct ChartLayoutInput {
    QRectF rect;
    QFont font;            // tick, legend and axis-title font
    double header = 26;    // title and badges row, when they share one row
    QStringList legend;    // unique legend keys, in order
    QString xLabel, yLabel;
    const Axis* x = nullptr;
    const Axis* y = nullptr;
    QString title;
    QFont titleFont;
    QStringList badges;    // APPROXIMATE, MODIFIED, ...: never elided for the title
    QFont badgeFont;
};

struct ChartLayout {
    struct Chip {
        QRectF rect;
        int key = -1;  // index into legend; -1 is the "+N more" chip
        QString text;  // possibly elided
    };
    struct Badge {
        QRectF rect;
        QString text;  // elided only when one badge is wider than the chart
    };
    QRectF title;      // empty without a title
    QString titleText; // elided only when it alone is wider than the chart
    std::vector<Badge> badges;
    QRectF plot;
    std::vector<Chip> legend;
    QStringList hiddenLegend;  // keys folded into "+N more"
    QRectF yTitle;             // empty without a y title
    bool yRotated = false;
    QRectF xTitle;
    int xTicks = 6, yTicks = 6;
};

ChartLayout chartLayout(const ChartLayoutInput& in);

} // namespace kapp
