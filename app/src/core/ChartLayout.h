#pragma once

// Where each part of a 2D chart goes. Every label owns its own rectangle:
// legend rows, the y title (rotated beside the axis when the plot is tall
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
    double header = 26;    // title and badges row
    QStringList legend;    // unique legend keys, in order
    QString xLabel, yLabel;
    const Axis* x = nullptr;
    const Axis* y = nullptr;
};

struct ChartLayout {
    struct Chip {
        QRectF rect;
        int key = -1;  // index into legend; -1 is the "+N more" chip
        QString text;  // possibly elided
    };
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
