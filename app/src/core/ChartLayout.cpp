#include "core/ChartLayout.h"

#include <QCoreApplication>
#include <QFontMetricsF>

#include <algorithm>

namespace kapp {

ChartLayout chartLayout(const ChartLayoutInput& in)
{
    ChartLayout l;
    const QRectF& r = in.rect;
    const QFontMetricsF fm(in.font);
    const double lineH = fm.height();
    const double left0 = r.left() + 10, right0 = r.right() - 10;

    // Legend chips wrap below the header. Long names are elided, and rows that
    // would squeeze the plot below a readable height fold into "+N more".
    const double chipH = lineH + 6;
    const double maxChipText = std::max(60.0, (right0 - left0) * 0.6 - 26);
    const double xAxisBlock = lineH + 6 + (in.xLabel.isEmpty() ? 2 : lineH + 4);
    const double minPlotH = 48;
    const int maxRows = std::max(1, int((r.height() - in.header - xAxisBlock - minPlotH) / chipH));
    double x = left0, y = r.top() + in.header;
    int row = 0;
    const auto more = [&](qsizetype n) { return QCoreApplication::translate("kapp::PlotWidget", "+%1 more").arg(n); };
    for (qsizetype k = 0; k < in.legend.size(); ++k) {
        const QString text = fm.elidedText(in.legend[k], Qt::ElideMiddle, maxChipText);
        const double w = fm.horizontalAdvance(text) + 26;
        if (x + w > right0 && x > left0) {
            x = left0;
            y += chipH;
            ++row;
        }
        // The last allowed row keeps room for the "+N more" chip.
        const bool overflow = row >= maxRows;
        const bool squeezed = row == maxRows - 1 && k + 1 < in.legend.size() && x + w + fm.horizontalAdvance(more(in.legend.size() - k)) + 14 > right0;
        if (overflow || squeezed) {
            if (overflow) {
                y -= chipH;
                row = maxRows - 1;
                // Start after the last chip on the final allowed row.
                x = l.legend.empty() ? left0 : l.legend.back().rect.right() + 6;
            }
            l.hiddenLegend = in.legend.mid(k);
            const QString t = more(l.hiddenLegend.size());
            const double mw = std::min(fm.horizontalAdvance(t) + 8, right0 - x);
            l.legend.push_back({QRectF(x, y, std::max(0.0, mw), chipH), -1, t});
            break;
        }
        l.legend.push_back({QRectF(x, y, std::min(w, right0 - x), chipH), int(k), text});
        x += w + 6;
    }
    const double legendBottom = l.legend.empty() ? r.top() + in.header : y + chipH;
    double top = legendBottom + 6;
    const double bottom = r.bottom() - xAxisBlock;

    double gutter = 0;
    if (!in.yLabel.isEmpty()) {
        if (bottom - top >= fm.horizontalAdvance(in.yLabel) + 8) {
            l.yRotated = true;
            gutter = lineH + 6;
        } else {
            l.yTitle = QRectF(left0, top, right0 - left0, lineH);
            top += lineH + 4;
        }
    }
    const double plotH = std::max(1.0, bottom - top);
    l.yTicks = std::clamp(int(plotH / (lineH * 2.2)), 2, 8);
    double yTickWidth = 0;
    if (in.y)
        for (const Tick& tk : in.y->ticks(l.yTicks))
            if (tk.major) yTickWidth = std::max(yTickWidth, fm.horizontalAdvance(tk.label));
    const double plotLeft = r.left() + 8 + gutter + yTickWidth + 6;
    const double plotRight = r.right() - 14;
    l.plot = QRectF(QPointF(plotLeft, top), QPointF(std::max(plotLeft + 1, plotRight), std::max(top + 1, bottom)));
    if (l.yRotated) l.yTitle = QRectF(r.left() + 6, l.plot.top(), lineH, l.plot.height());

    // As many x ticks as their labels allow without touching.
    double widest = 0;
    if (in.x)
        for (const Tick& tk : in.x->ticks(6)) widest = std::max(widest, fm.horizontalAdvance(tk.label));
    l.xTicks = std::clamp(int(l.plot.width() / (widest + 18)), 2, 10);
    if (!in.xLabel.isEmpty()) l.xTitle = QRectF(l.plot.left(), r.bottom() - lineH - 4, l.plot.width(), lineH);
    return l;
}

} // namespace kapp
