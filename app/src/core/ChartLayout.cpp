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

    // Header: the title, then its badges on the same row when all of them fit
    // beside it, else wrapped onto rows below it. The legend starts below
    // whichever ends lower, and never above in.header.
    const QFontMetricsF tfm(in.titleFont), bfm(in.badgeFont);
    const double avail = std::max(1.0, right0 - left0), top0 = r.top() + 8;
    const double badgeH = bfm.height() + 3;
    double headerBottom = top0;
    if (!in.title.isEmpty()) {
        l.titleText = tfm.elidedText(in.title, Qt::ElideRight, avail);
        l.title = QRectF(left0, top0, tfm.horizontalAdvance(l.titleText), tfm.height());
        headerBottom = l.title.bottom();
    }
    std::vector<double> widths;
    double badgesW = 0;
    for (const QString& b : in.badges) {
        widths.push_back(std::min(bfm.horizontalAdvance(b) + 12, avail));
        badgesW += widths.back() + (widths.size() > 1 ? 6 : 0);
    }
    const bool besideTitle = l.title.isEmpty() || l.title.width() + 10 + badgesW <= avail;
    double bx = l.title.isEmpty() ? left0 : besideTitle ? l.title.right() + 10 : left0;
    double by = besideTitle ? top0 : l.title.bottom() + 4;
    for (qsizetype i = 0; i < in.badges.size(); ++i) {
        const double w = widths[std::size_t(i)];
        if (bx + w > right0 + 0.5 && bx > left0) {
            bx = left0;
            by += badgeH + 4;
        }
        l.badges.push_back({QRectF(bx, by, w, badgeH), bfm.elidedText(in.badges[i], Qt::ElideRight, w - 12)});
        headerBottom = std::max(headerBottom, by + badgeH);
        bx += w + 6;
    }
    const double header = std::max(in.header, headerBottom + 2 - r.top());

    // Legend chips wrap below the header. Long names are elided, and rows that
    // would squeeze the plot below a readable height fold into "+N more".
    const double chipH = lineH + 6;
    const double maxChipText = std::max(60.0, (right0 - left0) * 0.6 - 26);
    const double xAxisBlock = lineH + 6 + (in.xLabel.isEmpty() ? 2 : lineH + 4);
    const double minPlotH = 48;
    const int maxRows = std::max(1, int((r.height() - header - xAxisBlock - minPlotH) / chipH));
    double x = left0, y = r.top() + header;
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
    const double legendBottom = l.legend.empty() ? r.top() + header : y + chipH;
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
