#include "ui/FlowLayout.h"

#include <QWidget>

#include <algorithm>

namespace kapp {

FlowLayout::FlowLayout(QWidget* parent, int hSpacing, int vSpacing) : QLayout(parent), m_h(hSpacing), m_v(vSpacing)
{
    setContentsMargins(0, 0, 0, 0);
}

FlowLayout::~FlowLayout()
{
    while (QLayoutItem* item = takeAt(0)) delete item;
}

void FlowLayout::addItem(QLayoutItem* item)
{
    // Stretches and spacers have no place in a wrapping row.
    if (item->spacerItem()) {
        delete item;
        return;
    }
    m_items.append(item);
}

QLayoutItem* FlowLayout::takeAt(int index)
{
    return index >= 0 && index < m_items.size() ? m_items.takeAt(index) : nullptr;
}

QSize FlowLayout::minimumSize() const
{
    // As narrow as the widest item: everything else wraps.
    QSize s;
    for (const QLayoutItem* item : m_items)
        if (!item->isEmpty()) s = s.expandedTo(item->minimumSize().expandedTo(item->sizeHint()));
    const QMargins m = contentsMargins();
    return s + QSize(m.left() + m.right(), m.top() + m.bottom());
}

QSize FlowLayout::sizeHint() const
{
    // One line.
    int w = 0, h = 0, n = 0;
    for (const QLayoutItem* item : m_items) {
        if (item->isEmpty()) continue;
        w += item->sizeHint().width() + (n++ ? m_h : 0);
        h = std::max(h, item->sizeHint().height());
    }
    const QMargins m = contentsMargins();
    return {w + m.left() + m.right(), h + m.top() + m.bottom()};
}

void FlowLayout::setGeometry(const QRect& rect)
{
    QLayout::setGeometry(rect);
    arrange(rect, true);
}

int FlowLayout::arrange(const QRect& rect, bool apply) const
{
    const QRect r = rect.marginsRemoved(contentsMargins());
    int x = r.x(), y = r.y(), line = 0;
    QList<std::pair<QLayoutItem*, QPoint>> pending;  // the current line, centred on it when it ends
    auto endLine = [&] {
        if (apply)
            for (const auto& [item, at] : pending) item->setGeometry(QRect(at + QPoint(0, (line - item->sizeHint().height()) / 2), item->sizeHint()));
        pending.clear();
    };
    for (QLayoutItem* item : m_items) {
        if (item->isEmpty()) continue;
        const QSize s = item->sizeHint();
        if (x > r.x() && x + s.width() > r.right() + 1) {
            endLine();
            x = r.x();
            y += line + m_v;
            line = 0;
        }
        pending.append({item, QPoint(x, y)});
        x += s.width() + m_h;
        line = std::max(line, s.height());
    }
    endLine();
    return y + line - rect.y() + contentsMargins().bottom();
}

} // namespace kapp
