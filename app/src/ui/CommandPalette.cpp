#include "ui/CommandPalette.h"

#include "core/TextMatch.h"
#include "ui/Theme.h"

#include <QApplication>
#include <QGraphicsDropShadowEffect>
#include <QKeyEvent>
#include <QLineEdit>
#include <QListWidget>
#include <QPainter>
#include <QStyledItemDelegate>
#include <QVBoxLayout>

#include <algorithm>

namespace kapp {

namespace {

class EntryDelegate : public QStyledItemDelegate {
public:
    using QStyledItemDelegate::QStyledItemDelegate;
    void paint(QPainter* p, const QStyleOptionViewItem& opt, const QModelIndex& idx) const override
    {
        const auto& t = tokens();
        p->save();
        p->setRenderHint(QPainter::Antialiasing);
        if (opt.state & QStyle::State_Selected) {
            p->setPen(Qt::NoPen);
            p->setBrush(t.accentSoft);
            p->drawRoundedRect(QRectF(opt.rect).adjusted(2, 1, -2, -1), 4, 4);
        }
        const QString kind = idx.data(Qt::UserRole).toString();
        const QString detail = idx.data(Qt::UserRole + 1).toString();
        QFont small = opt.font;
        small.setPointSizeF(small.pointSizeF() * 0.85);
        const QFontMetrics sfm(small);
        const QRect r = opt.rect.adjusted(10, 0, -10, 0);
        const int kindW = 64;
        p->setFont(small);
        p->setPen(t.faint);
        p->drawText(QRect(r.left(), r.top(), kindW, r.height()), Qt::AlignVCenter | Qt::AlignLeft, kind);
        const int detailW = detail.isEmpty() ? 0 : std::min(r.width() / 3, sfm.horizontalAdvance(detail) + 4);
        if (detailW) {
            p->setPen(t.muted);
            p->drawText(QRect(r.right() - detailW, r.top(), detailW, r.height()), Qt::AlignVCenter | Qt::AlignRight,
                        sfm.elidedText(detail, Qt::ElideLeft, detailW));
        }
        p->setFont(opt.font);
        p->setPen(t.text);
        const QRect tr(r.left() + kindW, r.top(), r.width() - kindW - detailW - 8, r.height());
        p->drawText(tr, Qt::AlignVCenter | Qt::AlignLeft, QFontMetrics(opt.font).elidedText(idx.data().toString(), Qt::ElideRight, tr.width()));
        p->restore();
    }
    QSize sizeHint(const QStyleOptionViewItem& opt, const QModelIndex& idx) const override
    {
        return {QStyledItemDelegate::sizeHint(opt, idx).width(), 30};
    }
};

} // namespace

CommandPalette::CommandPalette(QWidget* parent) : QFrame(parent)
{
    setObjectName("Palette");
    setAttribute(Qt::WA_StyledBackground);
    setAccessibleName(tr("Command palette"));
    auto* v = new QVBoxLayout(this);
    v->setContentsMargins(1, 1, 1, 6);
    v->setSpacing(0);
    m_edit = new QLineEdit(this);
    m_edit->setPlaceholderText(tr("Search commands, pages, settings and runs"));
    m_edit->setAccessibleName(tr("Command search"));
    m_list = new QListWidget(this);
    m_list->setItemDelegate(new EntryDelegate(m_list));
    m_list->setUniformItemSizes(true);
    m_list->setFocusPolicy(Qt::NoFocus);
    m_list->setAccessibleName(tr("Matches"));
    v->addWidget(m_edit);
    v->addWidget(m_list);
    auto* shadow = new QGraphicsDropShadowEffect(this);
    shadow->setBlurRadius(28);
    shadow->setOffset(0, 8);
    shadow->setColor(QColor(0, 0, 0, 110));
    setGraphicsEffect(shadow);
    connect(m_edit, &QLineEdit::textChanged, this, &CommandPalette::filter);
    connect(m_list, &QListWidget::itemClicked, this, [this](QListWidgetItem* it) { activate(m_list->row(it)); });
    m_edit->installEventFilter(this);
    parent->installEventFilter(this);
    hide();
}

void CommandPalette::place()
{
    QWidget* p = parentWidget();
    const int w = std::min(640, p->width() - 80);
    setGeometry((p->width() - w) / 2, 72, w, height());
}

void CommandPalette::open(const QString& initial)
{
    m_entries = m_source ? m_source() : std::vector<PaletteEntry>{};
    place();
    show();
    raise();
    m_edit->setText(initial);
    m_edit->selectAll();
    m_edit->setFocus();
    filter();
}

void CommandPalette::setQuery(const QString& q) { m_edit->setText(q); }

void CommandPalette::filter()
{
    const QString q = m_edit->text().trimmed();
    std::vector<std::pair<int, int>> scored;
    for (int i = 0; i < int(m_entries.size()); ++i) {
        const auto& e = m_entries[std::size_t(i)];
        int s = q.isEmpty() ? (e.kind == "Command" ? 1 : -1) : std::max(fuzzyScore(q, e.text), fuzzyScore(q, e.detail) - 200);
        if (s >= 0) scored.emplace_back(-s, i);
    }
    std::stable_sort(scored.begin(), scored.end(), [](auto a, auto b) { return a.first < b.first; });
    m_list->clear();
    m_shown.clear();
    for (std::size_t k = 0; k < scored.size() && k < 60; ++k) {
        const auto& e = m_entries[std::size_t(scored[k].second)];
        auto* it = new QListWidgetItem(e.text);
        it->setData(Qt::UserRole, e.kind);
        it->setData(Qt::UserRole + 1, e.detail);
        it->setData(Qt::AccessibleTextRole, e.kind + ": " + e.text);
        m_list->addItem(it);
        m_shown.push_back(scored[k].second);
    }
    if (m_list->count()) m_list->setCurrentRow(0);
    // Fit the list to its rows; the frame stays anchored at the top.
    const int rows = std::max(1, std::min(m_list->count(), 12));
    m_list->setFixedHeight(rows * 30 + 10);
    adjustSize();
    setFixedHeight(sizeHint().height());
}

void CommandPalette::activate(int row)
{
    if (row < 0 || row >= int(m_shown.size())) return;
    auto fn = m_entries[std::size_t(m_shown[std::size_t(row)])].run;
    hide();
    if (fn) fn();
}

bool CommandPalette::eventFilter(QObject* o, QEvent* e)
{
    if (o == parentWidget() && e->type() == QEvent::Resize && isVisible()) place();
    if (o == m_edit && e->type() == QEvent::KeyPress) {
        auto* k = static_cast<QKeyEvent*>(e);
        switch (k->key()) {
        case Qt::Key_Escape: hide(); return true;
        case Qt::Key_Down: m_list->setCurrentRow(std::min(m_list->count() - 1, m_list->currentRow() + 1)); return true;
        case Qt::Key_Up: m_list->setCurrentRow(std::max(0, m_list->currentRow() - 1)); return true;
        case Qt::Key_Return:
        case Qt::Key_Enter: activate(m_list->currentRow()); return true;
        default: break;
        }
    }
    if (o == m_edit && e->type() == QEvent::FocusOut && isVisible()) {
        // Clicking elsewhere closes the palette; clicks in the list keep it.
        QWidget* f = QApplication::focusWidget();
        if (!f || !isAncestorOf(f)) QMetaObject::invokeMethod(this, [this] { if (!m_list->underMouse()) hide(); }, Qt::QueuedConnection);
    }
    return QFrame::eventFilter(o, e);
}

} // namespace kapp
