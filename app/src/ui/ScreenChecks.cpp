#include "ui/ScreenChecks.h"

#include <QAbstractButton>
#include <QAbstractItemView>
#include <QAbstractScrollArea>
#include <QAbstractSpinBox>
#include <QComboBox>
#include <QDir>
#include <QLabel>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QTabBar>
#include <QTextEdit>

#include <algorithm>

namespace kapp {

namespace {

// Controls whose insides are their own business (a spin box's line edit, a
// combo box's button): only the outermost is checked.
bool isControl(const QWidget* w)
{
    return qobject_cast<const QLabel*>(w) || qobject_cast<const QAbstractButton*>(w) || qobject_cast<const QLineEdit*>(w) ||
           qobject_cast<const QComboBox*>(w) || qobject_cast<const QAbstractSpinBox*>(w);
}

QString describe(const QWidget* w)
{
    QString text;
    if (auto* l = qobject_cast<const QLabel*>(w)) text = l->text();
    else if (auto* b = qobject_cast<const QAbstractButton*>(w)) text = b->text();
    else if (auto* c = qobject_cast<const QComboBox*>(w)) text = c->currentText();
    if (text.isEmpty()) text = w->accessibleName();
    text.remove(QChar(0x200B));
    if (text.size() > 48) text = text.left(47) + "…";
    return QStringLiteral("%1 \"%2\"").arg(QString::fromLatin1(w->metaObject()->className()), text);
}

// The controls that are visible under root, outside scrolled content of
// their own (a table's cells scroll; the page itself is root).
QList<QWidget*> visibleControls(QWidget* root)
{
    QList<QWidget*> out;
    for (QWidget* w : root->findChildren<QWidget*>()) {
        if (!isControl(w) || !w->isVisibleTo(root)) continue;
        bool nested = false;
        for (QWidget* p = w->parentWidget(); p && p != root; p = p->parentWidget())
            if (isControl(p) || qobject_cast<QAbstractScrollArea*>(p)) nested = true;
        if (!nested) out << w;
    }
    return out;
}

} // namespace

QStringList layoutProblems(QWidget* root)
{
    QStringList problems;
    const QList<QWidget*> controls = visibleControls(root);
    struct Placed {
        QWidget* w;
        QRect r;  // in root's coordinates
    };
    std::vector<Placed> placed;
    for (QWidget* w : controls) {
        const QSize hint = w->minimumSizeHint();
        // An explicit minimum overrides the hint, as it does for layouts.
        const int minW = w->minimumWidth() > 0 ? w->minimumWidth() : hint.width();
        if (w->width() + 1 < minW) problems << QStringLiteral("%1 is %2 px wide, needs %3").arg(describe(w)).arg(w->width()).arg(minW);
        auto* l = qobject_cast<QLabel*>(w);
        if (l && l->wordWrap() && l->hasHeightForWidth()) {
            const int need = l->heightForWidth(l->width());
            if (l->height() + 1 < need) problems << QStringLiteral("%1 is %2 px tall, its text needs %3").arg(describe(w)).arg(l->height()).arg(need);
        } else if (w->height() + 1 < hint.height() && w->minimumHeight() == 0) {
            problems << QStringLiteral("%1 is %2 px tall, needs %3").arg(describe(w)).arg(w->height()).arg(hint.height());
        }
        // Every holder up to root must contain it, or part of it is cut off.
        QRect r = w->geometry();
        for (QWidget* p = w->parentWidget(); p; p = p->parentWidget()) {
            if (!p->rect().adjusted(-1, -1, 1, 1).contains(r)) {
                problems << QStringLiteral("%1 extends past %2 (%3 px wide)").arg(describe(w), QString::fromLatin1(p->metaObject()->className())).arg(p->width());
                break;
            }
            if (p == root) break;
            r.translate(p->pos());
        }
        placed.push_back({w, QRect(w->mapTo(root, QPoint(0, 0)), w->size())});
    }
    // Overlaps, by a sweep down the page.
    std::sort(placed.begin(), placed.end(), [](const Placed& a, const Placed& b) { return a.r.top() < b.r.top(); });
    for (std::size_t i = 0; i < placed.size(); ++i)
        for (std::size_t j = i + 1; j < placed.size() && placed[j].r.top() < placed[i].r.bottom(); ++j) {
            const QRect both = placed[i].r.intersected(placed[j].r);
            if (both.width() > 1 && both.height() > 1)
                problems << QStringLiteral("%1 overlaps %2").arg(describe(placed[i].w), describe(placed[j].w));
        }
    return problems;
}

QStringList homePathsOnScreen(QWidget* root)
{
    const QString home = QDir::homePath();
    const QString nativeHome = QDir::toNativeSeparators(home);
    QStringList found;
    auto check = [&](const QWidget* w, const QString& text) {
        if (home.size() > 1 && (text.contains(home) || text.contains(nativeHome))) {
            QString line = text.section(nativeHome, 0, 0).right(30) + nativeHome + text.section(nativeHome, 1).left(30);
            found << QStringLiteral("%1: …%2…").arg(QString::fromLatin1(w->metaObject()->className()), line.simplified());
        }
    };
    QList<QWidget*> widgets = root->findChildren<QWidget*>();
    widgets.prepend(root);
    for (QWidget* w : widgets) {
        if (!w->isVisible()) continue;
        check(w, w->windowTitle());
        if (auto* l = qobject_cast<QLabel*>(w)) check(w, l->text());
        else if (auto* b = qobject_cast<QAbstractButton*>(w)) check(w, b->text());
        else if (auto* e = qobject_cast<QLineEdit*>(w)) check(w, e->text());
        else if (auto* p = qobject_cast<QPlainTextEdit*>(w)) check(w, p->toPlainText());
        else if (auto* t = qobject_cast<QTextEdit*>(w)) check(w, t->toPlainText());
        else if (auto* c = qobject_cast<QComboBox*>(w)) check(w, c->currentText());
        else if (auto* tb = qobject_cast<QTabBar*>(w)) {
            for (int i = 0; i < tb->count(); ++i) check(w, tb->tabText(i));
        } else if (auto* v = qobject_cast<QAbstractItemView*>(w); v && v->model()) {
            const QAbstractItemModel* m = v->model();
            for (int r = 0; r < std::min(m->rowCount(), 2000); ++r)
                for (int col = 0; col < m->columnCount(); ++col) check(w, m->index(r, col).data().toString());
        }
    }
    return found;
}

} // namespace kapp
