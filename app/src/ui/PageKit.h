#pragma once

// Small building blocks shared by inspector pages.

#include "core/PathDisplay.h"
#include "ui/Theme.h"

#include <QDir>

#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QStyle>
#include <QVBoxLayout>

namespace kapp::kit {

// Page text. A path in it shows from "~", never the home folder (see
// core/PathDisplay.h); the screenshot check in main.cpp holds every page to that.
inline QLabel* label(const QString& text, const QString& role = {}, QWidget* parent = nullptr)
{
    auto* l = new QLabel(withHomeAsTilde(text), parent);
    if (!role.isEmpty()) l->setProperty("role", role);
    l->setWordWrap(true);
    l->setTextFormat(Qt::PlainText);
    return l;
}

// Break opportunities after "_" and "-" (zero-width spaces) in a name too
// long to sit on one line of a narrow panel (over 24 characters), so it wraps
// instead of widening the page. Shorter names keep their usual breaks.
inline QString breakable(const QString& name)
{
    if (name.size() <= 24) return name;
    QString out;
    out.reserve(name.size() + 8);
    for (const QChar ch : name) {
        out += ch;
        if (ch == '_' || ch == '-') out += QChar(0x200B);
    }
    return out;
}

// A file path in a page: compact (core/PathDisplay.h), the full path in the
// tooltip. It wraps at separators, and inside a long folder or file name.
inline QLabel* pathLabel(const QString& path, const QString& role = {}, QWidget* parent = nullptr)
{
    const QString shown = compactPath(path);
    const QChar sep = QDir::separator();
    QStringList parts = shown.split(sep);
    for (QString& part : parts) part = breakable(part);
    auto* l = label(parts.join(sep), role, parent);
    l->setToolTip(QDir::toNativeSeparators(path));
    l->setAccessibleName(shown);
    return l;
}

inline QLabel* badge(const QString& text, const QString& kind, QWidget* parent = nullptr)
{
    auto* l = new QLabel(text, parent);
    l->setProperty("badge", kind);
    l->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    l->setAccessibleName(text);
    return l;
}

inline QPushButton* button(const QString& text, const QString& icon = {}, const char* kind = nullptr, QWidget* parent = nullptr)
{
    auto* b = new QPushButton(text, parent);
    if (!icon.isEmpty()) b->setIcon(Theme::instance().icon(icon, kind && QByteArray(kind) == "primary" ? tokens().accentText : QColor()));
    if (kind) b->setProperty(kind, true);
    b->setAccessibleName(text);
    return b;
}

// A tinted box for notices: "amber" (approximation), "red" (blocked), "neutral".
inline QFrame* callout(const QString& kind, QWidget* parent = nullptr)
{
    auto* f = new QFrame(parent);
    f->setProperty("callout", kind);
    auto* v = new QVBoxLayout(f);
    v->setContentsMargins(10, 8, 10, 10);
    v->setSpacing(4);
    return f;
}

inline QFrame* rule(QWidget* parent = nullptr)
{
    auto* f = new QFrame(parent);
    f->setProperty("rule", "h");
    return f;
}

// A vertical page body with standard margins.
inline QVBoxLayout* body(QWidget* page)
{
    auto* v = new QVBoxLayout(page);
    v->setContentsMargins(14, 10, 14, 16);
    v->setSpacing(8);
    return v;
}

inline QLabel* section(const QString& title, QWidget* parent = nullptr)
{
    auto* l = new QLabel(title.toUpper(), parent);
    l->setProperty("role", "section");
    l->setAccessibleName(title);
    return l;
}

inline void repolish(QWidget* w)
{
    w->style()->unpolish(w);
    w->style()->polish(w);
}

} // namespace kapp::kit
