#pragma once

// Small building blocks shared by inspector pages.

#include "ui/Theme.h"

#include <QFrame>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QStyle>
#include <QVBoxLayout>

namespace kapp::kit {

inline QLabel* label(const QString& text, const QString& role = {}, QWidget* parent = nullptr)
{
    auto* l = new QLabel(text, parent);
    if (!role.isEmpty()) l->setProperty("role", role);
    l->setWordWrap(true);
    l->setTextFormat(Qt::PlainText);
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
