#include "ui/ConsoleView.h"

#include "core/PathDisplay.h"

#include "ui/Theme.h"

#include <QFontDatabase>
#include <QScrollBar>
#include <QTextBlock>

namespace kapp {

namespace {
constexpr int kMaxLines = 5000;

QTextCharFormat formatFor(ConsoleLevel level)
{
    const auto& t = tokens();
    QTextCharFormat f;
    switch (level) {
    case ConsoleLevel::Command: f.setForeground(t.muted); break;
    case ConsoleLevel::Info: f.setForeground(t.text); break;
    case ConsoleLevel::Table: f.setForeground(t.muted); break;
    case ConsoleLevel::Warning: f.setForeground(t.dark ? t.amber : t.amberText); break;
    case ConsoleLevel::Error: f.setForeground(t.red); break;
    case ConsoleLevel::Success: f.setForeground(t.dark ? t.accent : t.accent.darker(115)); break;
    }
    return f;
}

QString prefixed(ConsoleLevel level, const QString& text)
{
    switch (level) {
    case ConsoleLevel::Command: return "$ " + text;
    case ConsoleLevel::Warning: return "Warning: " + text;
    case ConsoleLevel::Error: return text.startsWith("  ") ? text : "Error: " + text;
    default: return text;
    }
}
} // namespace

ConsoleView::ConsoleView(QWidget* parent) : QPlainTextEdit(parent)
{
    setObjectName("Console");
    setReadOnly(true);
    setUndoRedoEnabled(false);
    // Messages are wrapped to the width when appended (see wrapped); the
    // residual table and commands never wrap, so their columns stay aligned.
    setLineWrapMode(QPlainTextEdit::NoWrap);
    setMaximumBlockCount(kMaxLines);
    QFont f = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    f.setPointSizeF(f.pointSizeF() * 0.95);
    setFont(f);
    setAccessibleName(tr("Console"));
    setPlaceholderText(tr("Commands Kutta runs and the solver's messages appear here."));
    connect(&Theme::instance(), &Theme::changed, this, &ConsoleView::restyle);
    m_rewrap.setSingleShot(true);
    m_rewrap.setInterval(150);
    connect(&m_rewrap, &QTimer::timeout, this, &ConsoleView::restyle);
}

void ConsoleView::resizeEvent(QResizeEvent* e)
{
    QPlainTextEdit::resizeEvent(e);
    const int columns = std::max(40, int(viewport()->width() / fontMetrics().horizontalAdvance(QLatin1Char('m'))) - 1);
    if (columns != m_columns) {
        m_columns = columns;
        m_rewrap.start();
    }
}

QString ConsoleView::wrapped(ConsoleLevel level, const QString& text) const
{
    if (level == ConsoleLevel::Table || level == ConsoleLevel::Command || m_columns <= 0 || text.size() <= m_columns) return text;
    // Break at spaces where possible; line separators keep one message one block.
    QString out;
    qsizetype start = 0;
    while (text.size() - start > m_columns) {
        qsizetype cut = text.lastIndexOf(' ', start + m_columns);
        if (cut <= start) cut = start + m_columns;
        out += text.mid(start, cut - start) + QChar::LineSeparator;
        start = text.at(cut) == ' ' ? cut + 1 : cut;
    }
    return out + text.mid(start);
}

void ConsoleView::append(ConsoleLevel level, const QString& rawText)
{
    // Commands and the solver's messages name files by absolute path; the
    // console shows them from "~" (still valid in a shell).
    const QString text = withHomeAsTilde(rawText);
    m_lines.emplace_back(level, text);
    if (m_lines.size() > std::size_t(kMaxLines)) m_lines.erase(m_lines.begin(), m_lines.begin() + 500);
    QScrollBar* sb = verticalScrollBar();
    const bool atEnd = sb->value() >= sb->maximum() - 4;
    QTextCursor c(document());
    c.movePosition(QTextCursor::End);
    if (!document()->isEmpty()) c.insertBlock();
    c.insertText(wrapped(level, prefixed(level, text)), formatFor(level));
    if (atEnd) sb->setValue(sb->maximum());
}

void ConsoleView::restyle()
{
    const auto lines = std::move(m_lines);
    m_lines.clear();
    clear();
    for (const auto& [level, text] : lines) append(level, text);
}

} // namespace kapp
