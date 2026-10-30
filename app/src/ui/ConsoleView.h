#pragma once

#include "core/Session.h"

#include <QPlainTextEdit>

namespace kapp {

// Fluent-style transcript: the commands the app ran, the CLI's messages and
// the residual table. Read-only; selectable and copyable.
class ConsoleView : public QPlainTextEdit {
    Q_OBJECT
public:
    explicit ConsoleView(QWidget* parent = nullptr);
    void append(ConsoleLevel level, const QString& text);
    QString plainText() const { return toPlainText(); }

private:
    void restyle();
    std::vector<std::pair<ConsoleLevel, QString>> m_lines;  // for re-colouring on theme change
};

} // namespace kapp
