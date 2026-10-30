#pragma once

#include "core/Session.h"

#include <QPlainTextEdit>
#include <QTimer>

namespace kapp {

// Fluent-style transcript: the commands the app ran, the CLI's messages and
// the residual table. Read-only; selectable and copyable.
class ConsoleView : public QPlainTextEdit {
    Q_OBJECT
public:
    explicit ConsoleView(QWidget* parent = nullptr);
    void append(ConsoleLevel level, const QString& text);
    QString plainText() const { return toPlainText(); }

protected:
    void resizeEvent(QResizeEvent* e) override;

private:
    void restyle();
    QString wrapped(ConsoleLevel level, const QString& text) const;
    std::vector<std::pair<ConsoleLevel, QString>> m_lines;  // for re-colouring and re-wrapping
    QTimer m_rewrap;
    int m_columns = 0;
};

} // namespace kapp
