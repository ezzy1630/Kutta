#pragma once

// Cmd/Ctrl-K: one search box over commands, outline pages, editable settings
// and runs. Enter runs the highlighted entry; Escape closes.

#include <QFrame>

#include <functional>
#include <vector>

class QLineEdit;
class QListWidget;

namespace kapp {

struct PaletteEntry {
    QString kind;      // "Command", "Go to", "Setting", "Run"
    QString text;      // searched and shown
    QString detail;    // shown faint on the right (shortcut, path)
    std::function<void()> run;
};

class CommandPalette : public QFrame {
    Q_OBJECT
public:
    explicit CommandPalette(QWidget* parent);
    void setSource(std::function<std::vector<PaletteEntry>()> source) { m_source = std::move(source); }
    void open(const QString& initial = {});
    void setQuery(const QString& q);

protected:
    bool eventFilter(QObject* o, QEvent* e) override;

private:
    void filter();
    void activate(int row);
    void place();

    std::function<std::vector<PaletteEntry>()> m_source;
    std::vector<PaletteEntry> m_entries;
    std::vector<int> m_shown;
    QLineEdit* m_edit;
    QListWidget* m_list;
};

} // namespace kapp
