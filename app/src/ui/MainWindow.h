#pragma once

#include "core/Session.h"

#include <QMainWindow>

class QComboBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QToolButton;

namespace ads {
class CDockManager;
class CDockWidget;
}

namespace kapp {

class ChartsPanel;
class CommandPalette;
class ConsoleView;
class GraphicsArea;
class Inspector;
class Outline;
struct PaletteEntry;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(Session* session, const QStringList& samples, QWidget* parent = nullptr);
    ~MainWindow() override;

    Session* session() const { return m_session; }
    GraphicsArea* graphics() const { return m_graphics; }
    ChartsPanel* charts() const { return m_charts; }
    Outline* outline() const { return m_outline; }
    Inspector* inspector() const { return m_inspector; }
    CommandPalette* palette() const { return m_palette; }
    ConsoleView* console() const { return m_console; }

    void openCaseDialog();
    void newCaseDialog();
    void createCaseDialog();
    void saveCaseAsDialog();
    void openCase(const QString& path, const QString& data = {});
    void navigate(const QString& id);
    void showPalette(const QString& query = {});
    void setTheme(bool dark);
    void exportActiveViewport();

protected:
    void closeEvent(QCloseEvent* e) override;
    void dragEnterEvent(QDragEnterEvent* e) override;
    void dropEvent(QDropEvent* e) override;

private:
    void buildMenus();
    QWidget* buildToolbar();
    void buildDocks();
    void refreshToolbar();
    void refreshTitle();
    std::vector<PaletteEntry> paletteEntries();

    Session* m_session;
    QStringList m_samples;
    ads::CDockManager* m_docks = nullptr;
    GraphicsArea* m_graphics = nullptr;
    ChartsPanel* m_charts = nullptr;
    Outline* m_outline = nullptr;
    Inspector* m_inspector = nullptr;
    ConsoleView* m_console = nullptr;
    CommandPalette* m_palette = nullptr;
    QList<QAction*> m_commands;
    QAction* m_saveAction = nullptr;
    QAction* m_saveAsAction = nullptr;

    // Toolbar
    QLabel* m_caseLabel = nullptr;
    QLabel* m_stateLabel = nullptr;
    QPushButton* m_runButton = nullptr;
    QLineEdit* m_iterations = nullptr;
    QComboBox* m_device = nullptr;
    QToolButton* m_layoutButtons[3] = {};
    QToolButton* m_themeButton = nullptr;
    QLabel* m_cliLabel = nullptr;
    QAction* m_runAction = nullptr;
    QAction* m_stopAction = nullptr;
    QAction* m_undoAction = nullptr;
    QAction* m_redoAction = nullptr;
    QMenu* m_recentMenu = nullptr;
};

} // namespace kapp
