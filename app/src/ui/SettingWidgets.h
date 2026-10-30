#pragma once

// Editors for inventory settings. One row per setting: label, value editor,
// unit, and state (modified with reset, or locked with the CLI's reason).

#include "core/SettingsModel.h"

#include <QGridLayout>
#include <QPointer>
#include <QWidget>

class QLabel;
class QToolButton;

namespace kapp {

// Display conventions layered over the inventory: Fluent's labels for terse
// leaf names, and percent for turbulent intensity (stored as a fraction).
QString settingLabel(const Setting& s);
struct DisplayUnit { double factor = 1; QString unit; };
DisplayUnit displayUnit(const Setting& s);
QString displayValue(const Setting& s, const QJsonValue& v);  // with unit

class SettingEditor : public QWidget {
    Q_OBJECT
public:
    SettingEditor(SettingsModel* model, const Setting& s, QWidget* parent = nullptr);
    QString path() const { return m_path; }
    void refresh();
    QWidget* focusTarget() const { return m_focus; }

private:
    void commit(const QJsonValue& v);
    void showError(const QString& message);

    SettingsModel* m_model;
    QString m_path;
    QWidget* m_focus = nullptr;
    std::function<void(const QJsonValue&)> m_set;  // writes a value into the editor widgets
    std::function<std::optional<QJsonValue>(QString* error)> m_get;
    QLabel* m_error = nullptr;
};

// A titled grid of settings rows.
class SettingsForm : public QWidget {
    Q_OBJECT
public:
    SettingsForm(SettingsModel* model, QWidget* parent = nullptr);
    void addSection(const QString& title, const QString& note = {});
    // Adds a row for the path if the inventory has it; returns false otherwise.
    bool addSetting(const QString& path, const QString& labelOverride = {});
    void addRow(const QString& label, QWidget* value, const QString& unit = {});
    void addNote(const QString& text, const QString& role = "note");
    void addWidget(QWidget* w);
    void finish();
    int rows() const { return m_rowCount; }
    // Rows the selected methods use; fields of other methods are hidden, as in Fluent.
    int activeRows() const;

private:
    void refreshRow(const QString& path);
    struct Row { QWidget* marker; QLabel* label; SettingEditor* editor; QLabel* unit; QWidget* status; QToolButton* reset; QLabel* lock; };
    SettingsModel* m_model;
    QGridLayout* m_grid;
    int m_row = 0, m_rowCount = 0;
    QHash<QString, Row> m_rows;
};

} // namespace kapp
