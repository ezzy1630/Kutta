#pragma once

// The case's settings inventory plus the user's overrides. The inventory is
// the CLI's; the model only records which editable values differ from the
// case, validates them with the inventory's own rules and produces the
// `solve --setup` file. Every change goes through an undo stack.

#include "core/Protocol.h"

#include <QHash>
#include <QJsonObject>
#include <QObject>
#include <QUndoStack>

#include <optional>

namespace kapp {

class SettingsModel : public QObject {
    Q_OBJECT
public:
    explicit SettingsModel(QObject* parent = nullptr);

    void reset(const std::vector<Setting>& inventory);
    const std::vector<Setting>& settings() const { return m_settings; }
    const Setting* find(const QString& path) const;

    // Current value: the override if any, else the case value.
    QJsonValue value(const QString& path) const;
    bool isModified(const QString& path) const { return m_overrides.contains(path); }
    QStringList modifiedPaths() const { return m_overrides.keys(); }
    int modifiedCount() const { return int(m_overrides.size()); }

    // Validation against the inventory entry; empty string means valid.
    QString validate(const QString& path, const QJsonValue& v) const;

    // Returns the validation error, or empty when the edit was applied (undoable).
    QString setValue(const QString& path, const QJsonValue& v);
    void resetValue(const QString& path);
    void resetAll();
    // Replace overrides wholesale (restoring a run's settings); undoable.
    QString applyOverrides(const QJsonObject& overrides, const QString& description);

    // Whether the case uses the setting with the current values: a field the
    // selected method does not use (k under Intensity and Hydraulic Diameter)
    // is inactive, as in Fluent's panels.
    bool isActive(const QString& path) const;
    // Settings whose activity follows this one (the fields of a method choice).
    QStringList dependents(const QString& path) const { return m_dependents.value(path); }

    QJsonObject overrides() const;
    // Overrides that change the case being solved. Run controls (iteration
    // count, monitor interval) do not make results MODIFIED, and values of
    // fields the selected method does not use change nothing.
    QJsonObject caseOverrides() const;
    int caseModifiedCount() const { return int(caseOverrides().size()); }
    static bool isRunControl(const QString& path)
    {
        return path == QLatin1String("run/iterations") || path == QLatin1String("run/monitor-interval");
    }
    // {"v":1,"overrides":{case overrides + run controls without a CLI flag}}.
    // The iteration count travels as --iterations instead.
    QJsonObject setupJson() const;
    bool needsSetupFile() const { return !setupJson().value("overrides").toObject().isEmpty(); }

    QUndoStack* undoStack() { return &m_undo; }

    // Plain display of a value with its unit, e.g. "1.2 m/s", "second-order-upwind".
    static QString formatValue(const Setting& s, const QJsonValue& v, bool withUnit = true);
    static bool sameValue(const QJsonValue& a, const QJsonValue& b);

signals:
    void valueChanged(const QString& path);
    void overridesChanged();

private:
    friend class SetOverrideCommand;
    void setOverrideRaw(const QString& path, const std::optional<QJsonValue>& v);

    std::vector<Setting> m_settings;
    QHash<QString, int> m_index;
    QHash<QString, QStringList> m_dependents;
    QHash<QString, QJsonValue> m_overrides;
    QUndoStack m_undo;
};

} // namespace kapp
