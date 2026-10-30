#include "core/SettingsModel.h"

#include <QJsonArray>
#include <QLocale>

#include <algorithm>
#include <cmath>

namespace kapp {

class SetOverrideCommand : public QUndoCommand {
public:
    SetOverrideCommand(SettingsModel* m, QString path, std::optional<QJsonValue> before,
                       std::optional<QJsonValue> after, const QString& text)
        : QUndoCommand(text), m_model(m), m_path(std::move(path)), m_before(std::move(before)), m_after(std::move(after))
    {
    }
    void undo() override { m_model->setOverrideRaw(m_path, m_before); }
    void redo() override { m_model->setOverrideRaw(m_path, m_after); }
    // Consecutive edits of one field (typing, spin arrows) undo as one step.
    int id() const override { return 1; }
    bool mergeWith(const QUndoCommand* other) override
    {
        const auto* o = dynamic_cast<const SetOverrideCommand*>(other);
        if (!o || o->m_path != m_path) return false;
        m_after = o->m_after;
        return true;
    }

private:
    SettingsModel* m_model;
    QString m_path;
    std::optional<QJsonValue> m_before, m_after;
};

SettingsModel::SettingsModel(QObject* parent) : QObject(parent) {}

void SettingsModel::reset(const std::vector<Setting>& inventory)
{
    m_settings = inventory;
    m_index.clear();
    m_dependents.clear();
    for (int i = 0; i < int(m_settings.size()); ++i) {
        const Setting& s = m_settings[std::size_t(i)];
        m_index.insert(s.path, i);
        if (!s.dependsOn.isEmpty()) m_dependents[s.dependsOn] << s.path;
    }
    m_overrides.clear();
    m_undo.clear();
    emit overridesChanged();
}

const Setting* SettingsModel::find(const QString& path) const
{
    auto it = m_index.find(path);
    return it == m_index.end() ? nullptr : &m_settings[std::size_t(*it)];
}

QJsonValue SettingsModel::value(const QString& path) const
{
    if (auto it = m_overrides.find(path); it != m_overrides.end()) return *it;
    const Setting* s = find(path);
    return s ? s->value : QJsonValue();
}

bool SettingsModel::isActive(const QString& path) const
{
    const Setting* s = find(path);
    if (!s) return false;
    if (s->dependsOn.isEmpty() || !find(s->dependsOn)) return s->active;
    const QJsonValue v = value(s->dependsOn);
    return std::any_of(s->activeWhen.begin(), s->activeWhen.end(), [&](const QJsonValue& w) { return sameValue(w, v); });
}

bool SettingsModel::sameValue(const QJsonValue& a, const QJsonValue& b)
{
    if (a.isDouble() && b.isDouble()) return a.toDouble() == b.toDouble();
    return a == b;
}

QString SettingsModel::validate(const QString& path, const QJsonValue& v) const
{
    const Setting* s = find(path);
    if (!s) return QStringLiteral("Kutta does not know the setting %1.").arg(path);
    if (!s->editable) return s->reason.isEmpty() ? QStringLiteral("This setting is read-only.") : s->reason;
    auto inRange = [s](double x) -> QString {
        if (!std::isfinite(x)) return QStringLiteral("Enter a finite number.");
        if (s->min && x < *s->min) return QStringLiteral("Must be at least %1.").arg(QLocale::c().toString(*s->min, 'g', 6));
        if (s->max && x > *s->max) return QStringLiteral("Must be at most %1.").arg(QLocale::c().toString(*s->max, 'g', 6));
        return {};
    };
    switch (s->kind) {
    case SettingKind::Number:
        if (!v.isDouble()) return QStringLiteral("Enter a number.");
        return inRange(v.toDouble());
    case SettingKind::Integer: {
        if (!v.isDouble()) return QStringLiteral("Enter a whole number.");
        const double x = v.toDouble();
        if (x != std::floor(x)) return QStringLiteral("Enter a whole number.");
        return inRange(x);
    }
    case SettingKind::Bool:
        return v.isBool() ? QString() : QStringLiteral("Choose on or off.");
    case SettingKind::Enum: {
        if (!v.isString() && !v.isDouble()) return QStringLiteral("Choose one of the listed options.");
        if (s->options.isEmpty() || std::any_of(s->options.begin(), s->options.end(), [&](const QJsonValue& o) { return sameValue(o, v); }))
            return {};
        QStringList names;
        for (qsizetype i = 0; i < s->options.size(); ++i)
            names << (i < s->optionLabels.size() ? s->optionLabels[i] : s->options[i].toVariant().toString());
        return QStringLiteral("Choose one of: %1.").arg(names.join(", "));
    }
    case SettingKind::Vec3: {
        const QJsonArray a = v.toArray();
        if (a.size() != 3) return QStringLiteral("Enter three components.");
        double norm2 = 0;
        for (const QJsonValue& c : a) {
            if (!c.isDouble()) return QStringLiteral("Enter three numbers.");
            if (QString e = inRange(c.toDouble()); !e.isEmpty()) return e;
            norm2 += c.toDouble() * c.toDouble();
        }
        // The CLI rejects zero directions (docs/app-protocol.md).
        if (s->leaf.endsWith("direction") && norm2 == 0) return QStringLiteral("A direction cannot be zero.");
        return {};
    }
    case SettingKind::String:
        if (!v.isString()) return QStringLiteral("Enter text.");
        if (!s->options.isEmpty() && !s->options.contains(v))
            return QStringLiteral("Choose one of the listed options.");
        return {};
    }
    return {};
}

QString SettingsModel::setValue(const QString& path, const QJsonValue& v)
{
    if (QString e = validate(path, v); !e.isEmpty()) return e;
    const Setting* s = find(path);
    std::optional<QJsonValue> before;
    if (auto it = m_overrides.find(path); it != m_overrides.end()) before = *it;
    // Setting a field back to the case value removes the override.
    std::optional<QJsonValue> after;
    if (!sameValue(v, s->value)) after = v;
    if (before == after || (before && after && sameValue(*before, *after))) return {};
    m_undo.push(new SetOverrideCommand(this, path, before, after,
                                       QStringLiteral("Change %1").arg(s->label)));
    return {};
}

void SettingsModel::resetValue(const QString& path)
{
    auto it = m_overrides.find(path);
    if (it == m_overrides.end()) return;
    const Setting* s = find(path);
    m_undo.push(new SetOverrideCommand(this, path, *it, std::nullopt,
                                       QStringLiteral("Reset %1").arg(s ? s->label : path)));
}

void SettingsModel::resetAll()
{
    if (m_overrides.isEmpty()) return;
    m_undo.beginMacro(QStringLiteral("Reset all settings to the case"));
    for (const QString& p : m_overrides.keys()) resetValue(p);
    m_undo.endMacro();
}

QString SettingsModel::applyOverrides(const QJsonObject& overrides, const QString& description)
{
    for (auto it = overrides.begin(); it != overrides.end(); ++it)
        if (QString e = validate(it.key(), it.value()); !e.isEmpty()) return it.key() + ": " + e;
    m_undo.beginMacro(description);
    for (const QString& p : m_overrides.keys())
        if (!overrides.contains(p)) resetValue(p);
    for (auto it = overrides.begin(); it != overrides.end(); ++it) setValue(it.key(), it.value());
    m_undo.endMacro();
    return {};
}

void SettingsModel::setOverrideRaw(const QString& path, const std::optional<QJsonValue>& v)
{
    if (v) m_overrides.insert(path, *v);
    else m_overrides.remove(path);
    emit valueChanged(path);
    // A method choice shows or hides its fields.
    for (const QString& d : m_dependents.value(path)) emit valueChanged(d);
    emit overridesChanged();
}

QJsonObject SettingsModel::overrides() const
{
    QJsonObject o;
    for (auto it = m_overrides.begin(); it != m_overrides.end(); ++it) o.insert(it.key(), *it);
    return o;
}

QJsonObject SettingsModel::caseOverrides() const
{
    QJsonObject o = overrides();
    for (const QString& k : o.keys())
        if (isRunControl(k) || !isActive(k)) o.remove(k);
    return o;
}

QJsonObject SettingsModel::runControls() const
{
    QJsonObject o;
    for (auto it = m_overrides.begin(); it != m_overrides.end(); ++it)
        if (isRunControl(it.key())) o.insert(it.key(), *it);
    return o;
}

QJsonObject SettingsModel::setupJsonFor(const QJsonObject& caseOverrides, const QJsonObject& runControls)
{
    QJsonObject o = caseOverrides;
    for (auto it = runControls.begin(); it != runControls.end(); ++it)
        if (it.key() != QLatin1String("run/iterations")) o.insert(it.key(), *it);
    return QJsonObject{{"v", 1}, {"overrides", o}};
}

QString SettingsModel::formatValue(const Setting& s, const QJsonValue& v, bool withUnit)
{
    QString text;
    const QLocale c = QLocale::c();
    auto fmt = [&](double x) { return c.toString(x, 'g', 8); };
    if (v.isNull() || v.isUndefined()) return QStringLiteral("—");
    switch (s.kind) {
    case SettingKind::Bool: return v.toBool() ? QStringLiteral("on") : QStringLiteral("off");
    case SettingKind::Vec3: {
        QStringList parts;
        for (const QJsonValue& x : v.toArray()) parts << (x.isDouble() ? fmt(x.toDouble()) : QStringLiteral("—"));
        text = "(" + parts.join(", ") + ")";
        break;
    }
    case SettingKind::Number:
    case SettingKind::Integer: text = v.isDouble() ? fmt(v.toDouble()) : v.toVariant().toString(); break;
    default: text = v.isString() ? v.toString() : v.toVariant().toString();
    }
    if (withUnit && !s.unit.isEmpty() && s.unit != "1" && (s.kind == SettingKind::Number || s.kind == SettingKind::Vec3 || s.kind == SettingKind::Integer))
        text += " " + s.unit;
    return text;
}

} // namespace kapp
