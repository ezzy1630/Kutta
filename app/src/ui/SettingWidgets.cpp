#include "ui/SettingWidgets.h"

#include "core/Fields.h"
#include "ui/Theme.h"

#include <QCheckBox>
#include <QComboBox>
#include <QHBoxLayout>
#include <QJsonArray>
#include <QLabel>
#include <QLineEdit>
#include <QToolButton>
#include <QVBoxLayout>

#include <cmath>

namespace kapp {

namespace {

const QHash<QString, QString>& leafLabels()
{
    static const QHash<QString, QString> m = {
        {"turb-intensity", "Turbulent intensity"},
        {"turb-hydraulic-diam", "Hydraulic diameter"},
        {"turb-length-scale", "Turbulent length scale"},
        {"turb-viscosity-ratio", "Turbulent viscosity ratio"},
        {"k", "Turbulent kinetic energy"},
        {"omega", "Specific dissipation rate"},
        {"epsilon", "Turbulent dissipation rate"},
        {"k-relaxation", "Turbulent kinetic energy"},
        {"omega-relaxation", "Specific dissipation rate"},
        {"epsilon-relaxation", "Turbulent dissipation rate"},
        {"intermit-relaxation", "Intermittency"},
        {"intermittency-scheme", "Intermittency"},
        {"intermittency", "Intermittency"},
        {"courant-number", "Flow Courant number"},
        {"mom-relaxation", "Momentum"},
        {"pressure-relaxation", "Pressure"},
        {"pseudo-time-scale-factor", "Pseudo time scale factor"},
    };
    return m;
}

std::optional<double> parseNumber(QString text, double factor)
{
    text = text.trimmed();
    if (text.isEmpty()) return std::nullopt;
    bool ok = false;
    double v = QLocale::c().toDouble(text, &ok);
    if (!ok) v = QLocale().toDouble(text, &ok);
    if (!ok || !std::isfinite(v)) return std::nullopt;
    return v / factor;
}

QString numberText(const QJsonValue& v, double factor)
{
    if (!v.isDouble()) return {};
    return QString::number(v.toDouble() * factor, 'g', 10);
}

QString enumText(const QString& v)
{
    QString s = v;
    s.replace('-', ' ');
    if (!s.isEmpty()) s[0] = s[0].toUpper();
    return s;
}

} // namespace

QString settingLabel(const Setting& s)
{
    const QString leaf = s.leaf.section('/', -1);
    if (auto it = leafLabels().find(leaf); it != leafLabels().end()) return *it;
    return s.label.isEmpty() ? leaf : s.label;
}

DisplayUnit displayUnit(const Setting& s)
{
    if (s.leaf == "turb-intensity") return {100.0, "%"};
    return {1.0, prettyUnit(s.unit)};
}

QString displayValue(const Setting& s, const QJsonValue& v)
{
    const DisplayUnit du = displayUnit(s);
    QString text;
    if (v.isNull() || v.isUndefined()) return QStringLiteral("—");
    switch (s.kind) {
    case SettingKind::Number: text = QLocale().toString(v.toDouble() * du.factor, 'g', 6); break;
    case SettingKind::Integer:
        if (QString l = s.labelFor(v); !l.isEmpty()) return l;
        text = QLocale().toString(qint64(v.toDouble()));
        break;
    case SettingKind::Bool: return v.toBool() ? QObject::tr("On") : QObject::tr("Off");
    case SettingKind::Enum: {
        // Fluent's words when the CLI gives them; never a bare code.
        if (QString l = s.labelFor(v); !l.isEmpty()) return l;
        return v.isString() ? enumText(v.toString()) : SettingsModel::formatValue(s, v, false);
    }
    case SettingKind::Vec3: {
        const QJsonArray a = v.toArray();
        bool any = false;
        QStringList parts;
        for (const QJsonValue& c : a) {
            any |= c.isDouble();
            parts << (c.isDouble() ? QLocale().toString(c.toDouble(), 'g', 5) : QStringLiteral("—"));
        }
        if (!any) return QStringLiteral("—");
        text = "(" + parts.join(", ") + ")";
        break;
    }
    case SettingKind::String: return v.toString();
    }
    return du.unit.isEmpty() ? text : text + " " + du.unit;
}

SettingEditor::SettingEditor(SettingsModel* model, const Setting& s, QWidget* parent)
    : QWidget(parent), m_model(model), m_path(s.path)
{
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(0, 0, 0, 0);
    outer->setSpacing(2);
    auto* row = new QHBoxLayout;
    row->setContentsMargins(0, 0, 0, 0);
    row->setSpacing(4);
    outer->addLayout(row);
    const QString name = settingLabel(s);
    const DisplayUnit du = displayUnit(s);

    if (!s.editable) {
        auto* l = new QLabel(this);
        l->setProperty("role", "muted");
        l->setTextInteractionFlags(Qt::TextSelectableByMouse);
        l->setAccessibleName(name);
        l->setAccessibleDescription(tr("Read only: %1").arg(s.reason));
        row->addWidget(l, 1);
        m_set = [l, s](const QJsonValue& v) {
            QString text = displayValue(s, v);
            const DisplayUnit u = displayUnit(s);
            if (!u.unit.isEmpty() && text.endsWith(" " + u.unit)) text.chop(u.unit.size() + 1);
            l->setText(text);
        };
        m_focus = l;
    } else if (s.kind == SettingKind::Bool) {
        auto* c = new QCheckBox(this);
        c->setAccessibleName(name);
        row->addWidget(c, 1);
        connect(c, &QCheckBox::toggled, this, [this](bool on) { commit(on); });
        m_set = [c](const QJsonValue& v) { const QSignalBlocker b(c); c->setChecked(v.toBool()); };
        m_focus = c;
    } else if (s.kind == SettingKind::Enum && !s.options.isEmpty()) {
        auto* c = new QComboBox(this);
        c->setAccessibleName(name);
        for (qsizetype i = 0; i < s.options.size(); ++i) {
            const QJsonValue o = s.options[i];
            const QString l = s.labelFor(o);
            c->addItem(!l.isEmpty() ? l : o.isString() ? enumText(o.toString()) : o.toVariant().toString());
        }
        row->addWidget(c, 1);
        const QJsonArray options = s.options;
        connect(c, &QComboBox::activated, this, [this, options](int i) { commit(options.at(i)); });
        m_set = [c, options](const QJsonValue& v) {
            const QSignalBlocker b(c);
            int at = -1;
            for (qsizetype i = 0; i < options.size(); ++i)
                if (SettingsModel::sameValue(options[i], v)) at = int(i);
            c->setCurrentIndex(at);
        };
        m_focus = c;
    } else if (s.kind == SettingKind::Vec3) {
        std::array<QLineEdit*, 3> edits{};
        for (int i = 0; i < 3; ++i) {
            edits[std::size_t(i)] = new QLineEdit(this);
            edits[std::size_t(i)]->setAccessibleName(QStringLiteral("%1 %2").arg(name, QString("XYZ"[i])));
            edits[std::size_t(i)]->setPlaceholderText(QString("xyz"[i]));
            row->addWidget(edits[std::size_t(i)], 1);
        }
        m_set = [edits](const QJsonValue& v) {
            const QJsonArray a = v.toArray();
            for (int i = 0; i < 3; ++i) edits[std::size_t(i)]->setText(numberText(a.at(i), 1));
        };
        m_get = [edits](QString* err) -> std::optional<QJsonValue> {
            QJsonArray a;
            for (auto* e : edits) {
                auto n = parseNumber(e->text(), 1);
                if (!n) { *err = tr("Enter three numbers."); return std::nullopt; }
                a.append(*n);
            }
            return a;
        };
        for (auto* e : edits)
            connect(e, &QLineEdit::editingFinished, this, [this] {
                QString err;
                if (auto v = m_get(&err)) commit(*v); else showError(err);
            });
        m_focus = edits[0];
    } else {
        auto* e = new QLineEdit(this);
        e->setAccessibleName(name);
        row->addWidget(e, 1);
        const bool numeric = s.kind == SettingKind::Number || s.kind == SettingKind::Integer;
        if (numeric) e->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
        const SettingKind kind = s.kind;
        m_set = [e, du](const QJsonValue& v) {
            e->setText(v.isDouble() ? numberText(v, du.factor) : v.toString());
            e->setCursorPosition(0);
        };
        m_get = [e, kind, du](QString* err) -> std::optional<QJsonValue> {
            if (kind == SettingKind::String) return QJsonValue(e->text());
            auto n = parseNumber(e->text(), du.factor);
            if (!n) { *err = tr("Enter a number."); return std::nullopt; }
            if (kind == SettingKind::Integer) {
                if (std::floor(*n) != *n) { *err = tr("Enter a whole number."); return std::nullopt; }
                return QJsonValue(qint64(*n));
            }
            return QJsonValue(*n);
        };
        connect(e, &QLineEdit::editingFinished, this, [this] {
            QString err;
            if (auto v = m_get(&err)) commit(*v); else showError(err);
        });
        m_focus = e;
    }
    m_error = new QLabel(this);
    m_error->setProperty("role", "error");
    m_error->setWordWrap(true);
    m_error->hide();
    outer->addWidget(m_error);
    refresh();
}

void SettingEditor::refresh()
{
    if (m_set) m_set(m_model->value(m_path));
    const bool mod = m_model->isModified(m_path);
    if (m_focus && m_focus->property("modified").toBool() != mod) {
        m_focus->setProperty("modified", mod);
        m_focus->style()->unpolish(m_focus);
        m_focus->style()->polish(m_focus);
    }
}

void SettingEditor::commit(const QJsonValue& v)
{
    if (SettingsModel::sameValue(v, m_model->value(m_path))) {
        showError({});
        return;
    }
    const QString err = m_model->setValue(m_path, v);
    showError(err);
    if (!err.isEmpty()) refresh();
}

void SettingEditor::showError(const QString& message)
{
    m_error->setText(message);
    m_error->setVisible(!message.isEmpty());
    if (m_focus) {
        m_focus->setProperty("invalid", !message.isEmpty());
        m_focus->style()->unpolish(m_focus);
        m_focus->style()->polish(m_focus);
    }
}

// ---- SettingsForm ---------------------------------------------------------------

SettingsForm::SettingsForm(SettingsModel* model, QWidget* parent) : QWidget(parent), m_model(model)
{
    m_grid = new QGridLayout(this);
    m_grid->setContentsMargins(14, 10, 12, 14);
    m_grid->setHorizontalSpacing(8);
    m_grid->setVerticalSpacing(6);
    m_grid->setColumnStretch(2, 1);
    m_grid->setColumnMinimumWidth(0, 3);
    connect(model, &SettingsModel::valueChanged, this, &SettingsForm::refreshRow);
}

void SettingsForm::addSection(const QString& title, const QString& note)
{
    if (m_row > 0) m_grid->setRowMinimumHeight(m_row++, 8);
    auto* l = new QLabel(title.toUpper(), this);
    l->setProperty("role", "section");
    l->setAccessibleName(title);
    m_grid->addWidget(l, m_row++, 0, 1, 5);
    if (!note.isEmpty()) addNote(note);
}

void SettingsForm::addNote(const QString& text, const QString& role)
{
    auto* l = new QLabel(text, this);
    l->setProperty("role", role);
    l->setWordWrap(true);
    l->setTextFormat(Qt::PlainText);
    m_grid->addWidget(l, m_row++, 1, 1, 4);
}

void SettingsForm::addWidget(QWidget* w)
{
    w->setParent(this);
    m_grid->addWidget(w, m_row++, 1, 1, 4);
}

void SettingsForm::addRow(const QString& label, QWidget* value, const QString& unit)
{
    auto* l = new QLabel(label, this);
    l->setBuddy(value);
    m_grid->addWidget(l, m_row, 1);
    value->setParent(this);
    m_grid->addWidget(value, m_row, 2);
    if (!unit.isEmpty()) {
        auto* u = new QLabel(unit, this);
        u->setProperty("role", "unit");
        m_grid->addWidget(u, m_row, 3);
    }
    ++m_row;
    ++m_rowCount;
}

bool SettingsForm::addSetting(const QString& path, const QString& labelOverride)
{
    const Setting* s = m_model->find(path);
    if (!s) return false;
    Row r{};
    r.marker = new QWidget(this);
    r.marker->setFixedWidth(3);
    r.marker->setAccessibleName(tr("Modified"));
    r.label = new QLabel(labelOverride.isEmpty() ? settingLabel(*s) : labelOverride, this);
    r.label->setToolTip(s->path);
    r.editor = new SettingEditor(m_model, *s, this);
    r.label->setBuddy(r.editor->focusTarget());
    // A method choice in Fluent's words ("Intensity and Hydraulic Diameter")
    // is too long beside its label: as in Fluent, it sits under it, full width.
    const bool stacked = !s->optionLabels.isEmpty() || !s->valueLabel.isEmpty();
    if (stacked) {
        m_grid->addWidget(r.label, m_row++, 1, 1, 3, Qt::AlignVCenter | Qt::AlignLeft);
        m_grid->addWidget(r.marker, m_row, 0);
        m_grid->addWidget(r.editor, m_row, 1, 1, 3);
    } else {
        m_grid->addWidget(r.marker, m_row, 0);
        m_grid->addWidget(r.label, m_row, 1, Qt::AlignVCenter | Qt::AlignLeft);
        m_grid->addWidget(r.editor, m_row, 2);
    }
    if (!stacked) {  // a method choice has no unit
        auto* unit = new QLabel(displayUnit(*s).unit, this);
        unit->setProperty("role", "unit");
        unit->setMinimumWidth(QFontMetrics(font()).horizontalAdvance("kg/m³") + 2);
        m_grid->addWidget(unit, m_row, 3);
        r.unit = unit;
    }
    auto* status = new QWidget(this);
    r.status = status;
    auto* sl = new QHBoxLayout(status);
    sl->setContentsMargins(0, 0, 0, 0);
    sl->setSpacing(0);
    status->setFixedWidth(24);
    if (s->editable) {
        r.reset = new QToolButton(status);
        r.reset->setIcon(Theme::instance().icon("rotate-ccw"));
        r.reset->setIconSize(QSize(14, 14));
        r.reset->setAccessibleName(tr("Reset %1 to the case value").arg(r.label->text()));
        const QString p = path;
        connect(r.reset, &QToolButton::clicked, this, [this, p] { m_model->resetValue(p); });
        sl->addWidget(r.reset);
    } else {
        r.lock = new QLabel(status);
        r.lock->setPixmap(Theme::instance().pixmap("lock", 13, tokens().faint, devicePixelRatioF()));
        r.lock->setToolTip(s->reason.isEmpty() ? tr("Read only") : s->reason);
        r.lock->setAccessibleName(tr("Read only: %1").arg(s->reason));
        sl->addWidget(r.lock);
    }
    m_grid->addWidget(status, m_row, 4);
    m_rows.insert(path, r);
    ++m_row;
    ++m_rowCount;
    refreshRow(path);
    return true;
}

void SettingsForm::refreshRow(const QString& path)
{
    auto it = m_rows.find(path);
    if (it == m_rows.end()) return;
    Row& r = *it;
    const bool mod = m_model->isModified(path);
    const Setting* s = m_model->find(path);
    r.editor->refresh();
    r.marker->setStyleSheet(mod ? QStringLiteral("background:%1; border-radius:1px;").arg(tokens().amber.name()) : QString());
    if (r.reset) {
        r.reset->setVisible(mod);
        r.reset->setToolTip(tr("Reset to the case value (%1)").arg(s ? displayValue(*s, s->value) : QString()));
    }
    r.label->setToolTip(mod && s ? tr("%1\nModified. Case value: %2").arg(path, displayValue(*s, s->value)) : path);
    const bool active = m_model->isActive(path);
    for (QWidget* w : {r.marker, static_cast<QWidget*>(r.label), static_cast<QWidget*>(r.editor), static_cast<QWidget*>(r.unit), r.status})
        if (w) w->setVisible(active);
}

int SettingsForm::activeRows() const
{
    int n = 0;
    for (auto it = m_rows.begin(); it != m_rows.end(); ++it) n += m_model->isActive(it.key());
    return n;
}

void SettingsForm::finish()
{
    m_grid->setRowStretch(m_row, 1);
}

} // namespace kapp
