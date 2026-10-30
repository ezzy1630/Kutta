#pragma once

// Field identity and Fluent vocabulary. The manifest names a field
// `quantity[@role]`; the app shows Fluent's own labels ("Static Pressure").

#include <QString>
#include <QStringList>

namespace kapp {

enum class FieldRole { Kutta, Fluent, Delta, Run };  // Run: another run, or run-vs-run delta

struct FieldRef {
    QString quantity;   // "pressure", "velocity", ...
    FieldRole role = FieldRole::Kutta;
    int component = -1; // -1 scalar or magnitude; 0..2 vector component

    static FieldRef parse(const QString& manifestName);
    QString manifestName() const;  // without component
    QString key() const;           // stable id including component
    static FieldRef fromKey(const QString& key);
    bool operator==(const FieldRef&) const = default;
};

// "Static Pressure", "Velocity Magnitude", "Wall Yplus", ...
QString fluentLabel(const QString& quantity);
// Label with role and component: "Static Pressure (Fluent)", "X Velocity", "Δ Static Pressure (Kutta − Fluent)"
QString displayLabel(const FieldRef& f);
QString roleLabel(FieldRole r);
// Pretty unit for display: "m^2/s^2" -> "m²/s²", "1" -> "" (dimensionless)
QString prettyUnit(const QString& unit);
// Which quantities a Fluent user expects first in pickers.
int quantityOrder(const QString& quantity);
bool isVectorQuantity(const QString& quantity);

} // namespace kapp
