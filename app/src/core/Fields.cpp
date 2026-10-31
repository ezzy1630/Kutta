#include "core/Fields.h"

#include <QHash>

namespace kapp {

FieldRef FieldRef::parse(const QString& name)
{
    FieldRef f;
    const qsizetype at = name.indexOf('@');
    f.quantity = at < 0 ? name : name.left(at);
    const QString role = at < 0 ? QString() : name.mid(at + 1);
    if (role == "fluent") f.role = FieldRole::Fluent;
    else if (role == "delta") f.role = FieldRole::Delta;
    else if (role == "run") f.role = FieldRole::Run;
    return f;
}

QString FieldRef::manifestName() const
{
    switch (role) {
    case FieldRole::Kutta: return quantity;
    case FieldRole::Fluent: return quantity + "@fluent";
    case FieldRole::Delta: return quantity + "@delta";
    case FieldRole::Run: return quantity + "@run";
    }
    return quantity;
}

QString FieldRef::key() const
{
    return component < 0 ? manifestName() : manifestName() + "#" + QString::number(component);
}

FieldRef FieldRef::fromKey(const QString& key)
{
    const qsizetype hash = key.indexOf('#');
    FieldRef f = parse(hash < 0 ? key : key.left(hash));
    if (hash >= 0) f.component = key.mid(hash + 1).toInt();
    return f;
}

QString fluentLabel(const QString& q)
{
    static const QHash<QString, QString> labels{
        {"pressure", "Static Pressure"},
        {"pressure_coefficient", "Pressure Coefficient"},
        {"total_pressure_coefficient", "Total Pressure Coefficient"},
        {"wall_shear", "Wall Shear Stress"},
        {"wall_shear_magnitude", "Wall Shear Stress Magnitude"},
        {"y_plus", "Wall Yplus"},
        {"y_star", "Wall Ystar"},
        {"epsilon", "Turbulent Dissipation Rate (Epsilon)"},
        {"velocity", "Velocity"},
        {"velocity_magnitude", "Velocity Magnitude"},
        {"k", "Turbulent Kinetic Energy"},
        {"omega", "Specific Dissipation Rate"},
        {"viscosity_ratio", "Turbulent Viscosity Ratio"},
        {"vorticity_magnitude", "Vorticity Magnitude"},
        {"vorticity", "Vorticity"},
        {"total_pressure", "Total Pressure"},
        {"velocity_x", "X Velocity"},
        {"velocity_y", "Y Velocity"},
        {"velocity_z", "Z Velocity"},
        {"q_criterion", "Q Criterion"},
        {"helicity", "Helicity"},
        {"skin_friction", "Skin Friction Coefficient"},
        {"skin_friction_vector", "Skin Friction"},
        {"arc_length", "Arc Length"},
    };
    if (auto it = labels.find(q); it != labels.end()) return *it;
    QString s = q;
    s.replace('_', ' ');
    if (!s.isEmpty()) s[0] = s[0].toUpper();
    return s;
}

QString roleLabel(FieldRole r)
{
    switch (r) {
    case FieldRole::Kutta: return QStringLiteral("Kutta");
    case FieldRole::Fluent: return QStringLiteral("Fluent");
    case FieldRole::Delta: return QStringLiteral("Kutta − Fluent");
    case FieldRole::Run: return QStringLiteral("Run A − Run B");
    }
    return {};
}

QString displayLabel(const FieldRef& f)
{
    QString base = fluentLabel(f.quantity);
    if (f.component >= 0) {
        static const char* axes[] = {"X", "Y", "Z"};
        if (f.quantity == "velocity") base = QStringLiteral("%1 Velocity").arg(axes[f.component]);
        else if (f.quantity == "wall_shear") base = QStringLiteral("%1 Wall Shear Stress").arg(axes[f.component]);
        else base = QStringLiteral("%1 (%2)").arg(base, axes[f.component]);
    } else if (isVectorQuantity(f.quantity)) {
        base += QStringLiteral(" Magnitude");
    }
    switch (f.role) {
    case FieldRole::Kutta: return base;
    case FieldRole::Fluent: return base + QStringLiteral(" · Fluent");
    case FieldRole::Delta: return QStringLiteral("Δ ") + base + QStringLiteral(" · Kutta − Fluent");
    case FieldRole::Run: return QStringLiteral("Δ ") + base + QStringLiteral(" · Run A − Run B");
    }
    return base;
}

QString prettyUnit(const QString& unit)
{
    if (unit == "1" || unit.isEmpty()) return {};
    QString u = unit;
    u.replace("^-1", "⁻¹").replace("^2", "²").replace("^3", "³").replace(" ", "·");
    return u;
}

int quantityOrder(const QString& q)
{
    static const QStringList order{"pressure", "pressure_coefficient", "total_pressure", "total_pressure_coefficient",
                                   "velocity_magnitude", "velocity", "velocity_x", "velocity_y", "velocity_z",
                                   "wall_shear_magnitude", "wall_shear", "skin_friction", "skin_friction_vector", "y_plus", "y_star",
                                   "k", "omega", "epsilon", "viscosity_ratio", "vorticity_magnitude", "vorticity", "q_criterion", "helicity"};
    const qsizetype i = order.indexOf(q);
    return i < 0 ? 1000 : int(i);
}

bool isVectorQuantity(const QString& q)
{
    return q == "velocity" || q == "wall_shear" || q == "vorticity" || q == "skin_friction_vector";
}

} // namespace kapp
