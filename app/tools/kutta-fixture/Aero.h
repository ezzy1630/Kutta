#pragma once

// The fixture CLI's `aero` and convergence checking, following
// docs/app-protocol.md so the app's pages, scenarios and tests run without the
// solver. The real CLI's implementation is the reference; app/tests checks
// this one against the same hand-computed cases as tests/ does for the CLI.

#include <QJsonArray>
#include <QJsonObject>
#include <QMap>
#include <QString>

#include <optional>
#include <vector>

namespace fixture {

// loads.json + spec -> aero.json object; empty error on success.
QJsonObject aeroFromLoads(const QJsonObject& loads, const QJsonObject& spec, QString* error);

struct ReportCondition { QString report; double criterion = 1e-3; int previous = 5, ignoreBefore = 100; };
struct ConvergenceConfig {
    bool enabled = false;
    bool any = false;
    QMap<QString, double> residuals;  // checked equations
    std::vector<ReportCondition> reports;
    QJsonObject toJson() const;       // the setup event's convergence object
};
// From the case inventory's values with the setup overrides applied.
ConvergenceConfig convergenceFrom(const QJsonArray& settings, const QJsonObject& overrides);

class ConvergenceCheck {
public:
    explicit ConvergenceCheck(ConvergenceConfig c) : m_config(std::move(c)) {}
    // One monitor sample; returns the iteration event's convergence object.
    QJsonObject sample(int iter, const QMap<QString, double>& residuals, const QMap<QString, double>& reportForces);
    bool met() const { return m_met; }
    QStringList metBy() const { return m_metBy; }

private:
    ConvergenceConfig m_config;
    QMap<QString, std::vector<double>> m_history;  // values after ignore-before
    bool m_met = false;
    QStringList m_metBy;
};

} // namespace fixture
