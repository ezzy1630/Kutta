#pragma once

// A design study: variants of one case that differ in one editable setting,
// run back to back on the one GPU. Each variant is an ordinary run (kept in
// Run History, MODIFIED when it changes the case); the study records which run
// belongs to which value. Studies live beside the case's runs, in
// <case dir>/studies/<id>.json.

#include "core/SettingsModel.h"

#include <QDateTime>
#include <QJsonObject>

#include <optional>
#include <vector>

namespace kapp {

enum class VariantStatus { Queued, Running, Done, Failed, Stopped, Cancelled };
QString variantStatusLabel(VariantStatus s);

struct StudyVariant {
    double value = 0;          // the parameter's value
    QJsonObject overrides;     // case overrides this variant adds to the study's base
    QString runId;             // empty until it starts
    VariantStatus status = VariantStatus::Queued;
};

struct Study {
    QString id;                // "study-0002"
    int number = 0;
    QString name;
    QString parameter;         // the varied setting's path
    QString parameterLabel, unit;
    QStringList followers;     // settings set to the same value (reference velocity following the inlet)
    QJsonObject base;          // case overrides every variant shares (the settings when it was made)
    QJsonObject runControls;   // iterations and convergence criteria, the same for every variant
    QStringList approximations;
    std::vector<StudyVariant> variants;
    QDateTime created;

    QString title() const { return name.isEmpty() ? QStringLiteral("Study %1").arg(number) : name; }
    // The case overrides a variant runs with: the base, then the variant's own.
    QJsonObject overridesFor(const StudyVariant& v) const;
    int nextQueued() const;    // index, or -1
    int indexOfRun(const QString& runId) const;
    bool finished() const { return nextQueued() < 0; }
    QJsonObject toJson() const;
    static std::optional<Study> fromJson(const QJsonObject& o);
};

// "16, 20, 24" or "16 to 28 step 4" (both ends included). Empty on error, with why.
std::vector<double> parseSweepValues(const QString& text, QString* error = nullptr);

// Settings a study may vary: editable numbers that change the case. Run
// controls (iterations, convergence) are the same for every variant.
bool studyParameterAllowed(const Setting& s);

// When the case takes its reference velocity from an inlet (equal values),
// a sweep of that inlet's speed must move the reference with it, or every
// coefficient would be scaled by (V/V_ref)^2. Returns the paths to follow.
QStringList referenceFollowers(const SettingsModel& model, const QString& parameter);

// A one-parameter sweep. Every value is validated against the inventory; a
// value equal to the case's leaves that setting at the case, so the variant
// is only MODIFIED by what really differs.
Result<Study> makeSweep(const SettingsModel& model, const QString& parameter, const std::vector<double>& values,
                        const QStringList& followers);

// Variant bookkeeping when its run starts and finishes.
void markStarted(Study& s, int index, const QString& runId);
void markFinished(Study& s, const QString& runId, VariantStatus status);
// Stop: the running variant is stopped by its run; queued ones are cancelled.
void cancelQueued(Study& s);
// Run the cancelled and failed variants again.
void requeue(Study& s);

class StudyStore {
public:
    void setDir(const QString& caseDir);  // <case dir>/studies
    std::vector<Study> load() const;
    bool save(const Study& s) const;
    bool remove(const QString& id) const;
    int nextNumber() const;

private:
    QString m_dir;
};

} // namespace kapp
