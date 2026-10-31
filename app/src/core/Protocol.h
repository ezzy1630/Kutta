#pragma once

// Typed view of docs/app-protocol.md (v1). Everything the CLI prints is parsed
// here once; the rest of the app never touches raw JSON from the solver.

#include "core/Result.h"
#include "core/ViewFormat.h"

#include <QJsonArray>
#include <QJsonObject>
#include <QJsonValue>
#include <QMap>
#include <QString>
#include <QStringList>

#include <optional>
#include <variant>
#include <vector>

namespace kapp {

struct Device {
    QString backend;  // "cuda" | "metal" | "cpu"
    QString name;
    qint64 memoryBytes = 0;
    int ordinal = 0;
    bool isDefault = false;
    QString label() const;  // "Apple M4 Pro · Metal · 17.8 GB"
};
Result<std::vector<Device>> parseDevices(const QJsonObject& json);

enum class IssueClass { Unsupported, Approximable, Informational };
enum class SupportStatus { Supported, Approximate, Unsupported };

struct SupportIssue {
    QString path, setting, label;
    IssueClass cls = IssueClass::Informational;
    QString approximation;  // name to pass to --approximate when approximable
    QString explanation, evidence;
    QString measuredEffect;  // approximable: the measured effect alone (newer CLIs)
    QString shortLabel;      // approximable: a few words for badges ("WFGC off"; newer CLIs)
    // The explanation without the measured effect the CLI appends to it.
    QString explanationOnly() const;
};

struct Zone {
    int id = -1;
    QString name, type;
    qint64 faces = 0;
    QJsonObject condition;
};

struct Material { QString name; std::optional<double> density, viscosity; };

struct ReportDef {
    QString name, kind;  // kind: "drag", "lift", "downforce", "force", ...
    QStringList zones;
    std::optional<QVector3D> direction;
};

enum class SettingKind { Number, Integer, Enum, Bool, Vec3, String };

struct Setting {
    QString path, label, group, unit, reason;
    SettingKind kind = SettingKind::String;
    QJsonValue value;     // null when unavailable
    QJsonArray options;   // enum choices (names or integer codes); empty when unrestricted
    QStringList optionLabels;  // Fluent's words for each option, parallel to options
    QString valueLabel;   // Fluent's words for the case value, when it is a code
    std::optional<double> min, max;
    bool editable = false;
    // Whether the case uses the setting. With dependsOn, it is used while that
    // enum's current value is one of activeWhen (Fluent shows only those fields).
    bool active = true;
    QString dependsOn;
    QJsonArray activeWhen;

    // Fluent's words for one value of this setting; empty when it has none.
    QString labelFor(const QJsonValue& v) const;

    // Derived from the path: "boundary/<zone>/velocity-magnitude" -> owner "<zone>".
    QString section;  // first segment: boundary, materials, controls, case, ...
    QString owner;    // decoded zone/material/report name, if any
    QString leaf;     // remaining path after owner
    bool raw = false; // retained Fluent setting without a Kutta control
};

// "%49nlet%20%41" -> "Inlet A"
QString decodePathSegment(const QString& segment);

struct CaseSummary {
    QString filePath, sha256, fluentVersion;
    qint64 cells = 0, faces = 0, nodes = 0;
    QMap<QString, qint64> cellTypes;
    Bounds bounds;
    std::vector<Zone> zones;
    QString viscousModel;
    QString transition;  // models.transition: "gamma-transport-eqn" etc.; empty when none
    std::optional<bool> energy, energyCoupled;
    std::vector<Material> materials;
    std::vector<ReportDef> reports;
    SupportStatus status = SupportStatus::Supported;
    std::vector<SupportIssue> issues;
    std::vector<Setting> settings;
    bool dataPresent = false;
    QStringList dataFields;
    std::optional<int> iterationsHint;
    bool synthetic = false;  // fixture data invented for demonstration
    std::optional<int> dimension;  // general/dimension: 2 for planar cases

    // A 2D case, which Kutta solves as one cell through the reference depth.
    bool planar() const { return dimension == 2; }
    QStringList approximations() const;  // names a run may opt into
    // Badge words for an approximation: the CLI's short label, else its name.
    QString approximationShort(const QString& name) const;
    std::vector<const SupportIssue*> blocking(const QStringList& optedIn) const;
    const Zone* zone(const QString& name) const;
};
Result<CaseSummary> parseCaseSummary(const QJsonObject& json);

// ---- solve --progress jsonl -------------------------------------------------------

struct ReportValue { std::optional<double> force, coefficient; };

// Convergence criteria as the solve applies them (setup event, done, summary).
struct ReportCriterion { QString report; double stopCriterion = 0; int previousValues = 0, ignoreBefore = 0; };
struct ConvergenceSetup {
    bool enabled = false;
    QString condition = "all";              // "all" | "any"
    QMap<QString, double> residuals;        // checked equations and their absolute criteria
    std::vector<ReportCriterion> reports;
    bool hasConditions() const { return !residuals.isEmpty() || !reports.empty(); }
};
// One condition's state at a monitor sample.
struct ConditionState {
    QString kind;                           // "residuals" | "report"
    QString report;                         // for reports
    bool met = false;
    std::optional<double> value, criterion; // reports: lag residual and stop criterion
    struct Equation { double value = 0, criterion = 0; bool met = false; };
    QMap<QString, Equation> equations;      // residuals
    QString id() const { return kind == "report" ? "report:" + report : kind; }
};
struct ConvergenceState { bool met = false; std::vector<ConditionState> conditions; };
ConvergenceSetup parseConvergenceSetup(const QJsonObject& o);
ConvergenceState parseConvergenceState(const QJsonObject& o);
// "drag changed less than 0.1 % over 5 values"; "every checked residual reached its criterion".
QString describeCondition(const QString& id, const ConvergenceSetup& setup);

struct SetupEvent {
    qint64 cells = 0;
    QString backend, device, init;
    QStringList approximations;
    QJsonObject overrides;
    bool modified = false;
    QString liveRunId;  // live_view.run_id: snapshots from any other run are stale
    ConvergenceSetup convergence;
    QString coupling;   // pressure_velocity_coupling that actually runs; empty from older CLIs
};
struct IterationEvent {
    int iter = 0;
    QMap<QString, double> residuals;  // ordered by name; display order is the app's
    QMap<QString, ReportValue> reports;
    std::optional<double> secPerIter, elapsed;
    std::optional<ConvergenceState> convergence;  // only when checking is on
};
struct WarningEvent { QString message; };
struct DoneEvent {
    bool converged = false;
    int iterations = 0;
    double elapsed = 0;
    QString out, stopReason;
    QStringList metBy;  // conditions met at the stopping sample ("residuals", "report:<name>")
};
struct ErrorEvent { QString message; std::vector<SupportIssue> issues; QStringList paths; };
// Proposed: a field snapshot written during the solve (docs/app-protocol-requests.md).
// A live surface view was published; `manifest` is an absolute path.
struct SnapshotEvent { int iter = 0; QString manifest; QString runId; int sequence = 0; };

using ProgressEvent = std::variant<SetupEvent, IterationEvent, WarningEvent, DoneEvent, ErrorEvent, SnapshotEvent>;
Result<ProgressEvent> parseProgressLine(const QByteArray& line);

// kutta's internal checks prefix messages with the source file that raised
// them ("/build/dir/src/kutta/post/request.cpp:30: unknown property"). The
// user reads the message first; the location stays, without the build
// machine's directories, for bug reports: "unknown property (request.cpp:30)".
QString cliMessage(const QString& text);

// Residual names in Fluent's display order.
QStringList residualOrder(const QStringList& present);
QString residualLabel(const QString& name);
// Residuals a Fluent user expects for the case: a 2D case has no z-velocity
// (Kutta solves it on the one-cell extrusion, where it stays at round-off).
bool residualShown(const QString& name, bool planar);
// Fluent's name for a pressure-velocity coupling id ("simplec" -> "SIMPLEC").
QString couplingLabel(const QString& id);

// ---- run outputs ----------------------------------------------------------------------

struct ForceRow {
    QString kind, name;  // kind: "wall" (per zone) or "report"
    int zoneId = -1;
    QVector3D pressure, viscous, total;
    std::optional<QVector3D> direction;
    std::optional<double> pressureProjection, viscousProjection, totalProjection, coefficient;
};
struct ForcesTable { QString notice; std::vector<ForceRow> rows; };
Result<ForcesTable> parseForcesCsv(const QByteArray& csv);

struct ResidualHistory {
    QString notice;
    QStringList names;
    std::vector<int> iterations;
    std::vector<std::vector<double>> values;  // per name, per sample
    std::vector<double> secPerIter;
};
Result<ResidualHistory> parseResidualsCsv(const QByteArray& csv);

struct FluxRow { QString name; int zoneId = -1; double massFlow = 0; };
Result<std::vector<FluxRow>> parseFluxesCsv(const QByteArray& csv, QString* notice = nullptr);

struct ApproximationRecord { QString name, setting, fluentValue, kuttaAction, measuredEffect, documentation, shortLabel; };
struct RunSummary {
    bool approximate = false, modified = false, converged = false;
    QString approximationNotice;
    std::vector<ApproximationRecord> approximations;
    QJsonObject overrides;
    QString caseFile, caseSha, dataFile, dataSha, kuttaVersion, kuttaCommit, backend, device, stopReason;
    qint64 cells = 0;
    int iterations = 0;
    double wallSeconds = 0, setupSeconds = 0;
    bool synthetic = false;
    ConvergenceSetup convergence;
    QStringList metBy;
    std::optional<int> convergedIteration;
    QString coupling;  // pressure_velocity_coupling that ran
};
Result<RunSummary> parseRunSummary(const QJsonObject& json);

// ---- kutta aero ------------------------------------------------------------------------

// Post-processing inputs, stored per case by the app and passed as --aero.
struct AeroSpec {
    struct Group { QString name; QStringList zones; };
    QVector3D drag{1, 0, 0}, lift{0, 0, 1};
    QVector3D momentCentre;
    QStringList carZones;  // empty: the CLI's default (report zones)
    std::vector<Group> groups;
    std::optional<QVector3D> frontAxle;
    std::optional<double> wheelbase;
    bool hasVehicle() const { return frontAxle && wheelbase && *wheelbase > 0; }
    QJsonObject toJson() const;
    static AeroSpec fromJson(const QJsonObject& o);
};
// Groups FSAE teams talk about, from zone names: front wing, rear wing,
// undertray/diffuser, wheels, body for the rest. Zones are the car's walls.
std::vector<AeroSpec::Group> suggestAeroGroups(const QStringList& zones);

struct AeroComponents { double drag = 0, lift = 0, side = 0; };
// Present when the run has positive reference density, velocity and area.
struct AeroCoefficients {
    double cd = 0, cl = 0, cs = 0, cdPressure = 0, cdViscous = 0, clPressure = 0, clViscous = 0;
    std::optional<double> cmPitch, cmRoll, cmYaw;  // also need the reference length
};
struct AeroRow {
    QString name, kind;  // kind: "zone" | "group" | "total"
    QStringList zones;
    AeroComponents pressure, viscous, total;
    double downforce = 0;
    double pitch = 0, roll = 0, yaw = 0;
    std::optional<AeroCoefficients> coefficients;
    std::optional<double> cla, cda;  // m², need only density and velocity
    std::optional<double> downforceShare;
};
struct AeroBalance {
    QString status;  // "defined" | "needs_vehicle" | "no_downforce"
    std::optional<double> wheelbase, frontLoad, rearLoad, frontPercent, copFromFrontAxle;
};
struct AeroResult {
    QString provenance;
    bool approximate = false, modified = false;
    int iterations = 0;
    std::optional<double> area, density, velocity, length, dynamicPressure;
    QVector3D drag, lift, side, momentCentre;
    std::vector<AeroRow> rows;
    QStringList ungrouped;
    std::optional<double> liftToDrag;
    AeroBalance balance;
    const AeroRow* total() const;
};
Result<AeroResult> parseAero(const QJsonObject& json);

// Split one CSV line honouring double quotes.
QStringList splitCsvLine(const QString& line);

} // namespace kapp
