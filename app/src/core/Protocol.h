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

    QStringList approximations() const;  // names a run may opt into
    std::vector<const SupportIssue*> blocking(const QStringList& optedIn) const;
    const Zone* zone(const QString& name) const;
};
Result<CaseSummary> parseCaseSummary(const QJsonObject& json);

// ---- solve --progress jsonl -------------------------------------------------------

struct ReportValue { std::optional<double> force, coefficient; };

struct SetupEvent {
    qint64 cells = 0;
    QString backend, device, init;
    QStringList approximations;
    QJsonObject overrides;
    bool modified = false;
    QString liveRunId;  // live_view.run_id: snapshots from any other run are stale
};
struct IterationEvent {
    int iter = 0;
    QMap<QString, double> residuals;  // ordered by name; display order is the app's
    QMap<QString, ReportValue> reports;
    std::optional<double> secPerIter, elapsed;
};
struct WarningEvent { QString message; };
struct DoneEvent {
    bool converged = false;
    int iterations = 0;
    double elapsed = 0;
    QString out, stopReason;
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

struct ApproximationRecord { QString name, setting, fluentValue, kuttaAction, measuredEffect, documentation; };
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
};
Result<RunSummary> parseRunSummary(const QJsonObject& json);

// Split one CSV line honouring double quotes.
QStringList splitCsvLine(const QString& line);

} // namespace kapp
