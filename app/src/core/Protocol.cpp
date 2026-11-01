#include "core/Protocol.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
#include <QRegularExpression>
#include <QUrl>

#include <cmath>

namespace kapp {

namespace {

std::optional<double> num(const QJsonValue& v)
{
    if (v.isDouble() && std::isfinite(v.toDouble())) return v.toDouble();
    return std::nullopt;
}

std::optional<QVector3D> vec3(const QJsonValue& v)
{
    const QJsonArray a = v.toArray();
    if (a.size() != 3 || !a[0].isDouble() || !a[1].isDouble() || !a[2].isDouble()) return std::nullopt;
    return QVector3D(float(a[0].toDouble()), float(a[1].toDouble()), float(a[2].toDouble()));
}

QStringList strings(const QJsonValue& v)
{
    QStringList out;
    for (const QJsonValue& s : v.toArray()) out << s.toString();
    return out;
}

bool versionOk(const QJsonObject& o) { return o.value("v").toInt() == 1; }

SupportIssue parseIssue(const QJsonObject& o)
{
    SupportIssue i;
    i.path = o.value("path").toString();
    i.setting = o.value("setting").toString(i.path);
    i.label = o.value("label").toString(i.setting);
    const QString c = o.value("class").toString();
    i.cls = c == "unsupported" ? IssueClass::Unsupported
          : c == "approximable" ? IssueClass::Approximable
                                : IssueClass::Informational;
    i.approximation = o.value("approximation").toString();
    i.explanation = o.value("explanation").toString();
    i.evidence = o.value("evidence").toString();
    i.measuredEffect = o.value("measured_effect").toString();
    i.shortLabel = o.value("short").toString();
    return i;
}

std::optional<SettingKind> parseKind(const QString& k)
{
    if (k == "number") return SettingKind::Number;
    if (k == "integer") return SettingKind::Integer;
    if (k == "enum") return SettingKind::Enum;
    if (k == "bool") return SettingKind::Bool;
    if (k == "vec3") return SettingKind::Vec3;
    if (k == "string") return SettingKind::String;
    return std::nullopt;
}

// Sections whose second segment names an owner (zone, material, report).
bool ownedSection(const QString& s)
{
    return s == "boundary" || s == "cell-zone" || s == "materials" || s == "reports" || s == "report";
}

Result<Setting> parseSetting(const QJsonObject& o)
{
    Setting s;
    s.path = o.value("path").toString();
    if (s.path.isEmpty()) return Result<Setting>::failure("setting without a path");
    s.label = o.value("label").toString(s.path);
    s.group = o.value("group").toString();
    s.unit = o.value("unit").toString();
    s.reason = o.value("reason").toString();
    const auto kind = parseKind(o.value("kind").toString());
    if (!kind) return Result<Setting>::failure(QStringLiteral("%1: unknown kind '%2'").arg(s.path, o.value("kind").toString()));
    s.kind = *kind;
    s.value = o.value("value");
    s.options = o.value("options").toArray();
    s.optionLabels = strings(o.value("option_labels"));
    if (s.optionLabels.size() != s.options.size()) s.optionLabels.clear();
    s.valueLabel = o.value("value_label").toString();
    s.min = num(o.value("min"));
    s.max = num(o.value("max"));
    s.editable = o.value("editable").toBool();
    // Older CLIs have no activity metadata: everything is in use.
    s.active = o.value("active").toBool(true);
    s.dependsOn = o.value("depends_on").toString();
    s.activeWhen = o.value("active_when").toArray();
    const QStringList seg = s.path.split('/');
    s.section = seg.value(0);
    int leafStart = 1;
    if (ownedSection(s.section) && seg.size() > 2) {
        s.owner = decodePathSegment(seg[1]);
        leafStart = 2;
    }
    // Retained Fluent state with no Kutta control: the whole "case" section,
    // per-owner "raw/" settings and the solver's case-config flags and sources.
    s.raw = s.section == "case" || seg.value(leafStart) == "raw"
            || (s.section == "general" && (seg.value(1) == "case-config" || seg.value(1) == "source"));
    QStringList leaf;
    for (int i = leafStart; i < seg.size(); ++i) leaf << decodePathSegment(seg[i]);
    s.leaf = leaf.join('/');
    return s;
}

} // namespace

QString Setting::labelFor(const QJsonValue& v) const
{
    const auto same = [](const QJsonValue& a, const QJsonValue& b) { return a.isDouble() && b.isDouble() ? a.toDouble() == b.toDouble() : a == b; };
    for (qsizetype i = 0; i < options.size() && i < optionLabels.size(); ++i)
        if (same(options[i], v) && !optionLabels[i].isEmpty()) return optionLabels[i];
    return same(v, value) ? valueLabel : QString();
}

QString SupportIssue::explanationOnly() const
{
    const QString tail = ' ' + measuredEffect + '.';
    if (!measuredEffect.isEmpty() && explanation.endsWith(tail)) return explanation.chopped(tail.size());
    return explanation;
}

QString decodePathSegment(const QString& segment)
{
    return QUrl::fromPercentEncoding(segment.toUtf8());
}

QString Device::label() const
{
    const QString api = backend == "cuda" ? "CUDA" : backend == "metal" ? "Metal" : "CPU";
    if (memoryBytes <= 0) return QStringLiteral("%1 · %2").arg(name, api);
    return QStringLiteral("%1 · %2 · %3 GB").arg(name, api, QString::number(memoryBytes / 1073741824.0, 'f', 1));
}

Result<std::vector<Device>> parseDevices(const QJsonObject& json)
{
    if (!versionOk(json)) return Result<std::vector<Device>>::failure("devices: unsupported protocol version");
    std::vector<Device> out;
    for (const QJsonValue& v : json.value("devices").toArray()) {
        const QJsonObject o = v.toObject();
        Device d;
        d.backend = o.value("backend").toString();
        d.name = o.value("name").toString();
        d.memoryBytes = qint64(o.value("memory_bytes").toDouble());
        d.ordinal = o.value("ordinal").toInt();
        d.isDefault = o.value("default").toBool();
        if (d.backend.isEmpty()) return Result<std::vector<Device>>::failure("devices: entry without a backend");
        out.push_back(d);
    }
    return out;
}

QStringList CaseSummary::approximations() const
{
    QStringList out;
    for (const auto& i : issues)
        if (i.cls == IssueClass::Approximable && !i.approximation.isEmpty() && !out.contains(i.approximation))
            out << i.approximation;
    return out;
}

QString CaseSummary::approximationShort(const QString& name) const
{
    for (const auto& i : issues)
        if (i.cls == IssueClass::Approximable && i.approximation == name && !i.shortLabel.isEmpty()) return i.shortLabel;
    return name;
}

std::vector<const SupportIssue*> CaseSummary::blocking(const QStringList& optedIn) const
{
    std::vector<const SupportIssue*> out;
    for (const auto& i : issues)
        if (i.cls == IssueClass::Unsupported || (i.cls == IssueClass::Approximable && !optedIn.contains(i.approximation)))
            out.push_back(&i);
    return out;
}

const Zone* CaseSummary::zone(const QString& name) const
{
    for (const auto& z : zones)
        if (z.name == name) return &z;
    return nullptr;
}

Result<CaseSummary> parseCaseSummary(const QJsonObject& json)
{
    using R = Result<CaseSummary>;
    if (!versionOk(json)) return R::failure("inspect: unsupported protocol version");
    CaseSummary c;
    c.synthetic = json.value("synthetic").toBool();
    const QJsonObject file = json.value("file").toObject();
    c.filePath = file.value("path").toString();
    c.sha256 = file.value("sha256").toString();
    c.fluentVersion = file.value("fluent_version").toString();
    if (file.value("format").isString()) c.format = file.value("format").toString();
    c.meshFile = file.value("mesh_file").toObject().value("path").toString();
    const QJsonObject mesh = json.value("mesh").toObject();
    c.cells = qint64(mesh.value("cells").toDouble());
    c.faces = qint64(mesh.value("faces").toDouble());
    c.nodes = qint64(mesh.value("nodes").toDouble());
    const QJsonObject types = mesh.value("cell_types").toObject();
    for (auto it = types.begin(); it != types.end(); ++it) c.cellTypes[it.key()] = qint64(it.value().toDouble());
    const QJsonArray b = mesh.value("bounds").toArray();
    if (b.size() == 2) {
        if (auto lo = vec3(b[0]), hi = vec3(b[1]); lo && hi) {
            c.bounds.min = *lo;
            c.bounds.max = *hi;
            c.bounds.valid = true;
        }
    }
    for (const QJsonValue& zv : json.value("zones").toArray()) {
        const QJsonObject o = zv.toObject();
        Zone z{o.value("id").toInt(-1), o.value("name").toString(), o.value("type").toString(),
               qint64(o.value("faces").toDouble()), o.value("condition").toObject()};
        if (z.name.isEmpty()) return R::failure("inspect: zone without a name");
        c.zones.push_back(z);
    }
    const QJsonObject models = json.value("models").toObject();
    c.viscousModel = models.value("viscous").toString();
    c.transition = models.value("transition").toString();
    if (models.value("energy").isBool()) c.energy = models.value("energy").toBool();
    if (models.value("energy_coupled").isBool()) c.energyCoupled = models.value("energy_coupled").toBool();
    for (const QJsonValue& mv : json.value("materials").toArray()) {
        const QJsonObject o = mv.toObject();
        c.materials.push_back({o.value("name").toString(), num(o.value("density")), num(o.value("viscosity"))});
    }
    for (const QJsonValue& rv : json.value("reports").toArray()) {
        const QJsonObject o = rv.toObject();
        c.reports.push_back({o.value("name").toString(), o.value("kind").toString(), strings(o.value("zones")),
                             vec3(o.value("direction"))});
    }
    const QJsonObject support = json.value("support").toObject();
    const QString st = support.value("status").toString();
    if (st == "supported") c.status = SupportStatus::Supported;
    else if (st == "approximate") c.status = SupportStatus::Approximate;
    else if (st == "unsupported") c.status = SupportStatus::Unsupported;
    else return R::failure(QStringLiteral("inspect: unknown support status '%1'").arg(st));
    for (const QJsonValue& iv : support.value("issues").toArray()) c.issues.push_back(parseIssue(iv.toObject()));
    const QJsonArray settings = json.value("settings").toArray();
    c.settings.reserve(std::size_t(settings.size()));
    for (const QJsonValue& sv : settings) {
        auto s = parseSetting(sv.toObject());
        if (!s) return R::failure("inspect: " + s.error());
        c.settings.push_back(std::move(s).value());
        if (c.settings.back().path == "general/dimension" && c.settings.back().value.isDouble())
            c.dimension = c.settings.back().value.toInt();
    }
    const QJsonObject data = json.value("data").toObject();
    c.dataPresent = data.value("present").toBool();
    c.dataFields = strings(data.value("fields"));
    if (data.value("iterations_hint").isDouble()) c.iterationsHint = data.value("iterations_hint").toInt();
    return c;
}

QString cliMessage(const QString& text)
{
    // Per line: an optional "kutta <command>: " lead, then [C:]path/file.ext:line: .
    static const QRegularExpression loc(QStringLiteral(R"(^((?:kutta [\w-]+: )?)(?:(?:[A-Za-z]:)?[^\s:]*[/\\])?([\w.-]+\.(?:cpp|hpp|cu|h)):(\d+): (.*)$)"));
    QStringList lines = text.split('\n');
    for (QString& l : lines)
        if (const auto m = loc.match(l); m.hasMatch()) l = QStringLiteral("%1%2 (%3:%4)").arg(m.captured(1), m.captured(4), m.captured(2), m.captured(3));
    return lines.join('\n');
}

Result<ProgressEvent> parseProgressLine(const QByteArray& line)
{
    using R = Result<ProgressEvent>;
    QJsonParseError err;
    const QJsonDocument doc = QJsonDocument::fromJson(line, &err);
    if (!doc.isObject()) return R::failure(QStringLiteral("not a JSON object: %1").arg(err.errorString()));
    const QJsonObject o = doc.object();
    if (!versionOk(o)) return R::failure("unsupported protocol version");
    const QString ev = o.value("event").toString();
    if (ev == "setup") {
        SetupEvent s;
        s.cells = qint64(o.value("cells").toDouble());
        s.backend = o.value("backend").toString();
        s.device = o.value("device").toString();
        s.init = o.value("init").toString();
        s.approximations = strings(o.value("approximations"));
        s.liveRunId = o.value("live_view").toObject().value("run_id").toString();
        s.overrides = o.value("overrides").toObject();
        s.modified = o.value("modified").toBool();
        s.convergence = parseConvergenceSetup(o.value("convergence").toObject());
        s.coupling = o.value("pressure_velocity_coupling").toString();
        return ProgressEvent{s};
    }
    if (ev == "iteration") {
        IterationEvent it;
        it.iter = o.value("iter").toInt();
        const QJsonObject res = o.value("residuals").toObject();
        for (auto r = res.begin(); r != res.end(); ++r)
            if (auto v = num(r.value())) it.residuals[r.key()] = *v;
        const QJsonObject reps = o.value("reports").toObject();
        for (auto r = reps.begin(); r != reps.end(); ++r) {
            const QJsonObject rv = r.value().toObject();
            it.reports[r.key()] = {num(rv.value("force")), num(rv.value("coefficient"))};
        }
        it.secPerIter = num(o.value("sec_per_iter"));
        it.elapsed = num(o.value("elapsed"));
        if (o.value("convergence").isObject()) it.convergence = parseConvergenceState(o.value("convergence").toObject());
        return ProgressEvent{it};
    }
    if (ev == "warning") return ProgressEvent{WarningEvent{o.value("message").toString()}};
    if (ev == "done") {
        DoneEvent d;
        d.converged = o.value("converged").toBool();
        d.iterations = o.value("iterations").toInt();
        d.elapsed = o.value("elapsed").toDouble();
        d.out = o.value("out").toString();
        d.stopReason = o.value("stop_reason").toString();
        d.metBy = strings(o.value("convergence").toObject().value("met_by"));
        return ProgressEvent{d};
    }
    if (ev == "error") {
        ErrorEvent e;
        e.message = cliMessage(o.value("message").toString());
        for (const QJsonValue& iv : o.value("issues").toArray()) {
            const QJsonObject io = iv.toObject();
            e.issues.push_back(parseIssue(io));
            if (io.contains("path")) e.paths << io.value("path").toString();
        }
        return ProgressEvent{e};
    }
    if (ev == "snapshot")
        return ProgressEvent{SnapshotEvent{o.value("iter").toInt(), o.value("manifest").toString(),
                                           o.value("run_id").toString(), o.value("sequence").toInt()}};
    return R::failure(QStringLiteral("unknown event '%1'").arg(ev));
}

QStringList residualOrder(const QStringList& present)
{
    static const QStringList fluent{"continuity", "x-velocity", "y-velocity", "z-velocity", "energy", "k", "omega", "epsilon", "intermit"};
    QStringList out;
    for (const QString& n : fluent)
        if (present.contains(n)) out << n;
    for (const QString& n : present)
        if (!out.contains(n)) out << n;
    return out;
}

QString residualLabel(const QString& name) { return name; }  // Fluent prints the same names

bool residualShown(const QString& name, bool planar) { return !(planar && name == QLatin1String("z-velocity")); }

QString couplingLabel(const QString& id)
{
    if (id == "simple") return QStringLiteral("SIMPLE");
    if (id == "simplec") return QStringLiteral("SIMPLEC");
    if (id == "coupled") return QStringLiteral("Coupled");
    if (id == "piso") return QStringLiteral("PISO");
    return id;
}

QStringList splitCsvLine(const QString& line)
{
    QStringList out;
    QString cur;
    bool quoted = false;
    for (qsizetype i = 0; i < line.size(); ++i) {
        const QChar ch = line[i];
        if (quoted) {
            if (ch == '"' && i + 1 < line.size() && line[i + 1] == '"') { cur += '"'; ++i; }
            else if (ch == '"') quoted = false;
            else cur += ch;
        } else if (ch == '"') quoted = true;
        else if (ch == ',') { out << cur; cur.clear(); }
        else cur += ch;
    }
    out << cur;
    return out;
}

namespace {

struct Csv {
    QString notice;
    QStringList header;
    std::vector<QStringList> rows;
};

Result<Csv> readCsv(const QByteArray& bytes)
{
    Csv csv;
    const QStringList lines = QString::fromUtf8(bytes).split('\n');
    for (QString line : lines) {
        line = line.trimmed();
        if (line.isEmpty()) continue;
        if (line.startsWith('#')) {
            if (!csv.notice.isEmpty()) csv.notice += '\n';
            csv.notice += line.mid(1).trimmed();
            continue;
        }
        if (csv.header.isEmpty()) csv.header = splitCsvLine(line);
        else csv.rows.push_back(splitCsvLine(line));
    }
    if (csv.header.isEmpty()) return Result<Csv>::failure("CSV without a header");
    return csv;
}

std::optional<double> cell(const QStringList& row, const QStringList& header, const QString& col)
{
    const qsizetype i = header.indexOf(col);
    if (i < 0 || i >= row.size() || row[i].isEmpty()) return std::nullopt;
    bool ok = false;
    const double v = row[i].toDouble(&ok);
    return ok && std::isfinite(v) ? std::optional(v) : std::nullopt;
}

QVector3D vecCols(const QStringList& row, const QStringList& h, const QString& prefix)
{
    return QVector3D(float(cell(row, h, prefix + "_x").value_or(0)), float(cell(row, h, prefix + "_y").value_or(0)),
                     float(cell(row, h, prefix + "_z").value_or(0)));
}

} // namespace

Result<ForcesTable> parseForcesCsv(const QByteArray& bytes)
{
    auto csv = readCsv(bytes);
    if (!csv) return Result<ForcesTable>::failure("forces.csv: " + csv.error());
    const QStringList& h = csv->header;
    if (!h.contains("kind") || !h.contains("name") || !h.contains("total_x"))
        return Result<ForcesTable>::failure("forces.csv: unexpected columns");
    ForcesTable t;
    t.notice = csv->notice;
    for (const QStringList& row : csv->rows) {
        ForceRow r;
        r.kind = row.value(h.indexOf("kind"));
        r.name = row.value(h.indexOf("name"));
        r.zoneId = int(cell(row, h, "zone_id").value_or(-1));
        r.pressure = vecCols(row, h, "pressure");
        r.viscous = vecCols(row, h, "viscous");
        r.total = vecCols(row, h, "total");
        if (cell(row, h, "direction_x")) r.direction = vecCols(row, h, "direction");
        r.pressureProjection = cell(row, h, "pressure_projection");
        r.viscousProjection = cell(row, h, "viscous_projection");
        r.totalProjection = cell(row, h, "total_projection");
        r.coefficient = cell(row, h, "coefficient");
        t.rows.push_back(r);
    }
    return t;
}

Result<ResidualHistory> parseResidualsCsv(const QByteArray& bytes)
{
    auto csv = readCsv(bytes);
    if (!csv) return Result<ResidualHistory>::failure("residuals.csv: " + csv.error());
    const QStringList& h = csv->header;
    if (h.value(0) != "iteration") return Result<ResidualHistory>::failure("residuals.csv: first column is not iteration");
    ResidualHistory r;
    r.notice = csv->notice;
    for (int i = 1; i < h.size(); ++i)
        if (h[i] != "seconds_per_iteration") r.names << h[i];
    r.values.resize(std::size_t(r.names.size()));
    for (const QStringList& row : csv->rows) {
        r.iterations.push_back(row.value(0).toInt());
        for (int n = 0; n < r.names.size(); ++n)
            r.values[std::size_t(n)].push_back(cell(row, h, r.names[n]).value_or(std::nan("")));
        r.secPerIter.push_back(cell(row, h, "seconds_per_iteration").value_or(std::nan("")));
    }
    return r;
}

Result<std::vector<FluxRow>> parseFluxesCsv(const QByteArray& bytes, QString* notice)
{
    auto csv = readCsv(bytes);
    if (!csv) return Result<std::vector<FluxRow>>::failure("fluxes.csv: " + csv.error());
    if (notice) *notice = csv->notice;
    std::vector<FluxRow> out;
    const QStringList& h = csv->header;
    for (const QStringList& row : csv->rows)
        out.push_back({row.value(h.indexOf("name")), int(cell(row, h, "zone_id").value_or(-1)),
                       cell(row, h, "mass_flow_outward").value_or(std::nan(""))});
    return out;
}

Result<RunSummary> parseRunSummary(const QJsonObject& o)
{
    if (!versionOk(o)) return Result<RunSummary>::failure("summary.json: unsupported version");
    RunSummary s;
    s.approximate = o.value("approximate").toBool();
    s.modified = o.value("modified").toBool();
    s.converged = o.value("converged").toBool();
    s.approximationNotice = o.value("approximation_notice").toString();
    for (const QJsonValue& av : o.value("approximations").toArray()) {
        const QJsonObject a = av.toObject();
        s.approximations.push_back({a.value("name").toString(), a.value("setting").toString(),
                                    a.value("fluent_value").toString(), a.value("kutta_action").toString(),
                                    a.value("measured_effect").toString(), a.value("effect_documentation").toString(),
                                    a.value("short").toString()});
    }
    s.overrides = o.value("overrides").toObject();
    s.caseFile = o.value("case_file").toString();
    s.caseSha = o.value("case_sha256").toString();
    s.dataFile = o.value("data_file").toString();
    s.dataSha = o.value("data_sha256").toString();
    s.kuttaVersion = o.value("kutta_version").toString();
    s.kuttaCommit = o.value("kutta_commit").toString();
    s.backend = o.value("backend").toString();
    s.device = o.value("device").toString();
    s.stopReason = o.value("stop_reason").toString();
    s.cells = qint64(o.value("cells").toDouble());
    s.poorQualityCells = qint64(o.value("poor_quality_cells").toDouble());
    s.poorCellQuality = num(o.value("poor_cell_quality"));
    s.iterations = o.value("iterations").toInt();
    s.wallSeconds = o.value("wall_time_seconds").toDouble();
    s.setupSeconds = o.value("setup_time_seconds").toDouble();
    s.synthetic = o.value("synthetic").toBool();
    const QJsonObject conv = o.value("convergence").toObject();
    s.convergence = parseConvergenceSetup(conv);
    s.metBy = strings(conv.value("met_by"));
    if (conv.value("iteration").isDouble()) s.convergedIteration = conv.value("iteration").toInt();
    s.coupling = o.value("pressure_velocity_coupling").toString();
    return s;
}

// ---- convergence ------------------------------------------------------------------------

ConvergenceSetup parseConvergenceSetup(const QJsonObject& o)
{
    ConvergenceSetup c;
    c.enabled = o.value("enabled").toBool();
    c.condition = o.value("condition").toString("all");
    const QJsonObject res = o.value("residuals").toObject();
    for (auto it = res.begin(); it != res.end(); ++it)
        if (auto v = num(it.value())) c.residuals[it.key()] = *v;
    for (const QJsonValue& rv : o.value("reports").toArray()) {
        const QJsonObject r = rv.toObject();
        c.reports.push_back({r.value("report").toString(), r.value("stop_criterion").toDouble(), r.value("previous_values").toInt(),
                             r.value("ignore_before").toInt()});
    }
    return c;
}

ConvergenceState parseConvergenceState(const QJsonObject& o)
{
    ConvergenceState s;
    s.met = o.value("met").toBool();
    for (const QJsonValue& cv : o.value("conditions").toArray()) {
        const QJsonObject c = cv.toObject();
        ConditionState st;
        st.kind = c.value("kind").toString();
        st.report = c.value("report").toString();
        st.met = c.value("met").toBool();
        st.value = num(c.value("value"));
        st.criterion = num(c.value("criterion"));
        const QJsonObject eq = c.value("equations").toObject();
        for (auto it = eq.begin(); it != eq.end(); ++it) {
            const QJsonObject e = it.value().toObject();
            st.equations[it.key()] = {e.value("value").toDouble(), e.value("criterion").toDouble(), e.value("met").toBool()};
        }
        s.conditions.push_back(std::move(st));
    }
    return s;
}

QString describeCondition(const QString& id, const ConvergenceSetup& setup)
{
    if (id == "residuals") return QStringLiteral("every checked residual reached its criterion");
    const QString name = id.startsWith("report:") ? id.mid(7) : id;
    for (const auto& r : setup.reports)
        if (r.report == name)
            return QStringLiteral("%1 changed less than %2 % between each of its last %3 values")
                .arg(name, QLocale::c().toString(r.stopCriterion * 100, 'g', 3)).arg(r.previousValues);
    return name;
}

// ---- aero -------------------------------------------------------------------------------

namespace {

QJsonArray vec3(const QVector3D& v) { return {double(v.x()), double(v.y()), double(v.z())}; }

} // namespace

QJsonObject AeroSpec::toJson() const
{
    QJsonObject o{{"v", 1}, {"axes", QJsonObject{{"drag", vec3(drag)}, {"lift", vec3(lift)}}}, {"moment_centre", vec3(momentCentre)}};
    if (!carZones.isEmpty()) o["car_zones"] = QJsonArray::fromStringList(carZones);
    QJsonArray g;
    for (const Group& x : groups)
        if (!x.zones.isEmpty()) g.append(QJsonObject{{"name", x.name}, {"zones", QJsonArray::fromStringList(x.zones)}});
    if (!g.isEmpty()) o["groups"] = g;
    if (hasVehicle()) o["vehicle"] = QJsonObject{{"front_axle", vec3(*frontAxle)}, {"wheelbase", *wheelbase}};
    return o;
}

AeroSpec AeroSpec::fromJson(const QJsonObject& o)
{
    AeroSpec s;
    const QJsonObject axes = o.value("axes").toObject();
    if (auto d = vec3(axes.value("drag"))) s.drag = *d;
    if (auto l = vec3(axes.value("lift"))) s.lift = *l;
    if (auto c = vec3(o.value("moment_centre"))) s.momentCentre = *c;
    s.carZones = strings(o.value("car_zones"));
    for (const QJsonValue& gv : o.value("groups").toArray())
        s.groups.push_back({gv.toObject().value("name").toString(), strings(gv.toObject().value("zones"))});
    const QJsonObject v = o.value("vehicle").toObject();
    s.frontAxle = vec3(v.value("front_axle"));
    s.wheelbase = num(v.value("wheelbase"));
    return s;
}

std::vector<AeroSpec::Group> suggestAeroGroups(const QStringList& zones)
{
    struct Rule { QString name; QRegularExpression pattern; };
    // First match wins: "front-wing-endplate" is front wing, "rear-tyre" a wheel.
    // A wing, flap or element with no front or rear in its name (a 2D section,
    // a lone wing) is a wing of its own, not body.
    static const std::vector<Rule> rules{
        {"Front wing", QRegularExpression(R"((^|[^a-z])(fw|front[-_ ]?wing|nose[-_ ]?wing))", QRegularExpression::CaseInsensitiveOption)},
        {"Rear wing", QRegularExpression(R"((^|[^a-z])(rw|rear[-_ ]?wing|drs))", QRegularExpression::CaseInsensitiveOption)},
        {"Wing", QRegularExpression(R"(wing|flap|element|main[-_ ]?plane|foil|slat|gurney|endplate)", QRegularExpression::CaseInsensitiveOption)},
        {"Undertray / diffuser", QRegularExpression(R"(undertray|under[-_ ]?tray|diffuser|floor|underbody|splitter)", QRegularExpression::CaseInsensitiveOption)},
        {"Wheels", QRegularExpression(R"(wheel|tyre|tire|rim)", QRegularExpression::CaseInsensitiveOption)},
    };
    std::vector<AeroSpec::Group> out;
    for (const Rule& r : rules) out.push_back({r.name, {}});
    out.push_back({"Body", {}});
    for (const QString& z : zones) {
        auto it = std::find_if(rules.begin(), rules.end(), [&](const Rule& r) { return r.pattern.match(z).hasMatch(); });
        out[std::size_t(it - rules.begin())].zones << z;
    }
    std::erase_if(out, [](const AeroSpec::Group& g) { return g.zones.isEmpty(); });
    return out;
}

const AeroRow* AeroResult::total() const
{
    auto it = std::find_if(rows.begin(), rows.end(), [](const AeroRow& r) { return r.kind == "total"; });
    return it == rows.end() ? nullptr : &*it;
}

Result<AeroResult> parseAero(const QJsonObject& o)
{
    using R = Result<AeroResult>;
    if (!versionOk(o)) return R::failure("aero: unsupported protocol version");
    AeroResult a;
    a.provenance = o.value("provenance").toString();
    a.approximate = o.value("approximate").toBool();
    a.modified = o.value("modified").toBool();
    a.iterations = o.value("iterations").toInt();
    const QJsonObject ref = o.value("reference").toObject();
    a.area = num(ref.value("area"));
    a.density = num(ref.value("density"));
    a.velocity = num(ref.value("velocity"));
    a.length = num(ref.value("length"));
    a.dynamicPressure = num(ref.value("dynamic_pressure"));
    const QJsonObject axes = o.value("axes").toObject();
    a.drag = vec3(axes.value("drag")).value_or(QVector3D(1, 0, 0));
    a.lift = vec3(axes.value("lift")).value_or(QVector3D(0, 0, 1));
    a.side = vec3(axes.value("side")).value_or(QVector3D::crossProduct(a.lift, a.drag));
    a.momentCentre = vec3(o.value("moment_centre")).value_or(QVector3D());
    auto comps = [](const QJsonValue& v) {
        const QJsonObject c = v.toObject();
        return AeroComponents{c.value("drag").toDouble(), c.value("lift").toDouble(), c.value("side").toDouble()};
    };
    for (const QJsonValue& rv : o.value("rows").toArray()) {
        const QJsonObject r = rv.toObject();
        AeroRow row;
        row.name = r.value("name").toString();
        row.kind = r.value("kind").toString();
        row.zones = strings(r.value("zones"));
        const QJsonObject f = r.value("force").toObject();
        row.pressure = comps(f.value("pressure"));
        row.viscous = comps(f.value("viscous"));
        row.total = comps(f.value("total"));
        row.downforce = r.value("downforce").toDouble(-row.total.lift);
        const QJsonObject m = r.value("moment").toObject();
        row.pitch = m.value("pitch").toDouble();
        row.roll = m.value("roll").toDouble();
        row.yaw = m.value("yaw").toDouble();
        // Every key is present; nulls mark the references the run lacks.
        const QJsonObject c = r.value("coefficients").toObject();
        if (c.value("cd").isDouble()) {
            AeroCoefficients k;
            k.cd = c.value("cd").toDouble();
            k.cl = c.value("cl").toDouble();
            k.cs = c.value("cs").toDouble();
            k.cdPressure = c.value("cd_pressure").toDouble();
            k.cdViscous = c.value("cd_viscous").toDouble();
            k.clPressure = c.value("cl_pressure").toDouble();
            k.clViscous = c.value("cl_viscous").toDouble();
            k.cmPitch = num(c.value("cm_pitch"));
            k.cmRoll = num(c.value("cm_roll"));
            k.cmYaw = num(c.value("cm_yaw"));
            row.coefficients = k;
        }
        row.cla = num(c.value("cla"));
        row.cda = num(c.value("cda"));
        row.downforceShare = num(r.value("downforce_share"));
        a.rows.push_back(std::move(row));
    }
    if (!a.total()) return R::failure("aero: no total row");
    a.ungrouped = strings(o.value("ungrouped"));
    a.liftToDrag = num(o.value("lift_to_drag"));
    const QJsonObject b = o.value("balance").toObject();
    a.balance.status = b.value("status").toString("needs_vehicle");
    a.balance.wheelbase = num(b.value("wheelbase"));
    a.balance.frontLoad = num(b.value("front_load"));
    a.balance.rearLoad = num(b.value("rear_load"));
    a.balance.frontPercent = num(b.value("front_percent"));
    a.balance.copFromFrontAxle = num(b.value("cop_from_front_axle"));
    return a;
}


std::vector<const NewCaseZone*> NewCaseDraft::boundaries() const
{
    std::vector<const NewCaseZone*> out;
    for (const auto& z : zones)
        if (z.kind == "boundary") out.push_back(&z);
    return out;
}

NewCaseZone* NewCaseDraft::zone(const QString& name)
{
    for (auto& z : zones)
        if (z.name == name) return &z;
    return nullptr;
}

Result<NewCaseDraft> parseNewCase(const QJsonObject& json)
{
    using R = Result<NewCaseDraft>;
    if (!versionOk(json)) return R::failure("new-case: unsupported protocol version");
    NewCaseDraft d;
    d.templateName = json.value("template").toString();
    d.speed = json.value("speed").toDouble();
    const QJsonObject mesh = json.value("mesh").toObject();
    d.meshPath = mesh.value("path").toString();
    d.dimension = mesh.value("dimension").toInt(3);
    d.cells = qint64(mesh.value("cells").toDouble());
    d.faces = qint64(mesh.value("faces").toDouble());
    d.nodes = qint64(mesh.value("nodes").toDouble());
    d.polyhedral = mesh.value("polyhedral").toBool();
    for (const QJsonValue& v : json.value("zones").toArray()) {
        const QJsonObject o = v.toObject();
        NewCaseZone z{o.value("id").toInt(-1), o.value("name").toString(), o.value("mesh_type").toString(), o.value("kind").toString(),
                      qint64(o.value("faces").toDouble()), o.value("role").toString(), o.value("suggested_role").toString(),
                      o.value("reason").toString(), o.value("group").toString()};
        if (z.name.isEmpty()) return R::failure("new-case: zone without a name");
        d.zones.push_back(z);
    }
    for (const QJsonValue& v : json.value("roles").toArray()) d.roles << v.toString();
    for (const QJsonValue& v : json.value("notes").toArray()) d.notes << v.toString();
    d.kuttaCase = json.value("case").toObject();
    d.out = json.value("out").toString();
    if (d.meshPath.isEmpty() || d.roles.isEmpty()) return R::failure("new-case: incomplete reply");
    return d;
}

Result<ConvertResult> parseConvert(const QJsonObject& json)
{
    using R = Result<ConvertResult>;
    if (!versionOk(json)) return R::failure("convert: unsupported protocol version");
    ConvertResult r;
    r.out = json.value("out").toString();
    for (const QJsonValue& v : json.value("omitted").toArray()) {
        const QJsonObject o = v.toObject();
        r.omitted.push_back({o.value("path").toString(), o.value("reason").toString(), o.value("class").toString() == "unsupported"});
    }
    for (const QJsonValue& v : json.value("run_controls_not_saved").toArray()) r.runControlsNotSaved << v.toString();
    if (r.out.isEmpty()) return R::failure("convert: no output path in the reply");
    return r;
}

} // namespace kapp
