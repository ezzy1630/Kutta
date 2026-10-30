#include "core/Protocol.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
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
    s.options = strings(o.value("options"));
    s.min = num(o.value("min"));
    s.max = num(o.value("max"));
    s.editable = o.value("editable").toBool();
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
    }
    const QJsonObject data = json.value("data").toObject();
    c.dataPresent = data.value("present").toBool();
    c.dataFields = strings(data.value("fields"));
    if (data.value("iterations_hint").isDouble()) c.iterationsHint = data.value("iterations_hint").toInt();
    return c;
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
        s.overrides = o.value("overrides").toObject();
        s.modified = o.value("modified").toBool();
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
        return ProgressEvent{d};
    }
    if (ev == "error") {
        ErrorEvent e;
        e.message = o.value("message").toString();
        for (const QJsonValue& iv : o.value("issues").toArray()) {
            const QJsonObject io = iv.toObject();
            e.issues.push_back(parseIssue(io));
            if (io.contains("path")) e.paths << io.value("path").toString();
        }
        return ProgressEvent{e};
    }
    if (ev == "snapshot") return ProgressEvent{SnapshotEvent{o.value("iter").toInt(), o.value("view").toString()}};
    return R::failure(QStringLiteral("unknown event '%1'").arg(ev));
}

QStringList residualOrder(const QStringList& present)
{
    static const QStringList fluent{"continuity", "x-velocity", "y-velocity", "z-velocity", "energy", "k", "omega"};
    QStringList out;
    for (const QString& n : fluent)
        if (present.contains(n)) out << n;
    for (const QString& n : present)
        if (!out.contains(n)) out << n;
    return out;
}

QString residualLabel(const QString& name) { return name; }  // Fluent prints the same names

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
                                    a.value("measured_effect").toString(), a.value("effect_documentation").toString()});
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
    s.iterations = o.value("iterations").toInt();
    s.wallSeconds = o.value("wall_time_seconds").toDouble();
    s.setupSeconds = o.value("setup_time_seconds").toDouble();
    s.synthetic = o.value("synthetic").toBool();
    return s;
}

} // namespace kapp
