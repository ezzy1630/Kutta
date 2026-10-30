#include "core/Study.h"

#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
#include <QRegularExpression>
#include <QSaveFile>

#include <cmath>

namespace kapp {

QString variantStatusLabel(VariantStatus s)
{
    switch (s) {
    case VariantStatus::Queued: return QStringLiteral("Queued");
    case VariantStatus::Running: return QStringLiteral("Running");
    case VariantStatus::Done: return QStringLiteral("Done");
    case VariantStatus::Failed: return QStringLiteral("Failed");
    case VariantStatus::Stopped: return QStringLiteral("Stopped");
    case VariantStatus::Cancelled: return QStringLiteral("Cancelled");
    }
    return {};
}

namespace {

QString statusKey(VariantStatus s) { return variantStatusLabel(s).toLower(); }

VariantStatus statusFrom(const QString& k)
{
    for (VariantStatus s : {VariantStatus::Queued, VariantStatus::Running, VariantStatus::Done, VariantStatus::Failed,
                            VariantStatus::Stopped, VariantStatus::Cancelled})
        if (statusKey(s) == k) return s;
    return VariantStatus::Failed;
}

} // namespace

QJsonObject Study::overridesFor(const StudyVariant& v) const
{
    QJsonObject o = base;
    o.remove(parameter);
    for (const QString& f : followers) o.remove(f);
    for (auto it = v.overrides.begin(); it != v.overrides.end(); ++it) o.insert(it.key(), it.value());
    return o;
}

int Study::nextQueued() const
{
    for (std::size_t i = 0; i < variants.size(); ++i)
        if (variants[i].status == VariantStatus::Queued) return int(i);
    return -1;
}

int Study::indexOfRun(const QString& runId) const
{
    for (std::size_t i = 0; i < variants.size(); ++i)
        if (!runId.isEmpty() && variants[i].runId == runId) return int(i);
    return -1;
}

QJsonObject Study::toJson() const
{
    QJsonArray vs;
    for (const StudyVariant& v : variants)
        vs.append(QJsonObject{{"value", v.value}, {"overrides", v.overrides}, {"run", v.runId}, {"status", statusKey(v.status)}});
    return {{"v", 1},
            {"id", id},
            {"number", number},
            {"name", name},
            {"parameter", parameter},
            {"parameter_label", parameterLabel},
            {"unit", unit},
            {"followers", QJsonArray::fromStringList(followers)},
            {"base", base},
            {"run_controls", runControls},
            {"approximations", QJsonArray::fromStringList(approximations)},
            {"variants", vs},
            {"created", created.toString(Qt::ISODateWithMs)}};
}

std::optional<Study> Study::fromJson(const QJsonObject& o)
{
    if (o.value("v").toInt() != 1 || o.value("id").toString().isEmpty() || o.value("parameter").toString().isEmpty()) return std::nullopt;
    Study s;
    s.id = o.value("id").toString();
    s.number = o.value("number").toInt();
    s.name = o.value("name").toString();
    s.parameter = o.value("parameter").toString();
    s.parameterLabel = o.value("parameter_label").toString();
    s.unit = o.value("unit").toString();
    for (const QJsonValue& f : o.value("followers").toArray()) s.followers << f.toString();
    s.base = o.value("base").toObject();
    s.runControls = o.value("run_controls").toObject();
    for (const QJsonValue& a : o.value("approximations").toArray()) s.approximations << a.toString();
    for (const QJsonValue& vv : o.value("variants").toArray()) {
        const QJsonObject v = vv.toObject();
        StudyVariant x{v.value("value").toDouble(), v.value("overrides").toObject(), v.value("run").toString(),
                       statusFrom(v.value("status").toString())};
        // A study saved while a variant ran (the app closed) did not finish it.
        if (x.status == VariantStatus::Running) x.status = VariantStatus::Stopped;
        s.variants.push_back(std::move(x));
    }
    s.created = QDateTime::fromString(o.value("created").toString(), Qt::ISODateWithMs);
    return s;
}

std::vector<double> parseSweepValues(const QString& text, QString* error)
{
    auto fail = [&](const QString& why) {
        if (error) *error = why;
        return std::vector<double>{};
    };
    const QString t = text.trimmed();
    if (t.isEmpty()) return fail(QStringLiteral("Enter the values to run, e.g. 16, 20, 24 or 16 to 28 step 4."));
    const QLocale c = QLocale::c();
    auto number = [&](const QString& s, double* out) {
        bool ok = false;
        *out = c.toDouble(s.trimmed(), &ok);
        return ok && std::isfinite(*out);
    };
    static const QRegularExpression range(QStringLiteral(R"(^(\S+)\s+to\s+(\S+)\s+step\s+(\S+)$)"), QRegularExpression::CaseInsensitiveOption);
    std::vector<double> out;
    if (const auto m = range.match(t); m.hasMatch()) {
        double a, b, step;
        if (!number(m.captured(1), &a) || !number(m.captured(2), &b) || !number(m.captured(3), &step))
            return fail(QStringLiteral("A range reads: start to end step size, e.g. 16 to 28 step 4."));
        if (step <= 0) return fail(QStringLiteral("The step must be positive."));
        const double n = std::floor((b - a) / step + 1e-9);
        if (b < a || n > 49) return fail(QStringLiteral("A range runs upward and has at most 50 values."));
        for (int i = 0; i <= int(n); ++i) out.push_back(a + i * step);
    } else {
        for (const QString& part : t.split(QRegularExpression(QStringLiteral(R"([,;\s]+)")), Qt::SkipEmptyParts)) {
            double v;
            if (!number(part, &v)) return fail(QStringLiteral("'%1' is not a number.").arg(part));
            out.push_back(v);
        }
    }
    for (std::size_t i = 0; i < out.size(); ++i)
        for (std::size_t j = 0; j < i; ++j)
            if (out[i] == out[j]) return fail(QStringLiteral("%1 is listed twice.").arg(c.toString(out[i], 'g', 8)));
    if (out.size() < 2) return fail(QStringLiteral("A study needs at least two values."));
    if (out.size() > 50) return fail(QStringLiteral("A study has at most 50 values."));
    if (error) error->clear();
    return out;
}

bool studyParameterAllowed(const Setting& s)
{
    return s.editable && !s.raw && !SettingsModel::isRunControl(s.path) && (s.kind == SettingKind::Number || s.kind == SettingKind::Integer);
}

QStringList referenceFollowers(const SettingsModel& model, const QString& parameter)
{
    if (!parameter.startsWith("boundary/") || !parameter.endsWith("/velocity-magnitude")) return {};
    const Setting* p = model.find(parameter);
    const Setting* ref = model.find("reference/velocity");
    if (!p || !ref || !ref->editable || !p->value.isDouble() || !ref->value.isDouble()) return {};
    return std::abs(p->value.toDouble() - ref->value.toDouble()) <= 1e-9 * std::max(1.0, std::abs(p->value.toDouble()))
               ? QStringList{"reference/velocity"}
               : QStringList{};
}

Result<Study> makeSweep(const SettingsModel& model, const QString& parameter, const std::vector<double>& values,
                        const QStringList& followers)
{
    using R = Result<Study>;
    const Setting* p = model.find(parameter);
    if (!p) return R::failure(QStringLiteral("The case has no setting %1.").arg(parameter));
    if (!studyParameterAllowed(*p))
        return R::failure(p->editable ? QStringLiteral("%1 is not a number that changes the case.").arg(p->label)
                                      : QStringLiteral("%1 cannot be changed in Kutta: %2").arg(p->label, p->reason));
    if (values.size() < 2) return R::failure(QStringLiteral("A study needs at least two values."));
    Study s;
    s.parameter = parameter;
    s.parameterLabel = p->label;
    s.unit = p->unit;
    s.followers = followers;
    s.base = model.caseOverrides();
    s.runControls = model.runControls();
    s.created = QDateTime::currentDateTime();
    for (double v : values) {
        StudyVariant var;
        var.value = v;
        for (const QString& path : QStringList{parameter} + followers) {
            const Setting* t = model.find(path);
            if (!t) return R::failure(QStringLiteral("The case has no setting %1.").arg(path));
            const QJsonValue jv = t->kind == SettingKind::Integer ? QJsonValue(qint64(std::llround(v))) : QJsonValue(v);
            if (const QString e = model.validate(path, jv); !e.isEmpty())
                return R::failure(QStringLiteral("%1 = %2: %3").arg(t->label, QLocale::c().toString(v, 'g', 8), e));
            if (!SettingsModel::sameValue(jv, t->value)) var.overrides.insert(path, jv);
        }
        s.variants.push_back(std::move(var));
    }
    return s;
}

void markStarted(Study& s, int index, const QString& runId)
{
    if (index < 0 || index >= int(s.variants.size())) return;
    s.variants[std::size_t(index)].runId = runId;
    s.variants[std::size_t(index)].status = VariantStatus::Running;
}

void markFinished(Study& s, const QString& runId, VariantStatus status)
{
    if (const int i = s.indexOfRun(runId); i >= 0) s.variants[std::size_t(i)].status = status;
}

void cancelQueued(Study& s)
{
    for (StudyVariant& v : s.variants)
        if (v.status == VariantStatus::Queued) v.status = VariantStatus::Cancelled;
}

void requeue(Study& s)
{
    for (StudyVariant& v : s.variants)
        if (v.status == VariantStatus::Cancelled || v.status == VariantStatus::Failed || v.status == VariantStatus::Stopped) {
            v.status = VariantStatus::Queued;
            v.runId.clear();
        }
}

void StudyStore::setDir(const QString& caseDir) { m_dir = caseDir.isEmpty() ? QString() : QDir(caseDir).filePath("studies"); }

std::vector<Study> StudyStore::load() const
{
    std::vector<Study> out;
    if (m_dir.isEmpty()) return out;
    for (const QString& f : QDir(m_dir).entryList({"study-*.json"}, QDir::Files, QDir::Name)) {
        QFile file(QDir(m_dir).filePath(f));
        if (!file.open(QIODevice::ReadOnly)) continue;
        if (auto s = Study::fromJson(QJsonDocument::fromJson(file.readAll()).object())) out.push_back(std::move(*s));
    }
    std::sort(out.begin(), out.end(), [](const Study& a, const Study& b) { return a.number < b.number; });
    return out;
}

bool StudyStore::save(const Study& s) const
{
    if (m_dir.isEmpty() || !QDir().mkpath(m_dir)) return false;
    QSaveFile f(QDir(m_dir).filePath(s.id + ".json"));
    return f.open(QIODevice::WriteOnly) && f.write(QJsonDocument(s.toJson()).toJson(QJsonDocument::Indented)) >= 0 && f.commit();
}

bool StudyStore::remove(const QString& id) const { return !m_dir.isEmpty() && QFile::remove(QDir(m_dir).filePath(id + ".json")); }

int StudyStore::nextNumber() const
{
    int n = 0;
    for (const Study& s : load()) n = std::max(n, s.number);
    return n + 1;
}

} // namespace kapp
