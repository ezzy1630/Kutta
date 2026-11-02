#include "core/CliClient.h"
#include "core/PathDisplay.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QLocale>
#include <QProcessEnvironment>
#include <QRegularExpression>
#include <QStandardPaths>

#include <optional>

#ifdef Q_OS_UNIX
#include <csignal>
#endif

namespace kapp {

QString CliInfo::describe() const
{
    return fixture ? QStringLiteral("Fixture CLI (recorded and synthetic data) · %1").arg(homeRelativePath(program))
                   : QStringLiteral("kutta CLI · %1").arg(homeRelativePath(program));
}

CliInfo locateCli(const QString& explicitPath, const QString& devRealCli, const QString& devFixture)
{
    auto make = [](const QString& p) {
        CliInfo i;
        i.program = QFileInfo(p).absoluteFilePath();
        i.fixture = QFileInfo(p).baseName().startsWith("kutta-fixture");
        return i;
    };
    if (!explicitPath.isEmpty()) return make(explicitPath);
    const QString env = qEnvironmentVariable("KUTTA_CLI");
    if (!env.isEmpty()) return make(env);
#ifdef Q_OS_WIN
    const QString suffix = ".exe";
#else
    const QString suffix;
#endif
    // Packaged builds ship the CLI in a folder of its own: the app is "Kutta"
    // and file systems on macOS and Windows ignore case, so "kutta" beside it
    // would be the app itself. macOS: Contents/Helpers; Windows: solver\.
    const QDir appDir(QCoreApplication::applicationDirPath());
    const QString self = QFileInfo(QCoreApplication::applicationFilePath()).canonicalFilePath();
    auto bundled = [&](const QString& name) -> std::optional<CliInfo> {
        for (const QString& rel : {QStringLiteral("../Helpers/"), QStringLiteral("solver/")}) {
            const QFileInfo fi(appDir.filePath(rel + name + suffix));
            if (fi.isExecutable() && fi.canonicalFilePath() != self) return make(fi.filePath());
        }
        return std::nullopt;
    };
    if (auto c = bundled("kutta")) return *c;
    if (!devRealCli.isEmpty() && QFileInfo(devRealCli).isExecutable()) return make(devRealCli);
    if (const QString onPath = QStandardPaths::findExecutable("kutta"); !onPath.isEmpty() && QFileInfo(onPath).canonicalFilePath() != self)
        return make(onPath);
    if (auto c = bundled("kutta-fixture")) return *c;
    if (!devFixture.isEmpty() && QFileInfo(devFixture).isExecutable()) {
        CliInfo f = make(devFixture);
        f.fixture = true;
        return f;
    }
    return {};
}

QString Plane::arg() const
{
    const QLocale c = QLocale::c();
    QStringList v;
    for (float x : {origin.x(), origin.y(), origin.z(), normal.x(), normal.y(), normal.z()}) v << c.toString(x, 'g', 9);
    return v.join(',');
}

QString shellQuote(const QStringList& argv)
{
    QStringList out;
    for (const QString& a : argv) {
        if (!a.isEmpty() && !a.contains(QRegularExpression("[\\s\"'$`\\\\]"))) out << a;
        else out << "\"" + QString(a).replace("\"", "\\\"") + "\"";
    }
    return out.join(' ');
}

// ---- SolveProcess -------------------------------------------------------------------

SolveProcess::SolveProcess(const CliInfo& cli, const SolveRequest& req, QObject* parent)
    : QObject(parent), m_cli(cli), m_req(req)
{
    m_proc.setProgram(cli.program);
    m_proc.setArguments(arguments());
    m_proc.setProcessChannelMode(QProcess::SeparateChannels);
    connect(&m_proc, &QProcess::readyReadStandardOutput, this, &SolveProcess::readStdout);
    connect(&m_proc, &QProcess::readyReadStandardError, this,
            [this] { emit stderrText(cliMessage(QString::fromUtf8(m_proc.readAllStandardError()))); });
    connect(&m_proc, &QProcess::finished, this, [this](int code, QProcess::ExitStatus st) {
        m_stopEscalation.stop();
        readStdout();
        if (!m_partial.trimmed().isEmpty()) {
            emit progressLine(m_partial.trimmed());
            auto ev = parseProgressLine(m_partial);
            if (ev) emit progress(ev.value());
            m_partial.clear();
        }
        emit finished(code, st == QProcess::CrashExit);
    });
    connect(&m_proc, &QProcess::errorOccurred, this, [this](QProcess::ProcessError e) {
        if (e == QProcess::FailedToStart)
            emit protocolError(QStringLiteral("Could not start %1: %2").arg(m_cli.program, m_proc.errorString()));
    });
    // A CLI that ignores stdin EOF gets an interrupt, then is killed. Only the
    // first path writes outputs; the console says which one happened.
    m_stopEscalation.setSingleShot(true);
    connect(&m_stopEscalation, &QTimer::timeout, this, [this] {
        if (!running()) return;
#ifdef Q_OS_UNIX
        emit stderrText(QStringLiteral("kutta did not stop after stdin closed; sending SIGINT\n"));
        ::kill(pid_t(m_proc.processId()), SIGINT);
#else
        emit stderrText(QStringLiteral("kutta did not stop after stdin closed; terminating\n"));
        m_proc.kill();
#endif
    });
}

SolveProcess::~SolveProcess()
{
    if (running()) {
        m_proc.closeWriteChannel();
        if (!m_proc.waitForFinished(3000)) m_proc.kill();
        m_proc.waitForFinished(1000);
    }
}

QStringList SolveProcess::arguments() const
{
    QStringList a{"solve", m_req.casePath};
    if (!m_req.dataPath.isEmpty()) a << "--data" << m_req.dataPath;
    if (!m_req.approximations.isEmpty()) a << "--approximate" << m_req.approximations.join(',');
    if (!m_req.backend.isEmpty()) a << "--backend" << m_req.backend;
    if (m_req.iterations > 0) a << "--iterations" << QString::number(m_req.iterations);
    if (!m_req.setupFile.isEmpty()) a << "--setup" << m_req.setupFile;
    if (m_req.liveViewInterval > 0)
        a << "--live-view" << m_req.outDir << "--live-every" << QString::number(m_req.liveViewInterval);
    a << "--out" << m_req.outDir << "--progress" << "jsonl";
    return a;
}

void SolveProcess::start()
{
    // stdin stays open for the life of the run: EOF is the stop request.
    m_proc.start(QIODevice::ReadWrite);
}

void SolveProcess::stop()
{
    if (!running()) return;
    m_proc.closeWriteChannel();
    m_stopEscalation.start(20000);
}

void SolveProcess::readStdout()
{
    m_partial += m_proc.readAllStandardOutput();
    qsizetype nl;
    while ((nl = m_partial.indexOf('\n')) >= 0) {
        const QByteArray line = m_partial.left(nl).trimmed();
        m_partial.remove(0, nl + 1);
        if (line.isEmpty()) continue;
        emit progressLine(line);
        auto ev = parseProgressLine(line);
        if (ev) emit progress(ev.value());
        else emit protocolError(QStringLiteral("Unreadable progress line (%1): %2")
                                    .arg(ev.error(), QString::fromUtf8(line.left(200))));
    }
}

// ---- CliClient ----------------------------------------------------------------------

CliClient::CliClient(CliInfo info, QObject* parent) : QObject(parent), m_info(std::move(info)) {}

void CliClient::run(const QStringList& args, RawDone done)
{
    auto* p = new QProcess(this);
    p->setProgram(m_info.program);
    p->setArguments(args);
    emit commandStarted(shellQuote(QStringList{QFileInfo(m_info.program).fileName()} + args));
    connect(p, &QProcess::finished, this, [this, p, done](int code, QProcess::ExitStatus st) {
        const QByteArray out = p->readAllStandardOutput(), err = p->readAllStandardError();
        if (!err.isEmpty()) emit commandOutput(QString::fromUtf8(err));
        p->deleteLater();
        done(st == QProcess::CrashExit ? -1 : code, out, err);
    });
    connect(p, &QProcess::errorOccurred, this, [this, p, done](QProcess::ProcessError e) {
        if (e != QProcess::FailedToStart) return;
        const QString msg = QStringLiteral("Could not start %1: %2").arg(m_info.program, p->errorString());
        emit commandOutput(msg + "\n");
        p->deleteLater();
        done(-1, {}, msg.toUtf8());
    });
    p->start(QIODevice::ReadOnly);
}

namespace {

QString failureText(int code, const QByteArray& out, const QByteArray& err)
{
    // Prefer the protocol's own error event over raw stderr.
    for (const QByteArray& line : out.split('\n')) {
        auto ev = parseProgressLine(line);
        if (ev && std::holds_alternative<ErrorEvent>(ev.value())) return std::get<ErrorEvent>(ev.value()).message;
    }
    const QString e = cliMessage(QString::fromUtf8(err).trimmed());
    return e.isEmpty() ? QStringLiteral("kutta exited with code %1").arg(code) : e;
}

Result<QJsonObject> jsonObject(const QByteArray& out)
{
    QJsonParseError pe;
    const QJsonDocument doc = QJsonDocument::fromJson(out, &pe);
    if (!doc.isObject()) return Result<QJsonObject>::failure("unreadable JSON from kutta: " + pe.errorString());
    return doc.object();
}

} // namespace

void CliClient::devices(std::function<void(Result<std::vector<Device>>)> done)
{
    run({"devices", "--json"}, [done](int code, const QByteArray& out, const QByteArray& err) {
        if (code != 0) return done(Result<std::vector<Device>>::failure(failureText(code, out, err)));
        auto obj = jsonObject(out);
        if (!obj) return done(Result<std::vector<Device>>::failure(obj.error()));
        done(parseDevices(obj.value()));
    });
}

void CliClient::inspect(const QString& casePath, const QString& dataPath, std::function<void(Result<CaseSummary>)> done)
{
    QStringList args{"inspect", casePath};
    if (!dataPath.isEmpty()) args << "--data" << dataPath;
    args << "--json";
    run(args, [done](int code, const QByteArray& out, const QByteArray& err) {
        if (code != 0) return done(Result<CaseSummary>::failure(failureText(code, out, err)));
        auto obj = jsonObject(out);
        if (!obj) return done(Result<CaseSummary>::failure(obj.error()));
        done(parseCaseSummary(obj.value()));
    });
}

void CliClient::exportView(const QString& runDir, const std::vector<Plane>& slices, const QString& postFile,
                           std::function<void(Result<ViewManifest>)> done)
{
    QStringList args{"export-view", runDir};
    if (!postFile.isEmpty()) args << "--post" << postFile;
    for (const Plane& p : slices) args << "--slice" << p.arg();
    run(args, [done, runDir](int code, const QByteArray& out, const QByteArray& err) {
        if (code != 0) return done(Result<ViewManifest>::failure(failureText(code, out, err)));
        done(loadViewManifest(QDir(runDir).filePath("view/manifest.json")));
    });
}

void CliClient::exportGeometry(const QString& casePath, const QString& outDir, std::function<void(Result<ViewManifest>)> done)
{
    run({"export-geometry", casePath, "--out", outDir, "--json"}, [done, outDir](int code, const QByteArray& out, const QByteArray& err) {
        // Exit 2 is a usage error: a CLI from before export-geometry.
        if (code == 2) return done(Result<ViewManifest>::failure("this kutta CLI cannot show a case's mesh before a run (no export-geometry)"));
        if (code != 0) return done(Result<ViewManifest>::failure(failureText(code, out, err)));
        done(loadViewManifest(QDir(outDir).filePath("geometry/manifest.json")));
    });
}

void CliClient::aero(const QString& runDir, const QString& specFile, std::function<void(Result<AeroResult>)> done)
{
    QStringList args{"aero", runDir};
    if (!specFile.isEmpty()) args << "--aero" << specFile;
    args << "--json";
    run(args, [done](int code, const QByteArray& out, const QByteArray& err) {
        // Exit 2 from a CLI without the command reads as a usage message.
        if (code != 0) return done(Result<AeroResult>::failure(failureText(code, out, err)));
        for (const QByteArray& line : out.split('\n')) {
            const QJsonObject o = QJsonDocument::fromJson(line).object();
            if (o.value("event").toString() == "aero") return done(parseAero(o));
        }
        done(Result<AeroResult>::failure("kutta aero printed no result"));
    });
}

void CliClient::newCase(const NewCaseRequest& req, std::function<void(Result<NewCaseDraft>)> done)
{
    QStringList args{"new-case", req.meshPath, "--json"};
    if (req.speed > 0) args << "--speed" << QString::number(req.speed, 'g', 17);
    for (auto it = req.roles.begin(); it != req.roles.end(); ++it) args << "--role" << it.key() + "=" + it.value();
    if (!req.out.isEmpty()) args << "--out" << req.out << "--force";
    run(args, [done](int code, const QByteArray& out, const QByteArray& err) {
        // Exit 2 is a usage error: a CLI from before new-case.
        if (code == 2 && !out.contains("\"event\"")) return done(Result<NewCaseDraft>::failure("this kutta CLI cannot set up a case from a mesh (no new-case)"));
        if (code != 0) return done(Result<NewCaseDraft>::failure(failureText(code, out, err)));
        auto obj = jsonObject(out);
        if (!obj) return done(Result<NewCaseDraft>::failure(obj.error()));
        done(parseNewCase(obj.value()));
    });
}

void CliClient::convert(const QString& casePath, const QString& setupFile, const QString& out, std::function<void(Result<ConvertResult>)> done)
{
    QStringList args{"convert", casePath, "--out", out, "--force", "--json"};
    if (!setupFile.isEmpty()) args << "--setup" << setupFile;
    run(args, [done](int code, const QByteArray& stdoutText, const QByteArray& err) {
        if (code == 2 && !stdoutText.contains("\"event\"")) return done(Result<ConvertResult>::failure("this kutta CLI cannot save Kutta cases (no convert)"));
        if (code != 0) return done(Result<ConvertResult>::failure(failureText(code, stdoutText, err)));
        auto obj = jsonObject(stdoutText);
        if (!obj) return done(Result<ConvertResult>::failure(obj.error()));
        done(parseConvert(obj.value()));
    });
}

SolveProcess* CliClient::solve(const SolveRequest& req, QObject* parent)
{
    auto* s = new SolveProcess(m_info, req, parent);
    emit commandStarted(shellQuote(QStringList{QFileInfo(m_info.program).fileName()} + s->arguments()));
    return s;
}

} // namespace kapp
