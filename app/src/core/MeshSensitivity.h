#pragma once

#include "core/MeshCli.h"

#include <QJsonArray>

namespace kapp {

struct MeshSetup;
class CliClient;
class SolveProcess;

class MeshSensitivity : public QObject {
    Q_OBJECT
public:
    explicit MeshSensitivity(CliClient* cli, QObject* parent);
    bool start(const QString& projectDir, const QStringList& inputs, const MeshSetup& mesh,
               double speed, const QString& casePath, const QJsonObject& overrides,
               const QJsonObject& aero, const QString& backend, int iterations, const QStringList& approximations,
               QString* error);
    void stop();
    bool running() const { return m_running; }
    QJsonObject latest(const QString& projectDir) const;

signals:
    void changed();

private:
    void nextGrid();
    void meshFinished(int exitCode, bool cancelled);
    void solveFinished(int exitCode);
    void update(const QString& state, const QString& message);
    void finish(const QString& state, const QString& message);

    CliClient* m_cli;
    MeshProcess* m_mesh = nullptr;
    SolveProcess* m_solve = nullptr;
    bool m_running = false, m_stopping = false, m_meshed = false;
    QString m_project, m_dir, m_gridDir, m_backend, m_error;
    QStringList m_inputs, m_approximations;
    QJsonObject m_case, m_record;
    QJsonArray m_bounds, m_rows;
    QStringList m_zones;
    int m_grid = 0, m_iterations = 0;
    double m_speed = 0;
};

} // namespace kapp
