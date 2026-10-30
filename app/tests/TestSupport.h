#pragma once

// Helpers shared by the app's unit tests: fixture data lives qCompress-ed.

#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QString>

inline QByteArray fixtureBytes(const QString& rel)
{
    const QString base = QStringLiteral(KUTTA_FIXTURE_DIR "/") + rel;
    QFile plain(base);
    if (plain.open(QIODevice::ReadOnly)) return plain.readAll();
    QFile packed(base + ".qz");
    if (packed.open(QIODevice::ReadOnly)) return qUncompress(packed.readAll());
    return {};
}

inline QJsonObject fixtureJson(const QString& rel) { return QJsonDocument::fromJson(fixtureBytes(rel)).object(); }
