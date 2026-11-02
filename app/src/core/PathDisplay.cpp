#include "core/PathDisplay.h"

#include <QDir>
#include <QRegularExpression>

#include <algorithm>

namespace kapp {

namespace {

// Home as both separator styles may appear in text: Qt's "/" and, on Windows,
// the native "\".
QString homeOrDefault(const QString& home)
{
    QString h = QDir::fromNativeSeparators(home.isEmpty() ? QDir::homePath() : home);
    while (h.size() > 1 && h.endsWith('/')) h.chop(1);
    return h;
}

bool usableHome(const QString& h) { return h.size() > 1 && h != QLatin1String("/"); }

Qt::CaseSensitivity pathCase()
{
#ifdef Q_OS_WIN
    return Qt::CaseInsensitive;
#else
    return Qt::CaseSensitive;
#endif
}

} // namespace

QString homeRelativePath(const QString& path, const QString& home)
{
    const QString h = homeOrDefault(home);
    const QString p = QDir::fromNativeSeparators(path);
    if (usableHome(h) && p.startsWith(h, pathCase()) && (p.size() == h.size() || p.at(h.size()) == '/'))
        return QDir::toNativeSeparators("~" + p.mid(h.size()));
    return QDir::toNativeSeparators(path);
}

QString compactPath(const QString& path, int keepDirs, const QString& home)
{
    const QString rel = QDir::fromNativeSeparators(homeRelativePath(path, home));
    const bool fromHome = rel == "~" || rel.startsWith("~/");
    QStringList parts = rel.split('/', Qt::SkipEmptyParts);
    if (fromHome) parts.removeFirst();
    // Folders before the file name, beyond keepDirs, fold into "…".
    if (parts.size() - 1 <= std::max(keepDirs, 0)) return QDir::toNativeSeparators(rel);
    const QStringList kept = parts.mid(parts.size() - 1 - std::max(keepDirs, 0));
    return QDir::toNativeSeparators((fromHome ? QStringLiteral("~/…/") : QStringLiteral("…/")) + kept.join('/'));
}

QString withHomeAsTilde(const QString& text, const QString& home)
{
    const QString h = homeOrDefault(home);
    if (!usableHome(h)) return text;
    // Pages put thousands of labels through this; most carry no path.
    const QString back = QString(h).replace('/', '\\');
    if (!text.contains(h, pathCase()) && !text.contains(back, pathCase())) return text;
    QString out = text;
    for (const QString& form : {h, back}) {
        // Only the whole folder: "/Users/ana" in "/Users/anabel" is not home.
        const QRegularExpression re(QRegularExpression::escape(form) + QStringLiteral(R"((?=[/\\]|$|[^\w.\-]))"),
                                    pathCase() == Qt::CaseInsensitive ? QRegularExpression::CaseInsensitiveOption : QRegularExpression::NoPatternOption);
        out.replace(re, QStringLiteral("~"));
    }
    return out;
}

} // namespace kapp
