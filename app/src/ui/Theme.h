#pragma once

// Two themes from one token set. Chrome stays neutral (graphite / paper);
// colour belongs to data plates and state. One instrument accent (light-sheet
// green) marks focus, selection and the cursor; amber is reserved for
// APPROXIMATE / MODIFIED, red for blocked and errors.

#include <QColor>
#include <QIcon>
#include <QObject>
#include <QString>

class QApplication;

namespace kapp {

struct ThemeTokens {
    bool dark = true;
    QColor window, panel, raised, input, viewport;
    QColor border, borderStrong;
    QColor text, muted, faint;
    QColor accent, accentText, accentSoft;  // accentSoft: selection fill
    QColor amber, amberText, amberSoft;
    QColor red, redSoft;
    QColor link;
    std::vector<QColor> series;             // categorical chart colours
};

class Theme : public QObject {
    Q_OBJECT
public:
    static Theme& instance();
    void apply(QApplication& app, bool dark);
    void toggle(QApplication& app) { apply(app, !m_tokens.dark); }
    const ThemeTokens& tokens() const { return m_tokens; }
    bool dark() const { return m_tokens.dark; }

    // Lucide icon tinted for the current theme (or an explicit colour).
    QIcon icon(const QString& name, QColor color = {}) const;
    QPixmap pixmap(const QString& name, int px, QColor color, qreal dpr) const;

signals:
    void changed();

private:
    Theme() = default;
    QString styleSheet() const;
    ThemeTokens m_tokens;
};

inline const ThemeTokens& tokens() { return Theme::instance().tokens(); }

} // namespace kapp
