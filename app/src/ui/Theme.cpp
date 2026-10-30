#include "ui/Theme.h"

#include <QApplication>
#include <QFile>
#include <QHash>
#include <QIconEngine>
#include <QPainter>
#include <QPalette>
#include <QStyleFactory>
#include <QSvgRenderer>

namespace kapp {

namespace {

ThemeTokens darkTokens()
{
    ThemeTokens t;
    t.dark = true;
    t.window = QColor(0x121416);
    t.panel = QColor(0x181a1d);
    t.raised = QColor(0x202327);
    t.input = QColor(0x101214);
    t.viewport = QColor(0x15171a);
    t.border = QColor(0x2a2e33);
    t.borderStrong = QColor(0x3b4047);
    t.text = QColor(0xe3e6ea);
    t.muted = QColor(0x9aa1a9);
    t.faint = QColor(0x6c727a);
    t.accent = QColor(0x3dd17a);
    t.accentText = QColor(0x06210f);
    t.accentSoft = QColor(61, 209, 122, 46);
    t.amber = QColor(0xe9a23b);
    t.amberText = QColor(0xf0b45a);
    t.amberSoft = QColor(233, 162, 59, 40);
    t.red = QColor(0xe5645a);
    t.redSoft = QColor(229, 100, 90, 40);
    t.link = QColor(0x8cc8ff);
    // Distinct in hue and lightness on graphite; checked for CVD separation.
    t.series = {QColor(0x6aa8ff), QColor(0xf2994a), QColor(0x4fd1a5), QColor(0xe86fa8),
                QColor(0xc6b6ff), QColor(0xe8d266), QColor(0x8fd3ff), QColor(0xff8a80)};
    return t;
}

ThemeTokens lightTokens()
{
    ThemeTokens t;
    t.dark = false;
    t.window = QColor(0xeef0f2);
    t.panel = QColor(0xfbfbfc);
    t.raised = QColor(0xffffff);
    t.input = QColor(0xffffff);
    t.viewport = QColor(0xffffff);
    t.border = QColor(0xdadde1);
    t.borderStrong = QColor(0xbfc4ca);
    t.text = QColor(0x1b1e22);
    t.muted = QColor(0x59606a);
    t.faint = QColor(0x868d96);
    t.accent = QColor(0x12924a);
    t.accentText = QColor(0xffffff);
    t.accentSoft = QColor(18, 146, 74, 34);
    t.amber = QColor(0xd98a12);
    t.amberText = QColor(0x8a5300);
    t.amberSoft = QColor(217, 138, 18, 38);
    t.red = QColor(0xc8402f);
    t.redSoft = QColor(200, 64, 47, 30);
    t.link = QColor(0x1f62c4);
    t.series = {QColor(0x2f6fdb), QColor(0xd9701c), QColor(0x0f9e74), QColor(0xc63d86),
                QColor(0x7357d6), QColor(0xa88a00), QColor(0x2b9ad1), QColor(0xd2453b)};
    return t;
}

QString hex(const QColor& c) { return c.alpha() == 255 ? c.name() : QStringLiteral("rgba(%1,%2,%3,%4)").arg(c.red()).arg(c.green()).arg(c.blue()).arg(c.alpha()); }

class SvgIconEngine : public QIconEngine {
public:
    SvgIconEngine(QString name, QColor color) : m_name(std::move(name)), m_color(color) {}
    void paint(QPainter* p, const QRect& rect, QIcon::Mode mode, QIcon::State) override
    {
        QColor c = m_color;
        if (mode == QIcon::Disabled) c.setAlphaF(0.35f);
        const QPixmap pm = Theme::instance().pixmap(m_name, std::max(rect.width(), rect.height()), c, p->device()->devicePixelRatioF());
        p->drawPixmap(rect, pm);
    }
    QPixmap pixmap(const QSize& size, QIcon::Mode mode, QIcon::State state) override
    {
        QPixmap pm(size);
        pm.fill(Qt::transparent);
        QPainter p(&pm);
        paint(&p, QRect(QPoint(0, 0), size), mode, state);
        return pm;
    }
    QIconEngine* clone() const override { return new SvgIconEngine(m_name, m_color); }

private:
    QString m_name;
    QColor m_color;
};

} // namespace

Theme& Theme::instance()
{
    static Theme t;
    return t;
}

QPixmap Theme::pixmap(const QString& name, int px, QColor color, qreal dpr) const
{
    static QHash<QString, QByteArray> sources;
    static QHash<QString, QPixmap> cache;
    const QString key = QStringLiteral("%1/%2/%3/%4").arg(name).arg(px).arg(color.name(QColor::HexArgb)).arg(dpr);
    if (auto it = cache.find(key); it != cache.end()) return *it;
    if (!sources.contains(name)) {
        QFile f(QStringLiteral(":/icons/%1.svg").arg(name));
        if (!f.open(QIODevice::ReadOnly)) qWarning("missing icon %s", qPrintable(name));
        sources.insert(name, f.readAll());
    }
    QByteArray svg = sources.value(name);
    svg.replace("currentColor", color.name(QColor::HexRgb).toUtf8());
    // Slightly lighter stroke than Lucide's default reads better at 16 px.
    svg.replace("stroke-width=\"2\"", "stroke-width=\"1.75\"");
    QSvgRenderer r(svg);
    QPixmap pm(QSize(px, px) * dpr);
    pm.fill(Qt::transparent);
    QPainter p(&pm);
    p.setOpacity(color.alphaF());
    r.render(&p);
    p.end();
    pm.setDevicePixelRatio(dpr);
    cache.insert(key, pm);
    return pm;
}

QIcon Theme::icon(const QString& name, QColor color) const
{
    return QIcon(new SvgIconEngine(name, color.isValid() ? color : m_tokens.muted));
}

void Theme::apply(QApplication& app, bool dark)
{
    m_tokens = dark ? darkTokens() : lightTokens();
    app.setStyle(QStyleFactory::create("Fusion"));
    QPalette p;
    const auto& t = m_tokens;
    p.setColor(QPalette::Window, t.window);
    p.setColor(QPalette::WindowText, t.text);
    p.setColor(QPalette::Base, t.input);
    p.setColor(QPalette::AlternateBase, t.panel);
    p.setColor(QPalette::Text, t.text);
    p.setColor(QPalette::PlaceholderText, t.faint);
    p.setColor(QPalette::Button, t.raised);
    p.setColor(QPalette::ButtonText, t.text);
    // Selection fill matches the QSS item fill so tree rows highlight evenly.
    p.setColor(QPalette::Highlight, t.accentSoft);
    p.setColor(QPalette::HighlightedText, t.text);
    p.setColor(QPalette::ToolTipBase, t.raised);
    p.setColor(QPalette::ToolTipText, t.text);
    p.setColor(QPalette::Link, t.link);
    p.setColor(QPalette::Mid, t.border);
    p.setColor(QPalette::Dark, t.borderStrong);
    p.setColor(QPalette::Light, t.raised);
    p.setColor(QPalette::Accent, t.accent);
    p.setColor(QPalette::Disabled, QPalette::Text, t.faint);
    p.setColor(QPalette::Disabled, QPalette::WindowText, t.faint);
    p.setColor(QPalette::Disabled, QPalette::ButtonText, t.faint);
    app.setPalette(p);
    app.setStyleSheet(styleSheet());
    emit changed();
}

QString Theme::styleSheet() const
{
    const auto& t = m_tokens;
    QString s = R"(
* { outline: 0; }
QMainWindow, QDialog { background: %window%; }
QWidget#Panel, QScrollArea#Panel > QWidget > QWidget { background: %panel%; }
QToolTip { background: %raised%; color: %text%; border: 1px solid %borderStrong%; padding: 5px 7px; border-radius: 4px; }

QMenuBar { background: %window%; color: %text%; }
QMenu { background: %raised%; color: %text%; border: 1px solid %borderStrong%; padding: 4px; border-radius: 6px; }
QMenu::item { padding: 5px 22px 5px 10px; border-radius: 4px; }
QMenu::item:selected { background: %accentSoft%; }
QMenu::item:disabled { color: %faint%; }
QMenu::separator { height: 1px; background: %border%; margin: 4px 6px; }

QTreeView, QListView, QTableView { background: %panel%; color: %text%; border: none; selection-background-color: %accentSoft%; selection-color: %text%; alternate-background-color: %panel%; }
QTreeView::item, QListView::item { padding: 3px 2px; border: none; }
QTreeView::item:hover, QListView::item:hover { background: %hover%; }
QTreeView::item:selected, QListView::item:selected, QTableView::item:selected { background: %accentSoft%; color: %text%; }
QTreeView { show-decoration-selected: 1; }
QTreeView::branch { background: transparent; }
QTreeView::branch:selected { background: %accentSoft%; }
QHeaderView { background: %panel%; border: none; }
QHeaderView::section { background: %panel%; color: %muted%; border: none; border-bottom: 1px solid %border%; padding: 4px 8px; font-weight: 500; }
QTableView { gridline-color: %border%; }
QTableCornerButton::section { background: %panel%; border: none; }

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QPlainTextEdit, QTextEdit {
  background: %input%; color: %text%; border: 1px solid %border%; border-radius: 4px; padding: 3px 6px; selection-background-color: %accentSel%; }
QLineEdit:hover, QSpinBox:hover, QComboBox:hover { border-color: %borderStrong%; }
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QPlainTextEdit:focus { border-color: %accent%; }
QLineEdit:read-only { background: transparent; border-color: transparent; color: %muted%; }
QLineEdit:disabled, QComboBox:disabled, QSpinBox:disabled { color: %faint%; background: transparent; }
QLineEdit[invalid="true"] { border-color: %red%; }
QLineEdit[modified="true"], QComboBox[modified="true"], QSpinBox[modified="true"] { border-color: %amber%; }
QComboBox::drop-down { border: none; width: 18px; }
QComboBox QAbstractItemView { background: %raised%; border: 1px solid %borderStrong%; selection-background-color: %accentSoft%; padding: 2px; }
QSpinBox::up-button, QSpinBox::down-button { width: 0; border: none; }

QPushButton { background: %raised%; color: %text%; border: 1px solid %borderStrong%; border-radius: 4px; padding: 4px 12px; }
QPushButton:hover { border-color: %muted%; }
QPushButton:pressed { background: %hover%; }
QPushButton:disabled { color: %faint%; border-color: %border%; }
QPushButton:focus { border-color: %accent%; }
QPushButton[primary="true"] { background: %accent%; color: %accentText%; border-color: %accent%; font-weight: 600; }
QPushButton[primary="true"]:hover { background: %accentHover%; }
QPushButton[primary="true"]:disabled { background: %raised%; color: %faint%; border-color: %border%; }
QPushButton[danger="true"] { background: transparent; color: %red%; border-color: %red%; font-weight: 600; }
QPushButton[flat="true"] { background: transparent; border-color: transparent; color: %muted%; padding: 2px 6px; }
QPushButton[flat="true"]:hover { color: %text%; background: %hover%; }

QToolButton { background: transparent; color: %muted%; border: 1px solid transparent; border-radius: 4px; padding: 3px; }
QToolButton:hover { background: %hover%; color: %text%; }
QToolButton:checked { background: %accentSoft%; color: %text%; border-color: transparent; }
QToolButton:focus { border-color: %accent%; }
QToolButton::menu-indicator { image: none; }

QCheckBox, QRadioButton { color: %text%; spacing: 6px; }
QCheckBox::indicator, QRadioButton::indicator { width: 14px; height: 14px; border: 1px solid %borderStrong%; background: %input%; }
QCheckBox::indicator { border-radius: 3px; }
QRadioButton::indicator { border-radius: 7px; }
QCheckBox::indicator:checked { background: %accent%; border-color: %accent%; image: url(:/icons/check-%checkVariant%.svg); }
QRadioButton::indicator:checked { background: %input%; border: 4px solid %accent%; width: 8px; height: 8px; }
QCheckBox::indicator:focus, QRadioButton::indicator:focus { border-color: %accent%; }
QCheckBox:disabled, QRadioButton:disabled { color: %faint%; }
QListView::indicator, QTableView::indicator, QTreeView::indicator { width: 13px; height: 13px; border: 1px solid %borderStrong%; background: %input%; border-radius: 3px; }
QListView::indicator:checked, QTableView::indicator:checked, QTreeView::indicator:checked { background: %accent%; border-color: %accent%; image: url(:/icons/check-%checkVariant%.svg); }

QScrollBar:vertical { background: transparent; width: 10px; margin: 0; }
QScrollBar:horizontal { background: transparent; height: 10px; margin: 0; }
QScrollBar::handle { background: %scroll%; border-radius: 3px; min-height: 24px; min-width: 24px; margin: 2px; }
QScrollBar::handle:hover { background: %muted%; }
QScrollBar::add-line, QScrollBar::sub-line, QScrollBar::add-page, QScrollBar::sub-page { background: none; border: none; height: 0; width: 0; }

QTabBar::tab { background: transparent; color: %muted%; padding: 5px 12px; border: none; border-bottom: 2px solid transparent; }
QTabBar::tab:selected { color: %text%; border-bottom-color: %accent%; }
QTabBar::tab:hover { color: %text%; }

QSplitter::handle { background: %border%; }
QSplitter::handle:horizontal { width: 1px; }
QSplitter::handle:vertical { height: 1px; }

QStatusBar { background: %window%; color: %muted%; border-top: 1px solid %border%; }
QStatusBar::item { border: none; }
QStatusBar QLabel { color: %muted%; padding: 0 6px; }

QLabel[role="muted"] { color: %muted%; }
QLabel[role="faint"] { color: %faint%; }
QLabel[role="title"] { color: %text%; font-size: 15px; font-weight: 600; }
QLabel[role="section"] { color: %muted%; font-weight: 600; font-size: 11px; letter-spacing: 0.4px; }
QLabel[role="unit"] { color: %muted%; }
QLabel[role="error"] { color: %red%; }
QLabel[role="warning"] { color: %amberInk2%; }
QLabel[role="note"] { color: %muted%; }
QLabel[badge="approximate"] { background: %amber%; color: %amberInk%; border-radius: 3px; padding: 1px 6px; font-weight: 700; font-size: 10px; letter-spacing: 0.6px; }
QLabel[badge="modified"] { color: %amberText%; border: 1px solid %amber%; border-radius: 3px; padding: 0px 5px; font-weight: 700; font-size: 10px; letter-spacing: 0.6px; }
QLabel[badge="supported"] { color: %accentInk%; border: 1px solid %accent%; border-radius: 3px; padding: 0px 5px; font-weight: 700; font-size: 10px; letter-spacing: 0.6px; }
QLabel[badge="unsupported"] { background: %red%; color: white; border-radius: 3px; padding: 1px 6px; font-weight: 700; font-size: 10px; letter-spacing: 0.6px; }
QLabel[badge="neutral"] { color: %muted%; border: 1px solid %borderStrong%; border-radius: 3px; padding: 0px 5px; font-weight: 700; font-size: 10px; letter-spacing: 0.6px; }
QLabel[badge="live"] { color: %accentInk%; border: 1px solid %accent%; border-radius: 3px; padding: 0px 5px; font-weight: 700; font-size: 10px; letter-spacing: 0.6px; }

QFrame[callout="amber"] { background: %amberSoft%; border: 1px solid %amber%; border-radius: 6px; }
QFrame[callout="red"] { background: %redSoft%; border: 1px solid %red%; border-radius: 6px; }
QFrame[callout="neutral"] { background: %raised%; border: 1px solid %border%; border-radius: 6px; }
QFrame[rule="h"] { background: %border%; max-height: 1px; min-height: 1px; border: none; }

/* Advanced Docking System */
ads--CDockContainerWidget { background: %border%; }
ads--CDockAreaWidget { background: %panel%; border: none; }
ads--CDockAreaTitleBar { background: %window%; border-bottom: 1px solid %border%; min-height: 28px; }
ads--CDockWidgetTab { background: transparent; border: none; border-bottom: 2px solid transparent; padding: 0 6px; }
ads--CDockWidgetTab[activeTab="true"] { border-bottom-color: %accent%; }
ads--CDockWidgetTab QLabel { color: %muted%; font-weight: 500; }
ads--CDockWidgetTab[activeTab="true"] QLabel { color: %text%; }
ads--CDockWidget { background: %panel%; border: none; }
ads--CTitleBarButton, #tabsMenuButton, #detachGroupButton, #dockAreaCloseButton, #dockAreaAutoHideButton { background: transparent; border: none; padding: 2px; }
ads--CTitleBarButton:hover, #tabsMenuButton:hover, #detachGroupButton:hover, #dockAreaCloseButton:hover, #dockAreaAutoHideButton:hover { background: %hover%; border-radius: 3px; }
#tabCloseButton { background: transparent; border: none; padding: 0; margin: 0; }
#tabCloseButton:hover { background: %hover%; border-radius: 3px; }
ads--CDockSplitter::handle { background: %border%; }
ads--CDockSplitter::handle:horizontal { width: 1px; }
ads--CDockSplitter::handle:vertical { height: 1px; }
ads--CFloatingDockContainer { background: %panel%; }
ads--CAutoHideTab { background: %window%; color: %muted%; padding: 4px; border: none; }
ads--CAutoHideTab[activeTab="true"] { color: %text%; }
ads--CAutoHideSideBar { background: %window%; border: none; }
ads--CDockAreaTitleBar QToolButton { color: %muted%; }

QWidget#Toolbar { background: %window%; border-bottom: 1px solid %border%; }
QWidget#ViewportHeader { background: %panel%; border-bottom: 1px solid %border%; }
QWidget#Palette { background: %raised%; border: 1px solid %borderStrong%; border-radius: 8px; }
QWidget#Palette QLineEdit { border: none; border-bottom: 1px solid %border%; border-radius: 0; padding: 10px 12px; font-size: 15px; background: transparent; }
QWidget#Palette QListView { background: transparent; padding: 4px; }
QWidget#Palette QListView::item { padding: 6px 10px; border-radius: 4px; }
QPlainTextEdit#Console { background: %input%; border: none; border-radius: 0; padding: 6px 8px; }
)";
    const QColor hover = t.dark ? QColor(255, 255, 255, 14) : QColor(0, 0, 0, 12);
    const QColor accentHover = t.dark ? t.accent.lighter(112) : t.accent.darker(110);
    const QColor scroll = t.dark ? QColor(255, 255, 255, 40) : QColor(0, 0, 0, 46);
    const QColor accentSel(t.accent.red(), t.accent.green(), t.accent.blue(), 90);
    const QColor accentInk = t.dark ? t.accent : t.accent.darker(115);
    const std::pair<const char*, QColor> vars[] = {
        {"%window%", t.window}, {"%panel%", t.panel}, {"%raised%", t.raised}, {"%input%", t.input},
        {"%borderStrong%", t.borderStrong}, {"%border%", t.border}, {"%text%", t.text}, {"%muted%", t.muted},
        {"%faint%", t.faint}, {"%accentSoft%", t.accentSoft}, {"%accentSel%", accentSel}, {"%accentText%", t.accentText},
        {"%accentHover%", accentHover}, {"%accentInk%", accentInk}, {"%accent%", t.accent},
        {"%amberSoft%", t.amberSoft}, {"%amberInk2%", t.dark ? t.amber : t.amberText}, {"%amberText%", t.amberText}, {"%amberInk%", QColor(0x2a1d05)}, {"%amber%", t.amber},
        {"%redSoft%", t.redSoft}, {"%red%", t.red}, {"%hover%", hover}, {"%scroll%", scroll}};
    for (const auto& [k, v] : vars) s.replace(QLatin1String(k), hex(v));
    s.replace(QLatin1String("%checkVariant%"), t.dark ? QLatin1String("dark") : QLatin1String("light"));
    return s;
}

} // namespace kapp
