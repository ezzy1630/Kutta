#include "ui/Outline.h"

#include "ui/Theme.h"

#include <QHeaderView>
#include <QPainter>
#include <QStyledItemDelegate>

namespace kapp {

namespace {

// Paints the zone type as a faint suffix and an amber dot on modified pages.
class OutlineDelegate : public QStyledItemDelegate {
public:
    using QStyledItemDelegate::QStyledItemDelegate;
    void paint(QPainter* p, const QStyleOptionViewItem& opt, const QModelIndex& idx) const override
    {
        QStyleOptionViewItem o = opt;
        initStyleOption(&o, idx);
        const QString suffix = idx.data(SuffixRole).toString();
        const int modified = idx.data(ModifiedRole).toInt();
        const int issue = idx.data(IssueRole).toInt();  // 1 approximable, 2 unsupported
        const auto& t = tokens();
        const int reserve = (modified ? 26 : 0) + (issue ? 18 : 0);
        o.rect.adjust(0, 0, -reserve, 0);
        QStyledItemDelegate::paint(p, o, idx);
        p->save();
        p->setRenderHint(QPainter::Antialiasing);
        const QFontMetrics fm(o.font);
        if (!suffix.isEmpty()) {
            const int textW = fm.horizontalAdvance(o.text);
            const int iconW = o.icon.isNull() ? 0 : o.decorationSize.width() + 6;
            QRect r = o.rect.adjusted(iconW + textW + 12, 0, 0, 0);
            if (r.width() > 20) {
                p->setPen(t.faint);
                p->drawText(r, Qt::AlignVCenter | Qt::AlignLeft, fm.elidedText(suffix, Qt::ElideRight, r.width()));
            }
        }
        int x = opt.rect.right() - 6;
        if (modified) {
            const QString n = QString::number(modified);
            const int w = std::max(16, fm.horizontalAdvance(n) + 8);
            const QRect r(x - w, opt.rect.center().y() - 8, w, 16);
            p->setPen(QPen(t.amber, 1));
            p->setBrush(Qt::NoBrush);
            p->drawRoundedRect(QRectF(r).adjusted(0.5, 0.5, -0.5, -0.5), 8, 8);
            p->setPen(t.dark ? t.amber : t.amberText);
            QFont f = o.font;
            f.setPointSizeF(f.pointSizeF() * 0.85);
            f.setBold(true);
            p->setFont(f);
            p->drawText(r, Qt::AlignCenter, n);
            x -= w + 4;
        }
        if (issue) {
            const QPixmap pm = Theme::instance().pixmap(issue == 2 ? "circle-alert" : "triangle-alert", 14,
                                                        issue == 2 ? t.red : t.amber, p->device()->devicePixelRatioF());
            p->drawPixmap(QPoint(x - 14, opt.rect.center().y() - 7), pm);
        }
        p->restore();
    }
    QSize sizeHint(const QStyleOptionViewItem& opt, const QModelIndex& idx) const override
    {
        QSize s = QStyledItemDelegate::sizeHint(opt, idx);
        s.setHeight(std::max(s.height(), 24));
        return s;
    }
};

const QHash<QString, QString>& titles()
{
    static const QHash<QString, QString> t = {
        {"setup", "Setup"},
        {"setup/general", "General"},
        {"setup/models", "Models"},
        {"setup/materials", "Materials"},
        {"setup/cell-zones", "Cell Zone Conditions"},
        {"setup/boundaries", "Boundary Conditions"},
        {"setup/reference", "Reference Values"},
        {"solution", "Solution"},
        {"solution/methods", "Methods"},
        {"solution/controls", "Controls"},
        {"solution/reports", "Report Definitions"},
        {"solution/monitors", "Monitors"},
        {"solution/initialization", "Initialization"},
        {"solution/run", "Run Calculation"},
        {"results", "Results"},
        {"results/graphics", "Graphics"},
        {"results/contours", "Contours"},
        {"results/vectors", "Vectors"},
        {"results/pathlines", "Pathlines"},
        {"results/iso", "Iso-Surfaces"},
        {"results/mesh", "Mesh"},
        {"results/planes", "Planes"},
        {"results/plots", "Plots"},
        {"results/xy", "XY Plot"},
        {"results/residuals", "Residuals"},
        {"results/report-plots", "Report Plots"},
        {"results/reports", "Reports"},
        {"results/forces", "Forces"},
        {"results/aero", "Aerodynamics"},
        {"results/integrals", "Surface Integrals"},
        {"runs", "Runs"},
        {"runs/history", "Run History"},
        {"runs/compare", "Compare Runs"},
        {"runs/studies", "Design Studies"},
        {"settings", "All Case Settings"},
    };
    return t;
}

} // namespace

QString Outline::titleOf(const QString& id)
{
    if (id.startsWith("setup/boundary/")) return id.mid(15);
    if (id.startsWith("setup/material/")) return id.mid(15);
    return titles().value(id, id);
}

Outline::Outline(Session* session, QWidget* parent) : QTreeWidget(parent), m_session(session)
{
    setHeaderHidden(true);
    setColumnCount(1);
    setRootIsDecorated(true);
    setIndentation(14);
    setUniformRowHeights(true);
    setItemDelegate(new OutlineDelegate(this));
    setAccessibleName(tr("Outline"));
    setIconSize(QSize(16, 16));
    setExpandsOnDoubleClick(true);
    connect(this, &QTreeWidget::currentItemChanged, this, [this](QTreeWidgetItem* cur) {
        if (!cur) return;
        QString id = cur->data(0, NodeIdRole).toString();
        // Group headers open their first page.
        if (cur->childCount() > 0 && !id.contains('/')) {
            if (id == "runs") id = "runs/history";
        }
        emit nodeSelected(id);
    });
    connect(session, &Session::caseOpened, this, &Outline::rebuild);
    connect(session->settings(), &SettingsModel::overridesChanged, this, &Outline::refreshMarkers);
    connect(session, &Session::approvalsChanged, this, &Outline::refreshMarkers);
    connect(session, &Session::historyChanged, this, &Outline::refreshMarkers);
    connect(&Theme::instance(), &Theme::changed, this, &Outline::rebuild);
    rebuild();
}

QTreeWidgetItem* Outline::add(QTreeWidgetItem* parent, const QString& id, const QString& text, const QString& icon, const QString& suffix)
{
    auto* it = parent ? new QTreeWidgetItem(parent) : new QTreeWidgetItem(this);
    it->setText(0, text);
    it->setData(0, NodeIdRole, id);
    it->setData(0, SuffixRole, suffix);
    if (!icon.isEmpty()) it->setIcon(0, Theme::instance().icon(icon));
    it->setData(0, Qt::AccessibleTextRole, suffix.isEmpty() ? text : text + ", " + suffix);
    m_items.insert(id, it);
    return it;
}

void Outline::rebuild()
{
    const QString cur = current();
    const QSignalBlocker b(this);
    clear();
    m_items.clear();
    const bool hasCase = m_session->hasCase();

    auto* setup = add(nullptr, "setup", tr("Setup"), "sliders-horizontal");
    add(setup, "setup/general", tr("General"));
    add(setup, "setup/models", tr("Models"));
    auto* mats = add(setup, "setup/materials", tr("Materials"));
    add(setup, "setup/cell-zones", tr("Cell Zone Conditions"));
    auto* bcs = add(setup, "setup/boundaries", tr("Boundary Conditions"));
    add(setup, "setup/reference", tr("Reference Values"));
    if (hasCase) {
        const auto& c = m_session->caseSummary();
        for (const auto& m : c.materials) add(mats, "setup/material/" + m.name, m.name);
        for (const auto& z : c.zones)
            if (z.type != "interior") add(bcs, "setup/boundary/" + z.name, z.name, {}, z.type);
    }
    auto* sol = add(nullptr, "solution", tr("Solution"), "activity");
    add(sol, "solution/methods", tr("Methods"));
    add(sol, "solution/controls", tr("Controls"));
    add(sol, "solution/reports", tr("Report Definitions"));
    add(sol, "solution/monitors", tr("Monitors"));
    add(sol, "solution/initialization", tr("Initialization"));
    add(sol, "solution/run", tr("Run Calculation"));
    auto* res = add(nullptr, "results", tr("Results"), "chart-spline");
    auto* gfx = add(res, "results/graphics", tr("Graphics"));
    add(gfx, "results/contours", tr("Contours"));
    add(gfx, "results/vectors", tr("Vectors"));
    add(gfx, "results/pathlines", tr("Pathlines"));
    add(gfx, "results/iso", tr("Iso-Surfaces"));
    add(gfx, "results/planes", tr("Planes"));
    add(gfx, "results/mesh", tr("Mesh"));
    auto* plots = add(res, "results/plots", tr("Plots"));
    add(plots, "results/xy", tr("XY Plot"));
    add(plots, "results/residuals", tr("Residuals"));
    add(plots, "results/report-plots", tr("Report Plots"));
    auto* reps = add(res, "results/reports", tr("Reports"));
    add(reps, "results/forces", tr("Forces"));
    add(reps, "results/aero", tr("Aerodynamics"));
    add(reps, "results/integrals", tr("Surface Integrals"));
    auto* runs = add(nullptr, "runs", tr("Runs"), "history");
    add(runs, "runs/history", tr("Run History"));
    add(runs, "runs/compare", tr("Compare Runs"));
    add(runs, "runs/studies", tr("Design Studies"));
    add(nullptr, "settings", tr("All Case Settings"), "list-filter");

    expandAll();
    for (const char* collapsed : {"setup/materials"})
        if (auto* it = m_items.value(collapsed)) it->setExpanded(false);
    refreshMarkers();
    if (auto* it = m_items.value(cur.isEmpty() ? QString("setup/general") : cur)) setCurrentItem(it);
}

void Outline::refreshMarkers()
{
    if (!m_session->hasCase()) return;
    const QStringList mod = m_session->settings()->modifiedPaths();
    auto count = [&](const QString& prefix) {
        int n = 0;
        for (const QString& p : mod) n += p.startsWith(prefix);
        return n;
    };
    for (auto it = m_items.begin(); it != m_items.end(); ++it) it.value()->setData(0, ModifiedRole, 0);
    auto set = [&](const QString& id, int n) {
        if (auto* it = m_items.value(id)) it->setData(0, ModifiedRole, n);
    };
    const auto& c = m_session->caseSummary();
    for (const auto& z : c.zones) set("setup/boundary/" + z.name, count("boundary/" + z.name + "/"));
    for (const auto& m : c.materials) set("setup/material/" + m.name, count("materials/" + m.name + "/"));
    set("setup/boundaries", count("boundary/"));
    set("setup/materials", count("materials/"));
    set("setup/reference", count("reference/"));
    set("setup/models", count("models/"));
    set("setup/cell-zones", count("cell-zone/"));
    set("solution/methods", count("methods/"));
    set("solution/controls", count("controls/"));
    set("solution/run", count("run/"));
    set("solution/monitors", count("monitors/"));
    set("solution/initialization", count("initialization/"));

    // Support issues: unsupported (red) or approximable and not yet approved (amber).
    int general = 0;
    for (const auto* i : c.blocking(m_session->approvedApproximations()))
        general = std::max(general, i->cls == IssueClass::Unsupported ? 2 : 1);
    if (auto* it = m_items.value("setup/general")) it->setData(0, IssueRole, general);
    viewport()->update();
}

void Outline::select(const QString& id)
{
    if (auto* it = m_items.value(id)) {
        setCurrentItem(it);
        scrollToItem(it);
    } else {
        emit nodeSelected(id);
    }
}

QString Outline::current() const
{
    return currentItem() ? currentItem()->data(0, NodeIdRole).toString() : QString();
}

std::vector<OutlineNode> Outline::nodes() const
{
    std::vector<OutlineNode> out;
    std::function<void(QTreeWidgetItem*, QString)> walk = [&](QTreeWidgetItem* it, QString prefix) {
        const QString path = prefix.isEmpty() ? it->text(0) : prefix + " ▸ " + it->text(0);
        out.push_back({it->data(0, NodeIdRole).toString(), path});
        for (int i = 0; i < it->childCount(); ++i) walk(it->child(i), path);
    };
    for (int i = 0; i < topLevelItemCount(); ++i) walk(topLevelItem(i), {});
    return out;
}

} // namespace kapp
