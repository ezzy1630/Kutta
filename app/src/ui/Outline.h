#pragma once

// The Fluent-order outline: Setup, Solution, Results, Runs. Each leaf is a
// page id the inspector knows how to show ("setup/boundary/<zone>").

#include "core/Session.h"

#include <QTreeWidget>

namespace kapp {

struct OutlineNode { QString id, path; };  // path: "Setup ▸ Boundary Conditions ▸ inlet"

class Outline : public QTreeWidget {
    Q_OBJECT
public:
    explicit Outline(Session* session, QWidget* parent = nullptr);
    void rebuild();
    void select(const QString& id);
    QString current() const;
    std::vector<OutlineNode> nodes() const;
    static QString titleOf(const QString& id);

signals:
    void nodeSelected(const QString& id);

private:
    void refreshMarkers();
    QTreeWidgetItem* add(QTreeWidgetItem* parent, const QString& id, const QString& text, const QString& icon = {},
                         const QString& suffix = {});
    Session* m_session;
    QHash<QString, QTreeWidgetItem*> m_items;
};

// Item data roles used by the outline's delegate.
enum OutlineRole { NodeIdRole = Qt::UserRole + 1, SuffixRole, ModifiedRole, IssueRole };

} // namespace kapp
