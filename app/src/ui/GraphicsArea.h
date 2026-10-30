#pragma once

// The graphics window: one, two or four viewports, each showing its own
// dataset and display. Picking in one viewport probes the same world point in
// the others (the linked cursor).

#include "core/Session.h"
#include "render/SceneView.h"

#include <QPointer>
#include <QWidget>

class QComboBox;
class QGridLayout;
class QLabel;
class QToolButton;

namespace kapp {

class ViewportPanel : public QWidget {
    Q_OBJECT
public:
    ViewportPanel(Session* session, int index, QWidget* parent = nullptr);

    SceneView* view() const { return m_view; }
    int index() const { return m_index; }
    QString datasetKey() const { return m_key; }
    void setDataset(const QString& key);
    std::shared_ptr<const SceneData> scene() const { return m_view->scene(); }

    const ViewDisplay& display() const { return m_view->display(); }
    void setDisplay(const ViewDisplay& d);
    void setActive(bool active);
    void refreshDatasets(const QList<QPair<QString, QString>>& datasets);  // key, label
    // Linked probe from another view; re-sampled here unless keepValue.
    void showProbe(const std::optional<PickResult>& r, bool keepValue = false);
    QImage exportImage() { return m_view->renderImage(true); }
    QString datasetLabel() const;

signals:
    void activated(kapp::ViewportPanel* panel);
    void displayChanged(kapp::ViewportPanel* panel);
    void probed(kapp::ViewportPanel* panel, const kapp::PickResult& r);
    void datasetChanged(kapp::ViewportPanel* panel);
    void exportRequested(kapp::ViewportPanel* panel);

private:
    void sceneArrived(std::shared_ptr<const SceneData> scene);
    void rebuildFieldCombo();
    void updateOverlay();
    QString probeText(const PickResult& r) const;

    Session* m_session;
    int m_index;
    QString m_key, m_pendingKey;
    SceneView* m_view;
    QComboBox* m_dataset;
    QComboBox* m_field;
    QToolButton* m_camera;
    std::optional<PickResult> m_probe;
    QStringList m_knownItems;
    bool m_active = false;
    bool m_fieldChosen = false;  // the user (or a page) picked the field
};

class GraphicsArea : public QWidget {
    Q_OBJECT
public:
    explicit GraphicsArea(Session* session, QWidget* parent = nullptr);

    enum class Layout { Single, SideBySide, Stacked, Quad };
    void setViewLayout(Layout l);
    Layout viewLayout() const { return m_layout; }
    int visibleCount() const;
    ViewportPanel* active() const { return m_active; }
    ViewportPanel* panel(int i) const { return m_panels[std::size_t(i)]; }
    void setActive(ViewportPanel* p);
    void setLinkedCameras(bool on) { m_linkCameras = on; }
    bool linkedCameras() const { return m_linkCameras; }

    // Datasets offered in every viewport besides the runs (comparisons).
    void addDataset(const QString& key, const QString& label);
    void refreshDatasets();
    // Show a dataset in the active viewport (or a given one).
    void show(const QString& key, int panelIndex = -1);

signals:
    void activeChanged(kapp::ViewportPanel* panel);
    void displayChanged(kapp::ViewportPanel* panel);

private:
    void relayout();
    Session* m_session;
    Layout m_layout = Layout::Single;
    std::array<ViewportPanel*, 4> m_panels{};
    ViewportPanel* m_active = nullptr;
    QGridLayout* m_grid;
    QList<QPair<QString, QString>> m_extra;
    bool m_linkCameras = false;
    bool m_syncing = false;
};

QStringList defaultVisibleItems(const SceneData& s);
FieldRef defaultField(const SceneData& s);

} // namespace kapp
