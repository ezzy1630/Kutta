#pragma once

#include "render/SceneRenderer.h"

#include <QElapsedTimer>
#include <QRhiWidget>

#include <optional>

namespace kapp {

class ViewOverlay;

// Everything the overlay paints besides the 3D image. Exported images paint
// the same state, so a saved PNG carries the same labels as the screen.
struct OverlayState {
    QString title;        // "Contours of Static Pressure [Pa]"
    QString subtitle;     // dataset: "Run 4 · mixing_elbow · iteration 600"
    QString unit;
    QStringList badges;   // "APPROXIMATE", "MODIFIED", "SYNTHETIC", "LIVE · iter 120"
    QString provenance;   // full notice, shown in exports
    bool legend = false;
    ColormapId colormap = ColormapId::Viridis;
    ColorRange range;
    int bands = 0;
    bool active = false;  // focused viewport: accent frame
    std::optional<PickResult> probe;
    QPointF probePixel;   // logical pixels
    QString probeText;
    QString emptyText;    // shown when there is nothing to draw
    // Mesh displays: one swatch per shown zone (invalid colour: no swatch).
    std::vector<std::pair<QString, QColor>> zoneKey;
    QMatrix4x4 view;      // for the axis triad
    bool triad = true;
    bool dark = true;
    // The viewport's clear colour. Labels sit on translucent plates of it, so
    // they read on a full-bleed 2D field and vanish into empty background.
    QColor backdrop;
};

class SceneView : public QRhiWidget {
    Q_OBJECT
public:
    explicit SceneView(QWidget* parent = nullptr);
    ~SceneView() override;

    void setScene(std::shared_ptr<const SceneData> scene);
    const std::shared_ptr<const SceneData>& scene() const { return m_renderer.scene(); }
    void setDisplay(const ViewDisplay& d);
    const ViewDisplay& display() const { return m_display; }
    void setOverlay(const OverlayState& o);
    const OverlayState& overlayState() const { return m_overlayState; }
    void setTheme(bool dark, QColor clear, QColor accent);

    void fitView();  // frames the visible items; the view then follows resizes until the user moves it
    Bounds visibleBounds() const;
    QPointF projectToLogical(const QVector3D& world) const;
    void applyPreset(Camera::Preset p);
    void setCamera(const Camera& c);

    // Image of the viewport including its overlay (legend, labels, provenance).
    QImage renderImage(bool withProvenanceFooter);
    ColorRange resolvedRange() const { return m_renderer.resolvedRange(); }
    bool fieldShown() const { return m_renderer.fieldShown(); }
    std::size_t trianglesDrawn() const { return m_renderer.trianglesDrawn(); }

    // Benchmark: render continuously for `frames` frames and report timings.
    struct FrameTimes { int frames = 0; double meanMs = 0, p95Ms = 0, maxMs = 0; std::size_t triangles = 0; };
    void measureFrames(int frames, std::function<void(FrameTimes)> done);
    int framesRendered() const { return m_framesRendered; }

signals:
    void activated(kapp::SceneView* view);
    void cameraChanged(const kapp::Camera& camera);
    void picked(const kapp::PickResult& result, QPointF logicalPos);

protected:
    void initialize(QRhiCommandBuffer* cb) override;
    void render(QRhiCommandBuffer* cb) override;
    void releaseResources() override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void wheelEvent(QWheelEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void resizeEvent(QResizeEvent* e) override;
    void focusInEvent(QFocusEvent* e) override;

private:
    QPointF toNdc(QPointF logical) const;
    std::vector<QVector3D> visiblePoints() const;
    void fitCamera();
    float aspect() const;
    void cameraEdited();
    void ensureCamera();

    SceneRenderer m_renderer;
    ViewDisplay m_display;
    OverlayState m_overlayState;
    ViewOverlay* m_overlay = nullptr;
    QPoint m_pressPos, m_lastPos;
    Qt::MouseButton m_dragButton = Qt::NoButton;
    bool m_dragged = false;
    bool m_dark = true;
    bool m_autoFit = false;  // framed automatically and not moved since
    int m_framesRendered = 0;
    // Benchmark state (continuous rendering only while measuring).
    int m_benchRemaining = 0;
    std::vector<double> m_benchTimes;
    QElapsedTimer m_benchClock;
    std::function<void(FrameTimes)> m_benchDone;
};

// Paints the overlay; shared by the on-screen overlay widget and exports.
void paintOverlay(QPainter& p, const QRect& r, const OverlayState& s, QColor accent, bool forExport);

} // namespace kapp
