#pragma once

// Draws one dataset with one ViewDisplay through QRhi (Metal / D3D11 / D3D12 /
// Vulkan). Used by the on-screen viewport and by the offscreen benchmark, so
// measured frame times are the viewport's own code path.

#include "core/SceneData.h"
#include "render/ViewDisplay.h"

#include <rhi/qrhi.h>

#include <QColor>
#include <QHash>

#include <functional>
#include <memory>

namespace kapp {

struct PickResult {
    bool hit = false;
    QVector3D world;
    QString itemId, itemLabel;
    std::optional<double> value;  // display field at the hit; nullopt when missing there
};

class SceneRenderer {
public:
    SceneRenderer();
    ~SceneRenderer();

    void setScene(std::shared_ptr<const SceneData> scene);
    const std::shared_ptr<const SceneData>& scene() const { return m_scene; }
    void setDisplay(const ViewDisplay& d);
    const ViewDisplay& display() const { return m_display; }
    void setClearColor(QColor c) { m_clear = c; }
    void setMissingGrey(float g) { m_missingGrey = g; }
    void setDevicePixelRatio(float dpr) { m_dpr = dpr > 0 ? dpr : 1.0f; }
    void setLineColor(QColor c) { m_lineColor = c; }

    // (Re)create pipelines for a render target; call when rhi or format changes.
    void initialize(QRhi* rhi, QRhiRenderPassDescriptor* rp, int sampleCount);
    void releaseResources();
    bool initializedFor(QRhi* rhi) const { return m_rhi == rhi && m_surfacePipe; }

    void render(QRhiCommandBuffer* cb, QRhiRenderTarget* rt, QSize pixelSize);

    // GPU pick at a framebuffer pixel; the callback runs when the readback completes.
    void requestPick(QPointF framebufferPixel, std::function<void(const PickResult&)> done);
    bool pickPending() const { return m_pickState != PickState::Idle; }

    // Colour range actually used for the current field (for the legend).
    ColorRange resolvedRange() const { return m_range; }
    bool fieldShown() const { return m_fieldShown; }
    std::size_t trianglesDrawn() const { return m_trianglesDrawn; }
    QMatrix4x4 viewProjection(QSize pixelSize) const;
    QVector3D projectToPixel(const QVector3D& world, QSize pixelSize) const;

private:
    struct GpuMesh;
    struct DrawCmd;
    GpuMesh* gpu(const PreparedMesh* mesh, QRhiResourceUpdateBatch* u);
    QRhiBuffer* valueBuffer(GpuMesh* g, const PreparedMesh* mesh, const FieldRef& f, QRhiResourceUpdateBatch* u);
    void ensureColormap(QRhiResourceUpdateBatch* u);
    void ensureUniformCapacity(int draws);
    void ensurePickTarget();
    QRhiGraphicsPipeline* makePipeline(const char* vs, const char* fs, QRhiRenderPassDescriptor* rp, int samples,
                                       QRhiGraphicsPipeline::Topology topo, bool depthWrite, bool bias, int vertexInputKind,
                                       int targets = 1);
    QStringList visibleItems() const;
    void resolveRange();
    std::vector<DrawCmd> buildDraws(QRhiResourceUpdateBatch* u, QSize px, bool forPick);

    QRhi* m_rhi = nullptr;
    QRhiRenderPassDescriptor* m_rp = nullptr;
    int m_samples = 1;
    std::shared_ptr<const SceneData> m_scene;
    ViewDisplay m_display;
    QColor m_clear{24, 26, 29};
    float m_missingGrey = 0.55f;
    float m_dpr = 1.0f;
    QColor m_lineColor{56, 204, 110};
    ColorRange m_range;
    bool m_fieldShown = false;
    std::size_t m_trianglesDrawn = 0;

    std::unique_ptr<QRhiBuffer> m_ubuf;
    int m_ubufDraws = 0;
    quint32 m_ubufStride = 0;
    std::unique_ptr<QRhiShaderResourceBindings> m_srb;
    std::unique_ptr<QRhiTexture> m_cmapTex;
    std::unique_ptr<QRhiSampler> m_sampler;
    int m_cmapUploaded = -1;
    std::unique_ptr<QRhiGraphicsPipeline> m_surfacePipe, m_surfaceBiasPipe, m_edgePipe, m_ribbonPipe, m_ribbonHaloPipe, m_glyphPipe;
    std::unique_ptr<QRhiBuffer> m_arrowVb, m_arrowIb;
    quint32 m_arrowIndexCount = 0;
    bool m_arrowUploaded = false;

    std::unordered_map<const PreparedMesh*, std::unique_ptr<GpuMesh>> m_gpu;
    std::shared_ptr<const SceneData> m_gpuScene;  // keeps cached meshes alive

    // Picking: 1x1 float targets rendered with a pixel-centred projection.
    enum class PickState { Idle, Requested, InFlight };
    PickState m_pickState = PickState::Idle;
    QPointF m_pickPixel;
    std::function<void(const PickResult&)> m_pickDone;
    std::unique_ptr<QRhiTexture> m_pickPos, m_pickVal;
    std::unique_ptr<QRhiRenderBuffer> m_pickDepth;
    std::unique_ptr<QRhiTextureRenderTarget> m_pickRt;
    std::unique_ptr<QRhiRenderPassDescriptor> m_pickRp;
    std::unique_ptr<QRhiGraphicsPipeline> m_pickPipe;
    std::unique_ptr<QRhiBuffer> m_pickUbuf;
    int m_pickUbufDraws = 0;
    std::unique_ptr<QRhiShaderResourceBindings> m_pickSrb;
    std::unique_ptr<QRhiReadbackResult> m_rbPos, m_rbVal;
    int m_rbPending = 0;
    QStringList m_pickIds;
};

} // namespace kapp
