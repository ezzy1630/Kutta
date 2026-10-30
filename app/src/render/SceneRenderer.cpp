#include "render/SceneRenderer.h"

#include <QFile>
#include <QtMath>

#include <cmath>
#include <cstring>

namespace kapp {

namespace {

struct DrawUniforms {
    float mvp[16];
    float view[16];
    float range[4];
    float solid[4];
    float params[4];
    float viewport[4];
};
static_assert(sizeof(DrawUniforms) == 192);

QShader loadShader(const char* name)
{
    QFile f(QStringLiteral(":/shaders/%1.qsb").arg(QLatin1String(name)));
    if (!f.open(QIODevice::ReadOnly)) qFatal("missing shader %s", name);
    return QShader::fromSerialized(f.readAll());
}

enum VertexInput { SurfaceInput = 0, LineInput = 1, RibbonInput = 2, GlyphInput = 3 };

// Arrow along +x from 0 to 1: hexagonal shaft and cone, radius in width units.
void arrowMesh(std::vector<float>* v, std::vector<quint32>* idx)
{
    constexpr int sides = 8;
    const float shaftR = 0.28f, coneR = 0.9f, shaftEnd = 0.62f;
    auto ring = [&](float x, float r, float nx, float nr) {
        const quint32 start = quint32(v->size() / 6);
        for (int i = 0; i < sides; ++i) {
            const float a = float(i) / sides * 2 * float(M_PI);
            const float c = std::cos(a), s = std::sin(a);
            v->insert(v->end(), {x, r * c, r * s, nx, nr * c, nr * s});
        }
        return start;
    };
    auto quads = [&](quint32 a, quint32 b) {
        for (int i = 0; i < sides; ++i) {
            const quint32 i0 = a + quint32(i), i1 = a + quint32((i + 1) % sides), j0 = b + quint32(i), j1 = b + quint32((i + 1) % sides);
            idx->insert(idx->end(), {i0, j0, i1, i1, j0, j1});
        }
    };
    const quint32 s0 = ring(0, shaftR, 0, 1), s1 = ring(shaftEnd, shaftR, 0, 1);
    quads(s0, s1);
    const float slope = coneR / (1 - shaftEnd);
    const float nlen = std::sqrt(1 + slope * slope);
    const quint32 c0 = ring(shaftEnd, coneR, slope / nlen, 1 / nlen), c1 = ring(1.0f, 0.0f, slope / nlen, 1 / nlen);
    quads(c0, c1);
    // Back face of the cone.
    const quint32 b0 = ring(shaftEnd, coneR, -1, 0);
    const quint32 centre = quint32(v->size() / 6);
    v->insert(v->end(), {shaftEnd, 0, 0, -1, 0, 0});
    for (int i = 0; i < sides; ++i) idx->insert(idx->end(), {centre, b0 + quint32((i + 1) % sides), b0 + quint32(i)});
}

} // namespace

struct SceneRenderer::GpuMesh {
    std::unique_ptr<QRhiBuffer> vb, ib, edgeIb;
    quint32 indexCount = 0, edgeCount = 0;
    QHash<QString, std::shared_ptr<QRhiBuffer>> values;  // by field key
    struct Ribbon { std::unique_ptr<QRhiBuffer> vb, ib; quint32 count = 0; };
    QHash<QString, std::shared_ptr<Ribbon>> ribbons;
    struct Glyphs { std::unique_ptr<QRhiBuffer> buf; quint32 count = 0; };
    QHash<QString, std::shared_ptr<Glyphs>> glyphs;
};

struct SceneRenderer::DrawCmd {
    enum Kind { Surface, Edges, Ribbon, Glyph } kind;
    QRhiGraphicsPipeline* pipe = nullptr;
    QRhiBuffer* vb = nullptr;
    QRhiBuffer* vb2 = nullptr;
    QRhiBuffer* ib = nullptr;
    quint32 count = 0, instances = 1;
    DrawUniforms u{};
    int objectId = -1;
};

SceneRenderer::SceneRenderer() = default;
SceneRenderer::~SceneRenderer() { releaseResources(); }

void SceneRenderer::setScene(std::shared_ptr<const SceneData> scene)
{
    if (scene == m_scene) return;
    m_scene = std::move(scene);
    resolveRange();
}

// The range is resolved eagerly so the legend never shows a previous field's range.
void SceneRenderer::setDisplay(const ViewDisplay& d)
{
    m_display = d;
    resolveRange();
}

void SceneRenderer::releaseResources()
{
    m_gpu.clear();
    m_gpuScene.reset();
    m_surfacePipe.reset();
    m_surfaceBiasPipe.reset();
    m_edgePipe.reset();
    m_ribbonPipe.reset();
    m_glyphPipe.reset();
    m_srb.reset();
    m_ubuf.reset();
    m_ubufDraws = 0;
    m_cmapTex.reset();
    m_sampler.reset();
    m_cmapUploaded = -1;
    m_arrowVb.reset();
    m_arrowIb.reset();
    m_arrowUploaded = false;
    m_pickPipe.reset();
    m_pickSrb.reset();
    m_pickUbuf.reset();
    m_pickUbufDraws = 0;
    m_pickRt.reset();
    m_pickRp.reset();
    m_pickPos.reset();
    m_pickVal.reset();
    m_pickDepth.reset();
    if (m_pickState != PickState::Idle && m_pickDone) m_pickDone(PickResult{});
    m_pickState = PickState::Idle;
    m_rhi = nullptr;
}

QRhiGraphicsPipeline* SceneRenderer::makePipeline(const char* vsName, const char* fsName, QRhiRenderPassDescriptor* rp,
                                                  int samples, QRhiGraphicsPipeline::Topology topo, bool depthWrite,
                                                  bool bias, int input, int targets)
{
    QRhiGraphicsPipeline* ps = m_rhi->newGraphicsPipeline();
    ps->setShaderStages({{QRhiShaderStage::Vertex, loadShader(vsName)}, {QRhiShaderStage::Fragment, loadShader(fsName)}});
    QRhiVertexInputLayout il;
    switch (input) {
    case SurfaceInput:
        // Value is a separate binding so one geometry buffer serves every field.
        il.setBindings({{6 * sizeof(float)}, {sizeof(float)}});
        il.setAttributes({{0, 0, QRhiVertexInputAttribute::Float3, 0},
                          {0, 1, QRhiVertexInputAttribute::Float3, 3 * sizeof(float)},
                          {1, 2, QRhiVertexInputAttribute::Float, 0}});
        break;
    case LineInput:
        il.setBindings({{6 * sizeof(float)}});
        il.setAttributes({{0, 0, QRhiVertexInputAttribute::Float3, 0}});
        break;
    case RibbonInput:
        il.setBindings({{8 * sizeof(float)}});
        il.setAttributes({{0, 0, QRhiVertexInputAttribute::Float3, 0},
                          {0, 1, QRhiVertexInputAttribute::Float3, 3 * sizeof(float)},
                          {0, 2, QRhiVertexInputAttribute::Float, 6 * sizeof(float)},
                          {0, 3, QRhiVertexInputAttribute::Float, 7 * sizeof(float)}});
        break;
    case GlyphInput:
        il.setBindings({{6 * sizeof(float)}, {7 * sizeof(float), QRhiVertexInputBinding::PerInstance}});
        il.setAttributes({{0, 0, QRhiVertexInputAttribute::Float3, 0},
                          {0, 1, QRhiVertexInputAttribute::Float3, 3 * sizeof(float)},
                          {1, 2, QRhiVertexInputAttribute::Float3, 0},
                          {1, 3, QRhiVertexInputAttribute::Float3, 3 * sizeof(float)},
                          {1, 4, QRhiVertexInputAttribute::Float, 6 * sizeof(float)}});
        break;
    }
    ps->setVertexInputLayout(il);
    ps->setShaderResourceBindings(targets == 1 ? m_srb.get() : m_pickSrb.get());
    ps->setRenderPassDescriptor(rp);
    ps->setSampleCount(samples);
    ps->setTopology(topo);
    ps->setDepthTest(true);
    ps->setDepthWrite(depthWrite);
    ps->setDepthOp(QRhiGraphicsPipeline::LessOrEqual);
    ps->setCullMode(QRhiGraphicsPipeline::None);
    if (bias) {
        // Push filled surfaces back so their own edges draw on top.
        ps->setDepthBias(1);
        ps->setSlopeScaledDepthBias(1.0f);
    }
    if (targets == 1) {
        QRhiGraphicsPipeline::TargetBlend blend;
        blend.enable = true;  // premultiplied alpha
        blend.srcColor = QRhiGraphicsPipeline::One;
        blend.dstColor = QRhiGraphicsPipeline::OneMinusSrcAlpha;
        blend.srcAlpha = QRhiGraphicsPipeline::One;
        blend.dstAlpha = QRhiGraphicsPipeline::OneMinusSrcAlpha;
        ps->setTargetBlends({blend});
    } else {
        ps->setTargetBlends({QRhiGraphicsPipeline::TargetBlend{}, QRhiGraphicsPipeline::TargetBlend{}});
    }
    if (!ps->create()) qWarning("pipeline %s/%s failed", vsName, fsName);
    return ps;
}

void SceneRenderer::initialize(QRhi* rhi, QRhiRenderPassDescriptor* rp, int sampleCount)
{
    if (m_rhi != rhi) releaseResources();
    m_rhi = rhi;
    m_rp = rp;
    m_samples = sampleCount;
    if (!m_sampler) {
        m_sampler.reset(rhi->newSampler(QRhiSampler::Linear, QRhiSampler::Linear, QRhiSampler::None,
                                        QRhiSampler::ClampToEdge, QRhiSampler::ClampToEdge));
        m_sampler->create();
        m_cmapTex.reset(rhi->newTexture(QRhiTexture::RGBA8, QSize(256, 1)));
        m_cmapTex->create();
    }
    ensureUniformCapacity(64);
    m_surfacePipe.reset(makePipeline("surface.vert", "surface.frag", rp, sampleCount, QRhiGraphicsPipeline::Triangles,
                                     true, false, SurfaceInput));
    m_surfaceBiasPipe.reset(makePipeline("surface.vert", "surface.frag", rp, sampleCount,
                                         QRhiGraphicsPipeline::Triangles, true, true, SurfaceInput));
    m_edgePipe.reset(makePipeline("lines.vert", "lines.frag", rp, sampleCount, QRhiGraphicsPipeline::Lines, false,
                                  false, LineInput));
    m_ribbonPipe.reset(makePipeline("ribbon.vert", "surface.frag", rp, sampleCount, QRhiGraphicsPipeline::Triangles,
                                    true, false, RibbonInput));
    m_glyphPipe.reset(makePipeline("glyph.vert", "surface.frag", rp, sampleCount, QRhiGraphicsPipeline::Triangles,
                                   true, false, GlyphInput));
}

void SceneRenderer::ensureUniformCapacity(int draws)
{
    if (m_ubuf && m_ubufDraws >= draws) return;
    const int n = std::max(draws, m_ubufDraws * 2);
    m_ubufStride = m_rhi->ubufAligned(sizeof(DrawUniforms));
    m_ubuf.reset(m_rhi->newBuffer(QRhiBuffer::Dynamic, QRhiBuffer::UniformBuffer, m_ubufStride * quint32(n)));
    m_ubuf->create();
    m_ubufDraws = n;
    if (!m_srb) m_srb.reset(m_rhi->newShaderResourceBindings());
    m_srb->setBindings({QRhiShaderResourceBinding::uniformBufferWithDynamicOffset(
                            0, QRhiShaderResourceBinding::VertexStage | QRhiShaderResourceBinding::FragmentStage,
                            m_ubuf.get(), sizeof(DrawUniforms)),
                        QRhiShaderResourceBinding::sampledTexture(1, QRhiShaderResourceBinding::FragmentStage,
                                                                  m_cmapTex.get(), m_sampler.get())});
    m_srb->create();
}

void SceneRenderer::ensureColormap(QRhiResourceUpdateBatch* u)
{
    const int id = int(m_display.effectiveColormap());
    if (id == m_cmapUploaded) return;
    u->uploadTexture(m_cmapTex.get(), colormap(ColormapId(id)).image());
    m_cmapUploaded = id;
}

SceneRenderer::GpuMesh* SceneRenderer::gpu(const PreparedMesh* mesh, QRhiResourceUpdateBatch* u)
{
    if (auto it = m_gpu.find(mesh); it != m_gpu.end()) return it->second.get();
    auto g = std::make_unique<GpuMesh>();
    const auto& pn = mesh->positionsNormals();
    if (!mesh->data().polylines && !pn.empty() && !mesh->data().indices.empty()) {
        g->vb.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::VertexBuffer, quint32(pn.size() * 4)));
        g->vb->create();
        u->uploadStaticBuffer(g->vb.get(), pn.data());
        const auto& idx = mesh->data().indices;
        g->ib.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::IndexBuffer, quint32(idx.size() * 4)));
        g->ib->create();
        u->uploadStaticBuffer(g->ib.get(), idx.data());
        g->indexCount = quint32(idx.size());
    }
    GpuMesh* raw = g.get();
    m_gpu.emplace(mesh, std::move(g));
    return raw;
}

QRhiBuffer* SceneRenderer::valueBuffer(GpuMesh* g, const PreparedMesh* mesh, const FieldRef& f, QRhiResourceUpdateBatch* u)
{
    const QString key = f.key();
    if (auto it = g->values.find(key); it != g->values.end()) return it->get();
    auto s = mesh->scalar(f);
    if (!s || s->empty()) return nullptr;
    auto buf = std::shared_ptr<QRhiBuffer>(
        m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::VertexBuffer, quint32(s->size() * 4)));
    buf->create();
    u->uploadStaticBuffer(buf.get(), s->data());
    g->values.insert(key, buf);
    return buf.get();
}

QStringList SceneRenderer::visibleItems() const
{
    if (m_display.visibleInitialized) return m_display.visible;
    QStringList out;
    if (m_scene)
        for (const auto& it : m_scene->items)
            if (it.defaultVisible) out << it.id;
    return out;
}

void SceneRenderer::resolveRange()
{
    m_fieldShown = false;
    if (!m_scene || !m_display.contours) return;
    if (m_display.rangeMode == RangeMode::Fixed && m_display.fixedMax > m_display.fixedMin) {
        m_range = {m_display.fixedMin, m_display.fixedMax};
        m_fieldShown = true;
        return;
    }
    auto r = m_scene->range(m_display.field, visibleItems());
    if (!r) return;
    m_range = autoRange(r->first, r->second, m_display.symmetricRange());
    m_fieldShown = true;
}

QMatrix4x4 SceneRenderer::viewProjection(QSize px) const
{
    const float aspect = px.height() > 0 ? float(px.width()) / float(px.height()) : 1.0f;
    const Camera& cam = m_display.camera;
    const float r = m_scene && m_scene->bounds.valid ? m_scene->bounds.radius() : cam.distance;
    const float dc = m_scene && m_scene->bounds.valid ? (cam.eye() - m_scene->bounds.center()).length() : cam.distance;
    const float nearZ = std::max(dc - r * 2.0f, std::max(cam.distance * 0.002f, r * 1e-4f));
    const float farZ = dc + r * 2.0f + cam.distance;
    QMatrix4x4 proj = cam.projection(aspect, nearZ, farZ);
    if (m_rhi) proj = m_rhi->clipSpaceCorrMatrix() * proj;
    return proj * cam.view();
}

QVector3D SceneRenderer::projectToPixel(const QVector3D& world, QSize px) const
{
    const QVector4D c = viewProjection(px) * QVector4D(world, 1);
    if (c.w() <= 0) return {-1, -1, -1};
    const float x = c.x() / c.w(), y = c.y() / c.w();
    const bool yUp = !m_rhi || m_rhi->isYUpInNDC();
    return {(x * 0.5f + 0.5f) * px.width(), (yUp ? (0.5f - y * 0.5f) : (y * 0.5f + 0.5f)) * px.height(), c.z() / c.w()};
}

std::vector<SceneRenderer::DrawCmd> SceneRenderer::buildDraws(QRhiResourceUpdateBatch* u, QSize px, bool forPick)
{
    std::vector<DrawCmd> draws;
    if (!m_scene) return draws;
    if (m_gpuScene != m_scene) {
        // A new dataset: drop GPU copies of meshes it no longer references.
        std::unordered_map<const PreparedMesh*, std::unique_ptr<GpuMesh>> keep;
        for (const auto& it : m_scene->items)
            if (auto f = m_gpu.find(it.mesh.get()); f != m_gpu.end()) keep.emplace(f->first, std::move(f->second));
        m_gpu = std::move(keep);
        m_gpuScene = m_scene;
    }
    const QMatrix4x4 mvp = viewProjection(px);
    const QMatrix4x4 view = m_display.camera.view();
    const QStringList visible = visibleItems();
    const float sceneR = m_scene->bounds.valid ? m_scene->bounds.radius() : 1.0f;
    m_trianglesDrawn = 0;
    if (forPick) m_pickIds.clear();

    auto base = [&](DrawUniforms& du) {
        std::memcpy(du.mvp, mvp.constData(), sizeof du.mvp);
        std::memcpy(du.view, view.constData(), sizeof du.view);
        du.viewport[0] = float(px.width());
        du.viewport[1] = float(px.height());
        du.viewport[2] = m_display.ribbonWidth * (px.width() > 0 ? 1.0f : 1.0f);
        du.params[2] = m_missingGrey;
    };
    auto fieldUniforms = [&](DrawUniforms& du, bool hasField) {
        du.range[0] = float(m_range.min);
        du.range[1] = float(m_range.max);
        du.range[2] = float(m_display.bands);
        du.range[3] = hasField ? 1.0f : 0.0f;
        du.params[1] = 0.85f;
        du.solid[0] = du.solid[1] = du.solid[2] = 0.74f;
        du.solid[3] = 1.0f;
    };

    for (const auto& item : m_scene->items) {
        if (!visible.contains(item.id) || !item.mesh) continue;
        const PreparedMesh* mesh = item.mesh.get();
        GpuMesh* g = gpu(mesh, u);
        if (item.kind == ItemKind::Streamlines) {
            if (forPick) continue;
            // Lines without the shown field are drawn neutral, never coloured by
            // another field against this legend.
            const FieldRef f = m_display.field;
            const bool hasField = m_fieldShown && mesh->hasField(f);
            const QString key = hasField ? f.key() : QStringLiteral("-");
            auto rit = g->ribbons.find(key);
            if (rit == g->ribbons.end()) {
                auto rib = std::make_shared<GpuMesh::Ribbon>();
                const auto& P = mesh->data().positions;
                const auto& off = mesh->data().indices;
                auto s = hasField ? mesh->scalar(f) : nullptr;
                std::vector<float> v;
                std::vector<quint32> idx;
                for (std::size_t l = 0; l + 1 < off.size(); ++l) {
                    for (std::uint32_t k = off[l]; k + 1 < off[l + 1]; ++k) {
                        const float* a = &P[3 * k];
                        const float* b = &P[3 * (k + 1)];
                        const float va = s ? (*s)[k] : 0.f, vb = s ? (*s)[k + 1] : 0.f;
                        const quint32 start = quint32(v.size() / 8);
                        v.insert(v.end(), {a[0], a[1], a[2], b[0], b[1], b[2], 1.f, va});
                        v.insert(v.end(), {a[0], a[1], a[2], b[0], b[1], b[2], -1.f, va});
                        v.insert(v.end(), {b[0], b[1], b[2], a[0], a[1], a[2], -1.f, vb});
                        v.insert(v.end(), {b[0], b[1], b[2], a[0], a[1], a[2], 1.f, vb});
                        idx.insert(idx.end(), {start, start + 1, start + 2, start + 1, start + 3, start + 2});
                    }
                }
                if (!v.empty()) {
                    rib->vb.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::VertexBuffer, quint32(v.size() * 4)));
                    rib->vb->create();
                    u->uploadStaticBuffer(rib->vb.get(), v.data());
                    rib->ib.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::IndexBuffer, quint32(idx.size() * 4)));
                    rib->ib->create();
                    u->uploadStaticBuffer(rib->ib.get(), idx.data());
                    rib->count = quint32(idx.size());
                }
                rit = g->ribbons.insert(key, rib);
            }
            if (!(*rit)->count) continue;
            DrawCmd d{DrawCmd::Ribbon};
            d.pipe = m_ribbonPipe.get();
            d.vb = (*rit)->vb.get();
            d.ib = (*rit)->ib.get();
            d.count = (*rit)->count;
            base(d.u);
            fieldUniforms(d.u, hasField);
            d.u.range[2] = 0;  // smooth colour along lines
            d.u.solid[0] = d.u.solid[1] = d.u.solid[2] = 0.62f;
            d.u.params[0] = 0;
            draws.push_back(d);
            continue;
        }
        if (!g->vb) continue;
        const bool hasField = m_fieldShown && mesh->hasField(m_display.field);
        QRhiBuffer* values = hasField ? valueBuffer(g, mesh, m_display.field, u) : nullptr;
        DrawCmd d{DrawCmd::Surface};
        d.pipe = forPick ? m_pickPipe.get() : (m_display.edges ? m_surfaceBiasPipe.get() : m_surfacePipe.get());
        d.vb = g->vb.get();
        // Solid-coloured draws still need a value stream; the geometry buffer
        // stands in (its floats are read but ignored when mode is solid).
        d.vb2 = values ? values : g->vb.get();
        d.ib = g->ib.get();
        d.count = g->indexCount;
        base(d.u);
        fieldUniforms(d.u, values != nullptr);
        d.u.params[0] = (item.kind == ItemKind::Slice || !m_display.lighting) ? 0.0f : 1.0f;
        if (forPick) {
            d.u.params[3] = float(m_pickIds.size());
            m_pickIds << item.id;
        }
        m_trianglesDrawn += g->indexCount / 3;
        draws.push_back(d);
        if (m_display.edges && !forPick) {
            if (!g->edgeIb) {
                const auto& e = mesh->edges();
                if (!e.empty()) {
                    g->edgeIb.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::IndexBuffer, quint32(e.size() * 4)));
                    g->edgeIb->create();
                    u->uploadStaticBuffer(g->edgeIb.get(), e.data());
                    g->edgeCount = quint32(e.size());
                }
            }
            if (g->edgeIb) {
                DrawCmd ed{DrawCmd::Edges};
                ed.pipe = m_edgePipe.get();
                ed.vb = g->vb.get();
                ed.ib = g->edgeIb.get();
                ed.count = g->edgeCount;
                base(ed.u);
                const bool darkBg = m_clear.lightnessF() < 0.5;
                const float ec = (hasField || !darkBg) ? 0.08f : 0.95f;
                ed.u.solid[0] = ed.u.solid[1] = ed.u.solid[2] = ec;
                ed.u.solid[3] = hasField ? 0.32f : 0.5f;
                draws.push_back(ed);
            }
        }
    }

    // Vector glyphs on the chosen item.
    if (!forPick && m_display.vectors) {
        const SceneItem* item = m_scene->item(m_display.vectorItem);
        if (!item) {
            for (const auto& it : m_scene->items)
                if (visible.contains(it.id) && it.mesh && (it.mesh->hasField({"velocity"}) || it.mesh->hasField({"wall_shear"}))) {
                    item = &it;
                    break;
                }
        }
        if (item && item->mesh) {
            const PreparedMesh* mesh = item->mesh.get();
            FieldRef vf{"velocity", m_display.field.role == FieldRole::Fluent ? FieldRole::Fluent : FieldRole::Kutta, -1};
            if (!mesh->hasField(vf)) vf.quantity = "wall_shear";
            const auto* fd = mesh->data().field(vf.manifestName());
            if (fd && fd->desc.components == 3) {
                GpuMesh* g = gpu(mesh, u);
                const int stride = std::max(1, m_display.vectorStride);
                const QString key = vf.key() + "/" + QString::number(stride);
                auto git = g->glyphs.find(key);
                double maxMag = 0;
                if (git == g->glyphs.end()) {
                    auto gl = std::make_shared<GpuMesh::Glyphs>();
                    std::vector<float> inst;
                    const auto& P = mesh->data().positions;
                    for (std::size_t i = 0; i < mesh->data().vertexCount(); i += std::size_t(stride)) {
                        const float vx = fd->values[3 * i], vy = fd->values[3 * i + 1], vz = fd->values[3 * i + 2];
                        if (!std::isfinite(vx) || !std::isfinite(vy) || !std::isfinite(vz)) continue;
                        const float mag = std::sqrt(vx * vx + vy * vy + vz * vz);
                        inst.insert(inst.end(), {P[3 * i], P[3 * i + 1], P[3 * i + 2], vx, vy, vz, mag});
                    }
                    if (!inst.empty()) {
                        gl->buf.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::VertexBuffer, quint32(inst.size() * 4)));
                        gl->buf->create();
                        u->uploadStaticBuffer(gl->buf.get(), inst.data());
                        gl->count = quint32(inst.size() / 7);
                    }
                    git = g->glyphs.insert(key, gl);
                }
                if (auto r = finiteRange(*fd)) maxMag = r->second;
                if ((*git)->count) {
                    if (!m_arrowUploaded) {
                        std::vector<float> av;
                        std::vector<quint32> ai;
                        arrowMesh(&av, &ai);
                        m_arrowVb.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::VertexBuffer, quint32(av.size() * 4)));
                        m_arrowVb->create();
                        u->uploadStaticBuffer(m_arrowVb.get(), av.data());
                        m_arrowIb.reset(m_rhi->newBuffer(QRhiBuffer::Immutable, QRhiBuffer::IndexBuffer, quint32(ai.size() * 4)));
                        m_arrowIb->create();
                        u->uploadStaticBuffer(m_arrowIb.get(), ai.data());
                        m_arrowIndexCount = quint32(ai.size());
                        m_arrowUploaded = true;
                    }
                    DrawCmd d{DrawCmd::Glyph};
                    d.pipe = m_glyphPipe.get();
                    d.vb = m_arrowVb.get();
                    d.vb2 = (*git)->buf.get();
                    d.ib = m_arrowIb.get();
                    d.count = m_arrowIndexCount;
                    d.instances = (*git)->count;
                    base(d.u);
                    fieldUniforms(d.u, true);
                    // Glyphs are coloured by their own magnitude over [0, max].
                    d.u.range[0] = 0;
                    d.u.range[1] = float(std::max(maxMag, 1e-30));
                    d.u.range[2] = 0;
                    d.u.params[0] = m_display.lighting ? 1.0f : 0.0f;
                    d.u.viewport[3] = m_display.vectorLength > 0 ? m_display.vectorLength : sceneR * 0.035f;
                    draws.push_back(d);
                }
            }
        }
    }
    return draws;
}

void SceneRenderer::ensurePickTarget()
{
    if (m_pickRt) return;
    const auto flags = QRhiTexture::RenderTarget | QRhiTexture::UsedAsTransferSource;
    m_pickPos.reset(m_rhi->newTexture(QRhiTexture::RGBA32F, QSize(1, 1), 1, flags));
    m_pickVal.reset(m_rhi->newTexture(QRhiTexture::RGBA32F, QSize(1, 1), 1, flags));
    m_pickPos->create();
    m_pickVal->create();
    m_pickDepth.reset(m_rhi->newRenderBuffer(QRhiRenderBuffer::DepthStencil, QSize(1, 1)));
    m_pickDepth->create();
    QRhiTextureRenderTargetDescription desc;
    desc.setColorAttachments({QRhiColorAttachment(m_pickPos.get()), QRhiColorAttachment(m_pickVal.get())});
    desc.setDepthStencilBuffer(m_pickDepth.get());
    m_pickRt.reset(m_rhi->newTextureRenderTarget(desc));
    m_pickRp.reset(m_pickRt->newCompatibleRenderPassDescriptor());
    m_pickRt->setRenderPassDescriptor(m_pickRp.get());
    m_pickRt->create();
    m_pickUbufDraws = 64;
    m_pickUbuf.reset(m_rhi->newBuffer(QRhiBuffer::Dynamic, QRhiBuffer::UniformBuffer, m_ubufStride * 64));
    m_pickUbuf->create();
    m_pickSrb.reset(m_rhi->newShaderResourceBindings());
    m_pickSrb->setBindings({QRhiShaderResourceBinding::uniformBufferWithDynamicOffset(
        0, QRhiShaderResourceBinding::VertexStage | QRhiShaderResourceBinding::FragmentStage, m_pickUbuf.get(),
        sizeof(DrawUniforms))});
    m_pickSrb->create();
    m_pickPipe.reset(makePipeline("pick.vert", "pick.frag", m_pickRp.get(), 1, QRhiGraphicsPipeline::Triangles, true,
                                  false, SurfaceInput, 2));
}

void SceneRenderer::requestPick(QPointF px, std::function<void(const PickResult&)> done)
{
    m_pickPixel = px;
    m_pickDone = std::move(done);
    m_pickState = PickState::Requested;
}

void SceneRenderer::render(QRhiCommandBuffer* cb, QRhiRenderTarget* rt, QSize px)
{
    if (!m_rhi || !m_surfacePipe) return;
    QRhiResourceUpdateBatch* u = m_rhi->nextResourceUpdateBatch();
    ensureColormap(u);
    resolveRange();

    auto record = [&](std::vector<DrawCmd>& draws, QRhiBuffer* ubuf, QRhiShaderResourceBindings* srb) {
        for (std::size_t i = 0; i < draws.size(); ++i)
            u->updateDynamicBuffer(ubuf, quint32(i) * m_ubufStride, sizeof(DrawUniforms), &draws[i].u);
        (void)srb;
    };
    auto submit = [&](std::vector<DrawCmd>& draws, QRhiShaderResourceBindings* srb, QSize vp) {
        cb->setViewport({0, 0, float(vp.width()), float(vp.height())});
        for (std::size_t i = 0; i < draws.size(); ++i) {
            const DrawCmd& d = draws[i];
            cb->setGraphicsPipeline(d.pipe);
            const QRhiCommandBuffer::DynamicOffset off(0, quint32(i) * m_ubufStride);
            cb->setShaderResources(srb, 1, &off);
            if (d.kind == DrawCmd::Surface || d.kind == DrawCmd::Glyph) {
                const QRhiCommandBuffer::VertexInput in[] = {{d.vb, 0}, {d.vb2, 0}};
                cb->setVertexInput(0, 2, in, d.ib, 0, QRhiCommandBuffer::IndexUInt32);
            } else {
                const QRhiCommandBuffer::VertexInput in[] = {{d.vb, 0}};
                cb->setVertexInput(0, 1, in, d.ib, 0, QRhiCommandBuffer::IndexUInt32);
            }
            cb->drawIndexed(d.count, d.instances);
        }
    };

    if (m_pickState == PickState::Requested) {
        ensurePickTarget();
        std::vector<DrawCmd> picks = buildDraws(u, px, true);
        // Map the requested pixel to the whole 1x1 target.
        const bool yUp = m_rhi->isYUpInNDC();
        const float cx = float(2.0 * (m_pickPixel.x() + 0.5) / px.width() - 1.0);
        float cy = float(2.0 * (m_pickPixel.y() + 0.5) / px.height() - 1.0);
        if (yUp) cy = -cy;
        QMatrix4x4 pickM;
        pickM.setRow(0, QVector4D(float(px.width()), 0, 0, -cx * float(px.width())));
        pickM.setRow(1, QVector4D(0, float(px.height()), 0, -cy * float(px.height())));
        const QMatrix4x4 full = pickM * viewProjection(px);
        if (int(picks.size()) > m_pickUbufDraws) picks.resize(std::size_t(m_pickUbufDraws));
        for (auto& d : picks) std::memcpy(d.u.mvp, full.constData(), sizeof d.u.mvp);
        record(picks, m_pickUbuf.get(), m_pickSrb.get());
        cb->beginPass(m_pickRt.get(), QColor(0, 0, 0, 0), {1.0f, 0}, u);
        submit(picks, m_pickSrb.get(), QSize(1, 1));
        QRhiResourceUpdateBatch* rb = m_rhi->nextResourceUpdateBatch();
        m_rbPos = std::make_unique<QRhiReadbackResult>();
        m_rbVal = std::make_unique<QRhiReadbackResult>();
        m_rbPending = 2;
        const QStringList ids = m_pickIds;
        auto finish = [this, ids] {
            if (--m_rbPending > 0) return;
            PickResult r;
            float pos[4] = {0, 0, 0, 0}, val[4] = {0, 0, 0, 0};
            if (m_rbPos->data.size() >= 16) std::memcpy(pos, m_rbPos->data.constData(), 16);
            if (m_rbVal->data.size() >= 16) std::memcpy(val, m_rbVal->data.constData(), 16);
            const int obj = int(std::lround(pos[3])) - 1;
            if (obj >= 0 && obj < ids.size()) {
                r.hit = true;
                r.world = QVector3D(pos[0], pos[1], pos[2]);
                r.itemId = ids[obj];
                if (const SceneItem* it = m_scene ? m_scene->item(r.itemId) : nullptr) {
                    r.itemLabel = it->label;
                    if (val[1] > 0.5f && val[0] < 1e37f && it->mesh->hasField(m_display.field) && m_fieldShown)
                        r.value = double(val[0]);
                }
            }
            m_pickState = PickState::Idle;
            if (auto done = std::move(m_pickDone)) done(r);
        };
        m_rbPos->completed = finish;
        m_rbVal->completed = finish;
        rb->readBackTexture(QRhiReadbackDescription(m_pickPos.get()), m_rbPos.get());
        rb->readBackTexture(QRhiReadbackDescription(m_pickVal.get()), m_rbVal.get());
        cb->endPass(rb);
        m_pickState = PickState::InFlight;
        u = m_rhi->nextResourceUpdateBatch();
    }

    std::vector<DrawCmd> draws = buildDraws(u, px, false);
    ensureUniformCapacity(int(draws.size()));
    record(draws, m_ubuf.get(), m_srb.get());
    cb->beginPass(rt, m_clear, {1.0f, 0}, u);
    submit(draws, m_srb.get(), px);
    cb->endPass();
}

} // namespace kapp
