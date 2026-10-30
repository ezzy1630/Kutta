#include "core/Colormap.h"

#include <algorithm>
#include <cmath>

namespace kapp {

namespace {

using Rgb = std::array<int, 3>;

// 33 samples per map from matplotlib 3.11 (viridis, cividis, inferno: BSD-style
// licence; coolwarm is Moreland's; RdBu_r from ColorBrewer, Apache-2.0).
constexpr Rgb kViridis[] = {{68,1,84},{71,13,96},{72,24,106},{72,35,116},{71,45,123},{69,55,129},{66,64,134},{62,73,137},{59,82,139},{55,91,141},{51,99,141},{47,107,142},{44,114,142},{41,122,142},{38,130,142},{35,137,142},{33,145,140},{31,152,139},{31,160,136},{34,167,133},{40,174,128},{50,182,122},{63,188,115},{78,195,107},{94,201,98},{112,207,87},{132,212,75},{152,216,62},{173,220,48},{194,223,35},{216,226,25},{236,229,27},{253,231,37}};
constexpr Rgb kCividis[] = {{0,34,78},{0,40,91},{0,46,106},{5,51,113},{26,56,111},{39,62,110},{50,67,109},{59,73,108},{67,78,108},{75,84,108},{83,90,109},{90,95,110},{97,101,111},{104,106,113},{111,112,115},{118,118,118},{125,124,120},{132,130,121},{140,136,120},{147,142,120},{155,148,118},{163,154,116},{171,160,114},{180,167,111},{188,174,108},{196,180,104},{205,187,99},{213,194,94},{222,201,88},{231,209,80},{240,216,70},{249,224,58},{254,232,56}};
constexpr Rgb kInferno[] = {{0,0,4},{4,3,18},{11,7,36},{21,11,55},{33,12,74},{47,10,91},{61,9,101},{74,12,107},{87,16,110},{100,21,110},{113,25,110},{125,30,109},{138,34,106},{151,39,102},{163,44,97},{176,49,91},{188,55,84},{199,62,76},{210,70,68},{219,80,59},{228,90,49},{235,102,40},{241,115,29},{246,128,19},{249,142,9},{251,157,7},{252,172,17},{251,188,33},{249,203,53},{245,219,76},{242,234,105},{243,246,138},{252,255,164}};
constexpr Rgb kCoolWarm[] = {{59,76,192},{68,90,204},{78,104,216},{88,117,225},{98,130,234},{108,143,241},{119,154,247},{130,166,251},{141,176,254},{152,185,255},{163,194,254},{174,201,252},{185,208,249},{195,213,244},{204,217,237},{213,219,229},{221,220,220},{229,216,209},{236,211,197},{241,204,184},{245,196,172},{247,186,159},{247,176,147},{246,165,134},{244,152,122},{240,139,110},{235,125,98},{228,110,86},{221,95,75},{212,78,65},{202,59,55},{190,36,46},{180,4,38}};
constexpr Rgb kBlueRed[] = {{5,48,97},{14,65,121},{23,82,144},{31,99,168},{42,113,178},{52,128,185},{63,142,192},{82,157,200},{107,172,209},{132,188,217},{155,201,224},{174,211,230},{194,221,236},{212,230,241},{224,236,243},{236,242,245},{247,246,246},{249,238,231},{251,229,216},{253,220,201},{251,204,180},{248,187,158},{245,170,137},{238,150,119},{228,128,102},{219,107,85},{208,85,72},{197,62,61},{186,40,50},{171,22,42},{147,14,38},{124,7,34},{103,0,31}};
// Fluent's default contour map: blue, cyan, green, yellow, red at equal steps.
constexpr Rgb kRainbow[] = {{0,0,255},{0,255,255},{0,255,0},{255,255,0},{255,0,0}};
constexpr Rgb kGray[] = {{20,20,20},{235,235,235}};

template <std::size_t N>
std::array<QRgb, 256> build(const Rgb (&pts)[N])
{
    std::array<QRgb, 256> lut{};
    for (int i = 0; i < 256; ++i) {
        const double x = i / 255.0 * double(N - 1);
        const std::size_t k = std::min<std::size_t>(std::size_t(x), N - 2);
        const double f = x - double(k);
        auto mix = [&](int c) { return int(std::lround(pts[k][c] * (1 - f) + pts[k + 1][c] * f)); };
        lut[std::size_t(i)] = qRgb(mix(0), mix(1), mix(2));
    }
    return lut;
}

const std::vector<Colormap>& table()
{
    static const std::vector<Colormap> maps{
        {ColormapId::Viridis, "Viridis", "viridis", false, true, build(kViridis)},
        {ColormapId::Cividis, "Cividis (colour-vision safe)", "cividis", false, true, build(kCividis)},
        {ColormapId::Inferno, "Inferno", "inferno", false, true, build(kInferno)},
        {ColormapId::CoolWarm, "Cool–warm (diverging)", "coolwarm", true, true, build(kCoolWarm)},
        {ColormapId::BlueRed, "Blue–red (diverging)", "bluered", true, false, build(kBlueRed)},
        {ColormapId::FluentRainbow, "Fluent rainbow", "rainbow", false, false, build(kRainbow)},
        {ColormapId::Grayscale, "Grayscale", "gray", false, true, build(kGray)},
    };
    return maps;
}

} // namespace

QColor Colormap::sample(double t) const
{
    if (!std::isfinite(t)) t = 0.5;
    const int i = std::clamp(int(std::lround(std::clamp(t, 0.0, 1.0) * 255.0)), 0, 255);
    return QColor(lut[std::size_t(i)]);
}

QImage Colormap::image(int bands) const
{
    QImage img(256, 1, QImage::Format_RGBA8888);
    for (int i = 0; i < 256; ++i) {
        double t = i / 255.0;
        if (bands > 0) t = (bandOf(t, bands) + 0.5) / bands;
        img.setPixelColor(i, 0, sample(t));
    }
    return img;
}

const Colormap& colormap(ColormapId id)
{
    for (const auto& m : table())
        if (m.id == id) return m;
    return table().front();
}

const std::vector<ColormapId>& allColormaps()
{
    static const std::vector<ColormapId> ids = [] {
        std::vector<ColormapId> v;
        for (const auto& m : table()) v.push_back(m.id);
        return v;
    }();
    return ids;
}

std::optional<ColormapId> colormapFromKey(const QString& key)
{
    for (const auto& m : table())
        if (m.key == key) return m.id;
    return std::nullopt;
}

ColormapId defaultColormap(bool delta) { return delta ? ColormapId::CoolWarm : ColormapId::Viridis; }

ColorRange autoRange(double lo, double hi, bool symmetric)
{
    if (!std::isfinite(lo) || !std::isfinite(hi)) return {0, 1};
    if (symmetric) {
        const double a = std::max(std::abs(lo), std::abs(hi));
        return a > 0 ? ColorRange{-a, a} : ColorRange{-1, 1};
    }
    if (hi <= lo) {
        const double pad = lo == 0 ? 1.0 : std::abs(lo) * 0.01;
        return {lo - pad, hi + pad};
    }
    return {lo, hi};
}

int bandOf(double t, int bands)
{
    if (bands <= 0) return 0;
    return std::clamp(int(std::floor(std::clamp(t, 0.0, 1.0) * bands)), 0, bands - 1);
}

} // namespace kapp
