#pragma once

#include <QColor>
#include <QImage>
#include <QString>

#include <array>
#include <optional>
#include <vector>

namespace kapp {

enum class ColormapId { Viridis, Cividis, Inferno, CoolWarm, BlueRed, FluentRainbow, Grayscale };

struct Colormap {
    ColormapId id;
    QString name;      // shown in pickers
    QString key;       // stable id for settings
    bool diverging;    // centre has meaning (deltas)
    bool perceptual;   // perceptually uniform lightness
    std::array<QRgb, 256> lut;

    QColor sample(double t) const;      // t in [0, 1], clamped
    QImage image(int bands = 0) const;  // 256 x 1; bands > 0 quantizes like the viewport
};

const Colormap& colormap(ColormapId id);
const std::vector<ColormapId>& allColormaps();
std::optional<ColormapId> colormapFromKey(const QString& key);
// Default map for a field: diverging for deltas, perceptual otherwise.
ColormapId defaultColormap(bool delta);

// Colour range for a field display.
enum class RangeMode { Auto, Fixed };
struct ColorRange {
    double min = 0, max = 1;
    bool valid() const { return max > min; }
    double normalize(double v) const { return valid() ? (v - min) / (max - min) : 0.5; }
};
// Auto range from data; deltas are made symmetric about zero so white means agreement.
ColorRange autoRange(double dataMin, double dataMax, bool symmetric);
// Band index a value falls into, matching the shader's quantization.
int bandOf(double t, int bands);

} // namespace kapp
