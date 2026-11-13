#pragma once

// STL geometry for New Case from Geometry: read straight into the app so the
// surfaces show before any CLI runs. ASCII files may hold several named solids
// (one wall zone each in `kutta mesh`); a binary file is one solid named after
// the file. Reading is a trust boundary: sizes and counts are checked first.

#include "core/Result.h"
#include "core/ViewFormat.h"

#include <QString>

#include <vector>

namespace kapp {

enum class LengthUnit { Metre, Millimetre };
double unitScale(LengthUnit u);  // to metres
QString unitName(LengthUnit u);  // "m", "mm"
LengthUnit unitFromName(const QString& s);

struct StlSolid {
    QString name;        // the solid's name in the file (binary: the file stem)
    MeshData mesh;       // welded triangles, in metres
    std::size_t triangles = 0;
    double area = 0;     // m^2
};

struct StlFile {
    QString path;
    bool binary = false;
    std::vector<StlSolid> solids;
};

// Reads an STL in the given unit; positions come out in metres with
// coincident vertices welded across smooth surfaces and split at sharp edges
// (over 40°) so the surfaces shade as they are.
Result<StlFile> readStl(const QString& path, LengthUnit unit);

// Scales positions (and areas) in place, e.g. after a unit change.
void rescale(StlFile& f, double factor);
Bounds stlBounds(const StlFile& f);

// Writes named solids as ASCII STL in `unit` (the inverse of readStl's scale).
// Used to hand the mesher renamed solids.
bool writeAsciiStl(const QString& path, const std::vector<std::pair<QString, const MeshData*>>& solids, LengthUnit unit, QString* error);

// Guess the unit from the size: a car-sized part above 20 "metres" is in mm.
LengthUnit guessUnit(const Bounds& rawBounds);

} // namespace kapp
