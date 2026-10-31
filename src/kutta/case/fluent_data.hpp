#pragma once
#include "kutta/case/fluent_case.hpp"
#include <span>

namespace kutta::fluent {
enum class FieldLocation { cell, face };
struct FieldSection {
  // Fluent global IDs, 1-based, inclusive; not zone-relative or PolyMesh indices.
  std::uint64_t first = 0, last = 0;
  std::size_t components = 0;
  // Row-major: values[(id-first)*components + component]. No padding for absent IDs.
  std::vector<double> values;
  std::span<const double> at(std::uint64_t fluent_id) const;
};
struct Field {
  std::string name;
  std::vector<FieldSection> sections;
  std::span<const double> at(std::uint64_t fluent_id) const;
  std::size_t size() const;
};
struct ZoneFieldView {
  std::uint64_t first = 0, last = 0;
  std::size_t components = 0;
  std::span<const double> values;
};
struct ResidualHistory {
  std::vector<std::uint64_t> iterations;
  // Four columns retained verbatim. Column 0 raw residual, column 1 scale;
  // columns 2/3 are undocumented, and NaN for legacy .dat files, which store
  // only the first two. Histories can be decimated, not per-iteration.
  std::vector<std::array<double, 4>> rows;
};
struct FluentData {
  std::optional<std::string> case_file;
  std::map<std::string, Field> cell_fields, face_fields;
  std::map<std::string, ResidualHistory> residuals;
  std::optional<std::uint64_t> iteration_count;
  // Distinct from iteration_count: last retained residual sample is not a save counter.
  std::optional<std::uint64_t> last_residual_iteration;
  RawSettings variables;
  std::map<std::string, std::string> raw_settings;
  MeshMetadata mesh;
  std::vector<SettingIssue> issues;
  // Views borrow this object's field storage. Empty means the field has no values
  // on that zone, not zero. Missing zone/field names throw.
  std::vector<ZoneFieldView> zone_field(FieldLocation, std::int64_t zone_id,
                                        std::string_view field) const;
};
// Legacy .dat or CFF .dat.h5, by extension, as read_case.
// Requires metadata from the SAME case: validates each field interval against
// the union of its cell/face zones, including gaps, rather than assuming ID==index.
// Present P/U/V/W fields must cover every fluid-cell zone. Counts and ranges
// cannot prove identity when unrelated meshes share the same topology sizes.
FluentData read_data(const std::filesystem::path &, const FluentCase &matching_case);
} // namespace kutta::fluent
