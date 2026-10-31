#include "kutta/case/fluent_data.hpp"
#include "kutta/case/hdf5.hpp"
#include "kutta/case/legacy.hpp"
#include <bit>
#include <charconv>
#include <sstream>
#include "kutta/core/parse.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace kutta::fluent {
std::span<const double> FieldSection::at(std::uint64_t id) const {
  if (id < first || id > last)
    return {};
  return std::span<const double>(values).subspan(static_cast<std::size_t>(id - first) * components,
                                                 components);
}
std::span<const double> Field::at(std::uint64_t id) const {
  auto it = std::upper_bound(sections.begin(), sections.end(), id,
                             [](auto n, const FieldSection &s) { return n < s.first; });
  return it == sections.begin() ? std::span<const double>{} : std::prev(it)->at(id);
}
std::size_t Field::size() const {
  std::size_t n = 0;
  for (auto &s : sections)
    n += s.values.size() / s.components;
  return n;
}
std::vector<ZoneFieldView> FluentData::zone_field(FieldLocation location, std::int64_t id,
                                                  std::string_view field) const {
  const auto &zones = location == FieldLocation::cell ? mesh.cell_zones : mesh.face_zones;
  auto z = std::find_if(zones.begin(), zones.end(), [id](auto &z) { return z.zone_id == id; });
  if (z == zones.end())
    throw std::out_of_range("unknown Fluent zone ID");
  const auto &fields = location == FieldLocation::cell ? cell_fields : face_fields;
  auto f = fields.find(std::string(field));
  if (f == fields.end())
    throw std::out_of_range("unknown Fluent field");
  std::vector<ZoneFieldView> views;
  for (auto &s : f->second.sections) {
    auto first = std::max(s.first, z->first), last = std::min(s.last, z->last);
    if (first <= last)
      views.push_back({first, last, s.components,
                       std::span<const double>(s.values).subspan(
                           static_cast<std::size_t>(first - s.first) * s.components,
                           static_cast<std::size_t>(last - first + 1) * s.components)});
  }
  return views;
}
namespace {
void contained(std::uint64_t first, std::uint64_t last, const std::vector<ZoneRange> &zones) {
  auto next = first;
  for (auto &z : zones) {
    if (z.last < next)
      continue;
    if (z.first > next)
      break;
    if (z.last >= last)
      return;
    next = z.last + 1;
  }
  throw std::runtime_error("data IDs lie outside matching case topology");
}
void check_components(const std::string &name, std::size_t components) {
  if (name == "SV_WALL_SHEAR" && components != 3)
    throw std::runtime_error("wall shear must have three row-major components");
  for (auto scalar : {"SV_P", "SV_U", "SV_V", "SV_W", "SV_K", "SV_O", "SV_D", "SV_MU_T",
                      "SV_DENSITY", "SV_WALL_YPLUS", "SV_FLUX"})
    if (name == scalar && components != 1)
      throw std::runtime_error("scalar field has multiple components: " + name);
}
void complete_primary_fields(const FluentData &out, const FluentCase &c) {
  for (auto field_name : {"SV_P", "SV_U", "SV_V", "SV_W"}) {
    auto field = out.cell_fields.find(field_name);
    if (field == out.cell_fields.end())
      continue;
    for (const auto &range : c.mesh.cell_zones) {
      if (!c.cell_zones.empty()) {
        auto zone = std::find_if(c.cell_zones.begin(), c.cell_zones.end(),
                                 [&](auto &z) { return z.id == range.zone_id; });
        if (zone == c.cell_zones.end())
          throw std::runtime_error("case topology has no matching cell-zone settings");
        if (zone->type != ZoneType::fluid)
          continue;
      }
      auto next = range.first;
      bool covered = false;
      for (const auto &chunk : field->second.sections) {
        if (chunk.last < next)
          continue;
        if (chunk.first > next)
          break;
        if (chunk.last >= range.last) {
          covered = true;
          break;
        }
        next = chunk.last + 1;
      }
      if (!covered)
        throw std::runtime_error("primary field " + std::string(field_name) +
                                 " does not cover expected fluid cell IDs in zone " +
                                 std::to_string(range.zone_id));
    }
  }
}
} // namespace
#if KUTTA_HAS_HDF5
namespace {
std::uint64_t required_id(hid_t f, const std::string &path, const char *name) {
  auto n = h5::attribute(f, path, name);
  if (!n || !*n)
    throw std::runtime_error("missing/zero ID attribute: " + path + "/" + name);
  return *n;
}
std::map<std::string, Field> fields(hid_t f, const std::string &base,
                                    const std::vector<ZoneRange> &zones) {
  std::map<std::string, Field> out;
  if (!h5::exists(f, base))
    return out;
  for (auto &name : h5::children(f, base)) {
    if (name == "fields")
      continue;
    const auto group = base + "/" + name;
    Field field;
    field.name = name;
    for (auto &part : h5::children(f, group)) {
      auto path = group + "/" + part;
      FieldSection s;
      s.first = required_id(f, path, "minId");
      s.last = required_id(f, path, "maxId");
      if (s.last < s.first)
        throw std::runtime_error("reversed field ID interval");
      contained(s.first, s.last, zones);
      std::vector<hsize_t> dims;
      s.values = h5::numbers<double>(f, path, &dims);
      if (dims.empty() || dims[0] != s.last - s.first + 1 || (dims.size() == 2 && dims[1] == 0))
        throw std::runtime_error("field shape disagrees with ID interval: " + path);
      s.components = dims.size() == 2 ? static_cast<std::size_t>(dims[1]) : 1;
      check_components(name, s.components);
      if (!field.sections.empty() && s.components != field.sections.front().components)
        throw std::runtime_error("inconsistent field component count");
      field.sections.push_back(std::move(s));
    }
    std::sort(field.sections.begin(), field.sections.end(),
              [](auto &a, auto &b) { return a.first < b.first; });
    for (std::size_t i = 1; i < field.sections.size(); ++i)
      if (field.sections[i - 1].last >= field.sections[i].first)
        throw std::runtime_error("overlapping field chunks: " + name);
    if (field.sections.empty())
      throw std::runtime_error("field group has no chunks: " + name);
    out.emplace(name, std::move(field));
  }
  return out;
}
void read_residuals(hid_t f, FluentData &out) {
  const std::string root = "/results/residuals/phase-1";
  if (!h5::exists(f, root)) {
    out.issues.push_back({"residual-history", "absent"});
    return;
  }
  for (auto &eq : h5::children(f, root)) {
    auto path = root + "/" + eq;
    ResidualHistory r;
    std::vector<hsize_t> shape;
    auto iterations = h5::numbers<double>(f, path + "/iterations", &shape);
    if (shape.size() != 1)
      throw std::runtime_error("residual iteration rank must be one");
    auto values = h5::numbers<double>(f, path + "/data", &shape);
    if (shape.size() != 2 || shape[0] != iterations.size() || shape[1] != 4)
      throw std::runtime_error("residual shape mismatch");
    for (std::size_t i = 0; i < iterations.size(); ++i) {
      auto n = iterations[i];
      if (!std::isfinite(n) || n < 0 ||
          n >= static_cast<double>(std::numeric_limits<std::uint64_t>::max()) || std::floor(n) != n)
        throw std::runtime_error("invalid residual iteration");
      auto it = static_cast<std::uint64_t>(n);
      if (!r.iterations.empty() && it <= r.iterations.back())
        throw std::runtime_error("non-increasing residual iterations");
      r.iterations.push_back(it);
      r.rows.push_back({values[4 * i], values[4 * i + 1], values[4 * i + 2], values[4 * i + 3]});
      out.last_residual_iteration = std::max(out.last_residual_iteration.value_or(0), it);
    }
    out.residuals.emplace(eq, std::move(r));
  }
}
} // namespace
#endif
namespace {
FluentData hdf5_data(const std::filesystem::path &path, const FluentCase &c) {
#if KUTTA_HAS_HDF5
  auto f = h5::open(path);
  FluentData out;
  out.mesh = c.mesh;
  for (auto &key : h5::children(f, "/settings"))
    out.raw_settings.emplace(key, h5::string(f, "/settings/" + key));
  if (auto saved = out.raw_settings.find("Case File"); saved != out.raw_settings.end()) {
    out.case_file = saved->second;
    auto basename = [](std::string name) {
      std::replace(name.begin(), name.end(), '\\', '/');
      return std::filesystem::path(name).filename().string();
    };
    if (c.source_path && basename(*out.case_file) != c.source_path->filename().string())
      out.issues.push_back({"data/Case File", "saved case name " + *out.case_file +
          " differs from loaded case " + c.source_path->filename().string(), IssueClass::informational});
  } else
    out.issues.push_back({"data/Case File", "absent: saved case name unavailable", IssueClass::informational});
  if (auto it = out.raw_settings.find("Data Variables"); it != out.raw_settings.end()) {
    auto s = parse_sexpr(it->second);
    if (s.elements.size() != 2 || s.elements[0].kind != SExpr::Kind::integer ||
        s.elements[0].integer != 37)
      throw std::runtime_error("malformed Data Variables");
    out.variables.alist = std::move(s.elements[1]);
    if (auto n = out.variables.value("number-of-iterations");
        n && n->kind == SExpr::Kind::integer && n->integer >= 0) {
      out.iteration_count = static_cast<std::uint64_t>(n->integer);
      out.variables.interpreted.insert("number-of-iterations");
    }
  }
  if (!out.iteration_count)
    out.issues.push_back(
        {"iteration-count", "absent or unparsed; last residual iteration exposed separately"});
  const std::string phase = "/results/1/phase-1";
  if (!h5::exists(f, phase))
    throw std::runtime_error("missing results/1/phase-1");
  for (auto &result : h5::children(f, "/results"))
    if (result != "1" && result != "residuals")
      out.issues.push_back({"results/" + result, "unparsed additional result set"});
  for (auto &p : h5::children(f, "/results/1"))
    if (p != "phase-1")
      out.issues.push_back({"results/1/" + p, "unparsed additional phase"});
  out.cell_fields = fields(f, phase + "/cells", c.mesh.cell_zones);
  out.face_fields = fields(f, phase + "/faces", c.mesh.face_zones);
  complete_primary_fields(out, c);
  read_residuals(f, out);
  return out;
#else
  (void)path;
  (void)c;
  throw std::runtime_error("HDF5 support is not built");
#endif
}
// Legacy data files (Fluent .dat): the same sections and byte-order flag as
// legacy case files. Field sections 300/2300/3300 have the decimal header
// (field-id zone-id size n-time-levels n-phases first-id last-id) and
// size values per ID, row major (SV_WALL_SHEAR is x,y,z per face; z is zero
// in 2D). The field name exists only in the preceding (0 "SV_NAME, ...")
// comment, so a field without one is kept under its numeric ID. Residual
// sections 302/2302/3302, (rows equation 1 1), store rows of iteration, raw
// residual and scale.
class LegacyData {
public:
  LegacyData(std::string bytes, const FluentCase &c) : bytes_(std::move(bytes)), s_(bytes_), case_(c) {}
  FluentData read() {
    out_.mesh = case_.mesh;
    std::size_t p = 0;
    for (;;) {
      p = legacy::skip_space(s_, p);
      if (p == s_.size())
        break;
      if (s_[p] != '(')
        throw std::runtime_error("expected legacy data section at byte " + std::to_string(p));
      const auto begin = p++;
      p = legacy::skip_space(s_, p);
      int section = 0;
      const auto conv = std::from_chars(s_.data() + p, s_.data() + s_.size(), section);
      if (conv.ec != std::errc{} || conv.ptr == s_.data() + p)
        throw std::runtime_error("invalid legacy data section number");
      p = static_cast<std::size_t>(conv.ptr - s_.data());
      const int code = section % 1000;
      const bool binary = section >= 2000 && section < 4000;
      if (code == 300 || code == 302) {
        p = code == 300 ? field(p, section) : residual(p, section);
        comment_.clear();
        continue;
      }
      if (binary) {
        // Data bodies close after the trailer (see values()).
        p = legacy::skip_space(s_, legacy::skip_binary(s_, p, section));
        if (p == s_.size() || s_[p++] != ')')
          throw std::runtime_error("unterminated legacy data section " + std::to_string(section));
        out_.issues.push_back({"data/section/" + std::to_string(section), "unparsed binary data section",
                               IssueClass::informational});
        comment_.clear();
        continue;
      }
      p = legacy::balanced_end(s_, begin);
      const auto e = parse_sexpr(s_.substr(begin, p - begin));
      if (section == 0) {
        if (e.elements.size() > 1 && e.elements[1].kind == SExpr::Kind::string)
          comment_ = e.elements[1].text;
        continue;
      }
      comment_.clear();
      if (section == 4) {
        if (e.elements.size() != 2 || e.elements[1].kind != SExpr::Kind::list ||
            e.elements[1].elements.empty() || e.elements[1].elements[0].kind != SExpr::Kind::integer)
          throw std::runtime_error("malformed legacy machine configuration");
        little_ = e.elements[1].elements[0].integer == 60;
      } else if (section == 33) {
        grid_size(e);
      } else if (section == 37) {
        if (e.elements.size() != 2)
          throw std::runtime_error("malformed legacy data variables");
        out_.variables.alist = e.elements[1];
        if (auto n = out_.variables.value("number-of-iterations");
            n && n->kind == SExpr::Kind::integer && n->integer >= 0) {
          out_.iteration_count = static_cast<std::uint64_t>(n->integer);
          out_.variables.interpreted.insert("number-of-iterations");
        }
      } else
        out_.issues.push_back({"data/section/" + std::to_string(section), "unparsed data section",
                               IssueClass::informational});
    }
    if (!grid_seen_)
      throw std::runtime_error("legacy data lacks its grid size section 33");
    for (auto *fields : {&out_.cell_fields, &out_.face_fields})
      for (auto &[name, f] : *fields) {
        std::sort(f.sections.begin(), f.sections.end(), [](auto &a, auto &b) { return a.first < b.first; });
        for (std::size_t i = 1; i < f.sections.size(); ++i)
          if (f.sections[i - 1].last >= f.sections[i].first)
            throw std::runtime_error("overlapping field chunks: " + name);
      }
    complete_primary_fields(out_, case_);
    name_residuals();
    out_.issues.push_back({"data/Case File", "absent: legacy data stores no saved case name", IssueClass::informational});
    if (!out_.iteration_count)
      out_.issues.push_back(
          {"iteration-count", "absent or unparsed; last residual iteration exposed separately"});
    return std::move(out_);
  }

private:
  std::string bytes_;
  std::string_view s_;
  const FluentCase &case_;
  FluentData out_;
  std::string comment_;
  bool little_ = true, grid_seen_ = false;
  std::map<int, std::string> names_;
  std::vector<std::pair<int, ResidualHistory>> residuals_;

  // Decimal integers between parentheses, e.g. "(1 11 1 0 1 1 25740)".
  std::vector<std::uint64_t> header(std::size_t &p) {
    p = legacy::skip_space(s_, p);
    if (p == s_.size() || s_[p] != '(')
      throw std::runtime_error("missing legacy data section header");
    const auto end = legacy::balanced_end(s_, p);
    std::vector<std::uint64_t> values;
    for (auto q = p + 1; (q = legacy::skip_space(s_, q)) < end - 1;) {
      std::uint64_t v = 0;
      const auto r = std::from_chars(s_.data() + q, s_.data() + end - 1, v);
      if (r.ec != std::errc{} || r.ptr == s_.data() + q)
        throw std::runtime_error("invalid legacy data section header");
      values.push_back(v);
      q = static_cast<std::size_t>(r.ptr - s_.data());
    }
    p = end;
    return values;
  }
  void grid_size(const SExpr &e) {
    if (e.elements.size() != 2 || e.elements[1].kind != SExpr::Kind::list || e.elements[1].elements.size() != 3)
      throw std::runtime_error("malformed legacy grid size");
    // Periodic-shadow faces belong to the case's face count but not, in the
    // corpus's periodic VKI pair, to the data's. There the shadows carry
    // type 1008 (0x3f0), no thread record, and their partners' face counts.
    std::uint64_t shadow = 0;
    for (const auto &z : case_.mesh.face_zones)
      if (z.type_code == 8 || z.type_code == 1008)
        shadow += z.size();
    const std::optional<std::uint64_t> expected[] = {case_.mesh.cells, case_.mesh.faces, case_.mesh.nodes};
    for (int k = 0; k < 3; ++k) {
      const auto &v = e.elements[1].elements[k];
      const auto n = static_cast<std::uint64_t>(v.integer);
      if (v.kind != SExpr::Kind::integer || v.integer < 0 ||
          (expected[k] && n != *expected[k] && !(k == 1 && n + shadow == *expected[k])))
        throw std::runtime_error("legacy data grid size disagrees with the matching case");
    }
    grid_seen_ = true;
  }
  // A body shorter than its header promises. Text-mode transfers (FTP ASCII
  // mode, git text normalisation) delete the 0x0D of every 0x0D 0x0A pair;
  // random binary data contains such pairs, so a file with none is suspect.
  std::string short_body(int section) const {
    std::string message = "binary data body length disagrees with its header in section " + std::to_string(section);
    if (s_.find("\r\n") == std::string_view::npos && s_.find('\r') != std::string_view::npos)
      message += "; the file has no CR LF byte pair, consistent with a text-mode transfer that stripped them";
    return message;
  }
  // Body values in the section's encoding, then the section's own closing.
  std::vector<double> values(std::size_t &p, int section, std::uint64_t count) {
    p = legacy::skip_space(s_, p);
    if (p == s_.size() || s_[p] != '(')
      throw std::runtime_error("missing legacy data section body");
    ++p;
    std::vector<double> out;
    if (section < 2000) {
      const auto end = legacy::balanced_end(s_, p - 1);
      std::istringstream text{std::string(s_.substr(p, end - 1 - p))};
      for (std::string token; text >> token;) {
        const auto v = parse_double(std::string_view(token).substr(token.starts_with('+') ? 1 : 0));
        if (!v)
          throw std::runtime_error("invalid ASCII data value '" + token + "'");
        out.push_back(*v);
      }
      p = legacy::skip_space(s_, end);
      if (p == s_.size() || s_[p++] != ')')
        throw std::runtime_error("unterminated legacy data section " + std::to_string(section));
    } else {
      const std::size_t width = section < 3000 ? 4 : 8;
      if (count > (s_.size() - p) / width)
        throw std::runtime_error("truncated binary data section " + std::to_string(section));
      out.resize(static_cast<std::size_t>(count));
      for (auto &v : out) {
        std::uint64_t bits = 0;
        for (std::size_t i = 0; i < width; ++i)
          bits |= static_cast<std::uint64_t>(static_cast<unsigned char>(s_[p + i]))
                  << (8 * (little_ ? i : width - 1 - i));
        v = width == 8 ? std::bit_cast<double>(bits)
                       : static_cast<double>(std::bit_cast<float>(static_cast<std::uint32_t>(bits)));
        p += width;
      }
      // Field bodies close after the trailer, "...\nEnd of Binary Section
      // 3300))"; residual bodies before it, "...)\n\nEnd of Binary Section
      // 3302)". Either way exactly one ')' closes the body. The separator is
      // matched exactly: skipping whitespace would accept a body one byte
      // short whose last read byte was the separator's LF.
      const bool closed = s_.substr(p, 3) == ")\n\n";
      const std::string_view separator = closed ? ")\n\nEnd of Binary Section" : "\nEnd of Binary Section";
      if (s_.substr(p, separator.size()) != separator)
        throw std::runtime_error(short_body(section));
      p += closed ? 3 : 1;
      p = legacy::binary_trailer(s_, p, section);
      if (!closed && (p == s_.size() || s_[p++] != ')'))
        throw std::runtime_error("unterminated legacy data section " + std::to_string(section));
    }
    if (out.size() != count)
      throw std::runtime_error("data body length disagrees with its header in section " + std::to_string(section));
    return out;
  }
  std::size_t field(std::size_t p, int section) {
    const auto h = header(p);
    if (h.size() != 7 || !h[2] || h[2] > 64 || !h[5] || h[6] < h[5] || h[1] > std::uint64_t(std::numeric_limits<int>::max()))
      throw std::runtime_error("invalid legacy field header");
    const auto id = static_cast<int>(h[0]);
    const auto zone = static_cast<std::int64_t>(h[1]);
    FieldSection chunk;
    chunk.first = h[5];
    chunk.last = h[6];
    chunk.components = static_cast<std::size_t>(h[2]);
    chunk.values = values(p, section, (chunk.last - chunk.first + 1) * chunk.components);
    std::string name = comment_.substr(0, comment_.find(','));
    if (!name.starts_with("SV_")) {
      name = "SV_ID_" + std::to_string(id);
      out_.issues.push_back({"data/field/" + name, "unnamed legacy field kept under its numeric ID",
                             IssueClass::informational});
    }
    if (auto [it, inserted] = names_.emplace(id, name); !inserted && it->second != name)
      throw std::runtime_error("legacy field ID " + std::to_string(id) + " names both " + it->second + " and " + name);
    // Time levels, phases other than the mixture and other domains are
    // not interpreted; dropping them must be visible.
    const auto domain = comment_.find(", domain ");
    const bool mixture = domain == std::string::npos || comment_.compare(domain, 11, ", domain 1,") == 0;
    if (h[3] != 0 || h[4] != 1 || !mixture) {
      out_.issues.push_back({"data/field/" + name, "unparsed time level or phase", IssueClass::informational});
      return p;
    }
    check_components(name, chunk.components);
    const auto in = [&](const std::vector<ZoneRange> &zones) {
      return std::find_if(zones.begin(), zones.end(), [&](auto &z) { return z.zone_id == zone; });
    };
    const auto cell = in(case_.mesh.cell_zones);
    const auto face = in(case_.mesh.face_zones);
    const bool is_cell = cell != case_.mesh.cell_zones.end();
    if (!is_cell && face == case_.mesh.face_zones.end())
      throw std::runtime_error("legacy field " + name + " names unknown zone " + std::to_string(zone));
    const auto &range = is_cell ? *cell : *face;
    if (chunk.first < range.first || chunk.last > range.last)
      throw std::runtime_error("data IDs lie outside matching case topology");
    auto &f = (is_cell ? out_.cell_fields : out_.face_fields)[name];
    f.name = name;
    if (!f.sections.empty() && f.sections.front().components != chunk.components)
      throw std::runtime_error("inconsistent field component count");
    for (const double v : chunk.values)
      if (!std::isfinite(v))
        throw std::runtime_error("nonfinite value in legacy field " + name);
    f.sections.push_back(std::move(chunk));
    return p;
  }
  std::size_t residual(std::size_t p, int section) {
    const auto h = header(p);
    if (h.size() != 4 || h[1] > std::uint64_t(std::numeric_limits<int>::max()) || h[0] > (1u << 26))
      throw std::runtime_error("invalid legacy residual header");
    const auto v = values(p, section, 3 * h[0]);
    ResidualHistory r;
    for (std::size_t i = 0; i < h[0]; ++i) {
      const double n = v[3 * i];
      if (!std::isfinite(n) || n < 0 || n >= 9.0e15 || std::floor(n) != n)
        throw std::runtime_error("invalid residual iteration");
      const auto it = static_cast<std::uint64_t>(n);
      if (!r.iterations.empty() && it <= r.iterations.back())
        throw std::runtime_error("non-increasing residual iterations");
      r.iterations.push_back(it);
      // Columns 2 and 3 exist only in CFF residual histories.
      const double absent = std::numeric_limits<double>::quiet_NaN();
      r.rows.push_back({v[3 * i + 1], v[3 * i + 2], absent, absent});
      out_.last_residual_iteration = std::max(out_.last_residual_iteration.value_or(0), it);
    }
    residuals_.emplace_back(static_cast<int>(h[1]), std::move(r));
    return p;
  }
  // Residual sections carry only an equation index. In every legacy pair of
  // the corpus (Fluent 17.2, 18.0 and 21.2) they appear in the order of the
  // case's residuals/settings names, which is how they are named here.
  void name_residuals() {
    std::vector<std::string> names;
    if (const auto *list = case_.variables.value("residuals/settings"); list && list->kind == SExpr::Kind::list)
      for (const auto &entry : list->elements)
        if (entry.kind == SExpr::Kind::list && !entry.elements.empty())
          if (auto t = entry.elements[0].atom_text())
            names.emplace_back(*t);
    const bool ordered = names.size() == residuals_.size();
    if (!residuals_.empty())
      out_.issues.push_back({"residual-history", ordered
          ? "legacy residual equations named in the order of the case's residuals/settings"
          : "legacy residual equations unnamed: their count differs from the case's residuals/settings",
          IssueClass::informational});
    for (std::size_t i = 0; i < residuals_.size(); ++i) {
      const auto name = ordered ? names[i] : "equation-" + std::to_string(residuals_[i].first);
      if (!out_.residuals.emplace(name, std::move(residuals_[i].second)).second)
        throw std::runtime_error("duplicate legacy residual equation " + name);
    }
    if (residuals_.empty())
      out_.issues.push_back({"residual-history", "absent"});
  }
};
FluentData legacy_data(const std::filesystem::path &path, const FluentCase &c) {
  return LegacyData(legacy::read_file(path, std::uintmax_t(2) << 30, "data"), c).read();
}
} // namespace
FluentData read_data(const std::filesystem::path &path, const FluentCase &c) {
  if (c.mesh.cell_zones.empty() || c.mesh.face_zones.empty())
    throw std::runtime_error("data reader needs matching case topology metadata");
  return path.extension() == ".h5" ? hdf5_data(path, c) : legacy_data(path, c);
}
} // namespace kutta::fluent
