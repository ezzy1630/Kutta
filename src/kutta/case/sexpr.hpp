#pragma once
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace kutta::fluent {
struct SExpr {
  enum class Kind { symbol, string, integer, real, boolean, list };
  Kind kind = Kind::list;
  std::string text;
  std::int64_t integer = 0;
  double real = 0;
  bool boolean = false;
  std::vector<SExpr> elements;
  std::shared_ptr<const SExpr> tail;
  static SExpr symbol(std::string text);
  std::optional<double> number() const;
  std::optional<std::string_view> atom_text() const;
  bool operator==(const SExpr &) const;
};
struct ParseLimits {
  std::size_t max_depth = 128;
  // Zero derives a budget from input bytes; explicit caller budgets remain hard limits.
  std::size_t max_nodes = 0;
  std::size_t max_bytes = 128 * 1024 * 1024;
};
struct ParseError : std::runtime_error {
  std::size_t offset;
  ParseError(std::size_t offset, std::string message);
};
// Data only: quoted forms and Scheme procedures are never evaluated.
std::vector<SExpr> parse_sexprs(std::string_view text, ParseLimits limits = {});
SExpr parse_sexpr(std::string_view text, ParseLimits limits = {});
std::string to_string(const SExpr &);
// Prefer an exact (key value) or (key . value) over rpvar registration metadata.
const SExpr *alist_entry(const SExpr &alist, std::string_view key);
const SExpr *entry_value(const SExpr &entry);
const SExpr *alist_value(const SExpr &alist, std::string_view key);
} // namespace kutta::fluent
