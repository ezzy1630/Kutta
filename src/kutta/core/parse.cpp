#include "kutta/core/parse.hpp"

#include <charconv>
#include <string>

#if !(defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L)
#include <clocale>
#include <cstdlib>
#if defined(_WIN32)
#include <locale.h>
#else
#include <locale.h>
#if defined(__APPLE__)
#include <xlocale.h>
#endif
#endif
#endif

namespace kutta {

std::optional<double> parse_double(std::string_view s) {
  if (s.empty() || s.front() == '+') return std::nullopt;
#if defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
  double value = 0;
  auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), value);
  if (error != std::errc{} || end != s.data() + s.size()) return std::nullopt;
  return value;
#else
  // strtod honours the global locale's decimal point, so parse with a private
  // "C" locale. Leading whitespace is accepted by strtod but not by us.
  if (s.front() == ' ' || s.front() == '\t' || s.front() == '\n' || s.front() == '\r') return std::nullopt;
  const std::string buffer(s);
  char* end = nullptr;
#if defined(_WIN32)
  static const _locale_t c_locale = _create_locale(LC_ALL, "C");
  const double value = _strtod_l(buffer.c_str(), &end, c_locale);
#else
  static const locale_t c_locale = newlocale(LC_ALL_MASK, "C", locale_t{});
  const double value = strtod_l(buffer.c_str(), &end, c_locale);
#endif
  if (end != buffer.c_str() + buffer.size()) return std::nullopt;
  return value;
#endif
}

}  // namespace kutta
