#pragma once

#include <optional>
#include <string_view>

namespace kutta {

// Parses all of `s` as a decimal or exponent-form double (no leading '+'),
// independent of the process locale. Returns nullopt unless every character
// was consumed. Exists because some shipping standard libraries (e.g. Xcode 16
// libc++) still lack floating-point std::from_chars.
std::optional<double> parse_double(std::string_view s);

}  // namespace kutta
