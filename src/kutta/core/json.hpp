#pragma once

#include <cmath>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string_view>

namespace kutta {
inline std::string json_string(std::string_view text) {
  std::ostringstream s; s.imbue(std::locale::classic()); s << '"';
  for (unsigned char ch : text) {
    if (ch == '"' || ch == '\\') s << '\\' << ch;
    else if (ch < 32) s << "\\u" << std::hex << std::setfill('0') << std::setw(4) << int(ch) << std::dec;
    else s << ch;
  }
  return s.str() + '"';
}
// JSON has no NaN/Infinity. Unknown optional metadata is null, never zero.
inline void json_number(std::ostream& out, double value) {
  if (std::isfinite(value)) out << value; else out << "null";
}
} // namespace kutta
