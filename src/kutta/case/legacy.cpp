#include "kutta/case/legacy.hpp"
#include <cctype>
#include <charconv>
#include <fstream>
#include <stdexcept>

namespace kutta::fluent::legacy {
std::string read_file(const std::filesystem::path &path, std::uintmax_t limit, std::string_view what) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  if (!f)
    throw std::runtime_error("cannot open " + std::string(what) + ": " + path.string());
  auto length = f.tellg();
  if (length < 0 || static_cast<std::uintmax_t>(length) > limit)
    throw std::runtime_error("legacy " + std::string(what) + " exceeds " +
                             std::to_string(limit >> 20) + " MiB limit");
  std::string bytes(static_cast<std::size_t>(length), '\0');
  f.seekg(0);
  if (!f.read(bytes.data(), static_cast<std::streamsize>(bytes.size())))
    throw std::runtime_error("cannot read " + std::string(what));
  return bytes;
}
std::size_t skip_space(std::string_view s, std::size_t p) {
  while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p])))
    ++p;
  return p;
}
std::size_t balanced_end(std::string_view s, std::size_t p) {
  std::size_t level = 0;
  bool quoted = false, escaped = false, comment = false;
  for (; p < s.size(); ++p) {
    char ch = s[p];
    if (comment) {
      if (ch == '\n')
        comment = false;
      continue;
    }
    if (quoted) {
      if (escaped)
        escaped = false;
      else if (ch == '\\')
        escaped = true;
      else if (ch == '"')
        quoted = false;
      continue;
    }
    if (ch == '"')
      quoted = true;
    else if (ch == ';')
      comment = true;
    else if (ch == '(')
      ++level;
    else if (ch == ')') {
      if (level == 0)
        throw std::runtime_error("unbalanced legacy section");
      if (--level == 0)
        return p + 1;
    }
  }
  throw std::runtime_error("truncated legacy section");
}
std::size_t binary_trailer(std::string_view s, std::size_t p, int section) {
  p = skip_space(s, p);
  const std::string_view marker = "End of Binary Section";
  if (s.substr(p, marker.size()) != marker)
    throw std::runtime_error("missing binary section terminator");
  p = skip_space(s, p + marker.size());
  int closing = 0;
  auto close = std::from_chars(s.data() + p, s.data() + s.size(), closing);
  if (close.ec != std::errc{} || closing != section)
    throw std::runtime_error("binary section terminator ID mismatch");
  p = skip_space(s, static_cast<std::size_t>(close.ptr - s.data()));
  if (p == s.size() || s[p] != ')')
    throw std::runtime_error("malformed binary section terminator");
  return p + 1;
}
std::size_t skip_binary(std::string_view s, std::size_t p, int section) {
  // Binary bodies are opaque here. Fluent's explicit terminator isolates
  // the following sections without interpreting the payload.
  const std::string_view marker = "End of Binary Section";
  auto end = s.find(marker, p);
  if (end == s.npos)
    throw std::runtime_error("missing binary section terminator");
  return binary_trailer(s, end, section);
}
} // namespace kutta::fluent::legacy
