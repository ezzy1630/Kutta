#pragma once
// Shared scanning for Fluent's legacy (non-HDF5) case and data files: a
// sequence of parenthesised sections "(id ...)". Ids 2000-2999 and 3000-3999
// carry single- and double-precision binary bodies closed by an explicit
// "End of Binary Section id" trailer.
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

namespace kutta::fluent::legacy {
std::string read_file(const std::filesystem::path &path, std::uintmax_t limit, std::string_view what);
// Index one past the ')' closing the S-expression opened at s[p] == '('.
// Quotes, escapes and ';' comments are honoured; binary bodies are not.
std::size_t balanced_end(std::string_view s, std::size_t p);
// Skip an opaque binary section whose id ended at p. Returns the index one
// past its closing ')'. Used only where the body is not interpreted.
std::size_t skip_binary(std::string_view s, std::size_t p, int section);
// From p: whitespace, then "End of Binary Section <section>)". Returns the
// index one past ')'.
std::size_t binary_trailer(std::string_view s, std::size_t p, int section);
std::size_t skip_space(std::string_view s, std::size_t p);
} // namespace kutta::fluent::legacy
