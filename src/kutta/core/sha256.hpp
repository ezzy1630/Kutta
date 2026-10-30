#pragma once
#include <filesystem>
#include <string>
#include <string_view>
namespace kutta {
std::string sha256_file(const std::filesystem::path& path);
std::string sha256(std::string_view text);
}
