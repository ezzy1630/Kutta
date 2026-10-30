#pragma once
#include <filesystem>
#include <string>
namespace kutta {
std::string sha256_file(const std::filesystem::path& path);
}
