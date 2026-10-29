#pragma once
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
struct MeshTemp {
  std::filesystem::path directory;
  MeshTemp() {
    static std::atomic<unsigned> counter{0};
    directory =
        std::filesystem::temp_directory_path() /
        ("kutta-mesh-" +
         std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count()) +
         "-" + std::to_string(counter++));
    std::filesystem::create_directory(directory);
  }
  ~MeshTemp() {
    std::error_code ec;
    std::filesystem::remove_all(directory, ec);
  }
  std::filesystem::path path(const std::string& name) const {
    return directory / name;
  }
};
inline std::string mesh_text(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), {}};
}
inline void mesh_text(const std::filesystem::path& path,
                      const std::string& text) {
  std::ofstream out(path, std::ios::binary);
  out << text;
}
