#include "live_view.hpp"
#include "kutta/core/path.hpp"
#include <atomic>
#include <chrono>
#include <fstream>
#include <system_error>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <thread>
#include <windows.h>
#endif

namespace kutta::view {
namespace {
void replace(const std::filesystem::path &source,
             const std::filesystem::path &target) {
#ifdef _WIN32
  // std::filesystem::rename cannot replace an existing file on Windows. A
  // polling reader's open handle makes the replace fail with a sharing
  // violation or, without FILE_SHARE_DELETE, access denied; retry until it
  // closes the manifest (about one second at most).
  for (int attempt = 0;; ++attempt) {
    if (MoveFileExW(source.c_str(), target.c_str(), MOVEFILE_REPLACE_EXISTING))
      return;
    const auto error = GetLastError();
    if ((error != ERROR_SHARING_VIOLATION && error != ERROR_ACCESS_DENIED) || attempt == 99)
      throw std::filesystem::filesystem_error(
          "publish live manifest", source, target,
          std::error_code(int(error), std::system_category()));
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }
#else
  std::filesystem::rename(source, target);
#endif
}
} // namespace
LiveView::LiveView(const std::filesystem::path &directory)
    : directory_(std::filesystem::absolute(directory) / "live") {
  std::filesystem::create_directories(directory_);
  static std::atomic<unsigned> next{0};
  do {
    run_id_ = "run-" +
              std::to_string(
                  std::chrono::steady_clock::now().time_since_epoch().count()) +
              "-" + std::to_string(next++);
  } while (!std::filesystem::create_directory(directory_ / run_id_));
}
Json::Object LiveView::publish(const PolyMesh &mesh, const Fields &fields,
                               int iteration, double elapsed) {
  const auto sequence = sequence_ + 1;
  const auto name = "snapshot-" + std::to_string(sequence);
  const auto stage = directory_ / run_id_ / ("." + name + ".tmp");
  const auto target = directory_ / run_id_ / name;
  const auto temporary = directory_ / run_id_ / "manifest.tmp";
  const auto relative = std::filesystem::path(run_id_) / name / "manifest.json";
  Json::Object pointer{{"v", 1},
                       {"run_id", run_id_},
                       {"sequence", double(sequence)},
                       {"iter", iteration},
                       {"elapsed", elapsed},
                       {"manifest", utf8_path(relative)}};
  auto event = pointer;
  event["manifest"] = utf8_path(directory_ / relative);
  bool published = false;
  try {
    write_surface_view(stage, mesh, fields);
    std::filesystem::rename(stage, target);
    std::ofstream file(temporary, std::ios::binary);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << dump_json(pointer) << '\n';
    file.close();
    replace(temporary, directory_ / "manifest.json");
    published = true;
    sequence_ = sequence;
    return event;
  } catch (...) {
    std::error_code ignored;
    std::filesystem::remove_all(stage, ignored);
    std::filesystem::remove(temporary, ignored);
    if (!published)
      std::filesystem::remove_all(target, ignored);
    throw;
  }
}
} // namespace kutta::view
