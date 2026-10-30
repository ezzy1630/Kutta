#pragma once

#include "kutta/core/json_value.hpp"
#include "view.hpp"

namespace kutta::view {
// One solver owns each live directory. Completed snapshots are immutable and
// retained, so a reader may keep using any previously published sequence.
class LiveView {
public:
  explicit LiveView(const std::filesystem::path &directory);
  LiveView(const LiveView &) = delete;
  LiveView &operator=(const LiveView &) = delete;
  const std::filesystem::path &directory() const { return directory_; }
  const std::string &run_id() const { return run_id_; }
  std::uint64_t sequence() const { return sequence_; }
  // Publish closed KUTV files first, then atomically replace the small pointer.
  // The returned manifest path is absolute UTF-8 for a child-process consumer.
  Json::Object publish(const PolyMesh &, const Fields &, int iteration,
                       double elapsed);

private:
  std::filesystem::path directory_;
  std::string run_id_;
  std::uint64_t sequence_ = 0;
};
} // namespace kutta::view
