#pragma once
#include <cstdlib>
#include <filesystem>

inline std::filesystem::path fluent_examples() {
  if(const auto* root=std::getenv("KUTTA_FLUENT_EXAMPLES")) return root;
  return KUTTA_SOURCE_DIR "/validation/fluent_examples";
}
