#pragma once
#include <filesystem>
#include <vector>
#include <fstream>
#include <set>
#include <sstream>
#include "kutta/core/types.hpp"

struct FluentReferenceFile {
  std::filesystem::path path;
  std::uintmax_t bytes;
  std::string sha256;
};
// The ledger covers parse-only research candidates, including known unsupported
// pairs. It does not promote them to the independently checked mesh corpus.
inline std::vector<FluentReferenceFile> fluent_reference_files() {
  std::ifstream in(KUTTA_SOURCE_DIR "/validation/fluent_reference_manifest.tsv");
  KUTTA_REQUIRE(in,"missing Fluent reference ledger");
  std::string line;std::getline(in,line);
  KUTTA_REQUIRE(line.starts_with("path\tbytes\tsha256\t"),"invalid Fluent reference ledger header");
  std::vector<FluentReferenceFile> files;std::set<std::string> seen;
  while(std::getline(in,line)) {
    std::istringstream row(line);std::string path,bytes,sha;
    KUTTA_REQUIRE(std::getline(row,path,'\t') && std::getline(row,bytes,'\t') && std::getline(row,sha,'\t'),"invalid Fluent reference ledger row");
    const std::filesystem::path relative(path);
    KUTTA_REQUIRE(!relative.is_absolute() && relative.lexically_normal()==relative && !path.starts_with("..") && seen.insert(path).second,"invalid or duplicate Fluent reference path");
    KUTTA_REQUIRE(sha.size()==64 && sha.find_first_not_of("0123456789abcdef")==std::string::npos,"invalid Fluent reference hash");
    KUTTA_REQUIRE(!bytes.empty() && bytes.find_first_not_of("0123456789")==std::string::npos,"invalid Fluent reference byte count");
    std::size_t end=0;const auto count=std::stoull(bytes,&end);
    KUTTA_REQUIRE(end==bytes.size(),"invalid Fluent reference byte count");
    files.push_back({relative,count,sha});
  }
  KUTTA_REQUIRE(!files.empty(),"empty Fluent reference ledger");
  return files;
}
