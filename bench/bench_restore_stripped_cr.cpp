// Write a CR-restored copy of a text-mode-damaged legacy Fluent .dat
// (fluent_restore_cr.hpp) and list every restored byte. The output belongs
// outside the corpus, e.g. under build/.
#include "fluent_restore_cr.hpp"
#include "kutta/core/sha256.hpp"
#include <fstream>
#include <iostream>
#include <iterator>

int main(int argc, char** argv) {
  try {
    if (argc != 3) {std::cerr << "Usage: bench_restore_stripped_cr damaged.dat restored.dat\n";return 2;}
    std::ifstream in(argv[1], std::ios::binary);
    if (!in) throw std::runtime_error(std::string("cannot read ") + argv[1]);
    const std::string bytes{std::istreambuf_iterator<char>(in), {}};
    const auto r = kutta::forensics::restore_stripped_cr(bytes);
    std::ofstream out(argv[2], std::ios::binary);
    out.write(r.bytes.data(), std::streamsize(r.bytes.size()));
    if (!out.flush()) throw std::runtime_error(std::string("cannot write ") + argv[2]);
    out.close();
    for (const auto& i : r.inserted)
      std::cout << "RESTORED_CR section=\"" << i.section << "\" value=" << i.value << " offset=" << i.offset << '\n';
    std::cout << "RESTORED bytes_in=" << bytes.size() << " bytes_out=" << r.bytes.size() << " restored=" << r.inserted.size()
              << " input_sha256=" << kutta::sha256_file(argv[1]) << " output_sha256=" << kutta::sha256_file(argv[2]) << '\n';
    return 0;
  } catch (const std::exception& e) {std::cerr << "restore: " << e.what() << '\n';return 1;}
}
