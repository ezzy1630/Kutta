#include "kutta/app/live.hpp"
#include "kutta/io/binary.hpp"
#include "kutta/io/live_view.hpp"
#include "test_support.hpp"
#include <doctest/doctest.h>
#include <future>
#include <thread>

namespace {
using namespace kutta;
const Json::Object &object(const Json &value) {
  return std::get<Json::Object>(value.value);
}
double number(const Json &value) { return std::get<double>(value.value); }
std::filesystem::path path(const Json &value) {
  const auto &s = std::get<std::string>(value.value);
  return std::filesystem::path(std::u8string(s.begin(), s.end()));
}
} // namespace
TEST_CASE("Live cadence reserves snapshot cost within the solve-time budget") {
  app::LiveCadence cadence(2);
  CHECK_FALSE(cadence.due(1, 1));
  CHECK(cadence.due(2, 1));
  cadence.record(2, .1);
  CHECK_FALSE(cadence.due(3, 10));
  CHECK_FALSE(cadence.due(4, 4));
  CHECK(cadence.due(4, 5));
  cadence.record(4, .1);
  CHECK_FALSE(cadence.due(6, 6));
  CHECK(cadence.due(6, 7));
  CHECK(cadence.count() == 2);
  CHECK(cadence.total_seconds() == doctest::Approx(.2));
  app::LiveCadence near_limit(INT32_MAX);
  near_limit.record(INT32_MAX, .1);
  CHECK_FALSE(near_limit.due(INT32_MAX, 100));
}
TEST_CASE("Live manifest swaps expose complete immutable views with monotonic "
          "sequences") {
  MeshTemp temporary;
  const auto mesh = make_box_hex(2, 1, 1, {0, 0, 0}, {2, 1, 1});
  view::Fields fields;
  fields.reference = true;
  fields.provenance = "atomic fixture";
  const auto nb = mesh.n_faces() - mesh.n_internal_faces();
  fields.boundary = {{"pressure", "Pa", 1, std::vector<float>(nb, 1)},
                     {"pressure@fluent", "Pa", 1, std::vector<float>(nb, 0)}};
  view::LiveView writer(temporary.directory);
  writer.publish(mesh, fields, 1, .1);
  std::atomic<bool> done{false};
  std::atomic<unsigned> reads{0};
  std::promise<void> ready;
  std::string error;
  std::thread reader([&] {
    bool announced = false;
    try {
      double previous = 0;
      while (!done || reads == 0) {
        const auto pointer =
            object(parse_json(mesh_text(writer.directory() / "manifest.json")));
        const auto sequence = number(pointer.at("sequence"));
        KUTTA_REQUIRE(sequence >= previous, "live sequence moved backwards");
        previous = sequence;
        const auto snapshot = writer.directory() / path(pointer.at("manifest"));
        const auto manifest = object(parse_json(mesh_text(snapshot)));
        for (const auto &entry :
             std::get<Json::Array>(manifest.at("zones").value)) {
          const auto &zone = object(entry);
          binary::Reader payload(snapshot.parent_path() /
                                 path(zone.at("file")));
          KUTTA_REQUIRE(payload.scalar<std::uint32_t>() == 0x5654554b &&
                            payload.scalar<std::uint32_t>() == 1 &&
                            payload.scalar<std::uint32_t>() == 32,
                        "incomplete live header");
          const auto nv = payload.scalar<std::uint32_t>(),
                     ni = payload.scalar<std::uint32_t>(),
                     nf = payload.scalar<std::uint32_t>();
          KUTTA_REQUIRE(nv == number(zone.at("vertices")) &&
                            ni == number(zone.at("indices")),
                        "inconsistent live geometry");
          (void)payload.scalar<std::uint32_t>();
          (void)payload.scalar<std::uint32_t>();
          std::uint64_t expected =
              std::uint64_t(nv) * 12 + std::uint64_t(ni) * 4;
          const auto &values = std::get<Json::Array>(zone.at("fields").value);
          KUTTA_REQUIRE(nf == values.size(), "incomplete live fields");
          for (const auto &field : values) {
            const auto &f = object(field);
            expected += std::uint64_t(nv) *
                        std::uint64_t(number(f.at("components"))) * 4;
            if (f.at("name") == Json("pressure"))
              KUTTA_REQUIRE(number(f.at("min")) == sequence &&
                                number(f.at("max")) == sequence,
                            "mixed live sequences");
          }
          KUTTA_REQUIRE(payload.remaining() == expected,
                        "live payload was not fully written");
        }
        ++reads;
        if (!announced) {
          ready.set_value();
          announced = true;
        }
      }
    } catch (const std::exception &e) {
      error = e.what();
      if (!announced)
        ready.set_value();
    }
  });
  ready.get_future().wait();
  try {
    for (int sequence = 2; sequence <= 12; ++sequence) {
      std::fill(fields.boundary[0].values.begin(),
                fields.boundary[0].values.end(), float(sequence));
      writer.publish(mesh, fields, sequence, sequence * .1);
    }
  } catch (...) {
    done = true;
    reader.join();
    throw;
  }
  done = true;
  reader.join();
  CHECK(error.empty());
  CHECK(reads > 0);
  CHECK(writer.sequence() == 12);
  const auto before = mesh_text(writer.directory() / "manifest.json");
  fields.boundary[0].values.pop_back();
  CHECK_THROWS(writer.publish(mesh, fields, 13, 1.3));
  CHECK(mesh_text(writer.directory() / "manifest.json") == before);
  CHECK(writer.sequence() == 12);
  // Fail after all payloads close, before replacing the public pointer.
  fields.boundary[0].values.push_back(12);
  std::filesystem::create_directory(writer.directory() / writer.run_id() /
                                    "manifest.tmp");
  CHECK_THROWS(writer.publish(mesh, fields, 13, 1.3));
  CHECK(mesh_text(writer.directory() / "manifest.json") == before);
  CHECK_FALSE(std::filesystem::exists(writer.directory() / writer.run_id() /
                                      "snapshot-13"));
  CHECK(writer.sequence() == 12);
  for (int sequence = 1; sequence <= 12; ++sequence)
    CHECK(std::filesystem::exists(writer.directory() / writer.run_id() /
                                  ("snapshot-" + std::to_string(sequence)) /
                                  "manifest.json"));
  for (const auto &entry : std::filesystem::directory_iterator(
           writer.directory() / writer.run_id()))
    CHECK(entry.path().extension() != ".tmp");
  view::LiveView next_run(temporary.directory);
  CHECK(next_run.run_id() != writer.run_id());
  CHECK(next_run.sequence() == 0);
  CHECK(mesh_text(writer.directory() / "manifest.json") == before);
}
