#include "kutta/core/json_value.hpp"
#include <doctest/doctest.h>
TEST_CASE("Setup JSON parser rejects ambiguous or non-JSON values at the boundary") {
  using namespace kutta;
  for(auto text:{"{\"a\":1,\"a\":2}","01","1.","+1","NaN","1e999","[true,]","{\"v\":1}junk","\"\\ud800\"","\"\\udfff\"","\"\\q\"","{unquoted:1}"})CHECK_THROWS(parse_json(text));
  const auto parsed=parse_json("{\"v\":1,\"overrides\":{\"run/iterations\":10},\"unicode\":\"\\ud83d\\ude00\",\"values\":[true,false,null,-1.25e2]}");
  CHECK(parse_json(dump_json(parsed))==parsed);
  CHECK_THROWS(parse_json(std::string(1024*1024+1,' ')));
}
