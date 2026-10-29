#include "kutta/case/sexpr.hpp"
#include <doctest/doctest.h>
#include <random>
using namespace kutta::fluent;
TEST_CASE("Scheme data handles atoms pairs quotation comments and escapes") {
  auto a = parse_sexpr(
      "; comment\n((x . #t) (n -12) (r +1.25e-3) (s \"a\\n\\\"\\\\\\x41;\") (q '(a . b)))");
  CHECK(alist_value(a, "x")->boolean);
  CHECK(alist_value(a, "n")->integer == -12);
  CHECK(*alist_value(a, "r")->number() == doctest::Approx(0.00125));
  CHECK(alist_value(a, "s")->text == "a\n\"\\A");
  CHECK(parse_sexpr(to_string(a)) == a);
  CHECK(to_string(parse_sexpr("#| outer #| inner |# |# #; (ignored) 'x")) == "(quote x)");
  CHECK(*alist_value(parse_sexpr("((flow/scheme 24 integer #f) (flow/scheme 20))"), "flow/scheme")
             ->number() == 20);
}
TEST_CASE("Scheme rejects malformed or excessive input without stack overflow") {
  for (auto s : {"(", ")", "'", "(. a)", "(a .)", "(a . b c)", "(a . b . c)", "\"bad", "#|", "#;",
                 "1e9999", "999999999999999999999999999"})
    CHECK_THROWS_AS(parse_sexpr(s), ParseError);
  CHECK_THROWS_AS(parse_sexpr(std::string(10000, '(')), ParseError);
  CHECK(parse_sexpr(R"scm("C:\Users\xyz")scm").text == R"(C:\Users\xyz)");
  CHECK_THROWS_AS(parse_sexpr(std::string(10000, '\'')), ParseError);
  CHECK_THROWS_AS(parse_sexpr("(a b c)", {128, 3, 100}), ParseError);
  CHECK_THROWS_AS(parse_sexpr("(a)", {128, 100, 2}), ParseError);
  std::mt19937 rng(42);
  std::string chars = "()'.;#tf012e+- abc\\\"\n|";
  for (int n = 0; n < 4000; ++n) {
    std::string s;
    for (unsigned i = 0, len = rng() % 200; i < len; ++i)
      s += chars[rng() % chars.size()];
    try {
      auto a = parse_sexprs(s, {40, 1000, 1000});
      for (auto &x : a)
        CHECK(parse_sexpr(to_string(x)) == x);
    } catch (const ParseError &) {
    }
  }
  std::string large = "(";
  for (int i = 0; i < 300000; ++i)
    large += "(key . 1) ";
  large += ')';
  REQUIRE(large.size() > 2 * 1024 * 1024);
  CHECK(parse_sexpr(large).elements.size() == 300000);
}
