#pragma once
#include <map>
#include <string>
#include <variant>
#include <vector>
#include <string_view>

namespace kutta {
struct Json {
  using Array=std::vector<Json>;
  using Object=std::map<std::string,Json>;
  using Value=std::variant<std::nullptr_t,bool,double,std::string,Array,Object>;
  Value value=nullptr;
  Json()=default;
  Json(std::nullptr_t):value(nullptr){}
  Json(bool v):value(v){}
  Json(double v):value(v){}
  Json(int v):value(double(v)){}
  Json(std::string v):value(std::move(v)){}
  Json(const char* v):value(std::string(v)){}
  Json(Array v):value(std::move(v)){}
  Json(Object v):value(std::move(v)){}
  bool operator==(const Json&) const=default;
};
// Strict JSON: duplicate keys, nonfinite numbers, invalid escapes and trailing
// text are rejected. Input is bounded by the setup-file boundary (1 MiB).
Json parse_json(std::string_view);
std::string dump_json(const Json&);
}
