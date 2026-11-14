#include "external_aero.hpp"
#include <cmath>
#include <algorithm>

#include <set>

namespace kutta::mesher {
namespace {
using Object = Json::Object;
[[noreturn]] void bad(const std::string& path, const std::string& why) { throw Error(path + ": " + why); }
const Object& object(const Json& j, const std::string& path, std::initializer_list<std::string_view> keys) {
  const auto* o = std::get_if<Object>(&j.value);
  if (!o) bad(path, "expected an object");
  for (const auto& [k, v] : *o) {
    (void)v;
    if (std::find(keys.begin(), keys.end(), k) == keys.end()) bad(path + "/" + k, "unknown property");
  }
  return *o;
}
double number(const Json& j, const std::string& path, double low, double high = 1e12) {
  const auto* d = std::get_if<double>(&j.value);
  if (!d || !std::isfinite(*d) || *d < low || *d > high) bad(path, "expected a finite number in [" + std::to_string(low) + ", " + std::to_string(high) + "]");
  return *d;
}
int integer(const Json& j, const std::string& path, int low, int high) {
  const auto d = number(j, path, low, high);
  if (std::floor(d) != d) bad(path, "expected an integer");
  return int(d);
}
std::string string(const Json& j, const std::string& path) {
  const auto* s = std::get_if<std::string>(&j.value);
  if (!s) bad(path, "expected a string");
  return *s;
}
bool boolean(const Json& j, const std::string& path) {
  const auto* b = std::get_if<bool>(&j.value);
  if (!b) bad(path, "expected a boolean");
  return *b;
}
Vec3d vector(const Json& j, const std::string& path) {
  const auto* a = std::get_if<Json::Array>(&j.value);
  if (!a || a->size() != 3) bad(path, "expected [x,y,z]");
  return {number((*a)[0], path + "/0", -1e9, 1e9), number((*a)[1], path + "/1", -1e9, 1e9), number((*a)[2], path + "/2", -1e9, 1e9)};
}
}
Settings parse_settings(std::string_view text) {
  Settings s;
  const auto j = parse_json(text);
  const auto& o = object(j, "mesh", {"v", "units", "flow", "up", "tunnel", "base_size", "surface", "refinement_boxes", "layers", "max_cells"});
  if (!o.contains("v")) bad("mesh/v", "required version 1");
  if (integer(o.at("v"), "mesh/v", 1, 1) != 1) bad("mesh/v", "expected 1");
  for (const auto& [k, v] : o) {
    const auto p = "mesh/" + k;
    if (k == "units") { s.units = string(v,p); if (s.units != "m" && s.units != "mm") bad(p,"expected m or mm"); }
    if (k == "flow" && string(v,p) != "+x") bad(p,"v1 supports +x");
    if (k == "up" && string(v,p) != "+z") bad(p,"v1 supports +z");
    if (k == "base_size") s.base_size = number(v,p,1e-7,1e5);
    if (k == "max_cells") s.max_cells = integer(v,p,1,100000000);
    if (k == "tunnel") {
      for (const auto& [a,b] : object(v,p,{"upstream","downstream","side","top","ground","symmetry"})) {
        const auto q=p+"/"+a;
        if (a=="upstream") s.upstream=number(b,q,0.01,100);
        if (a=="downstream") s.downstream=number(b,q,0.01,100);
        if (a=="side") s.side=number(b,q,0.01,100);
        if (a=="top") s.top=number(b,q,0.01,100);
        if (a=="ground") s.ground=boolean(b,q);
        if (a=="symmetry" && !std::holds_alternative<std::nullptr_t>(b.value)) { if(string(b,q)!="y") bad(q,"expected null or y"); s.symmetry_y=true; }
      }
    }
    if (k == "surface") {
      for (const auto& [a,b] : object(v,p,{"min_size","max_size","curvature_degrees","proximity_cells","proximity_min_size"})) {
        const auto q=p+"/"+a;
        if (a=="min_size") s.min_size=number(b,q,1e-7,1e5);
        if (a=="max_size") s.max_size=number(b,q,1e-7,1e5);
        if (a=="curvature_degrees") s.curvature_degrees=number(b,q,0.1,180);
        if (a=="proximity_cells") s.proximity_cells=integer(b,q,1,20);
        if (a=="proximity_min_size") s.proximity_min_size=number(b,q,0,1e5);  // 0: min_size/4
      }
    }
    if (k == "layers") {
      for (const auto& [a,b] : object(v,p,{"count","first_height","growth","zones"})) {
        const auto q=p+"/"+a;
        if(a=="count") s.layer_count=integer(b,q,0,100);
        if(a=="first_height") s.first_height=number(b,q,1e-9,1e3);
        if(a=="growth") s.growth=number(b,q,1.000001,3);
        if(a=="zones") {
          const auto* arr=std::get_if<Json::Array>(&b.value); if(!arr) bad(q,"expected an array of zone names");
          s.layer_zones.clear(); for(std::size_t i=0;i<arr->size();++i) s.layer_zones.push_back(string((*arr)[i],q+"/"+std::to_string(i)));
        }
      }
    }
    if (k == "refinement_boxes") {
      const auto* arr=std::get_if<Json::Array>(&v.value); if(!arr) bad(p,"expected an array");
      s.body_boxes=false;
      std::set<std::string> names;
      for(std::size_t i=0;i<arr->size();++i) {
        const auto q=p+"/"+std::to_string(i);
        const auto& b=object((*arr)[i],q,{"name","min","max","size"});
        for(const auto* key : {"name","min","max","size"}) if(!b.contains(key)) bad(q+"/"+key,"required");
        RefinementBox r{string(b.at("name"),q+"/name"),{vector(b.at("min"),q+"/min"),vector(b.at("max"),q+"/max")},number(b.at("size"),q+"/size",1e-7,1e5)};
        if(r.name.empty() || !names.insert(r.name).second) bad(q+"/name","expected a unique nonempty name");
        if(r.bounds.min.x>=r.bounds.max.x || r.bounds.min.y>=r.bounds.max.y || r.bounds.min.z>=r.bounds.max.z) bad(q,"min must be below max on every axis");
        s.boxes.push_back(r);
      }
    }
  }
  if (s.min_size>s.max_size) bad("mesh/surface/min_size","must not exceed max_size");
  if (s.max_size>s.base_size) bad("mesh/surface/max_size","must not exceed base_size");
  if (s.proximity_min_size>s.min_size) bad("mesh/surface/proximity_min_size","must not exceed min_size");
  return s;
}
}
