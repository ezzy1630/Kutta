#include "kutta/setup/settings.hpp"
#include "../fluent_examples.hpp"
#include "../io/test_support.hpp"
#include <doctest/doctest.h>
#include <algorithm>
#include <set>

TEST_CASE("Mixing elbow inventories every retained case and zone setting with Fluent groups") {
#if KUTTA_HAS_HDF5
  const auto path=fluent_examples()/"mixing_elbow/mixing_elbow.cas.h5";
  if(!std::filesystem::exists(path)){MESSAGE("SKIP: optional Fluent mixing elbow absent");return;}
  auto c=kutta::fluent::read_case(path);kutta::setup::RunControls run;
  const auto settings=kutta::setup::settings_inventory(c,run);
  std::set<std::string> paths,groups;
  for(const auto& s:settings) {
    CHECK(paths.insert(s.path).second);groups.insert(s.group);
    CHECK(!s.label.empty());CHECK(s.path.find(' ')==std::string::npos);
    CHECK(std::none_of(s.path.begin(),s.path.end(),[](char c){return c>='A' && c<='Z';}));
    CHECK(s.editable==bool(s.apply));if(!s.editable)CHECK(!s.reason.empty());
  }
  CHECK(settings.size()>=c.variables.alist.elements.size());
  for(auto g:{"General","Models","Materials","Cell Zone Conditions","Boundary Conditions","Reference Values","Methods","Controls","Monitors","Initialization","Run Calculation"})CHECK(groups.contains(g));
  for(const auto& z:c.boundaries) {
    // Fixture names are lowercase and path-safe; canonical identities are
    // asserted from the actual case rather than a guessed numeric zone ID.
    CHECK(paths.contains("boundary/"+z.name+"/type"));
    for(const auto& entry:z.raw.alist.elements)if(!entry.elements.empty() && entry.elements[0].atom_text()) {
      const std::string key(*entry.elements[0].atom_text());
      if(key=="vmag")CHECK(paths.contains("boundary/"+z.name+"/velocity-magnitude"));
    }
  }
  MeshTemp t;
  mesh_text(t.path("bad.json"),"{\"v\":2,\"overrides\":{\"missing/path\":1}}");
  CHECK_THROWS_AS(kutta::setup::apply_setup(t.path("bad.json"),settings),kutta::setup::RejectedSettings);
#else
  MESSAGE("SKIP: HDF5 disabled");
#endif
}
