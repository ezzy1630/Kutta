#include "kutta/core/path.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <numeric>

#include "../backend_helpers.hpp"
#include "support.hpp"
#include "../fluent_examples.hpp"

#if KUTTA_HAS_HDF5
#include <hdf5.h>

namespace {
using namespace sst_test;

// Private fixture reader, not a general Fluent solution-data API. Read numeric
// sections by minId/maxId and vectors as (n,components), never component-major.
struct H5 {
  hid_t id;
  herr_t (*close)(hid_t);
  H5(hid_t value, herr_t (*closer)(hid_t)) : id(value), close(closer) {
    KUTTA_REQUIRE(id >= 0,"SST fixture: HDF5 open failed");
  }
  ~H5() { close(id); }
  H5(const H5&) = delete;
  H5& operator=(const H5&) = delete;
  operator hid_t() const { return id; }
};
std::int64_t attribute(hid_t dataset, const char* name) {
  H5 a(H5Aopen(dataset,name,H5P_DEFAULT),H5Aclose);
  H5 space(H5Aget_space(a),H5Sclose);
  KUTTA_REQUIRE(H5Sget_simple_extent_npoints(space) == 1,"SST fixture: non-scalar id attribute");
  std::int64_t value = 0;
  KUTTA_REQUIRE(H5Aread(a,H5T_NATIVE_INT64,&value) >= 0,"SST fixture: attribute read failed");
  return value;
}
std::vector<double> numeric(hid_t dataset) {
  H5 space(H5Dget_space(dataset),H5Sclose);
  auto n = H5Sget_simple_extent_npoints(space);
  KUTTA_REQUIRE(n >= 0 && n <= 10000000,"SST fixture: unexpected dataset size");
  std::vector<double> out(static_cast<std::size_t>(n));
  KUTTA_REQUIRE(H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,out.data()) >= 0,
                "SST fixture: numeric read failed");
  return out;
}
std::vector<double> numeric(hid_t file, const std::string& path) {
  H5 dataset(H5Dopen2(file,path.c_str(),H5P_DEFAULT),H5Dclose);
  return numeric(dataset);
}
std::vector<std::string> children(hid_t group) {
  std::vector<std::string> names;
  auto visit = [](hid_t, const char* name, const H5L_info_t*, void* ptr) -> herr_t {
    static_cast<std::vector<std::string>*>(ptr)->push_back(name); return 0;
  };
  hsize_t position = 0;
  KUTTA_REQUIRE(H5Literate(group,H5_INDEX_NAME,H5_ITER_INC,&position,visit,&names) >= 0,
                "SST fixture: section iteration failed");
  return names;
}
std::vector<double> sections(hid_t file, const std::string& path, std::size_t n, int components = 1) {
  H5 group(H5Gopen2(file,path.c_str(),H5P_DEFAULT),H5Gclose);
  std::vector<double> out(n*components,std::numeric_limits<double>::quiet_NaN());
  std::vector<bool> covered(n);
  for (const auto& name : children(group)) {
    KUTTA_REQUIRE(!name.empty() && std::all_of(name.begin(),name.end(),[](char c) { return c >= '0' && c <= '9'; }),
                  "SST fixture: unexpected non-numeric section");
    H5 ds(H5Dopen2(group,name.c_str(),H5P_DEFAULT),H5Dclose);
    auto lo = attribute(ds,"minId"), hi = attribute(ds,"maxId");
    KUTTA_REQUIRE(lo > 0 && hi >= lo && std::uint64_t(hi) <= n,"SST fixture: invalid section ids");
    H5 space(H5Dget_space(ds),H5Sclose);
    hsize_t dims[2]{};
    const auto rank = H5Sget_simple_extent_ndims(space);
    KUTTA_REQUIRE(rank == 1 || rank == 2,"SST fixture: unsupported field rank");
    KUTTA_REQUIRE(H5Sget_simple_extent_dims(space,dims,nullptr) >= 0 &&
                  dims[0] == std::uint64_t(hi-lo+1) &&
                  (rank == 2 ? dims[1] == std::uint64_t(components) : components == 1),
                  "SST fixture: incorrect field shape");
    auto values = numeric(ds);
    for (auto id = lo-1; id < hi; ++id) {
      KUTTA_REQUIRE(!covered[id],"SST fixture: overlapping field sections");
      covered[id] = true;
      for (int c = 0; c < components; ++c) {
        const auto value = values[(id-lo+1)*components+c];
        KUTTA_REQUIRE(std::isfinite(value),"SST fixture: nonfinite saved field");
        out[id*components+c] = value;
      }
    }
  }
  return out;
}

// CFF stores boundaries before interiors; PolyMesh reorders them. Map by the
// source zone id and its explicit range, and cross-check each boundary owner.
std::vector<Index> boundary_source_ids(hid_t cas, const PolyMesh& m) {
  const std::string base = "/meshes/1/faces/zoneTopology/";
  auto ids = numeric(cas,base+"id"), low = numeric(cas,base+"minId"), high = numeric(cas,base+"maxId");
  KUTTA_REQUIRE(ids.size() == low.size() && ids.size() == high.size(),"SST fixture: zone shape mismatch");
  std::map<int,std::pair<Index,Index>> ranges;
  for (std::size_t j = 0; j < ids.size(); ++j) ranges.emplace(int(ids[j]),std::pair{Index(low[j])-1,Index(high[j])});
  auto c0 = sections(cas,"/meshes/1/faces/c0",m.n_faces());
  std::vector<Index> source(m.n_faces()-m.n_internal_faces(),-1);
  for (const auto& p : m.patches) {
    const auto [lo,hi] = ranges.at(p.source_zone_id);
    KUTTA_REQUIRE(hi-lo == p.size,"SST fixture: patch/source range mismatch");
    for (Index j = 0; j < p.size; ++j) {
      source[p.start-m.n_internal_faces()+j] = lo+j;
      KUTTA_REQUIRE(c0[lo+j] == m.owner[p.start+j]+1,"SST fixture: cell id mapping differs from saved fields");
    }
  }
  return source;
}
std::vector<double> boundary(const std::vector<double>& raw, const std::vector<Index>& ids, int components = 1) {
  std::vector<double> out;
  for (auto id : ids) for (int c = 0; c < components; ++c) {
    KUTTA_REQUIRE(id >= 0 && std::size_t(id*components+c) < raw.size() && std::isfinite(raw[id*components+c]),
                  "SST fixture: missing boundary field");
    out.push_back(raw[id*components+c]);
  }
  return out;
}
std::vector<Vec3d> xyz(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z) {
  KUTTA_REQUIRE(x.size() == y.size() && y.size() == z.size(),"SST fixture: vector size mismatch");
  std::vector<Vec3d> out;
  for (std::size_t i = 0; i < x.size(); ++i) out.push_back({x[i],y[i],z[i]});
  return out;
}

// tools/forensics/grad.py: inverse-distance-SQUARED LSQ. Values are at face
// centroids (model omega on walls, saved SV_* elsewhere).
// Solve the 3x3 normal equations in double with
// partial pivoting; upload component-row SoA floats to the actual SST kernels.
// This does not reproduce Fluent's undocumented WFGC / Gram-Schmidt weights.
std::vector<std::vector<Vec3d>> gradient(const PolyMesh& m, const MeshGeometry& g,
    const std::vector<std::vector<double>>& values, const std::vector<std::vector<double>>& bvalues) {
  const auto count = values.size();
  KUTTA_REQUIRE(count >= 1 && count <= 3 && bvalues.size() == count,"SST fixture: gradient component count");
  auto cf = build_cell_faces(m);
  std::vector<std::vector<Vec3d>> result(count,std::vector<Vec3d>(m.n_cells));
  for (Index c = 0; c < m.n_cells; ++c) {
    double a[3][6]{};
    for (Index j = cf.offsets[c]; j < cf.offsets[c+1]; ++j) {
      Index f = cf.faces[j]; const bool internal = f < m.n_internal_faces();
      Index other = internal ? (m.owner[f] == c ? m.neighbour[f] : m.owner[f]) : f-m.n_internal_faces();
      Vec3d delta = (internal ? g.cell_centroid[other] : g.face_centroid[f]) - g.cell_centroid[c];
      double dx[3]{delta.x,delta.y,delta.z}; const auto weight = 1/dot(delta,delta);
      for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) a[row][col] += weight*dx[row]*dx[col];
        for (std::size_t comp = 0; comp < count; ++comp)
          a[row][3+comp] += weight*dx[row]*((internal ? values[comp][other] : bvalues[comp][other])-values[comp][c]);
      }
    }
    for (int p = 0; p < 3; ++p) {
      int pivot = p;
      for (int row = p+1; row < 3; ++row) if (std::abs(a[row][p]) > std::abs(a[pivot][p])) pivot = row;
      KUTTA_REQUIRE(std::abs(a[pivot][p]) > 1e-15,"SST fixture: singular LSQ system");
      for (int col = 0; col < 3+int(count); ++col) std::swap(a[p][col],a[pivot][col]);
      auto diagonal = a[p][p];
      for (int col = p; col < 3+int(count); ++col) a[p][col] /= diagonal;
      for (int row = 0; row < 3; ++row) if (row != p) {
        auto factor = a[row][p];
        for (int col = p; col < 3+int(count); ++col) a[row][col] -= factor*a[p][col];
      }
    }
    for (std::size_t comp = 0; comp < count; ++comp)
      result[comp][c] = {a[0][3+comp],a[1][3+comp],a[2][3+comp]};
  }
  return result;
}
double quantile(std::vector<double> x, double q) {
  KUTTA_REQUIRE(!x.empty(),"SST fixture: empty comparison population");
  std::sort(x.begin(),x.end()); double p = q*(x.size()-1);
  auto j = std::size_t(p); auto k = std::min(j+1,x.size()-1);
  return x[j]+(p-j)*(x[k]-x[j]);
}
void report(const char* name, const std::vector<double>& errors) {
  std::cout << "  " << name << " n=" << errors.size() << " median=" << quantile(errors,.5)
    << " p95=" << quantile(errors,.95) << " max=" << quantile(errors,1) << '\n';
}
void check_errors(const std::vector<double>& e, double median, double p95, double maximum) {
  REQUIRE(!e.empty()); CHECK(quantile(e,.5) <= median); CHECK(quantile(e,.95) <= p95); CHECK(quantile(e,1) <= maximum);
}
}  // namespace
#endif

TEST_CASE("SST reproduces Fluent mixing-elbow saved viscosity and wall fields") {
#if !KUTTA_HAS_HDF5
  MESSAGE("SKIP: SST saved-field verification requires HDF5");
#else
  using namespace sst_test;
  const std::filesystem::path root = fluent_examples()/"mixing_elbow";
  const auto case_path = root/"mixing_elbow.cas.h5", data_path = root/"mixing_elbow.dat.h5";
  if (!std::filesystem::exists(case_path) || !std::filesystem::exists(data_path)) {
    MESSAGE("SKIP: external Fluent mixing_elbow .cas.h5/.dat.h5 pair is absent"); return;
  }
  auto poly = read_fluent_cff(case_path);
  // Production geometry must itself follow Fluent's face-pyramid convention;
  // the wall-omega/y+ tolerances below fail by ~1% otherwise.
  const auto geometry = compute_geometry(poly);
  REQUIRE(poly.n_cells == 17822); REQUIRE(poly.n_faces() == 91581);
  H5 cas(H5Fopen(kutta::utf8_path(case_path).c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
  H5 dat(H5Fopen(kutta::utf8_path(data_path).c_str(),H5F_ACC_RDONLY,H5P_DEFAULT),H5Fclose);
  const auto ids = boundary_source_ids(cas,poly);
  auto cell = [&](const char* name) {
    auto v = sections(dat,std::string("/results/1/phase-1/cells/")+name,poly.n_cells);
    KUTTA_REQUIRE(std::all_of(v.begin(),v.end(),[](double x){return std::isfinite(x);}),"SST fixture: incomplete cell field");
    return v;
  };
  const auto nb = ids.size();
  auto face = [&](const char* name) {
    return boundary(sections(dat,std::string("/results/1/phase-1/faces/")+name,poly.n_faces()),ids);
  };
  auto k = cell("SV_K"), w = cell("SV_O"), rho = cell("SV_DENSITY"), mu = cell("SV_MU_LAM");
  auto saved_mut = cell("SV_MU_T"), distance = cell("SV_WALL_DIST");
  auto u = cell("SV_U"), v = cell("SV_V"), velw = cell("SV_W");
  auto bu = face("SV_U"), bv = face("SV_V"), bw = face("SV_W"), bk = face("SV_K"), bo = face("SV_O");
  auto grad = gradient(poly,geometry,{u,v,velw},{bu,bv,bw});
  // k has zero normal gradient at walls; use the saved values here to reproduce
  // forensics, whose measured tiny face/cell discrepancy is separately checked.
  auto gk = gradient(poly,geometry,{k},{bk}), go = gradient(poly,geometry,{w},{bo});
  auto raw_yp = sections(dat,"/results/1/phase-1/faces/SV_WALL_YPLUS",poly.n_faces());
  auto raw_ypu = sections(dat,"/results/1/phase-1/faces/SV_WALL_YPLUS_UTAU",poly.n_faces());
  auto raw_shear = sections(dat,"/results/1/phase-1/faces/SV_WALL_SHEAR",poly.n_faces(),3);
  std::vector<bool> is_wall(nb);
  for (const auto& p : poly.patches) if (p.kind == PatchKind::Wall)
    for (Index f = p.start; f < p.start+p.size; ++f) is_wall[f-poly.n_internal_faces()] = true;
  REQUIRE(std::count(is_wall.begin(),is_wall.end(),true) == 2436);
  std::vector<bool> limited(poly.n_cells);
  for (Index c = 0; c < poly.n_cells; ++c) limited[c] = saved_mut[c]*w[c]/(rho[c]*k[c]) < 1-1e-9;
  REQUIRE(std::count(limited.begin(),limited.end(),true) == 6756);

  KUTTA_FOR_EACH_BACKEND(dev) {
    std::cout << "SST Fluent 22.2 mixing_elbow: " << backend_name(dev->backend()) << " / " << dev->name() << '\n';
    auto m = mesh(*dev,poly,geometry,distance); State s(*dev,m);
    s.velocity = vector(*dev,xyz(u,v,velw));
    s.rho = scalar(*dev,rho); s.mu = scalar(*dev,mu); s.k = scalar(*dev,k); s.omega = scalar(*dev,w);
    s.gk = vector(*dev,gk[0]); s.go = vector(*dev,go[0]);
    s.wall_velocity = vector(*dev,xyz(bu,bv,bw));
    SstModel model(*dev,m);
    // Gradient-free identity: turn off shear to test ALL saved inactive cells,
    // avoiding LSQ's known false limiter activations in ~10% of cells.
    model.update(s.inputs());
    auto no_shear = dev->download(model.fields().mu_t);
    std::vector<double> identity;
    for (Index c = 0; c < poly.n_cells; ++c) if (!limited[c])
      identity.push_back(std::abs(no_shear[c]/saved_mut[c]-1));
    REQUIRE(identity.size() == 11066); check_errors(identity,2e-7,3e-7,5e-7);
    report("mu_t inactive identity",identity);

    s.gu = vector(*dev,grad[0]); s.gv = vector(*dev,grad[1]); s.gw = vector(*dev,grad[2]);
    // Preserve the old saved-SV_O reconstruction as a comparison only. The
    // production-style path below replaces every wall value with the model's
    // reconstruction omega; it never uses saved wall SV_O for its gradients.
    model.update(s.inputs());
    const auto saved_face_mut = dev->download(model.fields().mu_t);
    const auto saved_face_f2 = dev->download(model.fields().f2);
    const auto saved_face_cdkw = dev->download(model.fields().cdkw);
    const auto model_face_w = dev->download(model.fields().wall_face_omega);
    auto model_bo = bo;
    std::vector<double> gradient_change;
    for (std::size_t b = 0; b < nb; ++b) if (is_wall[b]) model_bo[b] = model_face_w[b];
    const auto model_go = gradient(poly,geometry,{w},{model_bo});
    for (Index c = 0; c < poly.n_cells; ++c) if (norm(go[0][c]) > 0)
      gradient_change.push_back(norm(model_go[0][c]-go[0][c])/norm(go[0][c]));
    s.go = vector(*dev,model_go[0]);
    model.update(s.inputs()); const auto& f = model.fields();
    auto mut = dev->download(f.mu_t), strain = dev->download(f.strain), vort = dev->download(f.vorticity);
    auto f2 = dev->download(f.f2), wall_w = dev->download(f.fixed_omega);
    auto yp = dev->download(f.wall_y_plus), ypu = dev->download(f.wall_y_plus_utau);
    auto tau = dev->download(f.wall_tau), sx = dev->download(f.wall_shear.x);
    auto sy = dev->download(f.wall_shear.y), sz = dev->download(f.wall_shear.z);
    auto selected = dev->download(f.selected_wall_face);
    std::vector<double> active_errors, free_errors, discriminating_s, discriminating_o;
    int agreement = 0, limiter_changes = 0, cdkw_changes = 0;
    double max_mut_change = 0, max_production_ratio = 0;
    for (Index c = 0; c < poly.n_cells; ++c) {
      (limited[c] ? active_errors : free_errors).push_back(std::abs(mut[c]/saved_mut[c]-1));
      agreement += (strain[c]*f2[c] > .31*w[c]) == limited[c];
      limiter_changes += (strain[c]*f2[c] > .31*w[c]) != (strain[c]*saved_face_f2[c] > .31*w[c]);
      max_mut_change = std::max(max_mut_change,std::abs(double(mut[c])/saved_face_mut[c]-1));
      const auto exact_xf2 = rho[c]*.31*k[c]/saved_mut[c];
      if (limited[c] && std::abs(strain[c]-vort[c])/std::max(strain[c],vort[c]) > .2) {
        discriminating_s.push_back(strain[c]*f2[c]/exact_xf2);
        discriminating_o.push_back(vort[c]*f2[c]/exact_xf2);
      }
      max_production_ratio = std::max(max_production_ratio,saved_mut[c]*strain[c]*strain[c]/(10*rho[c]*.09*k[c]*w[c]));
    }
    const auto model_cdkw = dev->download(f.cdkw);
    for (Index c = 0; c < poly.n_cells; ++c) cdkw_changes += model_cdkw[c] != saved_face_cdkw[c];
    REQUIRE(cdkw_changes > 2000); // the reconstructed omega gradient really changed
    CHECK(max_mut_change == 0); CHECK(limiter_changes == 0);
    report("omega gradient change vs saved wall SV_O",gradient_change);
    std::cout << "  model wall omega: max mu_t change=" << max_mut_change
      << " limiter changes=" << limiter_changes << " CDkw changes=" << cdkw_changes << '\n';
    report("mu_t limited, inverse-square LSQ",active_errors);
    report("mu_t inactive, inverse-square LSQ",free_errors);
    std::cout << "  limiter agreement=" << double(agreement)/poly.n_cells
      << " discriminating n=" << discriminating_s.size() << " median S ratio=" << quantile(discriminating_s,.5)
      << " median Omega ratio=" << quantile(discriminating_o,.5)
      << " max production/cap=" << max_production_ratio << '\n';
    // Forensics §2: median 4.6%, p95 14.3% in limited cells; this tolerance
    // explicitly includes reconstruction scatter, not claimed SST model error.
    check_errors(active_errors,.055,.16,.32);
    CHECK(quantile(free_errors,.5) < 5e-7); CHECK(quantile(free_errors,.95) < .05);
    CHECK(double(agreement)/poly.n_cells > .895);
    REQUIRE(discriminating_s.size() > 2500);
    CHECK(std::abs(quantile(discriminating_s,.5)-1) < .02);
    CHECK(quantile(discriminating_o,.5) < .35);
    CHECK(max_production_ratio < 1); // inactive in this fixture; not evidence of its form

    const auto face_k_diffusion = dev->download(f.face_diffusivity_k);
    std::vector<double> face_omega_errors;
    std::vector<double> omega_errors, yp_errors, ypu_errors, shear_errors, shear_vector_errors, wall_k_errors;
    std::vector<double> shear_oracle_errors, direction_errors;
    for (std::size_t b = 0; b < nb; ++b) if (is_wall[b]) {
      const auto face_id = Index(b)+poly.n_internal_faces(), c = poly.owner[face_id], src = ids[b];
      REQUIRE(selected[c] == Index(b)); // this case has exactly one wall face/cell
      REQUIRE(bu[b] == 0); REQUIRE(bv[b] == 0); REQUIRE(bw[b] == 0);
      Vec3d force{raw_shear[3*src],raw_shear[3*src+1],raw_shear[3*src+2]};
      const auto area = norm(geometry.face_area[face_id]);
      REQUIRE(std::isfinite(norm(force))); REQUIRE(norm(force) > 0);
      face_omega_errors.push_back(std::abs(model_face_w[b]/bo[b]-1));
      CHECK(face_k_diffusion[face_id] == 0);
      omega_errors.push_back(std::abs(wall_w[c]/w[c]-1));
      yp_errors.push_back(std::abs(yp[b]/raw_yp[src]-1));
      ypu_errors.push_back(std::abs(ypu[b]/raw_ypu[src]-1));
      shear_errors.push_back(std::abs(tau[b]*area/norm(force)-1));
      shear_vector_errors.push_back(norm(Vec3d{sx[b]*area,sy[b]*area,sz[b]*area}-force)/norm(force));
      // Independently evaluate the full vector in double. Forensics only gives
      // magnitude-error quantiles and MEDIAN direction cosine. Do not apply its
      // magnitude tolerance to a different (full-vector) error metric.
      const auto normal = geometry.face_area[face_id]/area;
      const auto delta = geometry.face_centroid[face_id]-geometry.cell_centroid[c];
      const auto dy = dot(delta,normal);
      const Vec3d relative{u[c]-bu[b],v[c]-bv[b],velw[c]-bw[b]};
      const auto tangent = relative-normal*dot(relative,normal);
      const auto speed = norm(tangent), us = std::sqrt(mu[c]/rho[c]*speed/dy+.3*k[c]);
      const auto y = rho[c]/mu[c]*us*dy, log_up = std::log(9.793*std::max(y,.2))/.4187;
      const auto up = std::pow(std::pow(y,-4)+std::pow(log_up,-4),-.25);
      const auto oracle = tangent*(-rho[c]*us*area/up);
      shear_oracle_errors.push_back(norm(Vec3d{sx[b]*area,sy[b]*area,sz[b]*area}-oracle)/norm(oracle));
      direction_errors.push_back(1+dot(force,tangent)/(norm(force)*speed));
      wall_k_errors.push_back(std::abs(bk[b]/k[c]-1));
    }
    report("wall-face reconstruction omega",face_omega_errors);
    // Fixed constants retain the ~0.42% systematic offset in forensics 3.5.
    check_errors(face_omega_errors,.0043,.0044,.0045);
    report("wall-cell omega",omega_errors); report("wall y+ (u*)",yp_errors);
    report("wall y+ (u_tau)",ypu_errors); report("wall shear force magnitude",shear_errors);
    report("wall shear force vector",shear_vector_errors);
    report("wall shear vector kernel/double-oracle",shear_oracle_errors);
    report("saved shear directional scatter (1+cos)",direction_errors);
    report("wall k zero gradient",wall_k_errors);
    // Forensics §3 maxima: omega 2.93e-4, y+ 1.76e-3, shear 3.94e-3.
    // u_tau here uses OUR predicted shear; its error includes the shear fit.
    check_errors(omega_errors,1e-5,5e-5,3.2e-4);
    check_errors(yp_errors,1e-5,2.5e-4,2e-3);
    check_errors(ypu_errors,1.5e-5,2.5e-4,2.1e-3);
    check_errors(shear_errors,1.5e-5,4e-4,4e-3);
    check_errors(shear_oracle_errors,3e-7,8e-7,3e-6);
    CHECK(quantile(direction_errors,.5) < 1e-8);
    check_errors(wall_k_errors,1e-5,7e-4,5e-3);
  }
#endif
}
