#include "kutta/case/fluent_case.hpp"
#include "kutta/case/fluent_data.hpp"
#include "solve.hpp"
#include "protocol.hpp"
#include "export.hpp"
#include <algorithm>
#include <iostream>

using namespace kutta::fluent;
namespace {
template <class T> void value(const std::optional<T> &v) {
  if (v)
    std::cout << *v;
  else
    std::cout << "unparsed/absent";
}
void scalar(const std::optional<ScalarValue> &v) {
  if (!v) {
    std::cout << "unparsed/absent";
    return;
  }
  if (v->kind == ScalarValue::Kind::constant)
    value(v->constant);
  else {
    std::cout << "binding(" << v->binding << ')';
    if (v->resolved) {
      std::cout << " resolved=" << *v->resolved << " (resolved from expression";
      if (!v->expression_name.empty())
        std::cout << ' ' << v->expression_name;
      std::cout << ": " << v->definition << ')';
    }
  }
}
void vector(const ScalarVector &v) {
  std::cout << '(';
  for (int i = 0; i < 3; ++i) {
    if (i)
      std::cout << ',';
    scalar(v[i]);
  }
  std::cout << ')';
}
void mesh_summary(const FluentCase &c) {
  std::cout << "Fluent ";
  value(c.solver_version);
  std::cout << "; dimension=";
  value(c.dimension);
  std::cout << "; cells=";
  value(c.mesh.cells);
  std::cout << " faces=";
  value(c.mesh.faces);
  std::cout << " nodes=";
  value(c.mesh.nodes);
  std::cout << '\n';
  for (auto &z : c.mesh.cell_zones)
    std::cout << "  cell zone " << z.zone_id << ' ' << z.name << ": " << z.size() << " [" << z.first
              << ".." << z.last << "]\n";
  for (auto &z : c.mesh.face_zones)
    std::cout << "  face zone " << z.zone_id << ' ' << z.name << ": " << z.size() << " [" << z.first
              << ".." << z.last << "]\n";
}
void model_summary(const FluentCase &c) {
  std::cout << "Model: " << (c.viscous.model ? name(*c.viscous.model) : "unparsed")
            << "; production-limit=";
  value(c.viscous.production_limiter);
  std::cout << " factor=";
  value(c.viscous.production_factor);
  std::cout << " low-Re=";
  value(c.viscous.low_re);
  std::cout << " curvature=";
  value(c.viscous.curvature_correction);
  std::cout << " SST-strain=";
  value(c.viscous.sst_strain);
  std::cout << " wall-omega="
            << (c.viscous.wall_omega
                    ? (*c.viscous.wall_omega == WallOmega::correlation ? "correlation"
                                                                       : "tabulated")
                    : "unparsed/absent")
            << '\n';
}
void numerics_summary(const FluentCase &c) {
  auto scheme = [](const char *label, auto s) {
    std::cout << label << '=' << (s ? name(*s) : "unparsed") << ' ';
  };
  const bool epsilon = c.viscous.model == ViscousModel::ke_standard ||
                       c.viscous.model == ViscousModel::ke_rng ||
                       c.viscous.model == ViscousModel::ke_realizable;
  std::cout << "Schemes: ";
  scheme("p", c.numerics.pressure);
  scheme("mom", c.numerics.momentum);
  if (c.viscous.model != ViscousModel::laminar) {
    scheme("k", c.numerics.k);
    if (epsilon)
      scheme("epsilon", c.numerics.epsilon);
    else
      scheme("omega", c.numerics.omega);
  }
  std::cout << "coupling=" << (c.numerics.coupling ? name(*c.numerics.coupling) : "unparsed")
            << " gradient=";
  if (c.numerics.gradient)
    std::cout << (*c.numerics.gradient == Gradient::least_squares_cell
                      ? "LSQ-cell"
                      : (*c.numerics.gradient == Gradient::green_gauss_node ? "GG-node (unverified mapping)"
                                                                            : "GG-cell"));
  else
    std::cout << "unparsed";
  std::cout << " limiter=";
  if (c.numerics.limiter)
    std::cout << static_cast<int>(*c.numerics.limiter);
  else
    std::cout << "unparsed";
  std::cout << "\nControls: pseudo=";
  value(c.numerics.pseudo.formulation);
  std::cout << " dt-method=";
  value(c.numerics.coupling == Coupling::coupled ? c.numerics.pseudo.coupled_method
                                                 : c.numerics.pseudo.segregated_method);
  std::cout << " scale=";
  value(c.numerics.pseudo.scale_factor);
  std::cout << " CFL=";
  value(c.numerics.courant);
  std::cout << " length-scale-code=";
  value(c.numerics.pseudo.length_scale_method);
  const auto &n = c.numerics;
  const bool coupled = n.coupling == Coupling::coupled;
  const auto method = coupled ? n.pseudo.coupled_method : n.pseudo.segregated_method;
  if (!n.coupling || !n.pseudo.formulation ||
      (*n.pseudo.formulation && (!method || *method == 0))) {
    std::cout << "\nActive under-relaxation (unverified selection): unparsed";
  } else if (*n.pseudo.formulation) {
    std::cout << "\nActive under-relaxation ("
              << (coupled ? "coupled pseudo-transient" : "segregated pseudo-transient") << "):";
    std::vector<std::string> equations = {"pressure", "mom", "density", "body-force"};
    if (c.viscous.model != ViscousModel::laminar) {
      equations.push_back("k");
      equations.push_back(epsilon ? "epsilon" : "omega");
      equations.push_back("turb-viscosity");
    }
    for (const auto &eq : equations) {
      const auto key = coupled && (eq == "pressure" || eq == "mom")
          ? "pressure-coupled/" + eq + "/pseudo-explicit-relax" : eq + "/pseudo-relax";
      std::cout << ' ' << eq << '=';
      if (auto it = n.pseudo_relaxation.find(key); it != n.pseudo_relaxation.end())
        std::cout << it->second;
      else
        std::cout << "unparsed/absent";
    }
  } else {
    std::cout << "\nActive under-relaxation (segregated URFs):";
    for (auto &[eq, urf] : n.relaxation)
      std::cout << ' ' << eq << '=' << urf;
  }
  std::cout << '\n';
}
void material_summary(const FluentCase &c) {
  for (auto &m : c.materials) {
    std::cout << "Material " << m.name << " (" << m.type << "): rho=";
    scalar(m.density);
    std::cout << " mu=";
    scalar(m.viscosity);
    std::cout << '\n';
  }
  for (auto &z : c.cell_zones) {
    std::cout << "Cell zone " << z.id << ' ' << z.name << " (" << z.fluent_type << "): material=";
    value(z.material);
    std::cout << " MRF=";
    value(z.frame_motion);
    std::cout << " omega=";
    scalar(z.angular_velocity);
    std::cout << " porous=";
    value(z.porous);
    if (z.frame_motion == true) {
      std::cout << " axis=";
      vector(z.axis_direction);
      std::cout << " origin=";
      vector(z.axis_origin);
    }
    if (z.porous == true) {
      std::cout << " viscous=";
      vector(z.viscous_resistance);
      std::cout << " inertial=";
      vector(z.inertial_resistance);
    }
    std::cout << '\n';
  }
}
void boundary_summary(const FluentCase &c) {
  for (auto &z : c.boundaries) {
    std::cout << "BC " << z.id << ' ' << z.name << " (" << z.fluent_type << ')';
    if (z.type == ZoneType::velocity_inlet) {
      std::cout << " velocity-spec=";
      if (z.velocity_spec)
        std::cout << static_cast<int>(*z.velocity_spec);
      else
        std::cout << "unparsed";
      if (z.velocity_spec == VelocitySpec::components) {
        std::cout << " uvw=";
        vector(z.velocity);
      } else {
        std::cout << " speed=";
        scalar(z.velocity_magnitude);
        if (z.velocity_spec == VelocitySpec::magnitude_direction) {
          std::cout << " direction=";
          vector(z.direction);
        }
      }
    }
    if (z.type == ZoneType::pressure_outlet) {
      std::cout << " gauge-p=";
      scalar(z.gauge_pressure);
    }
    if (z.turbulence_spec) {
      std::cout << " turb-spec=" << static_cast<int>(*z.turbulence_spec);
      const bool omega = c.viscous.model == ViscousModel::kw_sst ||
                         c.viscous.model == ViscousModel::kw_bsl ||
                         c.viscous.model == ViscousModel::kw_standard;
      for (const auto &key : z.turbulence) {
        const auto &k = key.first;
        bool active = (*z.turbulence_spec == TurbulenceSpec::k_dissipation &&
                       (k == "k" || k == (omega ? "omega" : "epsilon"))) ||
                      (*z.turbulence_spec != TurbulenceSpec::k_dissipation &&
                       (k == "turb-intensity" ||
                        (*z.turbulence_spec == TurbulenceSpec::intensity_length &&
                         k == "turb-length-scale") ||
                        (*z.turbulence_spec == TurbulenceSpec::intensity_diameter &&
                         k == "turb-hydraulic-diam") ||
                        (*z.turbulence_spec == TurbulenceSpec::intensity_viscosity_ratio &&
                         k == "turb-viscosity-ratio")));
        if (active) {
          std::cout << ' ' << k << '=';
          scalar(key.second);
        }
      }
    }
    if (z.type == ZoneType::wall) {
      std::cout << " motion=";
      if (z.wall_motion)
        std::cout << (*z.wall_motion == WallMotion::stationary
                          ? "stationary"
                          : (*z.wall_motion == WallMotion::rotating ? "rotating" : "translating"));
      else
        std::cout << "unparsed";
      std::cout << " omega=";
      scalar(z.angular_velocity);
      std::cout << " roughness=";
      scalar(z.roughness_height);
      if (z.wall_motion == WallMotion::rotating) {
        std::cout << " axis-origin=";
        vector(z.axis_origin);
        std::cout << " axis=";
        vector(z.axis_direction);
      }
      if (z.wall_motion == WallMotion::translating) {
        if (z.wall_components == true) {
          std::cout << " components uvw=";
          vector(z.velocity);
        } else if (z.wall_components == false) {
          std::cout << " magnitude/direction speed=";
          scalar(z.velocity_magnitude);
          std::cout << " direction=";
          vector(z.direction);
        } else
          std::cout << " velocity representation=unparsed";
      }
      if (z.specified_shear == true) {
        std::cout << " specified-shear=";
        vector(z.shear);
      }
    }
    std::cout << '\n';
  }
}
void reference_summary(const FluentCase &c) {
  std::cout << "Operating pressure=";
  value(c.operating_pressure);
  std::cout << "\nReference: A=";
  value(c.reference.area);
  std::cout << " L=";
  value(c.reference.length);
  std::cout << " rho=";
  value(c.reference.density);
  std::cout << " U=";
  value(c.reference.velocity);
  std::cout << " p=";
  value(c.reference.pressure);
  std::cout << " mu=";
  value(c.reference.viscosity);
  std::cout << '\n';
}
void report_summary(const FluentCase &c) {
  std::cout << "Reports: " << c.reports.size() << '\n';
  for (auto &r : c.reports) {
    std::cout << "  " << r.name << " (" << r.fluent_type << ')';
    if (r.scaled)
      std::cout << (*r.scaled ? " coefficient" : " force (unscaled)");
    if (r.surface_ids) {
      std::cout << " surface-ids=";
      for (auto id : *r.surface_ids)
        std::cout << id << ',';
    }
    if (r.zone_ids) {
      std::cout << " zones=";
      for (auto id : *r.zone_ids)
        std::cout << id << ',';
    }
    if (r.zone_names) {
      std::cout << " names=";
      for (auto &s : *r.zone_names)
        std::cout << s << ',';
    }
    if (r.direction) {
      std::cout << " direction=";
      for (auto d : *r.direction)
        std::cout << d << ',';
    }
    std::cout << '\n';
  }
}
void unparsed_summary(const FluentCase &c) {
  const auto info = std::count_if(c.issues.begin(), c.issues.end(),
      [](const auto &i) { return i.classification == IssueClass::informational; });
  std::cout << "Issues: unsupported/required=" << c.issues.size() - info
            << " informational=" << info << '\n';
  for (auto &issue : c.issues)
    std::cout << "  [" << (issue.classification == IssueClass::informational ? "informational" : "unsupported")
              << "] " << issue.path << ": " << issue.reason << '\n';
  auto raw = c.variables.unparsed_names();
  std::cout << "Raw rpvars retained: " << raw.size() << " unparsed names";
  for (std::size_t i = 0; i < std::min<std::size_t>(8, raw.size()); ++i)
    std::cout << ' ' << raw[i];
  std::cout << " ...\n";
  for (auto &z : c.boundaries)
    std::cout << "  zone " << z.id << ": " << z.raw.unparsed_names().size()
              << " raw unparsed settings retained\n";
}
void summary(const FluentCase &c) {
  mesh_summary(c);
  model_summary(c);
  std::cout << "Physics: energy=";
  value(c.energy);
  std::cout << " gravity=";
  value(c.gravity);
  if (c.gravity == true) {
    std::cout << " vector=(";
    for (int i = 0; i < 3; ++i) {
      if (i)
        std::cout << ',';
      value(c.gravity_vector[i]);
    }
    std::cout << ')';
  }
  std::cout << '\n';
  numerics_summary(c);
  material_summary(c);
  boundary_summary(c);
  reference_summary(c);
  report_summary(c);
  unparsed_summary(c);
}
void data_summary(const FluentCase &c, const FluentData &d) {
  std::cout << "Data Case File=";
  value(d.case_file);
  std::cout << "; loaded case=";
  if (c.source_path)
    std::cout << c.source_path->filename().string();
  else
    std::cout << "unavailable";
  std::cout << '\n';
  std::cout << "Data: " << d.cell_fields.size() << " cell fields, " << d.face_fields.size()
            << " face fields; saved iteration=";
  value(d.iteration_count);
  std::cout << " last residual iteration=";
  value(d.last_residual_iteration);
  std::cout << '\n';
  for (auto &[eq, r] : d.residuals)
    std::cout << "  residual " << eq << ": " << r.rows.size() << " retained samples\n";
  for (auto &z : c.mesh.face_zones) {
    std::cout << "  face zone " << z.zone_id << ' ' << z.name << ':';
    for (auto name : {"SV_P", "SV_WALL_SHEAR", "SV_WALL_YPLUS", "SV_FLUX"})
      if (d.face_fields.contains(name)) {
        std::size_t n = 0;
        for (auto &v : d.zone_field(FieldLocation::face, z.zone_id, name))
          n += static_cast<std::size_t>(v.last - v.first + 1);
        std::cout << ' ' << name << '=' << n;
      }
    std::cout << '\n';
  }
  const auto info = std::count_if(d.issues.begin(), d.issues.end(),
      [](const auto &i) { return i.classification == IssueClass::informational; });
  std::cout << "Data issues: unsupported/required=" << d.issues.size() - info
            << " informational=" << info << '\n';
  for (auto &i : d.issues)
    std::cout << "  [" << (i.classification == IssueClass::informational ? "informational" : "unsupported")
              << "] " << i.path << ": " << i.reason << '\n';
}
} // namespace
int main(int argc, char **argv) {
  if(argc>=2 && std::string_view(argv[1])=="solve") return solve(argc,argv);
  if(argc>=2 && std::string_view(argv[1])=="export-view") return kutta::app::export_view(argc,argv);
  if(argc>=2 && std::string_view(argv[1])=="export-geometry") return kutta::app::export_geometry(argc,argv);
  if(argc>=2 && std::string_view(argv[1])=="aero") return kutta::app::aero(argc,argv);
  if(argc>=2 && std::string_view(argv[1])=="devices") return kutta::app::devices(argc,argv);
  if(argc>=2 && std::string_view(argv[1])=="inspect")
    for(int i=2;i<argc;++i) if(std::string_view(argv[i])=="--json") return kutta::app::inspect_json(argc,argv);
  if ((argc != 3 && argc != 5) || std::string_view(argv[1]) != "inspect" ||
      (argc == 5 && std::string_view(argv[3]) != "--data")) {
    std::cerr << "Usage: kutta solve <file.cas.h5> [--data file.dat.h5] [--backend cuda|metal|cpu|auto] [--iterations N] [--out DIR] [--approximate SETTING[,SETTING...]]\n       kutta inspect <file.cas.h5|file.cas> [--data file.dat.h5]\n";
    return 2;
  }
  try {
    auto c = read_case(argv[2]);
    summary(c);
    if (argc == 5) {
      auto d = read_data(argv[4], c);
      data_summary(c, d);
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "kutta inspect: " << e.what() << '\n';
    return 1;
  }
}
