#pragma once
#include "kutta/io/view.hpp"
#include "kutta/setup/fluent.hpp"
namespace kutta::app {
int export_view(int argc,char** argv);
int export_geometry(int argc,char** argv);
view::Fields collect_view_fields(Device&,const fv::SteadySolver&,const fluent::FluentCase&,
                     const fluent::FluentData*,const setup::FluentSetup&,std::string_view provenance,bool volume);
void save_view_input(const std::filesystem::path&,Device&,const fv::SteadySolver&,
                     const fluent::FluentCase&,const fluent::FluentData*,const setup::FluentSetup&,std::string_view provenance);
}
