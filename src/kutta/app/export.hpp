#pragma once
#include "kutta/io/view.hpp"
#include "kutta/setup/fluent.hpp"
namespace kutta::app {
int export_view(int argc,char** argv);
void save_view_input(const std::filesystem::path&,Device&,const fv::SteadySolver&,
                     const fluent::FluentCase&,const fluent::FluentData*,const setup::FluentSetup&,std::string_view provenance);
}
