#pragma once
// Run a test body once per backend that is compiled in and present on this machine.
#include <doctest/doctest.h>

#include <string>

#include "kutta/runtime/device.hpp"

#define KUTTA_FOR_EACH_BACKEND(dev)                                      \
  for (::kutta::Backend kutta_b_ : ::kutta::Device::available())         \
    for (auto dev = ::kutta::Device::create(kutta_b_); dev; dev.reset()) \
      if (INFO("backend: " << std::string(::kutta::backend_name(kutta_b_))); true)
