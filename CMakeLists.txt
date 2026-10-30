cmake_minimum_required(VERSION 3.26)
project(kutta LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(KUTTA_ENABLE_METAL "Build the Metal backend (Apple only)" ${APPLE})
option(KUTTA_ENABLE_CUDA "Build the CUDA backend" OFF)
option(KUTTA_BUILD_TESTS "Build tests" ON)

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# ---- Kernels: one source, generated entry points per backend -----------------
file(GLOB KUTTA_KERNEL_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/kutta/kernels/*.kern.h)
set(KUTTA_GEN_DIR ${CMAKE_BINARY_DIR}/generated)
set(KUTTA_GEN_OUTPUTS
  ${KUTTA_GEN_DIR}/kernels.hpp ${KUTTA_GEN_DIR}/kernels_host.cpp ${KUTTA_GEN_DIR}/kernels_cpu.cpp
  ${KUTTA_GEN_DIR}/kernels_cuda.cu ${KUTTA_GEN_DIR}/kernels_metal.cpp)
add_custom_command(
  OUTPUT ${KUTTA_GEN_OUTPUTS}
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/kernelgen.py
          --src-root ${CMAKE_SOURCE_DIR}/src --out ${KUTTA_GEN_DIR} ${KUTTA_KERNEL_SOURCES}
  DEPENDS ${CMAKE_SOURCE_DIR}/tools/kernelgen.py ${CMAKE_SOURCE_DIR}/src/kutta/kernels/prelude.h
          ${KUTTA_KERNEL_SOURCES}
  COMMENT "Generating kernel entry points")

# ---- Core library ---------------------------------------------------------------
file(GLOB_RECURSE KUTTA_CORE_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_SOURCE_DIR}/src/kutta/case/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/mesh/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/io/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/linalg/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/fv/*.cpp)
add_library(kutta
  ${KUTTA_CORE_SOURCES}
  src/kutta/runtime/device.cpp
  src/kutta/runtime/cpu/cpu_device.cpp
  ${KUTTA_GEN_DIR}/kernels_host.cpp
  ${KUTTA_GEN_DIR}/kernels_cpu.cpp)
target_include_directories(kutta PUBLIC ${CMAKE_SOURCE_DIR}/src ${KUTTA_GEN_DIR})
if(MSVC)
  target_compile_options(kutta PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4;/permissive->)
else()
  target_compile_options(kutta PRIVATE $<$<COMPILE_LANGUAGE:CXX,OBJCXX>:-Wall;-Wextra;-Wpedantic>)
endif()

# ---- Optional Fluent CFF reader (HDF5 C API) ---------------------------------
option(KUTTA_ENABLE_HDF5 "Enable Fluent CFF mesh import when HDF5 is available" ON)
set(KUTTA_HAS_HDF5 OFF)
if(KUTTA_ENABLE_HDF5)
  # FindHDF5's C component runs C probes, even though Kutta is C++.
  enable_language(C)
  find_package(HDF5 1.10.2 QUIET COMPONENTS C)
  if(HDF5_FOUND)
    set(KUTTA_HAS_HDF5 ON)
    find_package(ZLIB 1.2.9 REQUIRED)
    target_link_libraries(kutta PUBLIC HDF5::HDF5 ZLIB::ZLIB)
  endif()
endif()
target_compile_definitions(kutta PUBLIC KUTTA_HAS_HDF5=$<BOOL:${KUTTA_HAS_HDF5}>)
message(STATUS "Kutta HDF5 mesh reader: ${KUTTA_HAS_HDF5}")

if(KUTTA_ENABLE_METAL)
  enable_language(OBJCXX)
  set(CMAKE_OBJCXX_STANDARD 20)
  set(CMAKE_OBJCXX_STANDARD_REQUIRED ON)
  set_target_properties(kutta PROPERTIES OBJCXX_STANDARD 20 OBJCXX_STANDARD_REQUIRED ON)
  file(GLOB KUTTA_METAL_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/kutta/runtime/metal/*.mm)
  if(KUTTA_METAL_SOURCES)
    target_sources(kutta PRIVATE ${KUTTA_METAL_SOURCES} ${KUTTA_GEN_DIR}/kernels_metal.cpp)
    set_source_files_properties(${KUTTA_METAL_SOURCES} PROPERTIES COMPILE_FLAGS "-fobjc-arc")
    target_link_libraries(kutta PUBLIC "-framework Metal" "-framework Foundation")
    target_compile_definitions(kutta PUBLIC KUTTA_HAS_METAL=1)
  endif()
endif()

if(KUTTA_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75 86 89)
  endif()
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set_target_properties(kutta PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
  file(GLOB KUTTA_CUDA_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/kutta/runtime/cuda/*.cu)
  target_sources(kutta PRIVATE ${KUTTA_CUDA_SOURCES} ${KUTTA_GEN_DIR}/kernels_cuda.cu)
  target_link_libraries(kutta PUBLIC CUDA::cudart_static)
  target_compile_definitions(kutta PUBLIC KUTTA_HAS_CUDA=1)
endif()

# ---- Executables, tests, benchmarks ------------------------------------------------
file(GLOB KUTTA_APP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/kutta/app/*.cpp)
if(KUTTA_APP_SOURCES)
  add_executable(kutta-cli ${KUTTA_APP_SOURCES})
  set_target_properties(kutta-cli PROPERTIES OUTPUT_NAME kutta)
  target_link_libraries(kutta-cli PRIVATE kutta)
endif()

file(GLOB KUTTA_BENCH_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/bench/*.cpp)
foreach(src ${KUTTA_BENCH_SOURCES})
  get_filename_component(name ${src} NAME_WE)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE kutta)
  execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE KUTTA_BENCH_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
  target_compile_definitions(${name} PRIVATE KUTTA_BENCH_COMMIT="${KUTTA_BENCH_COMMIT}")
endforeach()

if(KUTTA_BUILD_TESTS)
  include(FetchContent)
  FetchContent_Declare(doctest
    URL https://github.com/doctest/doctest/archive/refs/tags/v2.4.12.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
  FetchContent_MakeAvailable(doctest)
  enable_testing()
  file(GLOB_RECURSE KUTTA_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/tests/*.cpp)
  add_executable(kutta-tests ${KUTTA_TEST_SOURCES})
  target_link_libraries(kutta-tests PRIVATE kutta doctest::doctest)
  add_dependencies(kutta-tests kutta-cli)
  target_compile_definitions(kutta-tests PRIVATE KUTTA_CLI_PATH="$<TARGET_FILE:kutta-cli>")
  target_compile_definitions(kutta-tests PRIVATE KUTTA_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  # Long numerical verification (doctest suite "slow"): `ctest -LE slow` skips it.
  add_test(NAME kutta-tests COMMAND kutta-tests --test-suite-exclude=slow)
  add_test(NAME kutta-tests-slow COMMAND kutta-tests --test-suite=slow)
  set_tests_properties(kutta-tests-slow PROPERTIES LABELS slow)
  add_test(NAME lint-portability COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/lint_portability.py)
endif()
