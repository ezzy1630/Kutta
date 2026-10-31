cmake_minimum_required(VERSION 3.26)
project(kutta VERSION 0.1.0 LANGUAGES C CXX)

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
  ${CMAKE_SOURCE_DIR}/src/kutta/setup/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/case/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/core/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/mesh/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/io/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/post/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/linalg/*.cpp
  ${CMAKE_SOURCE_DIR}/src/kutta/fv/*.cpp)
add_library(kutta
  ${KUTTA_CORE_SOURCES}
  src/kutta/runtime/device.cpp
  src/kutta/runtime/cpu/cpu_device.cpp
  ${KUTTA_GEN_DIR}/kernels_host.cpp
  ${KUTTA_GEN_DIR}/kernels_cpu.cpp)
target_include_directories(kutta PUBLIC ${CMAKE_SOURCE_DIR}/src ${KUTTA_GEN_DIR})
find_package(Threads REQUIRED)
target_link_libraries(kutta PUBLIC Threads::Threads)
# core/memory.cpp reads the process counters (GetProcessMemoryInfo).
if(WIN32)
  target_link_libraries(kutta PUBLIC psapi)
endif()
if(MSVC)
  target_compile_options(kutta PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/W4;/permissive->)
else()
  target_compile_options(kutta PRIVATE $<$<COMPILE_LANGUAGE:CXX,OBJCXX>:-Wall;-Wextra;-Wpedantic>)
endif()

# ---- Optional Fluent CFF reader (HDF5 C API) ---------------------------------
option(KUTTA_ENABLE_HDF5 "Enable Fluent CFF mesh import when HDF5 is available" ON)
option(KUTTA_REQUIRE_HDF5 "Fail configuration unless the Fluent HDF5 reader is enabled" OFF)
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
if(KUTTA_REQUIRE_HDF5 AND NOT KUTTA_HAS_HDF5)
  message(FATAL_ERROR "KUTTA_REQUIRE_HDF5 requires HDF5 >= 1.10.2 and zlib")
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
  # nvcc fuses a*b+c into FMA by default; the CPU oracle on x86 (GCC, MSVC)
  # does not, so CPU==CUDA bitwise parity fails. Measured on a T4: FV and AMG
  # throughput unchanged (src/kutta/runtime/cuda/GPU-RUN-2026-10-30.md).
  target_compile_options(kutta PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--fmad=false>)
  target_link_libraries(kutta PUBLIC CUDA::cudart_static)
  target_compile_definitions(kutta PUBLIC KUTTA_HAS_CUDA=1)
endif()

# ---- Executables, tests, benchmarks ------------------------------------------------
file(GLOB KUTTA_APP_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/src/kutta/app/*.cpp)
if(KUTTA_APP_SOURCES)
  add_executable(kutta-cli ${KUTTA_APP_SOURCES})
  # MSVC writes an import library next to kutta.exe; keep it from colliding
  # with the static library kutta.lib (LNK1149).
  set_target_properties(kutta-cli PROPERTIES OUTPUT_NAME kutta ARCHIVE_OUTPUT_NAME kutta-cli)
  target_link_libraries(kutta-cli PRIVATE kutta)
  execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE KUTTA_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT KUTTA_COMMIT)
    set(KUTTA_COMMIT "unknown")
  else()
    # Reconfigure after commits, including branches in linked worktrees, so an
    # incremental build cannot claim an older source commit in summary.json.
    execute_process(COMMAND git rev-parse --path-format=absolute --git-path HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE KUTTA_GIT_HEAD OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND git symbolic-ref -q HEAD WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
      OUTPUT_VARIABLE KUTTA_GIT_REF OUTPUT_STRIP_TRAILING_WHITESPACE)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${KUTTA_GIT_HEAD}")
    if(KUTTA_GIT_REF)
      execute_process(COMMAND git rev-parse --path-format=absolute --git-path ${KUTTA_GIT_REF} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE KUTTA_GIT_REF_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
      set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${KUTTA_GIT_REF_PATH}")
    endif()
  endif()
  target_compile_definitions(kutta-cli PRIVATE KUTTA_VERSION="${PROJECT_VERSION}" KUTTA_COMMIT="${KUTTA_COMMIT}")
endif()

file(GLOB KUTTA_BENCH_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/bench/*.cpp)
foreach(src ${KUTTA_BENCH_SOURCES})
  get_filename_component(name ${src} NAME_WE)
  # The memory probe reads POSIX rusage and backtraces; it measures Mac and Linux.
  if(WIN32 AND name STREQUAL "bench_memory")
    continue()
  endif()
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
  target_compile_definitions(kutta-tests PRIVATE KUTTA_PYTHON_PATH="${Python3_EXECUTABLE}")
  target_compile_definitions(kutta-tests PRIVATE KUTTA_CLI_PATH="$<TARGET_FILE:kutta-cli>")
  target_compile_definitions(kutta-tests PRIVATE KUTTA_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  # Long numerical verification (doctest suite "slow"): `ctest -LE slow` skips it.
  add_test(NAME kutta-tests COMMAND kutta-tests --test-suite-exclude=slow)
  add_test(NAME kutta-tests-slow COMMAND kutta-tests --test-suite=slow)
  set_tests_properties(kutta-tests-slow PROPERTIES LABELS slow)
  add_test(NAME lint-portability COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/lint_portability.py)
endif()

# The Qt 6 desktop app (app/). Off by default: it needs Qt and never links the
# solver; it runs the kutta CLI as a child process (docs/app-protocol.md).
option(KUTTA_BUILD_APP "Build the Qt desktop app" OFF)
if(KUTTA_BUILD_APP)
  add_subdirectory(app)
endif()
