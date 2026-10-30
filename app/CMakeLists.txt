# Kutta desktop app (Qt 6). Built from the root project with
# -DKUTTA_BUILD_APP=ON, or standalone: cmake -S app -B build-app.
# The solver is never linked in: the app runs the `kutta` CLI as a child
# process (docs/app-protocol.md).
cmake_minimum_required(VERSION 3.26)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(kutta-app VERSION 0.1.0 LANGUAGES CXX)
  set(CMAKE_CXX_STANDARD 20)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
  endif()
  set(KUTTA_APP_STANDALONE ON)
endif()

# QRhi (viewports) is semi-public API behind GuiPrivate; the app ships with the
# exact Qt it was built against, so the version tie is intended.
set(QT_NO_PRIVATE_MODULE_WARNING ON)
find_package(Qt6 6.8 REQUIRED COMPONENTS Core Gui GuiPrivate Widgets Svg ShaderTools Test)
qt_standard_project_setup()
set(CMAKE_AUTORCC ON)

# ---- Docking: Qt Advanced Docking System (shared library) ---------------------
include(FetchContent)
set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_STATIC OFF CACHE BOOL "" FORCE)
set(ADS_VERSION 5.1.1 CACHE STRING "" FORCE)
set(QT_VERSION_MAJOR 6 CACHE STRING "" FORCE)
FetchContent_Declare(qtads
  URL https://github.com/githubuser0xFFFF/Qt-Advanced-Docking-System/archive/refs/tags/v5.1.1.tar.gz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(qtads)
# Keep the library next to the executables so dev runs and tests find it.
set_target_properties(qtadvanceddocking-qt6 PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(MSVC)
  set(KUTTA_APP_WARNINGS /W4 /permissive-)
else()
  set(KUTTA_APP_WARNINGS -Wall -Wextra -Wpedantic)
endif()

# ---- Core: protocol, models and maths with no widgets (unit tested) ------------
file(GLOB KUTTA_APP_CORE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/core/*.cpp)
add_library(kutta_app_core STATIC ${KUTTA_APP_CORE_SOURCES})
target_include_directories(kutta_app_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(kutta_app_core PUBLIC Qt6::Core Qt6::Gui)
target_compile_options(kutta_app_core PRIVATE ${KUTTA_APP_WARNINGS})
set_target_properties(kutta_app_core PROPERTIES AUTOMOC ON)

# ---- Fixture CLI: speaks the protocol from recorded and synthetic data --------
file(GLOB KUTTA_FIXTURE_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tools/kutta-fixture/*.cpp)
qt_add_executable(kutta-fixture ${KUTTA_FIXTURE_SOURCES})
set_target_properties(kutta-fixture PROPERTIES MACOSX_BUNDLE OFF WIN32_EXECUTABLE OFF)
target_link_libraries(kutta-fixture PRIVATE kutta_app_core)
target_compile_options(kutta-fixture PRIVATE ${KUTTA_APP_WARNINGS})
target_compile_definitions(kutta-fixture PRIVATE
  KUTTA_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")

# ---- The app ----------------------------------------------------------------------
file(GLOB KUTTA_APP_UI_SOURCES CONFIGURE_DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/src/ui/*.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/render/*.cpp)
qt_add_executable(kutta-app ${KUTTA_APP_UI_SOURCES} src/main.cpp resources/app.qrc)
set_target_properties(kutta-app PROPERTIES
  MACOSX_BUNDLE ON
  WIN32_EXECUTABLE ON
  OUTPUT_NAME Kutta
  MACOSX_BUNDLE_BUNDLE_NAME Kutta
  MACOSX_BUNDLE_GUI_IDENTIFIER org.kutta.app
  MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
  MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION}
  MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/Info.plist.in)
target_link_libraries(kutta-app PRIVATE kutta_app_core Qt6::Widgets Qt6::Svg Qt6::GuiPrivate
  qtadvanceddocking-qt6)
target_compile_options(kutta-app PRIVATE ${KUTTA_APP_WARNINGS})
qt_add_shaders(kutta-app "kutta_shaders"
  PREFIX "/shaders"
  BASE ${CMAKE_CURRENT_SOURCE_DIR}/src/render/shaders
  FILES
    src/render/shaders/surface.vert src/render/shaders/surface.frag
    src/render/shaders/pick.vert src/render/shaders/pick.frag
    src/render/shaders/lines.vert src/render/shaders/lines.frag
    src/render/shaders/ribbon.vert
    src/render/shaders/glyph.vert)
# Dev runs find the fixture CLI without configuration; packaged builds look
# next to the executable first (see CliLocator).
target_compile_definitions(kutta-app PRIVATE
  KUTTA_APP_VERSION="${PROJECT_VERSION}"
  KUTTA_DEV_FIXTURE_CLI="$<TARGET_FILE:kutta-fixture>"
  KUTTA_DEV_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
if(TARGET kutta-cli)
  target_compile_definitions(kutta-app PRIVATE KUTTA_DEV_REAL_CLI="$<TARGET_FILE:kutta-cli>")
endif()
add_dependencies(kutta-app kutta-fixture)

# ---- Tests ----------------------------------------------------------------------
enable_testing()
file(GLOB KUTTA_APP_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/tst_*.cpp)
foreach(src ${KUTTA_APP_TESTS})
  get_filename_component(name ${src} NAME_WE)
  qt_add_executable(${name} ${src})
  set_target_properties(${name} PROPERTIES MACOSX_BUNDLE OFF WIN32_EXECUTABLE OFF)
  target_link_libraries(${name} PRIVATE kutta_app_core Qt6::Test)
  target_compile_definitions(${name} PRIVATE
    KUTTA_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures"
    KUTTA_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data"
    KUTTA_FIXTURE_CLI="$<TARGET_FILE:kutta-fixture>")
  add_dependencies(${name} kutta-fixture)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/tests)
  add_test(NAME app-${name} COMMAND ${name})
  set_tests_properties(app-${name} PROPERTIES LABELS app ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
endforeach()

# Screenshot smoke test: renders scenarios offscreen-in-a-window with the
# platform's GPU API (Metal / D3D; WARP on GPU-less Windows runners).
add_test(NAME app-screenshot-smoke
  COMMAND ${CMAKE_COMMAND} -DAPP=$<TARGET_FILE:kutta-app> -DOUT=${CMAKE_CURRENT_BINARY_DIR}/smoke
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/screenshot_smoke.cmake)
set_tests_properties(app-screenshot-smoke PROPERTIES LABELS "app;gui")

include(cmake/deploy.cmake)
