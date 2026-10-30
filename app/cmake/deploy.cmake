# Packaging: `cmake --build <build> --target package-app` writes
#   macOS:   Kutta-<version>-macos.zip    Kutta.app; the CLI in Contents/Helpers
#   Windows: Kutta-<version>-windows.zip  Kutta\Kutta.exe, Qt; the CLI in solver\
# The solver comes from KUTTA_APP_SOLVER_DIR (a folder holding kutta[.exe] and
# any libraries beside it, such as the windows-cuda CI package) or from the
# in-tree kutta-cli target. Without either, the fixture CLI and its recorded
# and synthetic data are bundled, and every view says FIXTURE or SYNTHETIC.
set(KUTTA_APP_SOLVER_DIR "" CACHE PATH "Folder containing a built kutta CLI to bundle")

get_target_property(_kutta_qmake Qt6::qmake IMPORTED_LOCATION)
get_filename_component(_kutta_qt_bin "${_kutta_qmake}" DIRECTORY)
if(APPLE)
  find_program(KUTTA_DEPLOYQT macdeployqt HINTS "${_kutta_qt_bin}" NO_DEFAULT_PATH)
elseif(WIN32)
  find_program(KUTTA_DEPLOYQT windeployqt HINTS "${_kutta_qt_bin}" NO_DEFAULT_PATH)
endif()

set(_kutta_solver "${KUTTA_APP_SOLVER_DIR}")
if(NOT _kutta_solver AND TARGET kutta-cli)
  set(_kutta_solver "$<TARGET_FILE:kutta-cli>")
endif()

if(KUTTA_DEPLOYQT)
  add_custom_target(package-app
    COMMAND ${CMAKE_COMMAND}
      -DAPP=$<TARGET_FILE:kutta-app>
      -DADS=$<TARGET_FILE:qtadvanceddocking-qt6>
      -DFIXTURE=$<TARGET_FILE:kutta-fixture>
      -DFIXTURE_DATA=${CMAKE_CURRENT_SOURCE_DIR}/fixtures
      -DSOLVER=${_kutta_solver}
      -DDEPLOYQT=${KUTTA_DEPLOYQT}
      -DLICENSES=${CMAKE_CURRENT_SOURCE_DIR}/../docs/app-licensing.md|${CMAKE_CURRENT_SOURCE_DIR}/resources/icons/LICENSE-lucide.txt|${qtads_SOURCE_DIR}/LICENSE|${qtads_SOURCE_DIR}/gnu-lgpl-v2.1.md
      -DVERSION=${PROJECT_VERSION}
      -DOUT=${CMAKE_CURRENT_BINARY_DIR}/package
      -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/package.cmake
    DEPENDS kutta-app kutta-fixture
    COMMENT "Packaging Kutta"
    VERBATIM)
endif()
