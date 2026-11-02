# Screenshot smoke test (cmake -P): render scenarios through the platform's GPU
# API and check each PNG has the expected size and is not blank. On GPU-less
# Windows runners Qt's D3D11 backend falls back to WARP.
#   -DAPP=<Kutta executable> -DOUT=<directory> [-DSCENARIOS=a;b] [-DCLI=<kutta>]
# Without CLI, scenarios use the build tree's fixture CLI when it exists.
cmake_minimum_required(VERSION 3.26)
# The recorded k-epsilon, laminar 2D and transition cases run too: each has
# its own equations, approximations and views.
if(NOT SCENARIOS)
  set(SCENARIOS start elbow-setup car-results param-run cylinder-run wing-aero new-case-car new-case-car-rotating new-case-car-boundary new-case-run new-case-save)
endif()
# Every scenario also checks its inspector page at each panel width and that
# no visible text shows the home folder (main.cpp); the saved case gets a long
# name, as users' cases have, which must wrap rather than widen the page.
set(ENV{KUTTA_SCENARIO_CASE_NAME} "fsae_front_wing_ground_effect_study_rev12_final")
file(REMOVE_RECURSE "${OUT}")
file(MAKE_DIRECTORY "${OUT}")
set(cli_args)
if(CLI)
  set(cli_args --cli "${CLI}")
endif()
foreach(s ${SCENARIOS})
  foreach(theme dark light)
    set(png "${OUT}/${s}-${theme}.png")
    execute_process(COMMAND "${APP}" --data-dir "${OUT}/data-${s}-${theme}" --scenario ${s} --theme ${theme}
                            --size 1200x800 --screenshot "${png}" ${cli_args}
                    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err TIMEOUT 300)
    if(NOT rc EQUAL 0 OR NOT EXISTS "${png}")
      message(FATAL_ERROR "${s} (${theme}) failed (${rc}):\n${out}${err}")
    endif()
    # The app exits 3 when a viewport showing data drew nothing; the window
    # itself must also be a real image rather than an empty grab.
    file(SIZE "${png}" bytes)
    if(bytes LESS 40000)
      message(FATAL_ERROR "${s} (${theme}) looks blank: ${bytes} bytes")
    endif()
    string(STRIP "${out}" out)
    message(STATUS "${s} (${theme}): ${bytes} bytes; ${out}")
  endforeach()
endforeach()
