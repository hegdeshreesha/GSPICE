cmake_minimum_required(VERSION 3.10)
project(GSPICE VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Find and enable OpenMP
find_package(OpenMP REQUIRED)

# Optional SuiteSparse/KLU backend. GSPICE still builds without it and falls
# back to the internal sparse engine.
find_path(SUITESPARSE_KLU_INCLUDE_DIR
    NAMES klu.h
    PATH_SUFFIXES suitesparse SuiteSparse include/suitesparse include/SuiteSparse
)
find_library(SUITESPARSE_KLU_LIBRARY NAMES klu libklu)
find_library(SUITESPARSE_AMD_LIBRARY NAMES amd libamd)
find_library(SUITESPARSE_COLAMD_LIBRARY NAMES colamd libcolamd)
find_library(SUITESPARSE_BTF_LIBRARY NAMES btf libbtf)
find_library(SUITESPARSE_CONFIG_LIBRARY NAMES suitesparseconfig SuiteSparseConfig libsuitesparseconfig)

set(GSPICE_HAVE_SUITESPARSE_KLU OFF)
if(SUITESPARSE_KLU_INCLUDE_DIR AND SUITESPARSE_KLU_LIBRARY AND
   SUITESPARSE_AMD_LIBRARY AND SUITESPARSE_COLAMD_LIBRARY AND
   SUITESPARSE_BTF_LIBRARY AND SUITESPARSE_CONFIG_LIBRARY)
    set(GSPICE_HAVE_SUITESPARSE_KLU ON)
endif()

# Core source files
set(GSPICE_CORE_SOURCES
    src/core/main.cpp
    src/core/parser.cpp
)

add_executable(gspice ${GSPICE_CORE_SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(gspice PUBLIC OpenMP::OpenMP_CXX)
endif()

if(GSPICE_HAVE_SUITESPARSE_KLU)
    target_compile_definitions(gspice PRIVATE GSPICE_HAVE_SUITESPARSE_KLU=1)
    target_include_directories(gspice PRIVATE ${SUITESPARSE_KLU_INCLUDE_DIR})
    target_link_libraries(gspice PRIVATE
        ${SUITESPARSE_KLU_LIBRARY}
        ${SUITESPARSE_BTF_LIBRARY}
        ${SUITESPARSE_AMD_LIBRARY}
        ${SUITESPARSE_COLAMD_LIBRARY}
        ${SUITESPARSE_CONFIG_LIBRARY}
    )
    message(STATUS "SuiteSparse/KLU backend enabled: ${SUITESPARSE_KLU_LIBRARY}")
else()
    target_compile_definitions(gspice PRIVATE GSPICE_HAVE_SUITESPARSE_KLU=0)
    message(STATUS "SuiteSparse/KLU backend not found; using internal sparse solver")
endif()

target_include_directories(gspice PRIVATE include)

enable_testing()

function(add_gspice_smoke_test test_name deck)
    add_test(
        NAME ${test_name}
        COMMAND gspice --threads 1 ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
    )
    set_tests_properties(${test_name} PROPERTIES
        PASS_REGULAR_EXPRESSION "Simulation Completed Successfully"
    )
endfunction()

function(add_gspice_error_test test_name deck expected_regex)
    add_test(
        NAME ${test_name}
        COMMAND powershell -NoProfile -ExecutionPolicy Bypass
            -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_expect_error.ps1
            $<TARGET_FILE:gspice>
            ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
            ${expected_regex}
    )
endfunction()

function(add_gspice_regex_test test_name deck expected_regex)
    add_test(
        NAME ${test_name}
        COMMAND gspice --threads 1 ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
    )
    set_tests_properties(${test_name} PROPERTIES
        PASS_REGULAR_EXPRESSION "${expected_regex}"
    )
endfunction()

add_gspice_smoke_test(smoke_rc_pulse_v tests/decks/rc_pulse_v.sp)
add_gspice_smoke_test(smoke_rc_pulse_i tests/decks/rc_pulse_i.sp)
add_gspice_smoke_test(smoke_options tests/decks/options.sp)
add_gspice_regex_test(smoke_solver_singleton tests/decks/solver_singleton.sp "singleton_elim=[1-9]")
if(GSPICE_HAVE_SUITESPARSE_KLU)
    add_gspice_regex_test(smoke_solver_klu tests/decks/solver_klu.sp "ext_klu_calls=[1-9]")
endif()
add_gspice_regex_test(smoke_tran_maxstep tests/decks/tran_maxstep.sp "max internal step=2\\.000000000e-11")
add_gspice_regex_test(smoke_tran_method_options tests/decks/tran_method_options.sp "method=Gear2/BDF2")
add_gspice_smoke_test(smoke_include tests/decks/include.sp)
add_gspice_smoke_test(smoke_lib_section tests/decks/lib_section.sp)
add_gspice_smoke_test(smoke_subckt_param tests/decks/subckt_param.sp)
add_gspice_smoke_test(smoke_ihp_lv_mos_wrapper tests/decks/ihp_lv_mos_wrapper.sp)
add_gspice_smoke_test(smoke_osdi_builtin tests/decks/osdi_builtin.sp)
add_gspice_error_test(error_missing_psp_osdi tests/decks/osdi_missing_psp.sp "OSDI model type 'psp103va' is not loaded")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/osdi/psp103.osdi")
    add_gspice_regex_test(smoke_osdi_psp_pmos_pullup tests/decks/osdi_psp_pmos_pullup.sp "Node 2=1\\.")
    add_gspice_regex_test(smoke_osdi_psp_inverter tests/decks/osdi_psp_inverter.sp "1\\.100000000e-06 \\| 2\\.000000000e\\+00 0\\.000000000e\\+00 1\\.")
    add_gspice_regex_test(smoke_osdi_psp_inverter_loaded tests/decks/osdi_psp_inverter_loaded.sp "e-11 \\| 2\\.000000000e\\+00")
endif()
