cmake_minimum_required(VERSION 3.10)
project(GSPICE VERSION 1.3.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include directories
include_directories(include)

# Find and enable OpenMP
find_package(OpenMP REQUIRED)

# SuiteSparse/KLU backend. Development builds may use the internal sparse
# fallback, while release presets can set GSPICE_REQUIRE_KLU=ON to prevent a
# missing or incomplete KLU installation from being silently accepted.
option(GSPICE_ENABLE_KLU "Enable the SuiteSparse/KLU sparse backend when available" ON)
option(GSPICE_REQUIRE_KLU "Fail configuration unless SuiteSparse/KLU is available" OFF)

if(GSPICE_REQUIRE_KLU)
    set(GSPICE_ENABLE_KLU ON)
endif()

set(GSPICE_HAVE_SUITESPARSE_KLU OFF)
set(GSPICE_KLU_LINK_TARGET "")

if(GSPICE_ENABLE_KLU)
    # SuiteSparse and vcpkg publish a transitive imported target containing
    # KLU plus AMD, COLAMD, BTF, and SuiteSparseConfig. Prefer it over
    # reconstructing the dependency graph by hand.
    find_package(KLU CONFIG QUIET)
    if(TARGET SuiteSparse::KLU)
        set(GSPICE_HAVE_SUITESPARSE_KLU ON)
        set(GSPICE_KLU_LINK_TARGET SuiteSparse::KLU)
    else()
        # Portable fallback for distributions that install only headers and
        # libraries without a KLUConfig.cmake package.
        find_path(SUITESPARSE_KLU_INCLUDE_DIR
            NAMES klu.h
            PATH_SUFFIXES suitesparse SuiteSparse include/suitesparse include/SuiteSparse
        )
        find_library(SUITESPARSE_KLU_LIBRARY NAMES klu libklu)
        find_library(SUITESPARSE_AMD_LIBRARY NAMES amd libamd)
        find_library(SUITESPARSE_COLAMD_LIBRARY NAMES colamd libcolamd)
        find_library(SUITESPARSE_BTF_LIBRARY NAMES btf libbtf)
        find_library(SUITESPARSE_CONFIG_LIBRARY NAMES suitesparseconfig SuiteSparseConfig libsuitesparseconfig)

        if(SUITESPARSE_KLU_INCLUDE_DIR AND SUITESPARSE_KLU_LIBRARY AND
           SUITESPARSE_AMD_LIBRARY AND SUITESPARSE_COLAMD_LIBRARY AND
           SUITESPARSE_BTF_LIBRARY AND SUITESPARSE_CONFIG_LIBRARY)
            set(GSPICE_HAVE_SUITESPARSE_KLU ON)
        endif()
    endif()
endif()

if(GSPICE_REQUIRE_KLU AND NOT GSPICE_HAVE_SUITESPARSE_KLU)
    message(FATAL_ERROR
        "GSPICE_REQUIRE_KLU=ON, but SuiteSparse/KLU was not found. "
        "Configure with the vcpkg toolchain and install the suitesparse-klu manifest dependency.")
endif()

# Core source files
set(GSPICE_CORE_SOURCES
    src/core/main.cpp
    src/core/parser.cpp
)

add_executable(gspice ${GSPICE_CORE_SOURCES})
target_compile_definitions(gspice PRIVATE GSPICE_VERSION="${PROJECT_VERSION}-academic-beta")

if(OpenMP_CXX_FOUND)
    target_link_libraries(gspice PUBLIC OpenMP::OpenMP_CXX)
endif()

if(GSPICE_HAVE_SUITESPARSE_KLU)
    target_compile_definitions(gspice PRIVATE GSPICE_HAVE_SUITESPARSE_KLU=1)
    if(GSPICE_KLU_LINK_TARGET)
        target_link_libraries(gspice PRIVATE ${GSPICE_KLU_LINK_TARGET})
        message(STATUS "SuiteSparse/KLU backend enabled through ${GSPICE_KLU_LINK_TARGET}")
    else()
        target_include_directories(gspice PRIVATE ${SUITESPARSE_KLU_INCLUDE_DIR})
        target_link_libraries(gspice PRIVATE
            ${SUITESPARSE_KLU_LIBRARY}
            ${SUITESPARSE_BTF_LIBRARY}
            ${SUITESPARSE_AMD_LIBRARY}
            ${SUITESPARSE_COLAMD_LIBRARY}
            ${SUITESPARSE_CONFIG_LIBRARY}
        )
        message(STATUS "SuiteSparse/KLU backend enabled: ${SUITESPARSE_KLU_LIBRARY}")
    endif()
else()
    target_compile_definitions(gspice PRIVATE GSPICE_HAVE_SUITESPARSE_KLU=0)
    if(GSPICE_ENABLE_KLU)
        message(STATUS "SuiteSparse/KLU backend not found; using internal sparse solver")
    else()
        message(STATUS "SuiteSparse/KLU backend disabled; using internal sparse solver")
    endif()
endif()

target_include_directories(gspice PRIVATE include)

enable_testing()

add_executable(gspice_dae_state_test tests/dae_state_test.cpp)
target_include_directories(gspice_dae_state_test PRIVATE include)
if(OpenMP_CXX_FOUND)
    target_link_libraries(gspice_dae_state_test PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME unit_dae_state COMMAND gspice_dae_state_test)

find_program(GSPICE_POWERSHELL NAMES pwsh powershell)
find_program(GSPICE_PYTHON NAMES python3 python)
find_program(GSPICE_XYCE NAMES Xyce xyce)
find_program(GSPICE_VACASK NAMES vacask)
if(NOT GSPICE_POWERSHELL)
    message(WARNING "PowerShell was not found; file/error regression helpers will be skipped")
endif()

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
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_expect_error.ps1
                $<TARGET_FILE:gspice>
                ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                ${expected_regex}
        )
    endif()
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

function(add_gspice_env_regex_test test_name deck env_assignment expected_regex)
    add_test(
        NAME ${test_name}
        COMMAND gspice --threads 1 ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
    )
    set_tests_properties(${test_name} PROPERTIES
        ENVIRONMENT "${env_assignment}"
        PASS_REGULAR_EXPRESSION "${expected_regex}"
    )
endfunction()

function(add_gspice_raw_expect_test test_name deck signal expected tolerance)
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_raw_expect.ps1
                $<TARGET_FILE:gspice>
                ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                ${signal}
                ${expected}
                ${tolerance}
        )
    endif()
endfunction()

function(add_gspice_raw_header_expect_test test_name deck expected_variables expected_signal absent_signal)
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_raw_header_expect.ps1
                -Exe $<TARGET_FILE:gspice>
                -Deck ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                -ExpectedVariables ${expected_variables}
                -ExpectedSignals ${expected_signal}
                -AbsentSignals ${absent_signal}
        )
    endif()
endfunction()

function(add_gspice_csv_expect_test test_name deck)
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_csv_expect.ps1
                -Exe $<TARGET_FILE:gspice>
                -Deck ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                -Output ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.csv
        )
    endif()
endfunction()

add_gspice_smoke_test(smoke_rc_pulse_v tests/decks/rc_pulse_v.sp)
add_gspice_smoke_test(smoke_rc_pulse_i tests/decks/rc_pulse_i.sp)
add_gspice_smoke_test(smoke_options tests/decks/options.sp)
add_gspice_regex_test(smoke_options_policy tests/decks/options_policy.sp "Numerical policy:.*source_stepping=on.*gmin_stepping=on.*line_search=on.*residual_check=on")
add_gspice_regex_test(smoke_ic_temp_global tests/decks/ic_temp_global.sp "Initial conditions: applied 1 node voltage seed")
add_gspice_regex_test(smoke_solver_singleton tests/decks/solver_singleton.sp "singleton_elim=[1-9]")
if(GSPICE_HAVE_SUITESPARSE_KLU)
    add_gspice_regex_test(smoke_solver_klu tests/decks/solver_klu.sp "ext_klu_calls=[1-9]")
endif()
add_gspice_regex_test(smoke_tran_maxstep tests/decks/tran_maxstep.sp "max internal step=2\\.000000000e-11")
add_gspice_regex_test(smoke_tran_method_options tests/decks/tran_method_options.sp "method=Gear2/BDF2")
add_gspice_regex_test(smoke_tran_bdf5_options tests/decks/tran_bdf5_options.sp "method=Variable-step BDF.*max_order_used=5")
add_gspice_regex_test(smoke_tran_adams5_options tests/decks/tran_adams5_options.sp "method=Variable-step Adams-Moulton.*max_order_used=5")
add_gspice_regex_test(smoke_tran_predictor_options tests/decks/tran_predictor_options.sp "predictor=on")
add_gspice_regex_test(smoke_tran_predictor_lte tests/decks/rc_pulse_v.sp "pc_lte_steps=[1-9]")
add_gspice_regex_test(smoke_tran_adaptive_order tests/decks/tran_adaptive_order.sp "max_order_used=[3-5]")
add_gspice_regex_test(smoke_nodeset_solver tests/decks/nodeset_solver.sp "Nodesets: applied 1 temporary nonlinear constraint")
add_gspice_regex_test(smoke_dae_audit tests/decks/dae_audit.sp "DAE audit: checked=[1-9].*PASS")
add_gspice_smoke_test(smoke_include tests/decks/include.sp)
add_gspice_smoke_test(smoke_lib_section tests/decks/lib_section.sp)
add_gspice_smoke_test(smoke_subckt_param tests/decks/subckt_param.sp)
add_gspice_regex_test(smoke_param_expression tests/decks/param_expression.sp "Node 1=1\\.333333V \\[out\\]")
add_gspice_regex_test(smoke_subckt_conditional tests/decks/subckt_conditional.sp "Node 1=0\\.500000V \\[out\\]")
add_gspice_smoke_test(smoke_ihp_lv_mos_wrapper tests/decks/ihp_lv_mos_wrapper.sp)
set_tests_properties(smoke_ihp_lv_mos_wrapper PROPERTIES
    ENVIRONMENT "GSPICE_ALLOW_PRIMITIVE_IHP_FALLBACK=1"
)
add_gspice_smoke_test(smoke_osdi_builtin tests/decks/osdi_builtin.sp)
add_gspice_regex_test(smoke_dc_sweep_model tests/decks/dc_sweep_model.sp "Starting DC Sweep Analysis")
add_gspice_regex_test(smoke_dc_nested tests/decks/dc_nested.sp "sweep\\(V1\\).*sweep\\(V2\\)")
add_gspice_regex_test(smoke_controlled_sources tests/decks/controlled_sources.sp "Node 1=1\\.000000V \\[out\\].*Node 2=2\\.000000V \\[buf\\]")
add_gspice_smoke_test(smoke_current_controlled_sources tests/decks/current_controlled_sources.sp)
add_gspice_regex_test(smoke_tf_divider tests/decks/tf_divider.sp "magnitude=4\\.99999999[0-9]e-01")
add_gspice_smoke_test(smoke_bjt_npn tests/decks/bjt_npn.sp)
add_gspice_regex_test(smoke_noise_resistor tests/decks/noise_resistor.sp "Starting Noise Analysis")
add_gspice_regex_test(smoke_measure_tran tests/decks/measure_tran.sp "MEASURE vfinal = 6\\.[0-9]+e-01")
add_gspice_regex_test(smoke_step_op tests/decks/step_op.sp "Starting STEP Analysis")
add_gspice_regex_test(smoke_sens_source tests/decks/sens_source.sp "Sensitivity: dV\\(out\\)/dV1")
add_gspice_regex_test(smoke_pz_rc tests/decks/pz_rc.sp "Dominant pole estimate")
add_gspice_regex_test(smoke_mc_source tests/decks/mc_source.sp "MC summary")
add_gspice_regex_test(smoke_mc_uniform_lhs tests/decks/mc_uniform_lhs.sp "MC sampling: distribution=uniform sampling=latin-hypercube seed=17")
add_gspice_regex_test(smoke_behavioral_current tests/decks/behavioral_current.sp "Node 1=0\\.500000V \\[out\\]")
add_gspice_regex_test(smoke_behavioral_voltage tests/decks/behavioral_voltage.sp "Node 1=1\\.400000V \\[out\\]")
add_gspice_regex_test(smoke_corners_spec tests/decks/corners_spec.sp "corner\\(low\\).*SPEC\\(vout\\)=PASS")
add_gspice_regex_test(smoke_mc_spec tests/decks/mc_spec.sp "Production summary: spec=vout")
add_gspice_raw_expect_test(regression_rc_step tests/decks/regression_rc_step.sp "V(out)" 0.630 0.03)
add_gspice_raw_header_expect_test(regression_save_voltage_selection tests/decks/save_voltage_selection.sp 2 "V(out)" "V(in)")
add_gspice_csv_expect_test(regression_csv_output tests/decks/regression_rc_step.sp)
add_gspice_raw_expect_test(regression_lumen_ihp_inverter tests/decks/lumen_ihp_inverter_short.sp "V(net2)" 0.0 0.2)
set_tests_properties(regression_lumen_ihp_inverter PROPERTIES
    ENVIRONMENT "GSPICE_OSDI_DIR=${CMAKE_CURRENT_SOURCE_DIR}/osdi"
)
if(GSPICE_PYTHON AND GSPICE_XYCE)
    add_test(
        NAME reference_transient_xyce
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_transient.py
            --gspice $<TARGET_FILE:gspice>
            --xyce ${GSPICE_XYCE}
            --source ${CMAKE_CURRENT_SOURCE_DIR}
    )
    set_tests_properties(reference_transient_xyce PROPERTIES TIMEOUT 180)
endif()
if(GSPICE_PYTHON AND GSPICE_VACASK)
    if(GSPICE_XYCE)
        add_test(
            NAME reference_transient_vacask
            COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_transient.py
                --gspice $<TARGET_FILE:gspice>
                --xyce ${GSPICE_XYCE}
                --vacask ${GSPICE_VACASK}
                --source ${CMAKE_CURRENT_SOURCE_DIR}
        )
    else()
        add_test(
            NAME reference_transient_vacask
            COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_transient.py
                --gspice $<TARGET_FILE:gspice>
                --vacask ${GSPICE_VACASK}
                --source ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()
    set_tests_properties(reference_transient_vacask PROPERTIES TIMEOUT 180)
endif()
add_gspice_error_test(error_missing_psp_osdi tests/decks/osdi_missing_psp.sp "OSDI model type 'psp103va' is not loaded")
add_gspice_error_test(error_unsupported_element tests/decks/error_unsupported_element.sp "unsupported element")
add_gspice_error_test(error_ihp_primitive_fallback tests/decks/error_ihp_primitive_fallback.sp "Refusing primitive MOS fallback")
add_gspice_error_test(error_compact_mos_fallback tests/decks/error_compact_mos_fallback.sp "Refusing primitive Level-1 fallback")
add_gspice_error_test(error_pss_unimplemented tests/decks/error_pss_unimplemented.sp "has no validated execution engine")
add_gspice_error_test(error_trtol_min_step tests/decks/error_trtol_min_step.sp "Transient progress starved after 1024 LTE violations")

add_test(NAME cli_version COMMAND gspice --version)
set_tests_properties(cli_version PROPERTIES PASS_REGULAR_EXPRESSION "GSPICE 1\\.3\\.0-academic-beta")
add_test(NAME cli_capabilities COMMAND gspice --capabilities)
set_tests_properties(cli_capabilities PROPERTIES PASS_REGULAR_EXPRESSION "\"pss\": \"unsupported\"")
if(GSPICE_HAVE_SUITESPARSE_KLU)
    add_test(NAME cli_capabilities_klu COMMAND gspice --capabilities)
    set_tests_properties(cli_capabilities_klu PROPERTIES
        PASS_REGULAR_EXPRESSION "\"sparse_backend\": \"SuiteSparse-KLU\""
    )
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/osdi/psp103.osdi")
    add_gspice_regex_test(smoke_osdi_model_status_loaded tests/decks/osdi_psp_inverter.sp "MODEL_STATUS: OSDI_LOADED")
    add_gspice_regex_test(smoke_osdi_model_status_device tests/decks/osdi_psp_inverter.sp "MODEL_STATUS: OSDI_DEVICE")
    add_gspice_env_regex_test(smoke_osdi_limiting_rhs_option tests/decks/osdi_psp_inverter.sp "GSPICE_USE_OSDI_LIMITING_RHS=1" "OSDI advanced options: limiting_rhs=requested")
    add_gspice_regex_test(smoke_osdi_options_deck tests/decks/osdi_psp_options.sp "OSDI advanced options: limiting_rhs=requested.*full_model_params=off.*spice_rhs=requested")
    add_gspice_regex_test(smoke_osdi_internal_nodes tests/decks/osdi_internal_nodes.sp "internal_nodes=expanded internal_unknowns=[1-9]")
    add_gspice_regex_test(smoke_osdi_internal_solution tests/decks/osdi_internal_nodes.sp "Node 2=1\\.")
    add_gspice_regex_test(smoke_osdi_adaptive_homotopy tests/decks/osdi_internal_nodes.sp "adaptive continuation step\\(s\\)")
    add_gspice_regex_test(smoke_osdi_psp_pmos_pullup tests/decks/osdi_psp_pmos_pullup.sp "Node 2=1\\.")
    add_gspice_regex_test(smoke_osdi_psp_dae_audit tests/decks/osdi_psp_dae_audit.sp "DAE audit: checked=[1-9].*PASS")
    add_gspice_regex_test(smoke_osdi_psp_inverter tests/decks/osdi_psp_inverter.sp "1\\.100000000e-06 \\| 2\\.000000000e\\+00 0\\.000000000e\\+00 1\\.")
    add_gspice_regex_test(smoke_osdi_psp_inverter_loaded tests/decks/osdi_psp_inverter_loaded.sp "e-11 \\| 2\\.000000000e\\+00")
endif()
