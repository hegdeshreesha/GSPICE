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
option(GSPICE_ENABLE_NOISE_VALIDATION "Enable external CSV-based stochastic noise validation tests" OFF)
option(GSPICE_ENABLE_RF_NOISE_REGRESSION "Enable experimental RF/noise/PSS regression tests" OFF)
option(GSPICE_ENABLE_EXTERNAL_ORACLE_VALIDATION "Enable optional EXTERNAL_ORACLE oracle parity tests" OFF)
option(GSPICE_ENABLE_GMC_GENERATED "Build generated GMC Verilog-A model smoke targets" ON)
find_program(GSPICE_NGSPICE ngspice ngspice.exe)
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

if(GSPICE_HAVE_SUITESPARSE_KLU)
    target_compile_definitions(gspice PRIVATE GSPICE_HAVE_SUITESPARSE_KLU=1)
    if(SUITESPARSE_KLU_INCLUDE_DIR)
        target_include_directories(gspice PRIVATE ${SUITESPARSE_KLU_INCLUDE_DIR})
    endif()
    if(GSPICE_KLU_LINK_TARGET)
        target_link_libraries(gspice PRIVATE ${GSPICE_KLU_LINK_TARGET})
    else()
        target_link_libraries(gspice PRIVATE
            ${SUITESPARSE_KLU_LIBRARY}
            ${SUITESPARSE_AMD_LIBRARY}
            ${SUITESPARSE_COLAMD_LIBRARY}
            ${SUITESPARSE_BTF_LIBRARY}
            ${SUITESPARSE_CONFIG_LIBRARY}
        )
    endif()
endif()

add_executable(psp103_parameters_test tests/psp103_parameters_test.cpp)
add_executable(gmc_dual_test tests/gmc_dual_test.cpp)
add_executable(psp103_core_test tests/psp103_core_test.cpp)
add_executable(psp103_model_dae_test tests/psp103_model_dae_test.cpp)
add_executable(psp103_temperature_test tests/psp103_temperature_test.cpp)
add_executable(psp103_self_heating_test tests/psp103_self_heating_test.cpp)
add_executable(psp103_leakage_test tests/psp103_leakage_test.cpp)
add_executable(psp103_gsdi_descriptor_test tests/psp103_gsdi_descriptor_test.cpp)
add_executable(psp103_noise_test tests/psp103_noise_test.cpp)
add_executable(psp103_idvg_probe tests/psp103_idvg_probe.cpp)
add_executable(psp103_charges_probe tests/psp103_charges_probe.cpp)
add_executable(psp103_internal_dump tests/psp103_internal_dump.cpp)
add_executable(juncap_express_test tests/juncap_express_test.cpp)
add_executable(gmc_juncap_express_test tests/gmc_juncap_express_test.cpp)
add_executable(gsdi_gdi_adapter_test tests/gsdi_gdi_adapter_test.cpp)
add_executable(gsdi_metadata_test tests/gsdi_metadata_test.cpp)
add_executable(gsdi_collapse_test tests/gsdi_collapse_test.cpp)
add_executable(juncap2_core_test tests/juncap2_core_test.cpp)
add_executable(feature_registry_test tests/feature_registry_test.cpp)
add_executable(bsim3_dual_jacobian_test tests/bsim3_dual_jacobian_test.cpp)
add_executable(bsim3_charge_conservation_test tests/bsim3_charge_conservation_test.cpp)
add_executable(bsim3_parameters_test tests/bsim3_parameters_test.cpp)
add_executable(bsim_translator_test tests/bsim_translator_test.cpp)
add_executable(bsim3_dc_core_test tests/bsim3_dc_core_test.cpp)
add_executable(bsim3_dynamic_stamp_test tests/bsim3_dynamic_stamp_test.cpp)
add_executable(bsim3_noise_test tests/bsim3_noise_test.cpp)
add_executable(bsim3_reference_probe tests/bsim3_reference_probe.cpp)
add_executable(bsim4_parameters_test tests/bsim4_parameters_test.cpp)
add_executable(bsim4_dc_core_test tests/bsim4_dc_core_test.cpp)
add_executable(bsim4_reference_probe tests/bsim4_reference_probe.cpp)
add_executable(bsim4_diagnostics_probe tests/bsim4_diagnostics_probe.cpp)
add_executable(bsim4_abulk_test tests/bsim4_abulk_test.cpp)
add_executable(bsim4_model_dae_test tests/bsim4_model_dae_test.cpp)
add_executable(bsim4_leakage_noise_test tests/bsim4_leakage_noise_test.cpp)
add_executable(primitive_model_param_probe tests/primitive_model_param_probe.cpp)

if(OpenMP_CXX_FOUND)
    target_link_libraries(gspice PUBLIC OpenMP::OpenMP_CXX)
endif()
if(CMAKE_DL_LIBS)
    target_link_libraries(gspice PRIVATE ${CMAKE_DL_LIBS})
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
add_test(NAME psp103_parameters COMMAND psp103_parameters_test)
add_test(NAME gmc_dual COMMAND gmc_dual_test)
add_test(NAME psp103_core COMMAND psp103_core_test)
add_test(NAME psp103_model_dae COMMAND psp103_model_dae_test)
add_test(NAME psp103_temperature COMMAND psp103_temperature_test)
add_test(NAME psp103_self_heating COMMAND psp103_self_heating_test)
add_test(NAME psp103_leakage COMMAND psp103_leakage_test)
add_test(NAME psp103_gsdi_descriptor COMMAND psp103_gsdi_descriptor_test)
add_test(NAME psp103_noise COMMAND psp103_noise_test)
add_test(NAME psp103_idvg_parity COMMAND psp103_idvg_probe)
add_test(NAME psp103_charges_parity COMMAND psp103_charges_probe)
add_test(NAME psp103_internal_dump COMMAND psp103_internal_dump)
add_test(NAME juncap_express COMMAND juncap_express_test)
add_test(NAME gmc_juncap_express COMMAND gmc_juncap_express_test)
add_test(NAME gsdi_gdi_adapter COMMAND gsdi_gdi_adapter_test)
add_test(NAME gsdi_metadata COMMAND gsdi_metadata_test)
add_test(NAME gsdi_collapse COMMAND gsdi_collapse_test)
add_test(NAME juncap2_core COMMAND juncap2_core_test)
add_test(NAME feature_registry COMMAND feature_registry_test)
add_test(NAME bsim3_dual_jacobian COMMAND bsim3_dual_jacobian_test)
add_test(NAME bsim3_charge_conservation COMMAND bsim3_charge_conservation_test)
add_test(NAME bsim3_parameters COMMAND bsim3_parameters_test)
add_test(NAME bsim_translator COMMAND bsim_translator_test)
add_test(NAME bsim3_dc_core COMMAND bsim3_dc_core_test)
add_test(NAME bsim3_dynamic_stamp COMMAND bsim3_dynamic_stamp_test)
add_test(NAME bsim3_noise COMMAND bsim3_noise_test)
add_test(NAME bsim4_parameters COMMAND bsim4_parameters_test)
add_test(NAME bsim4_dc_core COMMAND bsim4_dc_core_test)
add_test(NAME bsim4_abulk COMMAND bsim4_abulk_test)
add_test(NAME bsim4_model_dae COMMAND bsim4_model_dae_test)
add_test(NAME bsim4_leakage_noise COMMAND bsim4_leakage_noise_test)
add_test(NAME primitive_model_param_probe COMMAND primitive_model_param_probe)

find_program(GSPICE_POWERSHELL NAMES pwsh powershell)
find_program(GSPICE_PYTHON NAMES python3 python)
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

function(add_gspice_no_warning_test test_name deck expected_regex)
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_no_warning.ps1
                -Exe $<TARGET_FILE:gspice>
                -Deck ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                -ExpectedRegex ${expected_regex}
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

function(add_gspice_raw_bounds_test test_name deck signal min_value max_value)
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_raw_bounds.ps1
                -Exe $<TARGET_FILE:gspice>
                -Deck ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                -Signal ${signal}
                -Min ${min_value}
                -Max ${max_value}
        )
    endif()
endfunction()

function(add_gspice_raw_spacing_test test_name deck max_min_dt min_points)
    if(GSPICE_POWERSHELL)
        add_test(
            NAME ${test_name}
            COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
                -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_raw_spacing.ps1
                -Exe $<TARGET_FILE:gspice>
                -Deck ${CMAKE_CURRENT_SOURCE_DIR}/${deck}
                -MaxMinDt ${max_min_dt}
                -MinPoints ${min_points}
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
add_gspice_regex_test(smoke_solver_singleton tests/decks/solver_singleton.sp "singleton_elim=[0-9]")
if(GSPICE_HAVE_SUITESPARSE_KLU)
    add_gspice_regex_test(smoke_solver_klu tests/decks/solver_klu.sp "ext_klu_calls=[1-9]")
endif()
add_gspice_regex_test(smoke_tran_maxstep tests/decks/tran_maxstep.sp "max internal step=2\\.000000000e-11.*adaptive=on")
add_gspice_regex_test(smoke_tran_maxstep_applied tests/decks/tran_maxstep.sp "Transient summary:.*max_step=2\\.000000000e-11")
add_gspice_regex_test(smoke_tran_maxstep_auto tests/decks/tran_maxstep_auto.sp "max internal step=1\\.000000000e-09.*adaptive=on")
add_gspice_regex_test(smoke_tran_stamp_cache_default tests/decks/tran_maxstep_auto.sp "stamp_cache=on")
add_gspice_regex_test(smoke_tran_accuracy_high_preset tests/decks/tran_accuracy_high_preset.sp "lte_reference=HISTORY.*lte_audit_interval=0.*trtol=3\\.500000000e\\+00")
add_gspice_regex_test(smoke_tran_method_options tests/decks/tran_method_options.sp "method=Gear2/BDF2")
add_gspice_regex_test(smoke_tran_bdf5_options tests/decks/tran_bdf5_options.sp "method=Variable-step BDF.*max_order_used=[2-5]")
add_gspice_regex_test(smoke_tran_adams5_options tests/decks/tran_adams5_options.sp "method=Variable-step Adams-Moulton.*max_order_used=[2-5]")
add_gspice_regex_test(smoke_tran_predictor_options tests/decks/tran_predictor_options.sp "predictor=on")
add_gspice_regex_test(smoke_tran_predictor_lte tests/decks/rc_pulse_v.sp "lte_rejected=[1-9]")
add_gspice_regex_test(smoke_tran_adaptive_order tests/decks/tran_adaptive_order.sp "max_order_used=[3-5]")
add_gspice_regex_test(smoke_nodeset_solver tests/decks/nodeset_solver.sp "Nodesets: applied 1 temporary nonlinear constraint")
add_gspice_regex_test(smoke_dae_audit tests/decks/dae_audit.sp "DAE audit: checked=[1-9].*PASS")
add_gspice_regex_test(smoke_gmc_primitive_mos tests/decks/dae_audit.sp "MODEL_STATUS: GMC primitive MOS: M1 type=NMOS")
if(GSPICE_POWERSHELL)
    add_test(
        NAME psp103_idvg_e2e_parity
        COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
            -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_psp103_idvg_parity.ps1
            $<TARGET_FILE:gspice>
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/psp103_idvg_e2e.sp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/psp103_idvg_ref.txt
    )
    add_test(
        NAME psp103_pmos_idvg_e2e_parity
        COMMAND ${GSPICE_POWERSHELL} -NoProfile -ExecutionPolicy Bypass
            -File ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_psp103_pmos_idvg_parity.ps1
            $<TARGET_FILE:gspice>
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/psp103_pmos_idvg_e2e.sp
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/reference/psp103_pmos_idvg_ref.txt
    )
endif()
add_gspice_regex_test(smoke_gmc_juncap_express_parser tests/decks/juncap_express_parser.sp "MODEL_STATUS: GMC JUNCAP Express: D1")
add_gspice_regex_test(smoke_gmc_juncap2_parser tests/decks/juncap2_parser.sp "MODEL_STATUS: GMC JUNCAP2: D1")
add_gspice_regex_test(smoke_gmc_juncap2_tat_avalanche_parser tests/decks/juncap2_tat_avalanche_parser.sp "MODEL_STATUS: GMC JUNCAP2: D1")
add_gspice_regex_test(smoke_bsim3_gsdi_factory tests/decks/bsim3_reference_dc.sp "MODEL_STATUS: GMC BSIM3: M1")
add_gspice_regex_test(smoke_bsim4_gsdi_factory tests/decks/bsim4_reference_dc.sp "MODEL_STATUS: GMC BSIM4: M1")
add_gspice_no_warning_test(no_warning_bsim3_gsdi tests/decks/bsim3_reference_dc.sp "MODEL_STATUS: GMC BSIM3: M1")
add_gspice_no_warning_test(no_warning_bsim4_gsdi tests/decks/bsim4_reference_dc.sp "MODEL_STATUS: GMC BSIM4: M1")
add_gspice_smoke_test(smoke_include tests/decks/include.sp)
add_gspice_smoke_test(smoke_lib_section tests/decks/lib_section.sp)
add_gspice_smoke_test(smoke_subckt_param tests/decks/subckt_param.sp)
add_gspice_regex_test(smoke_param_expression tests/decks/param_expression.sp "Node 1=1\\.333333V \\[out\\]")
add_gspice_regex_test(smoke_subckt_conditional tests/decks/subckt_conditional.sp "Node 1=0\\.500000V \\[out\\]")
add_gspice_smoke_test(smoke_ihp_lv_mos_wrapper tests/decks/ihp_lv_mos_wrapper.sp)
add_gspice_regex_test(smoke_ihp_corner_moslv_native tests/decks/ihp_corner_moslv_native.sp "GSDI_BUILTIN: psp103va")
add_gspice_regex_test(smoke_ihp_lv_mos_op_native tests/decks/ihp_lv_mos_op_native.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_ihp_lv_mos_dc_native tests/decks/ihp_lv_mos_dc_native.sp "Starting DC Sweep Analysis")
add_gspice_regex_test(smoke_ihp_lv_mos_ac_native tests/decks/ihp_lv_mos_ac_native.sp "AC summary: points=5 sweep=LIN")
add_gspice_regex_test(smoke_ihp_lv_inverter_tran_native tests/decks/ihp_lv_inverter_tran_native.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_ihp_lv_mos_noise_native tests/decks/ihp_lv_mos_noise_native.sp "Starting Noise Analysis")
add_gspice_regex_test(smoke_ihp_lv_mos_corners_native tests/decks/ihp_lv_mos_corners_native.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_ihp_lv_mos_temp_native tests/decks/ihp_lv_mos_temp_native.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_ihp_lv_mos_geometry_native tests/decks/ihp_lv_mos_geometry_native.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_ihp_lv_mos_mismatch_native tests/decks/ihp_lv_mos_mismatch_native.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_ihp_lv_mos_rf_native tests/decks/ihp_lv_mos_rf_native.sp "AC summary: points=3 sweep=LIN")
add_gspice_regex_test(smoke_ihp_lv_mos_rf_gate_resistance_native tests/decks/ihp_lv_mos_rf_native.sp "MODEL_STATUS: GMC PSP RF gate resistance: M1")
add_gspice_regex_test(smoke_ihp_lv_mos_rf_ac_charge_native tests/decks/ihp_lv_mos_rf_native.sp "1\\.00e\\+06 \\| .*\\([-+0-9.eE]+,[1-9][0-9.]*e-0[0-9]\\)")
add_gspice_no_warning_test(no_warning_ihp_lv_op tests/decks/ihp_lv_mos_op_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_dc tests/decks/ihp_lv_mos_dc_native.sp "Starting DC Sweep Analysis")
add_gspice_no_warning_test(no_warning_ihp_lv_ac tests/decks/ihp_lv_mos_ac_native.sp "AC summary: points=5 sweep=LIN")
add_gspice_no_warning_test(no_warning_ihp_lv_tran tests/decks/ihp_lv_inverter_tran_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_noise tests/decks/ihp_lv_mos_noise_native.sp "Starting Noise Analysis")
add_gspice_no_warning_test(no_warning_lumen_ihp_inverter tests/decks/lumen_ihp_inverter_short.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_corners tests/decks/ihp_lv_mos_corners_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_temp tests/decks/ihp_lv_mos_temp_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_geometry tests/decks/ihp_lv_mos_geometry_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_mismatch tests/decks/ihp_lv_mos_mismatch_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_lv_rf tests/decks/ihp_lv_mos_rf_native.sp "AC summary: points=3 sweep=LIN")
add_gspice_regex_test(smoke_ihp_hv_mos_native tests/decks/ihp_hv_mos_native.sp "GSDI_BUILTIN: psp103va")
add_gspice_regex_test(smoke_ihp_hv_mos_mismatch_native tests/decks/ihp_hv_mos_mismatch_native.sp "GSDI_BUILTIN: psp103va")
add_gspice_env_regex_test(smoke_ihp_hv_svaricap_native tests/decks/ihp_hv_svaricap_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "native voltage-dependent MOSVAR capacitance")
add_gspice_regex_test(smoke_ihp_hv_svaricap_ac_native tests/decks/ihp_hv_svaricap_ac_native.sp "AC summary: points=4 sweep=LIN")
add_gspice_smoke_test(smoke_ihp_tap_resistor_native tests/decks/ihp_tap_resistor_native.sp)
add_gspice_env_regex_test(smoke_ihp_r3cmc_resistor_native tests/decks/ihp_r3cmc_resistor_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "RSH\\*LEFF/WEFF with foundry wrapper width/contact/corner/temperature terms")
add_gspice_env_regex_test(smoke_ihp_r3cmc_rhigh_native tests/decks/ihp_r3cmc_rhigh_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "RSH\\*LEFF/WEFF with foundry wrapper width/contact/corner/temperature terms")
add_gspice_env_regex_test(smoke_ihp_r3cmc_rppd_native tests/decks/ihp_r3cmc_rppd_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "RSH\\*LEFF/WEFF with foundry wrapper width/contact/corner/temperature terms")
add_gspice_smoke_test(smoke_ihp_r3cmc_resistor_temp_native tests/decks/ihp_r3cmc_resistor_temp_native.sp)
add_gspice_smoke_test(smoke_ihp_cparasitic_native tests/decks/ihp_cparasitic_native.sp)
add_gspice_env_regex_test(smoke_ihp_cmim_capacitor_native tests/decks/ihp_cmim_capacitor_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "routed as \\(CJ\\*L\\*W \\+ CJSW\\*perimeter\\)\\*SCALE with TC1/TC2 temperature scaling")
add_gspice_env_regex_test(smoke_ihp_cmim_mismatch_native tests/decks/ihp_cmim_mismatch_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "routed as \\(CJ\\*L\\*W \\+ CJSW\\*perimeter\\)\\*SCALE with TC1/TC2 temperature scaling")
add_gspice_env_regex_test(smoke_ihp_cmim_temp_native tests/decks/ihp_cmim_temp_native.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "TC1/TC2 temperature scaling")
add_gspice_no_warning_test(no_warning_ihp_cmim tests/decks/ihp_cmim_capacitor_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_cmim_mismatch tests/decks/ihp_cmim_mismatch_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_cmim_temp tests/decks/ihp_cmim_temp_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_r3cmc tests/decks/ihp_r3cmc_resistor_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_r3cmc_rhigh tests/decks/ihp_r3cmc_rhigh_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_r3cmc_rppd tests/decks/ihp_r3cmc_rppd_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_r3cmc_temp tests/decks/ihp_r3cmc_resistor_temp_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_hv_svaricap tests/decks/ihp_hv_svaricap_native.sp "Simulation Completed Successfully")
add_gspice_no_warning_test(no_warning_ihp_hv_svaricap_ac tests/decks/ihp_hv_svaricap_ac_native.sp "AC summary: points=4 sweep=LIN")
add_gspice_smoke_test(smoke_ihp_antenna_diode_native tests/decks/ihp_antenna_diode_native.sp)
add_gspice_smoke_test(smoke_ihp_schottky_native tests/decks/ihp_schottky_native.sp)
add_gspice_smoke_test(smoke_ihp_esd_native tests/decks/ihp_esd_native.sp)
add_gspice_smoke_test(smoke_ihp_bondpad_native tests/decks/ihp_bondpad_native.sp)
add_gspice_smoke_test(smoke_ihp_hbt_native tests/decks/ihp_hbt_native.sp)
add_gspice_regex_test(smoke_dc_sweep_model tests/decks/dc_sweep_model.sp "Starting DC Sweep Analysis")
add_gspice_smoke_test(smoke_primitive_model_params tests/decks/primitive_model_params.sp)
if(GSPICE_PYTHON)
    foreach(CURRENT_PROBE_CASE diode mos bjt psp103)
        if(CURRENT_PROBE_CASE STREQUAL "diode")
            set(CURRENT_PROBE_DECK tests/decks/current_probe_diode.sp)
            set(CURRENT_PROBE_SIGNAL "I(D1)")
            set(CURRENT_PROBE_MEASURE id1)
            set(CURRENT_PROBE_MIN 1e-5)
        elseif(CURRENT_PROBE_CASE STREQUAL "mos")
            set(CURRENT_PROBE_DECK tests/decks/current_probe_mos.sp)
            set(CURRENT_PROBE_SIGNAL "I(M1)")
            set(CURRENT_PROBE_MEASURE im1)
            set(CURRENT_PROBE_MIN 1e-7)
        elseif(CURRENT_PROBE_CASE STREQUAL "bjt")
            set(CURRENT_PROBE_DECK tests/decks/current_probe_bjt.sp)
            set(CURRENT_PROBE_SIGNAL "I(Q1)")
            set(CURRENT_PROBE_MEASURE iq1)
            set(CURRENT_PROBE_MIN 1e-5)
        else()
            set(CURRENT_PROBE_DECK tests/decks/current_probe_psp103.sp)
            set(CURRENT_PROBE_SIGNAL "I(M1)")
            set(CURRENT_PROBE_MEASURE im1)
            set(CURRENT_PROBE_MIN 1e-6)
        endif()
        add_test(NAME current_probe_${CURRENT_PROBE_CASE}
            COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_current_probe.py
                --exe $<TARGET_FILE:gspice>
                --deck ${CMAKE_CURRENT_SOURCE_DIR}/${CURRENT_PROBE_DECK}
                --signal ${CURRENT_PROBE_SIGNAL}
                --measure ${CURRENT_PROBE_MEASURE}
                --min ${CURRENT_PROBE_MIN})
    endforeach()
endif()
add_gspice_regex_test(smoke_dc_nested tests/decks/dc_nested.sp "sweep\\(V1\\).*sweep\\(V2\\)")
add_gspice_regex_test(smoke_ac_lin tests/decks/ac_rc.sp "AC summary: points=5 sweep=LIN")
add_gspice_regex_test(smoke_ac_values_phase tests/decks/ac_values_phase.sp "1\\.00e\\+03 \\| \\(6\\.[0-9]+e-17,1\\.00e\\+00\\).*AC summary: points=3 sweep=VALUES")
add_gspice_regex_test(smoke_ac_current_phase tests/decks/ac_current_phase.sp "1\\.00e\\+03 \\| \\(-6\\.[0-9]+e-14,-1\\.00e\\+03\\).*AC summary: points=1 sweep=LIN")
add_gspice_regex_test(smoke_mutual_inductor_ac tests/decks/mutual_inductor_ac.sp "1\\.00e\\+03 \\| .*\\(4\\.[0-9]+e-01,.*AC summary: points=1 sweep=LIN")
add_gspice_regex_test(smoke_tline_grounded_op tests/decks/tline_grounded_op.sp "Node [0-9]+=0\\.500000V \\[out\\]")
add_gspice_regex_test(smoke_tline_floating_ac tests/decks/tline_floating_ac.sp "AC summary: points=1 sweep=LIN")
add_gspice_regex_test(smoke_controlled_sources tests/decks/controlled_sources.sp "Node 1=1\\.000000V \\[out\\].*Node 2=2\\.000000V \\[buf\\]")
add_gspice_smoke_test(smoke_current_controlled_sources tests/decks/current_controlled_sources.sp)
add_gspice_regex_test(smoke_jfet_n_channel tests/decks/jfet_n_channel.sp "Node [0-9]+=3\\.[0-9]+V \\[drain\\]")
add_gspice_regex_test(smoke_mesfet_alias tests/decks/mesfet_alias.sp "Node [0-9]+=3\\.[0-9]+V \\[drain\\]")
add_gspice_regex_test(smoke_tf_divider tests/decks/tf_divider.sp "magnitude=4\\.99999999[0-9]e-01")
add_gspice_regex_test(smoke_acxf_alias tests/decks/acxf_alias.sp "Starting ACXF Small-Signal Analysis")
add_gspice_regex_test(smoke_dcxf_alias tests/decks/dcxf_alias.sp "Starting DCXF Small-Signal Analysis")
add_gspice_regex_test(smoke_dcinc_alias tests/decks/dcinc_alias.sp "Starting DCINC Small-Signal Analysis")
add_gspice_regex_test(smoke_spectre_errpreset_conservative tests/decks/spectre_errpreset_conservative.sp "Numerical policy: reltol=1\\.000000000e-04")
add_gspice_regex_test(smoke_hbac_alias tests/decks/hbac_alias.sp "Native HB Converged:.*HBAC conversion matrix: sidebands=")
add_gspice_regex_test(smoke_hbnoise_alias tests/decks/hbnoise_alias.sp "Native HB Converged:.*HBNOISE conversion matrix: sidebands=")
add_gspice_regex_test(smoke_hbsp_alias tests/decks/hbsp_alias.sp "Native HB Converged:.*HBSP conversion matrix: sidebands=")
add_gspice_regex_test(smoke_hbstb_alias tests/decks/hbstb_alias.sp "Native HB Converged:.*HBSTB conversion matrix: sidebands=")
add_gspice_regex_test(smoke_hbxf_alias tests/decks/hbxf_alias.sp "Native HB Converged:.*HBAC conversion matrix: sidebands=")
add_gspice_regex_test(smoke_hb_alias tests/decks/hb_alias.sp "Native HB Converged:.*HB Converged: f0=")
add_gspice_regex_test(smoke_hb_multitone_alias tests/decks/hb_multitone_alias.sp "Native HB Converged:.*HB Converged: f0=5\\.000000000e\\+02.*tones=2")
add_gspice_regex_test(smoke_hb_signoff_options tests/decks/hb_signoff_options.sp "Native HB Converged:.*native_hb=yes")
add_gspice_regex_test(smoke_hb_diode_nonlinear tests/decks/hb_diode_nonlinear.sp "Native HB Converged:.*native_hb=yes.*2\\.000000000e\\+03 \\| .*\\| [1-9]\\.[0-9]+e-03")
add_gspice_regex_test(smoke_hb_psp103_nonlinear tests/decks/hb_psp103_nonlinear.sp "Native HB Converged:.*native_hb=yes.*2\\.000000000e\\+03 \\| .*\\| .*\\| [2-3]\\.[0-9]+e-04")
add_gspice_regex_test(smoke_hbac_diode_nonlinear tests/decks/hbac_diode_nonlinear.sp "Native HB Converged:.*HBAC conversion matrix: sidebands=5")
add_gspice_regex_test(smoke_hbnoise_diode_nonlinear tests/decks/hbnoise_diode_nonlinear.sp "Native HB Converged:.*HBNOISE conversion matrix: sidebands=5.*PNOISE summary:")
add_gspice_regex_test(smoke_pxf_alias tests/decks/pxf_alias.sp "Starting PAC from converged PSS operating point")
add_gspice_regex_test(smoke_psssp_alias tests/decks/psssp_alias.sp "PSSSP conversion matrix: sidebands=")
add_gspice_regex_test(smoke_pnoise_phase_jitter tests/decks/pnoise_phase_jitter.sp "phase_noise_dbc_per_hz.*jitter_s_per_rtHz")
if(GSPICE_PYTHON)
    set(RF_SIGNOFF_GOLDEN ${CMAKE_CURRENT_SOURCE_DIR}/tests/rf_signoff_golden.json)
    add_test(NAME signoff_hb_linear_rc
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hb_alias.sp
            --case hb_linear_rc
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_hb_diode_nonlinear
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hb_diode_nonlinear.sp
            --case hb_diode
            --golden ${RF_SIGNOFF_GOLDEN})
    if(GSPICE_NGSPICE)
        add_test(NAME validation_rf_ngspice_hb_diode
            COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_rf_ngspice_hb.py
                --gspice $<TARGET_FILE:gspice>
                --ngspice ${GSPICE_NGSPICE}
                --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hb_diode_nonlinear.sp
                --relative-tolerance 1e-2)
    endif()
    add_test(NAME signoff_hb_multitone
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hb_multitone_alias.sp
            --case hb_multitone
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_hb_psp103_nonlinear
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hb_psp103_nonlinear.sp
            --case hb_psp103
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_hbac_diode_nonlinear
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hbac_diode_nonlinear.sp
            --case hbac_diode
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_hbnoise_diode_nonlinear
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hbnoise_diode_nonlinear.sp
            --case hbnoise_diode
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_pnoise_phase_jitter
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/pnoise_phase_jitter.sp
            --case pnoise_phase_jitter
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_hbsp_matched
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hbsp_alias.sp
            --case hbsp
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_psssp_matched
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/psssp_alias.sp
            --case psssp
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_hbstb_loop
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/hbstb_alias.sp
            --case hbstb
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_pssstb_loop
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/pssstb_rc_feedback.sp
            --case pssstb
            --golden ${RF_SIGNOFF_GOLDEN})
    add_test(NAME signoff_psspac_linear
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tests/validate_rf_signoff.py
            --exe $<TARGET_FILE:gspice>
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/psspac_linear.sp
            --case psspac_linear
            --golden ${RF_SIGNOFF_GOLDEN})
endif()
add_gspice_smoke_test(smoke_bjt_npn tests/decks/bjt_npn.sp)
add_gspice_regex_test(smoke_noise_resistor tests/decks/noise_resistor.sp "Starting Noise Analysis")
add_test(
    NAME regression_pnoise_default_sweep
    COMMAND gspice --threads 1 --output ${CMAKE_CURRENT_BINARY_DIR}/pnoise_default_sweep.raw
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/pnoise_default_sweep.sp
)
set_tests_properties(regression_pnoise_default_sweep PROPERTIES
    PASS_REGULAR_EXPRESSION "PNOISE summary:.*offsets=251"
)
if(GSPICE_ENABLE_RF_NOISE_REGRESSION)
    add_gspice_regex_test(smoke_sp_oneport tests/decks/sp_oneport_matched.sp "1\\.000000000e\\+03 \\| -2\\.[0-9]+e-11 0\\.000000000e\\+00.*SP summary: points=1 ports=1 method=multiport-y-to-s")
    add_gspice_regex_test(smoke_sp_twoport tests/decks/sp_twoport_resistor.sp "1\\.000000000e\\+03 \\| 6\\.66666666[0-9]e-02 0\\.000000000e\\+00 2\\.666666667e-01 0\\.000000000e\\+00 2\\.666666667e-01 0\\.000000000e\\+00 6\\.66666666[0-9]e-02 0\\.000000000e\\+00.*SP summary: points=3 ports=2 method=multiport-y-to-s")
    add_gspice_regex_test(smoke_noise_integrated tests/decks/noise_resistor.sp "Integrated output noise: [1-9].* Vrms")
    add_gspice_regex_test(smoke_noise_input_referred tests/decks/noise_resistor.sp "Noise summary: points=2 output=V\\(out\\) input=V1")
    add_gspice_smoke_test(smoke_transient_noise tests/decks/transient_noise.sp)
    add_gspice_smoke_test(smoke_transient_noise_sde tests/decks/transient_noise_sde.sp)
    add_gspice_regex_test(smoke_stb_return_ratio tests/decks/stb_rc_feedback.sp "STB summary: points=[1-9].*method=twoport-return-ratio")
    add_gspice_regex_test(smoke_pss_linearized tests/decks/pss_linear.sp "PSS summary:.*shooting_iterations=[1-9].*residual=")
    add_gspice_regex_test(smoke_pss_shooting_residual tests/decks/pss_newton_rc.sp "PSS iteration [1-9]+: residual=.*PSS summary:.*shooting_iterations=[1-9]")
    add_gspice_regex_test(smoke_pss_variational_newton tests/decks/pss_newton_long_rc.sp "PSS Newton correction:.*PSS summary:.*newton_iterations=[1-9].*monodromy_columns=[1-9].*variational_steps=[1-9]")
    add_gspice_regex_test(smoke_pss_oscillator_option tests/decks/pss_oscillator_option.sp "Autonomous PSS phase condition:.*PSS summary: oscillator=yes autonomous=yes.*direct_period_columns=[1-9].*phase_unknown=[0-9]")
    add_gspice_regex_test(smoke_pss_use_initial_conditions tests/decks/pss_use_initial_conditions.sp "Initial conditions: applied 1 node voltage seed\\(s\\) with UIC\\.")
    add_gspice_regex_test(smoke_pss_tstab_adaptive tests/decks/pss_tstab_adaptive.sp "PSS summary:.*tstab_seconds=[1-9].*adaptive_orbit=yes")
    add_gspice_error_test(error_pnoise_small_signal tests/decks/pnoise_resistor.sp ".PNOISE requires a preceding .PSS")
    add_gspice_regex_test(smoke_pnoise_phase_jitter_rf tests/decks/pnoise_phase_jitter.sp "PNOISE summary:.*carrier=")
    add_gspice_regex_test(smoke_psspac_sidebands tests/decks/psspac_linear.sp "Starting PAC from converged PSS operating point")
    add_gspice_regex_test(smoke_pssstb_twoport tests/decks/pssstb_rc_feedback.sp "Starting PSSSTB from converged PSS operating point")
    add_gspice_error_test(error_pnoise_sidebands tests/decks/pnoise_sidebands.sp ".PNOISE requires a preceding .PSS")
    add_gspice_regex_test(smoke_psspac_diode_lptv tests/decks/psspac_diode_lptv.sp "Starting PAC from converged PSS operating point")
    add_gspice_error_test(error_pnoise_diode_lptv tests/decks/pnoise_diode_lptv.sp ".PNOISE requires a preceding .PSS")
endif()
add_gspice_regex_test(smoke_measure_tran tests/decks/measure_tran.sp "MEASURE vfinal = 6\\.[0-9]+e-01")
add_gspice_regex_test(smoke_measure_tran_ops tests/decks/measure_tran_ops.sp "MEASURE area = 5\\.000000000e-04.*MEASURE slope = 1\\.000000000e\\+03.*MEASURE slope_at = 1\\.000000000e\\+03")
add_gspice_regex_test(smoke_measure_tran_when tests/decks/measure_tran_when.sp "MEASURE t50 = 5\\.000000000e-04")
add_gspice_regex_test(smoke_measure_tran_delay tests/decks/measure_tran_delay.sp "MEASURE delay = 2\\.000000000e-04")
add_gspice_regex_test(smoke_four_sine tests/decks/four_sine.sp "FOUR V\\(out\\).*fundamental=1\\.000000000e\\+03.*thd_percent=.*harmonic frequency magnitude phase_deg normalized")
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
add_gspice_raw_expect_test(regression_op_raw_divider tests/decks/op_raw_divider.sp "V(out)" 0.900 0.001)
add_gspice_raw_header_expect_test(regression_save_voltage_selection tests/decks/save_voltage_selection.sp 2 "V(out)" "V(in)")
add_gspice_raw_header_expect_test(regression_save_none_options tests/decks/save_none_options.sp 1 "time" "V(in)")
add_gspice_raw_header_expect_test(regression_save_selected_options tests/decks/save_selected_options.sp 2 "V(out)" "V(in)")
add_gspice_raw_header_expect_test(regression_save_all_voltage_only tests/decks/save_all_voltage_only.sp 3 "V(in)" "I(V1)")
add_gspice_raw_header_expect_test(regression_save_explicit_current tests/decks/save_explicit_current.sp 5 "I(V1)" "__NO_SUCH_SIGNAL__")
add_gspice_raw_spacing_test(regression_saveadaptive_raw_spacing tests/decks/tran_saveadaptive.sp 2e-10 50)
add_gspice_csv_expect_test(regression_csv_output tests/decks/regression_rc_step.sp)
add_gspice_raw_expect_test(regression_lumen_ihp_inverter tests/decks/lumen_ihp_inverter_short.sp "V(net2)" 0.0 0.2)
set_tests_properties(regression_lumen_ihp_inverter PROPERTIES
    ENVIRONMENT "GSPICE_ALLOW_PRIMITIVE_IHP_FALLBACK=1"
)
add_gspice_raw_bounds_test(regression_ihp_ring_unloaded_outnet_bounds tests/decks/ihp_ring_unloaded_outnet_tran.sp "V(OUTNET)" -0.05 1.05)
set_tests_properties(regression_ihp_ring_unloaded_outnet_bounds PROPERTIES
    TIMEOUT 120
)
if(GSPICE_PYTHON AND GSPICE_ENABLE_NOISE_VALIDATION)
    add_test(
        NAME validation_noise_rc_ktc
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_noise.py
            --gspice $<TARGET_FILE:gspice>
            --source ${CMAKE_CURRENT_SOURCE_DIR}
            --case rc
    )
    set_tests_properties(validation_noise_rc_ktc PROPERTIES TIMEOUT 180)

    add_test(
        NAME validation_transient_noise_seed
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_noise.py
            --gspice $<TARGET_FILE:gspice>
            --source ${CMAKE_CURRENT_SOURCE_DIR}
            --case seed
    )
    set_tests_properties(validation_transient_noise_seed PROPERTIES TIMEOUT 180)
endif()
if(GSPICE_PYTHON)
    add_test(
        NAME validation_bsim4_reference
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_bsim4_reference.py
            --python ${GSPICE_PYTHON}
            --gspice-probe $<TARGET_FILE:bsim4_reference_probe>
            --source ${CMAKE_CURRENT_SOURCE_DIR}
            --relative-tolerance 0.01
    )
    if(GSPICE_NGSPICE)
        add_test(
            NAME validation_bsim3_ngspice_full_region
            COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/validate_bsim3_reference.py
                --ngspice ${GSPICE_NGSPICE}
                --gspice $<TARGET_FILE:bsim3_reference_probe>
                --source ${CMAKE_CURRENT_SOURCE_DIR}
                --relative-tolerance 0.10
                --absolute-tolerance 1e-12
        )
    endif()
    add_test(
        NAME validation_external_oracle_harness_cases
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/external_oracle_parity.py
            --list-cases
    )
    if(GSPICE_NGSPICE)
        add_test(
            NAME validation_external_oracle_ngspice_diode_dc
            COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/external_oracle_parity.py
                --source ${CMAKE_CURRENT_SOURCE_DIR}
                --gspice $<TARGET_FILE:gspice>
                --case ngspice_diode_dc
                --external_oracle-command "${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/run_ngspice_oracle.py --ngspice ${GSPICE_NGSPICE} --case {case} --deck {deck} --csv {csv} --workdir {workdir}"
                --require-external_oracle
                --relative-tolerance 1e-2
                --absolute-tolerance 1e-6
        )
    endif()
    add_test(
        NAME scanner_psp103_approximate_param_enforced
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/rank_psp103_ignored.py
            --deck ${CMAKE_CURRENT_SOURCE_DIR}/tests/decks/psp103_approximate_param.sp
            --fail-on-approximate
    )
    set_tests_properties(scanner_psp103_approximate_param_enforced PROPERTIES WILL_FAIL TRUE)
    if(GSPICE_ENABLE_EXTERNAL_ORACLE_VALIDATION)
        foreach(EXTERNAL_ORACLE_CASE ihp_lv_dc ihp_lv_ac ihp_lv_tran ihp_lv_noise)
            add_test(
                NAME validation_external_oracle_${EXTERNAL_ORACLE_CASE}
                COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/external_oracle_parity.py
                    --source ${CMAKE_CURRENT_SOURCE_DIR}
                    --gspice $<TARGET_FILE:gspice>
                    --case ${EXTERNAL_ORACLE_CASE}
            )
            set_tests_properties(validation_external_oracle_${EXTERNAL_ORACLE_CASE} PROPERTIES
                SKIP_RETURN_CODE 77
                TIMEOUT 180
            )
        endforeach()
    endif()
endif()
if(GSPICE_PYTHON AND GSPICE_ENABLE_GMC_GENERATED)
    set(GMC_GENERATED_VA_HEADER ${CMAKE_CURRENT_BINARY_DIR}/gmc_va_probe.hpp)
    set(GMC_GENERATED_VA_GSDI ${CMAKE_CURRENT_BINARY_DIR}/gmc_va_probe.gsdi)
    add_custom_command(
        OUTPUT ${GMC_GENERATED_VA_HEADER}
        BYPRODUCTS ${GMC_GENERATED_VA_GSDI}
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_compile.py
            --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/gmc_va_probe.va
            --output ${GMC_GENERATED_VA_HEADER}
        DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/gmc_va_probe.va
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_compile.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_parser.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_cpp_emitter.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_ast.py
    )
    add_executable(gmc_generated_va_test tests/gmc_generated_va_test.cpp ${GMC_GENERATED_VA_HEADER})
    target_include_directories(gmc_generated_va_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME gmc_generated_va COMMAND gmc_generated_va_test)

    set(GMC_GENERATED_HIDDEN_HEADER ${CMAKE_CURRENT_BINARY_DIR}/gmc_hidden.hpp)
    set(GMC_GENERATED_HIDDEN_GSDI ${CMAKE_CURRENT_BINARY_DIR}/gmc_hidden.gsdi)
    add_custom_command(
        OUTPUT ${GMC_GENERATED_HIDDEN_HEADER}
        BYPRODUCTS ${GMC_GENERATED_HIDDEN_GSDI}
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_compile.py
            --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/gmc_hidden.va
            --output ${GMC_GENERATED_HIDDEN_HEADER}
        DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/gmc_hidden.va
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_compile.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_parser.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_cpp_emitter.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_ast.py
    )
    add_executable(gmc_hidden_nodes_test tests/gmc_hidden_nodes_test.cpp ${GMC_GENERATED_HIDDEN_HEADER})
    target_include_directories(gmc_hidden_nodes_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME gmc_hidden_nodes COMMAND gmc_hidden_nodes_test)

    set(GMC_GENERATED_DEEP_HEADER ${CMAKE_CURRENT_BINARY_DIR}/gmc_deep.hpp)
    set(GMC_GENERATED_DEEP_GSDI ${CMAKE_CURRENT_BINARY_DIR}/gmc_deep.gsdi)
    add_custom_command(
        OUTPUT ${GMC_GENERATED_DEEP_HEADER}
        BYPRODUCTS ${GMC_GENERATED_DEEP_GSDI}
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_compile.py
            --input ${CMAKE_CURRENT_SOURCE_DIR}/tests/gmc_deep.va
            --output ${GMC_GENERATED_DEEP_HEADER}
        DEPENDS
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/gmc_deep.va
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_compile.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_parser.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_cpp_emitter.py
            ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/gmc_ast.py
    )
    add_executable(gmc_deep_va_test tests/gmc_deep_va_test.cpp ${GMC_GENERATED_DEEP_HEADER})
    target_include_directories(gmc_deep_va_test PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    add_test(NAME gmc_deep_va COMMAND gmc_deep_va_test)

    # Generated GSDI models are part of the simulator binary: the elaborator
    # registers them in GmcRegistry so N devices can instantiate them.
    target_sources(gspice PRIVATE ${GMC_GENERATED_VA_HEADER} ${GMC_GENERATED_HIDDEN_HEADER})
    target_include_directories(gspice PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_compile_definitions(gspice PRIVATE GSPICE_HAVE_GMC_GENERATED=1)
    add_gspice_regex_test(
        gmc_hidden_internal_elaboration tests/decks/gmc_hidden_internal.sp
        "Node 2=4\\.800000V"
    )

    add_test(
        NAME gmc_veriloga_parser
        COMMAND ${GSPICE_PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/tools/gmc/test_gmc_parser.py
    )
endif()
add_gspice_error_test(error_unsupported_compact_model tests/decks/compact_model_unsupported.sp ".GSDI artifact for model type")
add_gspice_error_test(error_osdi_directive_rejected tests/decks/osdi_directive_rejected.sp "external compact-model runtime loading is disabled")
add_gspice_error_test(error_gsdi_bad_envelope tests/decks/gsdi_bad_envelope.sp "invalid .GSDI artifact envelope")
if(GSPICE_ENABLE_GMC_GENERATED)
    add_gspice_error_test(error_gsdi_terminal_mismatch tests/decks/gsdi_terminal_mismatch.sp "declares 2")
    add_gspice_error_test(error_gsdi_unknown_parameter tests/decks/gsdi_unknown_parameter.sp "not declared by .GSDI model type")
    add_gspice_regex_test(smoke_gsdi_artifact_load tests/decks/gsdi_artifact_load.sp "MODEL_STATUS: GSDI_LOADED: va_probe")
    add_gspice_regex_test(smoke_gsdi_va_probe_instance tests/decks/gsdi_va_probe_instance.sp "MODEL_STATUS: GSDI_LOADED: va_probe")
endif()
add_gspice_regex_test(smoke_psp103_native_op tests/decks/psp103_native_op.sp "MODEL_STATUS: GSDI_LOADED: psp103va")
add_gspice_regex_test(smoke_psp103_native_m_op tests/decks/psp103_native_m_op.sp "MODEL_STATUS: GSDI_LOADED: psp103va")
add_gspice_regex_test(smoke_psp103_generic_no_gsdi_m_op tests/decks/psp103_generic_no_gsdi_m_op.sp "MODEL_STATUS: GSDI_BUILTIN: psp103va")
add_gspice_regex_test(smoke_psp103_generic_alias_no_gsdi_n_op tests/decks/psp103_generic_alias_no_gsdi_n_op.sp "MODEL_STATUS: GSDI_BUILTIN: psp")
add_gspice_regex_test(smoke_psp103_generic_foundry_alias_m_op tests/decks/psp103_generic_foundry_alias_m_op.sp "MODEL_STATUS: GSDI_BUILTIN: psp103e")
add_gspice_regex_test(smoke_psp103_generic_nqs_alias_n_op tests/decks/psp103_generic_nqs_alias_n_op.sp "MODEL_STATUS: GSDI_BUILTIN: pspnqs103e")
add_gspice_env_regex_test(smoke_psp103_unsupported_param_warning tests/decks/psp103_unsupported_param_warning.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "unsupported PSP parameter\\(s\\) ignored by native evaluator: XJ")
add_gspice_env_regex_test(smoke_bsim3_unsupported_param_warning tests/decks/bsim3_unsupported_param_warning.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "unsupported BSIM3 model parameter\\(s\\) ignored by native evaluator: PHANTOMVAR")
add_gspice_env_regex_test(smoke_bsim4_unsupported_param_warning tests/decks/bsim4_unsupported_param_warning.sp "GSPICE_VERBOSE_COMPAT_WARNINGS=1" "unsupported BSIM4 model parameter\\(s\\) ignored by native evaluator: PHANTOMVAR")
add_gspice_regex_test(smoke_psp103_level103_nmos_m_op tests/decks/psp103_level103_nmos_m_op.sp "MODEL_STATUS: GSDI_BUILTIN: PSP103VA")
add_gspice_regex_test(smoke_psp103_level103_pmos_m_op tests/decks/psp103_level103_pmos_m_op.sp "MODEL_STATUS: GSDI_BUILTIN: PSP103VA")
add_gspice_regex_test(smoke_psp103_common_aliases_m_op tests/decks/psp103_common_aliases_m_op.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_psp103_tox_mobility_m_op tests/decks/psp103_tox_mobility_m_op.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_psp103_gate_leakage_op tests/decks/psp103_gate_leakage_op.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_psp103_gidl_op tests/decks/psp103_gidl_op.sp "Simulation Completed Successfully")
add_gspice_regex_test(smoke_psp103_tran tests/decks/psp103_tran_smoke.sp "Simulation Completed Successfully")
if(GSPICE_ENABLE_GMC_GENERATED)
    set_tests_properties(error_gsdi_terminal_mismatch error_gsdi_unknown_parameter smoke_gsdi_artifact_load smoke_gsdi_va_probe_instance PROPERTIES DEPENDS gmc_generated_va)
endif()
add_gspice_error_test(error_unsupported_element tests/decks/error_unsupported_element.sp "unsupported element")
add_gspice_error_test(error_ihp_primitive_fallback tests/decks/error_ihp_primitive_fallback.sp "model card not found")
add_gspice_error_test(error_compact_mos_fallback tests/decks/error_compact_mos_fallback.sp "Refusing primitive Level-1 fallback")
add_gspice_error_test(error_trtol_min_step tests/decks/error_trtol_min_step.sp "Transient progress starved after 1024 LTE violations")

add_test(NAME cli_version COMMAND gspice --version)
set_tests_properties(cli_version PROPERTIES PASS_REGULAR_EXPRESSION "GSPICE 1\\.3\\.0-academic-beta")
add_test(NAME cli_capabilities COMMAND gspice --capabilities)
set_tests_properties(cli_capabilities PROPERTIES PASS_REGULAR_EXPRESSION "\"pss\": \"prototype\"")
add_test(NAME cli_capabilities_psp103_native COMMAND gspice --capabilities)
set_tests_properties(cli_capabilities_psp103_native PROPERTIES PASS_REGULAR_EXPRESSION "\"psp103_native\": \\{\"maturity\": \"tested\", \"available\": true\\}")
add_test(NAME cli_capabilities_psp103_full_reference_gate COMMAND gspice --capabilities)
set_tests_properties(cli_capabilities_psp103_full_reference_gate PROPERTIES PASS_REGULAR_EXPRESSION "\"psp103_full_reference\": \\{\"maturity\": \"prototype\", \"available\": false\\}")
if(GSPICE_HAVE_SUITESPARSE_KLU)
    add_test(NAME cli_capabilities_klu COMMAND gspice --capabilities)
    set_tests_properties(cli_capabilities_klu PROPERTIES
        PASS_REGULAR_EXPRESSION "\"sparse_backend\": \"SuiteSparse-KLU\""
    )
endif()
