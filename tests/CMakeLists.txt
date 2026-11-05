add_test(NAME lamp_help COMMAND LAMP --help)
set_tests_properties(lamp_help PROPERTIES PASS_REGULAR_EXPRESSION "Usage:" TIMEOUT 30)
include(architecture/check.cmake)
add_test(NAME repository_layout COMMAND ${CMAKE_COMMAND}
    "-DPROJECT_ROOT=${PROJECT_SOURCE_DIR}" -P "${CMAKE_CURRENT_SOURCE_DIR}/repository_layout.cmake")
if(TARGET LAMP_hydro_csf_static_benchmark)
    add_test(NAME hydro_csf_static_smoke COMMAND LAMP_hydro_csf_static_benchmark
        --smoke --output "${PROJECT_SOURCE_DIR}/runs/checks/hydro_csf_static_smoke/series")
    add_test(NAME hydro_csf_static_input_contract COMMAND LAMP_hydro_csf_static_benchmark --reject-input)
    add_test(NAME hydro_csf_near_axis_reference_smoke COMMAND LAMP_hydro_csf_static_benchmark
        --reference-smoke --output "${PROJECT_SOURCE_DIR}/runs/checks/hydro_csf_near_axis_reference_smoke/series")
    add_test(NAME hydro_csf_azimuthal_projection_smoke COMMAND LAMP_hydro_csf_static_benchmark
        --azimuthal-smoke --output "${PROJECT_SOURCE_DIR}/runs/checks/hydro_csf_azimuthal_projection_smoke/series")
    add_test(NAME hydro_csf_face_average_projection_smoke COMMAND LAMP_hydro_csf_static_benchmark
        --face-average-smoke --output "${PROJECT_SOURCE_DIR}/runs/checks/hydro_csf_face_average_projection_smoke/series")
    add_test(NAME hydro_csf_cell_projection_smoke COMMAND LAMP_hydro_csf_static_benchmark
        --cell-projection-smoke --output "${PROJECT_SOURCE_DIR}/runs/checks/hydro_csf_cell_projection_smoke/series")
endif()

function(lamp_add_module_contract name library)
    file(GLOB_RECURSE _source "${CMAKE_CURRENT_SOURCE_DIR}/contracts/*/${name}.cpp")
    list(LENGTH _source _count)
    if(NOT _count EQUAL 1)
        message(FATAL_ERROR "Contract ${name} must have exactly one registered source")
    endif()
    add_executable(lamp_${name} ${_source})
    target_include_directories(lamp_${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/support")
    target_link_libraries(lamp_${name} PRIVATE lamp_compile_options ${library})
    lamp_copy_linalg_runtime(lamp_${name})
    add_test(NAME ${name} COMMAND lamp_${name})
endfunction()

# VOF production and retained module reference contracts.
lamp_add_module_contract(quadratic_surface_axis_contract lamp_core)
add_executable(lamp_runtime_registry_contract contracts/runtime/runtime_registry_contract.cpp
    "${PROJECT_SOURCE_DIR}/src/Apps/AppRegistry.cpp")
target_include_directories(lamp_runtime_registry_contract PRIVATE "${PROJECT_SOURCE_DIR}/include/Apps")
# D1b: the existing runtime contract consumes the actual eight-document reader.
target_sources(lamp_runtime_registry_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/src/Apps/IO/CaseDocuments.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/IO/TwoPhaseInputComponents.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/HydrodynamicsInput.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/VofInputComponents.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.cpp")
# D3b: actual passive writer + joint transaction, on the original Hydro/module chain.
target_sources(lamp_runtime_registry_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/src/Apps/IO/HydroOutputComponents.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/HydroVofOutputComponents.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/HydrodynamicsVofRegrid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/HydrodynamicsRun.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/VofCheckpoint.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp")
# D3b compiles the same exact initial-Q implementation in this directory scope.
# Parent-directory per-source properties do not propagate to child targets.
if(MSVC)
    set_source_files_properties("${PROJECT_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp"
        "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/HydrodynamicsRun.cpp"
        "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/VofCheckpoint.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties("${PROJECT_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp"
        "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/HydrodynamicsRun.cpp"
        "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/VofCheckpoint.cpp"
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()
target_link_libraries(lamp_runtime_registry_contract PRIVATE lamp_compile_options lamp_hydrodynamics)
lamp_copy_linalg_runtime(lamp_runtime_registry_contract)
add_test(NAME vof_ehd_input_contract COMMAND lamp_runtime_registry_contract --vof-ehd-input-only
    "${PROJECT_SOURCE_DIR}/runs/checks/vof_ehd_input_contract")
# Charge-FCT scope: this registration file is a necessary dependency exception
# beyond the two existing contract sources; no new contract executable is added.
add_test(NAME vof_charge_advection_input_contract COMMAND lamp_runtime_registry_contract
    --vof-charge-advection-input-only "${PROJECT_SOURCE_DIR}/runs/checks/vof_charge_advection_input_contract")
add_test(NAME vof_charge_regrid_contract COMMAND lamp_runtime_registry_contract
    --vof-charge-regrid-only "${PROJECT_SOURCE_DIR}/runs/checks/vof_charge_regrid_contract")
add_test(NAME runtime_registry_contract COMMAND lamp_runtime_registry_contract)
# Actual public controller + original raw CSG asset and physical boundary IDs.
add_test(NAME vof_csg_initial_app_contract COMMAND lamp_runtime_registry_contract
    --vof-csg-initial-only "D:/Chen/Codex/temp/lamp_electrospray_csg_initial_app_20261019/ctest")
add_test(NAME vof_csg_electric_boundary_contract COMMAND lamp_runtime_registry_contract
    --vof-csg-electric-boundary-only "D:/Chen/Codex/temp/lamp_electrospray_csg_electric_boundary_20261019/ctest")
set_tests_properties(vof_csg_initial_app_contract vof_csg_electric_boundary_contract PROPERTIES TIMEOUT 180)
add_test(NAME vof_pressure_driven_reader_contract COMMAND lamp_runtime_registry_contract
    --vof-pressure-driven-reader-only "${PROJECT_SOURCE_DIR}/runs/checks/vof_pressure_driven_reader_contract")
add_test(NAME vof_pressure_driven_app_contract COMMAND lamp_runtime_registry_contract
    --vof-pressure-driven-app-only "${PROJECT_SOURCE_DIR}/runs/checks/vof_pressure_driven_app_contract")
set_tests_properties(vof_pressure_driven_reader_contract vof_pressure_driven_app_contract PROPERTIES TIMEOUT 180)

add_test(NAME vof_b1_zero_wiring_contract COMMAND lamp_runtime_registry_contract
    --vof-b1-zero-wiring-only "${PROJECT_SOURCE_DIR}/tests/fixtures/vof_b1_zero_wiring/system"
    "${PROJECT_SOURCE_DIR}/runs/checks/vof_b1_zero_wiring_contract")
set_tests_properties(vof_b1_zero_wiring_contract PROPERTIES TIMEOUT 180)
lamp_add_module_contract(interface_hermite_gals_state_contract lamp_interface)
lamp_add_module_contract(geometric_vof_publication_contract lamp_interface)
add_test(NAME geometric_vof_pinned_corner_contract
    COMMAND lamp_geometric_vof_publication_contract --pinned-corner-only)
set_tests_properties(geometric_vof_pinned_corner_contract PROPERTIES TIMEOUT 60)
lamp_add_module_contract(geometric_vof_checkpoint_contract lamp_interface)
if(MSVC)
    target_compile_options(lamp_geometric_vof_checkpoint_contract PRIVATE /fp:strict /GL-)
else()
    target_compile_options(lamp_geometric_vof_checkpoint_contract PRIVATE -fno-fast-math -ffp-contract=off -frounding-math -fno-lto)
endif()
set_target_properties(lamp_geometric_vof_checkpoint_contract PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE UNITY_BUILD OFF)
set_tests_properties(geometric_vof_checkpoint_contract PROPERTIES
    TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
# Selected fluid domains retain dense native publication rows and exclude solids.
add_test(NAME geometric_vof_selected_fluid_domain_contract
    COMMAND lamp_geometric_vof_publication_contract --selected-domain-only)
set_tests_properties(geometric_vof_selected_fluid_domain_contract PROPERTIES TIMEOUT 60)
lamp_add_module_contract(geometric_vof_mesh_transfer_contract lamp_interface)
add_test(NAME geometric_vof_selected_domain_amr_contract
    COMMAND lamp_geometric_vof_mesh_transfer_contract --selected-domain-amr-only)
set_tests_properties(geometric_vof_selected_domain_amr_contract PROPERTIES TIMEOUT 120)
lamp_add_module_contract(geometric_vof_transport_contract lamp_hydrodynamics)
add_test(NAME geometric_vof_pinned_corner_transport_contract
    COMMAND lamp_geometric_vof_transport_contract --pinned-corner-transport-only)
set_tests_properties(geometric_vof_pinned_corner_transport_contract PROPERTIES TIMEOUT 60)
add_test(NAME geometric_vof_selected_fluid_transport_contract
    COMMAND lamp_geometric_vof_transport_contract --selected-fluid-transport-only)
set_tests_properties(geometric_vof_selected_fluid_transport_contract PROPERTIES TIMEOUT 60)
add_test(NAME geometric_vof_selected_fluid_hydro_contract
    COMMAND lamp_geometric_vof_transport_contract --selected-fluid-hydro-only)
set_tests_properties(geometric_vof_selected_fluid_hydro_contract PROPERTIES TIMEOUT 120)
add_test(NAME geometric_vof_selected_coupled_step_contract
    COMMAND lamp_geometric_vof_transport_contract --selected-fluid-ehd-step-only)
set_tests_properties(geometric_vof_selected_coupled_step_contract PROPERTIES TIMEOUT 120)
# Conservative VOF momentum: independent P0/Coons continuity, actual stored
# endpoint masses and no-uR closed/open load and energy identities.
add_test(NAME geometric_vof_conservative_momentum_contract
    COMMAND lamp_geometric_vof_transport_contract --local-plic-conservative-momentum-only)
set_tests_properties(geometric_vof_conservative_momentum_contract PROPERTIES TIMEOUT 120)

add_test(NAME geometric_vof_transport_local_smooth_contract
    COMMAND lamp_geometric_vof_transport_contract --local-smooth-only)
lamp_add_module_contract(scalar_expression_contract lamp_core)
lamp_add_module_contract(tecplot_csv_history_contract lamp_core)
lamp_add_module_contract(planar_cell_cut_contract lamp_core)
add_test(NAME planar_corner_inverse_contract
    COMMAND lamp_planar_cell_cut_contract --corner-inverse-only)
set_tests_properties(planar_corner_inverse_contract PROPERTIES TIMEOUT 60)
lamp_add_module_contract(mesh_adaptation_contract lamp_core)
lamp_add_module_contract(interface_hermite_gals_trace_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_vertex_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_edge_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_mixed_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_amr_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_volume_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_flux_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_reflux_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_publication_contract lamp_interface)
lamp_add_module_contract(interface_temporal_flux_interval_contract lamp_interface)
lamp_add_module_contract(interface_polynomial_velocity_contract lamp_interface)
lamp_add_module_contract(interface_characteristic_enclosure_contract lamp_interface)
lamp_add_module_contract(interface_characteristic_phase_contract lamp_interface)
lamp_add_module_contract(interface_bounded_spacetime_flux_contract lamp_interface)
lamp_add_module_contract(interface_planar_temporal_flux_contract lamp_interface)
lamp_add_module_contract(interface_integrated_inventory_contract lamp_interface)
lamp_add_module_contract(interface_integrated_reflux_contract lamp_interface)
lamp_add_module_contract(interface_phase_region_contract lamp_interface)
lamp_add_module_contract(interface_inventory_region_contract lamp_interface)
lamp_add_module_contract(interface_inventory_sweep_contract lamp_interface)
lamp_add_module_contract(interface_curved_sweep_contract lamp_interface)
lamp_add_module_contract(interface_constant_inflow_contract lamp_interface)
lamp_add_module_contract(interface_exact_refinement_contract lamp_interface)
lamp_add_module_contract(interface_refinement_inventory_contract lamp_interface)
lamp_add_module_contract(interface_hermite_gals_skeleton_contract lamp_core)
lamp_add_module_contract(interface_affine_characteristic_contract lamp_core)
lamp_add_module_contract(interface_affine_phase_contract lamp_core)
lamp_add_module_contract(interface_shifted_green_contract lamp_core)
lamp_add_module_contract(interface_root_enclosure_contract lamp_core)
target_sources(lamp_interface_hermite_gals_mixed_contract PRIVATE contracts/interface/interface_gals_contract.cpp)
add_executable(lamp_first_jet_contract contracts/interface/first_jet_contract.cpp)
target_include_directories(lamp_first_jet_contract PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/support")
target_link_libraries(lamp_first_jet_contract PRIVATE lamp_compile_options lamp_interface)
lamp_copy_linalg_runtime(lamp_first_jet_contract)
foreach(_mode publication step mesh-transfer crossing)
    add_test(NAME interface_first_jet_${_mode} COMMAND lamp_first_jet_contract ${_mode})
endforeach()
# H02a exercises the actual HYDRO library, with no legacy Hydro/Interface link.
lamp_add_module_contract(hydro_dof_layout_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_integral_divergence_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_pressure_load_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_constant_normal_stress_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_capillary_traction_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_single_phase_saddle_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_physical_divergence_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_energy_viscosity_contract lamp_hydrodynamics)
add_test(NAME hydro_stationary_no_slip_strain_contract COMMAND lamp_hydro_energy_viscosity_contract --no-slip-strain-only)
set_tests_properties(hydro_stationary_no_slip_strain_contract PROPERTIES TIMEOUT 60)
lamp_add_module_contract(hydro_phase_geometry_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_regularized_color_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_balanced_csf_load_contract lamp_hydrodynamics)
add_test(NAME hydro_pinned_corner_hf_contract
    COMMAND lamp_hydro_balanced_csf_load_contract --pinned-corner-hf-only)
set_tests_properties(hydro_pinned_corner_hf_contract PROPERTIES TIMEOUT 120)

lamp_add_module_contract(hydro_electric_interface_load_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_sharp_maxwell_load_contract lamp_hydrodynamics)
target_link_libraries(lamp_hydro_electric_interface_load_contract PRIVATE lamp_electric)
add_test(NAME hydro_balanced_csf_domains_contract
    COMMAND lamp_hydro_balanced_csf_load_contract --domains)
lamp_add_module_contract(hydro_csf_curvature_contract lamp_hydrodynamics)
add_test(NAME hydro_csf_subarc_selected_contract
    COMMAND lamp_hydro_csf_curvature_contract --extra-provenance)
lamp_add_module_contract(hydro_phase_viscosity_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_material_transport_contract lamp_hydrodynamics)
add_test(NAME hydro_spacetime_material_endpoints_contract
    COMMAND lamp_hydro_material_transport_contract --spacetime-endpoints)
foreach(_rz 0 1)
    foreach(_amr 0 1)
        add_test(NAME hydro_spacetime_material_${_rz}_${_amr}_contract
            COMMAND lamp_hydro_material_transport_contract --spacetime ${_rz} ${_amr})
    endforeach()
    add_test(NAME hydro_spacetime_material_amr_crossing_${_rz}_contract
        COMMAND lamp_hydro_material_transport_contract --spacetime-amr-crossing ${_rz})
endforeach()
lamp_add_module_contract(hydro_momentum_geometry_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_compatible_source_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_characteristic_contract lamp_hydrodynamics)
target_compile_definitions(lamp_hydro_characteristic_contract PRIVATE
    LAMP_CHARACTERISTIC_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures/hydrodynamics/characteristic/")
lamp_add_module_contract(hydro_dual_material_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_kep_convection_contract lamp_hydrodynamics)
lamp_add_module_contract(hydro_two_phase_step_contract lamp_hydrodynamics)
add_test(NAME hydro_vof_checkpoint_contract COMMAND lamp_hydro_two_phase_step_contract --vof-checkpoint-only)
set_tests_properties(hydro_vof_checkpoint_contract PROPERTIES
    TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=2;OPENBLAS_NUM_THREADS=1")
add_test(NAME vof_pressure_driven_reservoir_contract
    COMMAND lamp_hydro_two_phase_step_contract --vof-pressure-driven-only)
set_tests_properties(vof_pressure_driven_reservoir_contract PROPERTIES TIMEOUT 120)
add_test(NAME vof_stationary_no_slip_pressure_contract COMMAND lamp_hydro_two_phase_step_contract --vof-no-slip-pressure-only)
set_tests_properties(vof_stationary_no_slip_pressure_contract PROPERTIES TIMEOUT 120)
add_test(NAME hydro_vof_transport_provider_contract
    COMMAND lamp_hydro_two_phase_step_contract --vof-transport-provider-only)
lamp_add_module_contract(hydro_vof_refinement_contract lamp_hydrodynamics)
target_sources(lamp_hydro_vof_refinement_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/HydrodynamicsVofRegrid.cpp")
target_include_directories(lamp_hydro_vof_refinement_contract PRIVATE "${PROJECT_SOURCE_DIR}/include/Apps")
lamp_add_module_contract(hydro_vof_mesh_transfer_contract lamp_hydrodynamics)
add_test(NAME hydro_selected_pressure_amr_contract COMMAND lamp_hydro_vof_mesh_transfer_contract --selected-domain-amr-pressure-only)
set_tests_properties(hydro_selected_pressure_amr_contract PROPERTIES TIMEOUT 120)
add_test(NAME hydro_stationary_no_slip_amr_contract COMMAND lamp_hydro_vof_mesh_transfer_contract --selected-domain-no-slip-amr-only)
set_tests_properties(hydro_stationary_no_slip_amr_contract PROPERTIES TIMEOUT 120)
# D2: consume the actual passive VOF-EHD writer in the existing coupled contract.
target_include_directories(lamp_hydro_vof_mesh_transfer_contract PRIVATE "${PROJECT_SOURCE_DIR}/include/Apps")
target_sources(lamp_hydro_vof_mesh_transfer_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Hydrodynamics/IO/HydroVofOutputComponents.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/IO/HydroOutputComponents.cpp")
add_test(NAME vof_ehd_output_contract COMMAND lamp_hydro_vof_mesh_transfer_contract --vof-ehd-output-only
    "${PROJECT_SOURCE_DIR}/runs/checks/vof_ehd_output_contract")
lamp_add_module_contract(electric_solver_contract lamp_electric)
add_test(NAME electric_vof_selected_charge_amr_contract COMMAND lamp_electric_solver_contract --electric-vof-selected-amr-only)
set_tests_properties(electric_vof_selected_charge_amr_contract PROPERTIES TIMEOUT 60)
lamp_add_module_contract(electric_vof_refinement_contract lamp_electric)
add_test(NAME electric_vof_selected_fluid_domain_contract
    COMMAND lamp_electric_solver_contract --electric-vof-selected-domain-only)
set_tests_properties(electric_vof_selected_fluid_domain_contract PROPERTIES TIMEOUT 60)
add_test(NAME electric_vof_selected_chain_contract
    COMMAND lamp_electric_solver_contract --electric-vof-selected-chain-only)
set_tests_properties(electric_vof_selected_chain_contract PROPERTIES TIMEOUT 60)
add_test(NAME electric_vof_charge_second_order_contract
    COMMAND lamp_electric_solver_contract --electric-vof-charge-second-order-only)
add_test(NAME electric_vof_charge_remap_second_order_contract
    COMMAND lamp_electric_solver_contract --electric-vof-charge-remap-second-order-only)
set_tests_properties(electric_vof_charge_remap_second_order_contract PROPERTIES TIMEOUT 60)
lamp_add_module_contract(vof_smooth_geometry_contract lamp_interface)
lamp_add_module_contract(vof_sharp_field_contract lamp_electric)
lamp_add_module_contract(vof_sharp_maxwell_contract lamp_electric)
# Keep the complete original Q1 source available for explicit diagnostics.
# Only these supported submodes are registered; the full mode's historical
# hairpin geometry failure is retained and is not counted as a successful test.
add_executable(lamp_electric_surface_charge_transport_contract
    contracts/electric/electric_surface_charge_transport_contract.cpp)
target_include_directories(lamp_electric_surface_charge_transport_contract PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/support")
target_link_libraries(lamp_electric_surface_charge_transport_contract PRIVATE lamp_compile_options lamp_electric)
lamp_copy_linalg_runtime(lamp_electric_surface_charge_transport_contract)
add_test(NAME electric_surface_charge_transport_stationary_contract
    COMMAND lamp_electric_surface_charge_transport_contract --stationary-only)
add_test(NAME electric_surface_charge_transport_closed_projection_contract
    COMMAND lamp_electric_surface_charge_transport_contract --closed-projection-cost)
add_test(NAME electric_surface_charge_projection_certificate_contract
    COMMAND lamp_electric_surface_charge_transport_contract --projection-certificate-only)
# B1 tests share strict IEEE semantics with the actual production operators.
foreach(_b1_contract vof_smooth_geometry_contract vof_sharp_field_contract vof_sharp_maxwell_contract electric_surface_charge_transport_contract hydro_sharp_maxwell_load_contract)
    if(MSVC)
        target_compile_options(lamp_${_b1_contract} PRIVATE /fp:strict /GL-)
    else()
        target_compile_options(lamp_${_b1_contract} PRIVATE -fno-fast-math -ffp-contract=off -frounding-math -fno-lto)
    endif()
    set_target_properties(lamp_${_b1_contract} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE UNITY_BUILD OFF)
endforeach()
unset(_b1_contract)
lamp_add_module_contract(dictionary_index_contract lamp_core)
add_test(NAME csg_fitted_solid_contract COMMAND lamp_dictionary_index_contract --csg-fitted-solid-only)
set_tests_properties(csg_fitted_solid_contract PROPERTIES TIMEOUT 60 WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}")
lamp_add_module_contract(pic_collision_contract lamp_plasma)
lamp_add_module_contract(plasma_particle_contract lamp_plasma)
lamp_add_module_contract(movingeb_input_contract lamp_plasma)
lamp_add_module_contract(grid_motion_contract lamp_plasma)
lamp_add_module_contract(pic_checkpoint_contract lamp_plasma)
target_compile_definitions(lamp_pic_checkpoint_contract PRIVATE
    LAMP_MESH_COEFFICIENT_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures/mesh/legacy_coefficients/")
lamp_add_module_contract(pic_field_output_contract lamp_plasma)
lamp_add_module_contract(pic_step_contract lamp_plasma)
lamp_add_module_contract(pic_app_config_contract lamp_plasma)
target_sources(lamp_pic_app_config_contract PRIVATE
    "${PROJECT_SOURCE_DIR}/src/Apps/Pic/PicApp.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Pic/Setup/PicRunSetup.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Pic/IO/PicInput.cpp"
    "${PROJECT_SOURCE_DIR}/src/Apps/Pic/IO/PicOutput.cpp")
target_include_directories(lamp_pic_app_config_contract PRIVATE "${PROJECT_SOURCE_DIR}/include/Apps")
add_executable(lamp_superlu_contract contracts/solvers/superlu_contract.cpp)
target_link_libraries(lamp_superlu_contract PRIVATE Eigen3::Eigen lamp_linalg)
lamp_copy_linalg_runtime(lamp_superlu_contract)
add_test(NAME superlu_contract COMMAND lamp_superlu_contract)

# Formal VOF HYDRO App schema and streamed output contract.
if(NOT DEFINED Python3_EXECUTABLE AND EXISTS "${PROJECT_SOURCE_DIR}/.venv/Scripts/python.exe")
    set(Python3_EXECUTABLE "${PROJECT_SOURCE_DIR}/.venv/Scripts/python.exe")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
if(WIN32)
    add_test(NAME vof_restart_launcher_contract COMMAND "${Python3_EXECUTABLE}" -B
        "${CMAKE_CURRENT_SOURCE_DIR}/contracts/runtime/vof_restart_launcher_contract.py"
        "${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/runs/checks/vof_restart_launcher_contract")
    set_tests_properties(vof_restart_launcher_contract PROPERTIES TIMEOUT 120)
    add_test(NAME vof_restart_contract COMMAND "${Python3_EXECUTABLE}" -B
        "${CMAKE_CURRENT_SOURCE_DIR}/contracts/runtime/vof_restart_contract.py"
        "$<TARGET_FILE:LAMP>" "${PROJECT_SOURCE_DIR}/runs/checks/vof_restart_contract")
    set_tests_properties(vof_restart_contract PROPERTIES
        LABELS "app;runtime;vof;restart" TIMEOUT 600)
endif()
add_test(NAME hydro_app_contract COMMAND "${Python3_EXECUTABLE}" -B
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/hydrodynamics/hydro_app_contract.py"
    "$<TARGET_FILE:LAMP>" "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runs/checks/hydro_app_contract")
set_tests_properties(hydro_app_contract PROPERTIES LABELS "app;hydrodynamics;HYDRO")

# Formal HHO with surface charge Electric App schema, clock/charge transaction and Tecplot output.
# Test-only registration: uses the production LAMP executable, no new numerical target.
add_test(NAME electric_app_contract COMMAND "${Python3_EXECUTABLE}" -B
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/electric/electric_app_contract.py"
    "$<TARGET_FILE:LAMP>" "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runs/checks/electric_app_contract")
set_tests_properties(electric_app_contract PROPERTIES LABELS "app;electric;HHO")

# Retired EHD App-owner contracts are archived and no longer applicable.


# Fixture-dependent historical contracts remain visible and are never counted
# as passing if the original archived input is absent from this checkout.
foreach(_test pic_checkpoint_contract)
    set_tests_properties(${_test} PROPERTIES LABELS "compatibility;missing_fixture" DISABLED TRUE)
endforeach()
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/mesh/legacy_coefficients/qmc/mesh_cache.bin" AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/fixtures/mesh/legacy_coefficients/bmc/cache/mesh_cache.bin")
    set_tests_properties(pic_checkpoint_contract PROPERTIES DISABLED FALSE LABELS "compatibility")
endif()

# Strict directed-rounding tests share the production arithmetic boundary.
file(GLOB _strict "${CMAKE_CURRENT_SOURCE_DIR}/contracts/interface/interface_*characteristic_contract.cpp")
list(APPEND _strict
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/interface/interface_polynomial_velocity_contract.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/interface/interface_characteristic_enclosure_contract.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/interface/interface_characteristic_phase_contract.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/interface/interface_bounded_spacetime_flux_contract.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/interface/interface_affine_phase_contract.cpp")
if(MSVC)
    set_source_files_properties(${_strict} PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties(${_strict} PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()
# The historical EHD joint App owner is archived; its contract is not a PASS.

get_property(_tests DIRECTORY PROPERTY TESTS)
foreach(_test IN LISTS _tests)
    set(_work "${PROJECT_SOURCE_DIR}/runs/checks/${_test}")
    file(MAKE_DIRECTORY "${_work}")
    set_tests_properties(${_test} PROPERTIES WORKING_DIRECTORY "${_work}"
        TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=2;OMP_THREAD_LIMIT=2;OMP_DYNAMIC=false;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endforeach()
# This contract deliberately compares real 1/4-thread executions; the shared
# two-thread cap above must not silently reduce the requested four-thread team.
set_tests_properties(electric_vof_charge_second_order_contract PROPERTIES TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=4;OMP_THREAD_LIMIT=4;OMP_DYNAMIC=false;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
set_tests_properties(vof_charge_advection_input_contract PROPERTIES TIMEOUT 120)
set_tests_properties(dictionary_index_contract planar_cell_cut_contract geometric_vof_publication_contract PROPERTIES TIMEOUT 5)

set_tests_properties(geometric_vof_mesh_transfer_contract PROPERTIES TIMEOUT 10)

set_tests_properties(geometric_vof_transport_contract PROPERTIES TIMEOUT 30)







# Shared Core quadratic geometry reference contracts.
lamp_add_module_contract(quadratic_exact_graph_contract lamp_core)
lamp_add_module_contract(quadratic_moment_plan_contract lamp_core)
lamp_add_module_contract(quadratic_axis_volume_graph_contract lamp_core)
lamp_add_module_contract(quadratic_fixed_port_boundary_jet_contract lamp_core)
lamp_add_module_contract(quadratic_three_port_graph_contract lamp_core)
target_include_directories(lamp_quadratic_three_port_graph_contract PRIVATE
    "${LAMP_MULTIPRECISION_INCLUDE_DIR}"
    "${LAMP_BOOST_CONFIG_INCLUDE_DIR}")
target_include_directories(lamp_quadratic_fixed_port_boundary_jet_contract PRIVATE
    "${LAMP_MULTIPRECISION_INCLUDE_DIR}"
    "${LAMP_BOOST_CONFIG_INCLUDE_DIR}")
lamp_add_module_contract(quadratic_curve_endpoint_contract lamp_core)




# Actual VOF native quad DAT and direct numeric tables from the formal App.
add_test(NAME hydro_tecplot_output_contract COMMAND "${Python3_EXECUTABLE}" -B
    "${CMAKE_CURRENT_SOURCE_DIR}/contracts/hydrodynamics/hydro_tecplot_output_contract.py"
    "$<TARGET_FILE:LAMP>" "${PROJECT_SOURCE_DIR}"
    "${PROJECT_SOURCE_DIR}/runs/checks/hydro_tecplot_output_contract")
file(MAKE_DIRECTORY "${PROJECT_SOURCE_DIR}/runs/checks/hydro_tecplot_output_contract")
set_tests_properties(hydro_tecplot_output_contract PROPERTIES
    LABELS "app;hydrodynamics;tecplot" TIMEOUT 180
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/runs/checks/hydro_tecplot_output_contract"
    ENVIRONMENT "OMP_NUM_THREADS=2;OMP_THREAD_LIMIT=2;OMP_DYNAMIC=false;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

# Original CSG nozzle physical boundary roles and selected target-tree workflow.
add_test(NAME vof_csg_fluid_boundary_contract COMMAND lamp_runtime_registry_contract
    --vof-csg-fluid-boundary-only "D:/Chen/Codex/temp/lamp_electrospray_csg_fluid_boundary_20261019/ctest")
add_test(NAME vof_csg_target_amr_contract COMMAND lamp_runtime_registry_contract
    --vof-csg-target-amr-only "D:/Chen/Codex/temp/lamp_electrospray_csg_target_amr_20261019/ctest")
set_tests_properties(vof_csg_fluid_boundary_contract vof_csg_target_amr_contract PROPERTIES TIMEOUT 180)

add_test(NAME electrospray_profile_input_contract COMMAND lamp_runtime_registry_contract
    --electrospray-profile-input-only "D:/Chen/Codex/temp/lamp_electrospray_app_profile_20261019/ctest")
set_tests_properties(electrospray_profile_input_contract PROPERTIES TIMEOUT 180)

add_test(NAME vof_csg_native_leaf_overlap_contract COMMAND lamp_geometric_vof_mesh_transfer_contract --native-leaf-overlap-only)
add_test(NAME electric_vof_weak_field_gauss_contract COMMAND lamp_electric_solver_contract --electric-vof-weak-field-gauss-only)
set_tests_properties(vof_csg_native_leaf_overlap_contract electric_vof_weak_field_gauss_contract PROPERTIES TIMEOUT 180)
add_test(NAME hydro_hf_dyadic_width_contract COMMAND lamp_hydro_balanced_csf_load_contract --hf-dyadic-width-only)
