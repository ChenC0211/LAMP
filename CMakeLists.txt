cmake_minimum_required(VERSION 3.20)

project(LAMP LANGUAGES CXX)

# MSVC: 并行编译写 PDB 可能冲突（C1041），/FS 用文件锁序列化写入
if(MSVC)
    # Serialize PDB writes, support large objects, and read sources as UTF-8.
    add_compile_options(/FS /bigobj /utf-8)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Keep shared load/charge arithmetic and local electric QR in the same
# precise non-LTO boundary as the measured sharp path. Do not propagate this
# choice to unrelated modules or alter the strict Hermite engine object.
set(_quantum_precise_interface_sources
    "${CMAKE_SOURCE_DIR}/src/Core/Base/ScalarExpression.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/InterfaceGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PlanarCellCut.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/QuadraticCellCut.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/CubicGraphCellMeasure.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/AxisymmetricNeckCap.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ArcScalarRegularization.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceConservativePhase.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceMeshTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceConservativeGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceConservativeEventGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceHeightCurvature.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceCurvature.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceColumnMode.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceNeckCap.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfacePlicTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceSurfaceTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceTransportInput.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfacePlicSpaceTime.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroInterfaceLoadGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/SurfaceChargeTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/SurfaceChargeField.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/SurfaceChargeMeshTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/SurfaceChargePreparedTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/SurfaceChargeRelaxation.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetInterfacePanels.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetNativeProjection.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetNativePanels.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/AxisymmetricBoundaryIntegral.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveHermiteVolumeSensitivity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteReconstruction.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteGals.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/InterfaceCellMeasure.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/InterfaceTopology.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedDomain.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/NativeMeshIdentity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/NativeLeafDomain.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedMomentumDuals.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedMomentumLayout.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedBoundaryMomentum.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedBoundaryMomentumAssembly.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedInterfaceLoadGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFittedInterfaceTractionLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroUniformLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroConservativeTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroConservativePressure.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroNativeVelocityTrace.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroConservativeStress.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydrodynamicsConservativeNative.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydrodynamicsConservativeTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydrodynamicsConservative.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Field/StaggeredVelocitySampler.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/SharpDielectricElectrostatics.cpp")
if(MSVC)
    set_source_files_properties(${_quantum_precise_interface_sources}
        PROPERTIES COMPILE_OPTIONS "/fp:precise;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties(${_quantum_precise_interface_sources}
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()
unset(_quantum_precise_interface_sources)

# Dependency locations and optional toolchains are supplied by presets or
# standard CMAKE_PREFIX_PATH/CMAKE_TOOLCHAIN_FILE cache entries.
find_package(Eigen3 CONFIG QUIET)
if(NOT Eigen3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core PATH_SUFFIXES eigen3 include/eigen3 REQUIRED)
endif()

# Older CMake FindEigen3 modes may not define Eigen3::Eigen.
if(NOT TARGET Eigen3::Eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    if(DEFINED EIGEN3_INCLUDE_DIRS)
        target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIRS})
    elseif(DEFINED EIGEN3_INCLUDE_DIR)
        target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "Eigen3 found but include directory variables are missing")
    endif()
endif()
find_package(OpenMP)
# Shared compiler settings belong to the build, not to the main executable.
# Core/Solvers retain their existing narrower compiler contract below.
add_library(quantum_compile_options INTERFACE)
if(MSVC)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(_quantum_openmp_option -openmp:llvm)
    else()
        set(_quantum_openmp_option /openmp:llvm)
    endif()
    target_compile_options(quantum_compile_options INTERFACE /W3 ${_quantum_openmp_option})
    target_compile_definitions(quantum_compile_options INTERFACE _CRT_SECURE_NO_WARNINGS)
elseif(OpenMP_CXX_FOUND)
    target_link_libraries(quantum_compile_options INTERFACE OpenMP::OpenMP_CXX)
endif()

# Explicit existing libraries also support source-only copies without vcpkg.
# Keep machine paths in CMakeUserPresets.json, never in the project targets.
if(LAMP_SUPERLU_INCLUDE_DIR AND LAMP_SUPERLU_LIBRARY AND LAMP_BLAS_LIBRARY)
    foreach(_path LAMP_SUPERLU_INCLUDE_DIR LAMP_SUPERLU_LIBRARY LAMP_BLAS_LIBRARY)
        if(NOT EXISTS "${${_path}}")
            message(FATAL_ERROR "${_path} does not exist: ${${_path}}")
        endif()
    endforeach()
    add_library(SuperLU::SuperLU INTERFACE IMPORTED)
    set_target_properties(SuperLU::SuperLU PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${LAMP_SUPERLU_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${LAMP_SUPERLU_LIBRARY};${LAMP_BLAS_LIBRARY}")
elseif(LAMP_SUPERLU_INCLUDE_DIR OR LAMP_SUPERLU_LIBRARY OR LAMP_BLAS_LIBRARY)
    message(FATAL_ERROR "Specify all three LAMP_SUPERLU_INCLUDE_DIR, LAMP_SUPERLU_LIBRARY and LAMP_BLAS_LIBRARY")
endif()
if(NOT TARGET SuperLU::SuperLU)
    find_package(SuperLU QUIET)
    if(NOT SuperLU_FOUND AND NOT TARGET SuperLU::SuperLU AND NOT TARGET superlu::superlu)
        # Some installations/vcpkg ports may export a lower-case package name.
        find_package(superlu QUIET)
    endif()
    # Package builds may need METIS and an OpenBLAS runtime. Explicit local
    # libraries already name their BLAS dependency and need no package probes.
    find_package(metis CONFIG QUIET)
    find_package(METIS QUIET)
    if(WIN32)
        find_package(OpenBLAS CONFIG QUIET)
    endif()
endif()

# One public dependency contract for production solvers and standalone Hydro
# tests. The macro also changes Eigen-containing public class layouts, so it
# must propagate with the include paths and link libraries to every consumer.
add_library(quantum_linalg INTERFACE)
target_compile_definitions(quantum_linalg INTERFACE EIGEN_SUPERLU_SUPPORT)
set(_quantum_superlu_target "")
if(TARGET SuperLU::SuperLU)
    set(_quantum_superlu_target SuperLU::SuperLU)
elseif(TARGET superlu::superlu)
    set(_quantum_superlu_target superlu::superlu)
elseif(SuperLU_FOUND AND SuperLU_LIBRARIES)
    target_include_directories(quantum_linalg INTERFACE ${SuperLU_INCLUDE_DIRS} ${SuperLU_INCLUDE_DIR})
    target_link_libraries(quantum_linalg INTERFACE ${SuperLU_LIBRARIES})
    message(STATUS "SuperLU: required backend, legacy libraries ${SuperLU_LIBRARIES}")
else()
    message(FATAL_ERROR
        "LAMP requires SuperLU for Hydro/EHD; SparseLU fallback is disabled. "
        "Set LAMP_SUPERLU_INCLUDE_DIR, LAMP_SUPERLU_LIBRARY and LAMP_BLAS_LIBRARY, "
        "or set CMAKE_PREFIX_PATH to an existing SuperLU installation and its dependencies.")
endif()
if(_quantum_superlu_target)
    target_link_libraries(quantum_linalg INTERFACE ${_quantum_superlu_target})
    message(STATUS "SuperLU: required backend, target ${_quantum_superlu_target}")
endif()

set(_quantum_metis_target "")
foreach(_candidate metis::metis METIS::metis metis)
    if(TARGET ${_candidate})
        set(_quantum_metis_target ${_candidate})
        target_link_libraries(quantum_linalg INTERFACE ${_candidate})
        break()
    endif()
endforeach()
unset(_candidate)

# /openmp:llvm needs its runtime beside each executable on a clean shell.
if(MSVC)
    string(REGEX REPLACE "/Tools/MSVC/.*" "" _quantum_vc_root "${CMAKE_CXX_COMPILER}")
    file(GLOB _quantum_openmp_runtimes
        "${_quantum_vc_root}/Redist/MSVC/*/debug_nonredist/x64/Microsoft.VC143.OpenMP.LLVM/libomp140.x86_64.dll")
    list(SORT _quantum_openmp_runtimes COMPARE NATURAL ORDER DESCENDING)
    if(_quantum_openmp_runtimes)
        list(GET _quantum_openmp_runtimes 0 _quantum_openmp_default)
    endif()
    set(LAMP_OPENMP_RUNTIME "${_quantum_openmp_default}" CACHE FILEPATH "MSVC LLVM OpenMP runtime DLL")
    if(NOT EXISTS "${LAMP_OPENMP_RUNTIME}")
        message(FATAL_ERROR "Set LAMP_OPENMP_RUNTIME to the installed libomp140.x86_64.dll")
    endif()
endif()
function(quantum_copy_linalg_runtime target)
    if(WIN32)
        if(MSVC)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${LAMP_OPENMP_RUNTIME}" "$<TARGET_FILE_DIR:${target}>" VERBATIM)
        endif()
        foreach(_dependency ${_quantum_superlu_target} OpenBLAS::OpenBLAS ${_quantum_metis_target})
            if(TARGET ${_dependency})
                get_target_property(_dependency_type ${_dependency} TYPE)
                if(_dependency_type STREQUAL "SHARED_LIBRARY")
                    add_custom_command(TARGET ${target} POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                            "$<TARGET_FILE:${_dependency}>" "$<TARGET_FILE_DIR:${target}>"
                        VERBATIM)
                endif()
            endif()
        endforeach()
    endif()
endfunction()

# Core sources (no module/runtime business logic)
# Cartesian coordinates have one out-of-line arithmetic path, shared by
# construction/refinement/cache validation. Do not contract min + index*h or
# re-inline it into callers built with different floating-point options.
# Bounded circulation likewise owns its exact residual arithmetic in one
# strict TU; consumers receive only checked binary64 flows/residual bounds.
if(MSVC)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Core/Mesh/DyadicGridCoordinates.cpp"
        "${CMAKE_SOURCE_DIR}/src/Core/Discrete/BoundedCirculation.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Core/Mesh/DyadicGridCoordinates.cpp"
        "${CMAKE_SOURCE_DIR}/src/Core/Discrete/BoundedCirculation.cpp"
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()
file(GLOB_RECURSE LAMP_CORE_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Core/Base/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Discrete/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Field/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Lagrangian/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Units/*.cpp"
    # Minimal IO in Core (generic dictionary + mesh/geometry schema + filesystem)
    "${CMAKE_SOURCE_DIR}/src/Core/IO/Dictionary.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/DictionaryBlockText.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/DictionaryMaterials.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/DictionaryMesh.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/DictionaryGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/DictionaryUnits.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/FileManagement.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/OMPPlacementDiag.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/TecplotWriter.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/IO/TecplotMeshOutput.cpp"
)

# Geometry's error-free arithmetic is implementation-only and TU-local. Keep
# its FP/environment contract out of unrelated Core/App compilation flags.
set(_quantum_hermite_sources
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PolynomialVelocitySlab.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/NativeLeafCoverage.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetHermiteGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteSkeleton.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveHermiteTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteTransfer.cpp")
list(REMOVE_ITEM LAMP_CORE_SOURCES ${_quantum_hermite_sources})
add_library(quantum_hermite_geometry OBJECT ${_quantum_hermite_sources})
target_include_directories(quantum_hermite_geometry PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(quantum_hermite_geometry PRIVATE Eigen3::Eigen)
target_compile_options(quantum_hermite_geometry PRIVATE $<TARGET_PROPERTY:quantum_core,COMPILE_OPTIONS>)
target_compile_definitions(quantum_hermite_geometry PRIVATE $<TARGET_PROPERTY:quantum_core,COMPILE_DEFINITIONS>)
set_target_properties(quantum_hermite_geometry PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION FALSE UNITY_BUILD OFF)
foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(quantum_hermite_geometry PROPERTIES "INTERPROCEDURAL_OPTIMIZATION_${_cfg}" FALSE)
endforeach()
if(MSVC)
    target_compile_options(quantum_hermite_geometry PRIVATE /fp:strict /GL-)
else()
    target_compile_options(quantum_hermite_geometry PRIVATE
        -fno-fast-math -ffp-contract=off -frounding-math -fno-lto)
endif()
add_library(quantum_core STATIC ${LAMP_CORE_SOURCES} $<TARGET_OBJECTS:quantum_hermite_geometry>)
unset(_quantum_hermite_sources)

target_include_directories(quantum_core PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(quantum_core PUBLIC Eigen3::Eigen)

# -------------------------------
# Solvers layer (top-level FV/EB solvers; depends on Core, no module dependency)
# -------------------------------

file(GLOB_RECURSE LAMP_SOLVERS_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Solvers/*.cpp"
)

add_library(quantum_solvers STATIC ${LAMP_SOLVERS_SOURCES})

target_include_directories(quantum_solvers PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(quantum_solvers PUBLIC quantum_core Eigen3::Eigen quantum_linalg)

# Runtime infrastructure and named application workflows. The main entrypoint
# remains explicit below; moved files retain the same target/compiler options.
file(GLOB_RECURSE LAMP_WORKFLOW_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Apps/*.cpp"
)
list(REMOVE_ITEM LAMP_WORKFLOW_SOURCES "${CMAKE_SOURCE_DIR}/src/Apps/lamp_main.cpp")

# The interface state is a real module library. Its independent link check
# below has neither App sources nor the include/Apps search path. Future
# transport services belong to this target, not to a second executable copy.
file(GLOB_RECURSE _quantum_interface_sources CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/*.cpp")
list(FILTER _quantum_interface_sources EXCLUDE REGEX "/Compatibility/")
add_library(quantum_interface STATIC ${_quantum_interface_sources})
target_link_libraries(quantum_interface PRIVATE quantum_compile_options PUBLIC quantum_core PRIVATE quantum_linalg)
unset(_quantum_interface_sources)

# Existing Hydro/EHD adapters have an explicit owner, separate from FirstJet.
# This target is not an alternative Interface App and admits no new consumers.
file(GLOB_RECURSE _lamp_interface_compatibility CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/Compatibility/*.cpp")
add_library(lamp_interface_compatibility STATIC ${_lamp_interface_compatibility})
target_link_libraries(lamp_interface_compatibility PRIVATE quantum_compile_options
    PUBLIC quantum_interface PRIVATE quantum_linalg)
unset(_lamp_interface_compatibility)

# All electric services compile/link without Hydro or App runtime.
file(GLOB_RECURSE LAMP_ELECTRIC_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/*.cpp")
# The old vertex/particle model is an explicit compatibility target below.
set(_quantum_legacy_electricity_source
    "${CMAKE_SOURCE_DIR}/src/Modules/Electricity/Legacy/Electricity.cpp")
list(REMOVE_ITEM LAMP_ELECTRIC_SOURCES ${_quantum_legacy_electricity_source})
add_library(quantum_electricity STATIC ${LAMP_ELECTRIC_SOURCES})
target_link_libraries(quantum_electricity PRIVATE quantum_compile_options PUBLIC quantum_solvers)

# Native CCH has its own implementation owner, with no legacy Hydro/Interface
# compatibility link. Keep one explicit source list for both ownership and
# exclusion from the frozen Hydro glob below.
set(_lamp_cch_sources
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchDofLayout.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchIntegralDivergence.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchPressureLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchInterfaceNormalLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchSinglePhaseSaddle.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchPhysicalDivergence.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchEnergyViscosity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchPhaseGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchRegularizedColor.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchBalancedCsfLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchPhaseViscosity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchMaterialTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchMomentumGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchSpatialAdvector.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchRegularAdvectorGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchCharacteristic.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchDualMaterial.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchKepConvection.cpp")
add_library(lamp_hydrodynamics_cch STATIC ${_lamp_cch_sources})
target_link_libraries(lamp_hydrodynamics_cch PRIVATE quantum_compile_options PUBLIC quantum_core quantum_interface)

# BF-CSF visibility reuses the exact Core expansion engine on original
# selected-leaf ray endpoints. Keep this one consumer's arithmetic strict.
if(MSVC)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchBalancedCsfLoad.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/CchBalancedCsfLoad.cpp"
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()

# Frozen Hydro owns the existing NS implementation and its existing Interface
# compatibility dependency; CCH sources must never also compile into this lib.
file(GLOB_RECURSE LAMP_HYDRO_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/*.cpp")
list(REMOVE_ITEM LAMP_HYDRO_SOURCES ${_lamp_cch_sources})
unset(_lamp_cch_sources)
add_library(quantum_hydrodynamics STATIC ${LAMP_HYDRO_SOURCES})
target_link_libraries(quantum_hydrodynamics PRIVATE quantum_compile_options PUBLIC lamp_interface_compatibility quantum_solvers)

# Plasma compiles against physical inputs, with no Runtime/App include path.
file(GLOB_RECURSE LAMP_PLASMA_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Plasma/*.cpp")
add_library(quantum_plasma STATIC ${LAMP_PLASMA_SOURCES})
target_link_libraries(quantum_plasma PRIVATE quantum_compile_options PUBLIC quantum_solvers)

# Only the old init-only workflow uses this mixed vertex/particle model.
# Formal Electricity remains independently linked without Plasma.
add_library(quantum_legacy_electricity STATIC ${_quantum_legacy_electricity_source})
target_link_libraries(quantum_legacy_electricity PRIVATE quantum_compile_options PUBLIC quantum_electricity quantum_plasma quantum_hydrodynamics)
unset(_quantum_legacy_electricity_source)

# Share the actual spray workflow between the main executable and its coupled
# state checks. Physics stays in the module libraries linked by this App.
file(GLOB_RECURSE _quantum_spray_sources CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Apps/Electrospray/*.cpp")
list(REMOVE_ITEM LAMP_WORKFLOW_SOURCES ${_quantum_spray_sources})
add_library(quantum_electrospray STATIC ${_quantum_spray_sources})
target_include_directories(quantum_electrospray PRIVATE
    "${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/include/Apps")
target_link_libraries(quantum_electrospray PRIVATE quantum_compile_options PUBLIC quantum_hydrodynamics quantum_electricity)
unset(_quantum_spray_sources)

add_executable(LAMP ${LAMP_WORKFLOW_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/Apps/lamp_main.cpp")

target_include_directories(LAMP PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    # Runtime / Apps layer headers live under include/Apps
    ${CMAKE_SOURCE_DIR}/include/Apps
)

# Keep the executable and its runtime DLLs together in bin/.
set_target_properties(LAMP PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    if(_cfg STREQUAL "RELEASE")
        set(_quantum_runtime_dir "${CMAKE_SOURCE_DIR}/bin")
    else()
        string(TOLOWER "${_cfg}" _quantum_runtime_config)
        set(_quantum_runtime_dir "${CMAKE_SOURCE_DIR}/bin/${_quantum_runtime_config}")
    endif()
    set_target_properties(LAMP PROPERTIES "RUNTIME_OUTPUT_DIRECTORY_${_cfg}" "${_quantum_runtime_dir}")
endforeach()
unset(_quantum_runtime_dir)
unset(_quantum_runtime_config)

target_link_libraries(LAMP PRIVATE quantum_compile_options quantum_electrospray quantum_plasma quantum_legacy_electricity)

quantum_copy_linalg_runtime(LAMP)

if(MSVC)
    # Emit a linker map so Windows APPCRASH offsets (e.g. 0x335b0) can be resolved quickly.
    target_link_options(LAMP PRIVATE "/MAP:${CMAKE_BINARY_DIR}/LAMP.map")

    # 代码中使用了 OpenMP reduction(max:...) 等特性；MSVC 传统 /openmp 仅支持较旧规范。
    # 使用 LLVM OpenMP 前端以支持更多 OpenMP 语法（VS 2022 17.x）。
    target_compile_options(quantum_core PRIVATE ${_quantum_openmp_option})
    target_compile_options(quantum_solvers PRIVATE ${_quantum_openmp_option})
endif()
unset(_quantum_openmp_option)

include(CTest)
option(LAMP_BUILD_BENCHMARKS "Build complete research/benchmark drivers under cases/benchmarks" OFF)
if(LAMP_BUILD_BENCHMARKS)
    add_subdirectory(cases/benchmarks)
endif()
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
