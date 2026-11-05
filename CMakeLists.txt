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
set(_lamp_precise_interface_sources
    "${CMAKE_SOURCE_DIR}/src/Core/Base/ScalarExpression.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/InterfaceGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PlanarCellCut.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PlanarFaceTimeIntegral.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofPublication.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofMeshTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofPhaseBudget.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroPhaseGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMomentumGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/QuadraticCellCut.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/CubicGraphCellMeasure.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/AxisymmetricNeckCap.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ArcScalarRegularization.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/SurfaceChargeTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/SurfaceChargeField.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/SurfaceChargeMeshTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/SurfaceChargePreparedTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetInterfacePanels.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetNativeProjection.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetNativePanels.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/AxisymmetricBoundaryIntegral.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveHermiteVolumeSensitivity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteReconstruction.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteGals.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/InterfaceCellMeasure.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/NativeMeshIdentity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/NativeLeafDomain.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Field/StaggeredVelocitySampler.cpp")
if(MSVC)
    set_source_files_properties(${_lamp_precise_interface_sources}
        PROPERTIES COMPILE_OPTIONS "/fp:precise;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties(${_lamp_precise_interface_sources}
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()
unset(_lamp_precise_interface_sources)
# Actual-point polygon predicates/measures, the pure VOF budget and its
# direct exact-rounding contract and VOF Hydro material consumers reuse
# Core error-free transforms.
# Keep these strict and outside LTO.
if(MSVC)
    set_source_files_properties(
        "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PlanarCellCut.cpp"
        "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PlanarFaceTimeIntegral.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofPhaseBudget.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroPhaseGeometry.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMomentumGeometry.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
endif()

# Direct VOF fraction/initial-inventory arithmetic uses Core exact cuts.
# Keep these standalone translation units outside contraction/LTO and unity.
if(MSVC)
    set_source_files_properties(
        "${CMAKE_SOURCE_DIR}/src/Modules/Interface/VofPhaseFraction.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofInitialization.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties(
        "${CMAKE_SOURCE_DIR}/src/Modules/Interface/VofPhaseFraction.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Interface/GeometricVofInitialization.cpp"
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto"
        SKIP_UNITY_BUILD_INCLUSION ON)
endif()

# B1 candidate operators retain the exact-arithmetic environment of their
# public VOF source and derived tensor contractions; keep the scope local.
set(_lamp_b1_sources
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/VofSmoothGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofSharpField.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofSharpMaxwellStress.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricSurfaceChargeTransport.cpp")
if(MSVC)
    set_source_files_properties(${_lamp_b1_sources}
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties(${_lamp_b1_sources}
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto"
        SKIP_UNITY_BUILD_INCLUSION ON)
endif()
unset(_lamp_b1_sources)

# Dependency locations and optional toolchains are supplied by presets or
# standard CMAKE_PREFIX_PATH/CMAKE_TOOLCHAIN_FILE cache entries.
find_package(Eigen3 CONFIG QUIET)
if(NOT Eigen3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core PATH_SUFFIXES eigen3 include/eigen3 REQUIRED)
endif()

# The exact framed-strip backend is header-only. Machine locations stay in
# CMakeUserPresets.json or the normal cache/search paths, never in this file.
find_path(LAMP_MULTIPRECISION_INCLUDE_DIR boost/multiprecision/cpp_int.hpp REQUIRED)
find_path(LAMP_BOOST_CONFIG_INCLUDE_DIR boost/config.hpp REQUIRED)

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
add_library(lamp_compile_options INTERFACE)
if(MSVC)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(_lamp_openmp_option -openmp:llvm)
    else()
        set(_lamp_openmp_option /openmp:llvm)
    endif()
    target_compile_options(lamp_compile_options INTERFACE /W3 ${_lamp_openmp_option})
    target_compile_definitions(lamp_compile_options INTERFACE _CRT_SECURE_NO_WARNINGS)
elseif(OpenMP_CXX_FOUND)
    target_link_libraries(lamp_compile_options INTERFACE OpenMP::OpenMP_CXX)
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
add_library(lamp_linalg INTERFACE)
target_compile_definitions(lamp_linalg INTERFACE EIGEN_SUPERLU_SUPPORT)
set(_lamp_superlu_target "")
if(TARGET SuperLU::SuperLU)
    set(_lamp_superlu_target SuperLU::SuperLU)
elseif(TARGET superlu::superlu)
    set(_lamp_superlu_target superlu::superlu)
elseif(SuperLU_FOUND AND SuperLU_LIBRARIES)
    target_include_directories(lamp_linalg INTERFACE ${SuperLU_INCLUDE_DIRS} ${SuperLU_INCLUDE_DIR})
    target_link_libraries(lamp_linalg INTERFACE ${SuperLU_LIBRARIES})
    message(STATUS "SuperLU: required backend, legacy libraries ${SuperLU_LIBRARIES}")
else()
    message(FATAL_ERROR
        "LAMP requires SuperLU for Hydro/EHD; SparseLU fallback is disabled. "
        "Set LAMP_SUPERLU_INCLUDE_DIR, LAMP_SUPERLU_LIBRARY and LAMP_BLAS_LIBRARY, "
        "or set CMAKE_PREFIX_PATH to an existing SuperLU installation and its dependencies.")
endif()
if(_lamp_superlu_target)
    target_link_libraries(lamp_linalg INTERFACE ${_lamp_superlu_target})
    message(STATUS "SuperLU: required backend, target ${_lamp_superlu_target}")
endif()

set(_lamp_metis_target "")
foreach(_candidate metis::metis METIS::metis metis)
    if(TARGET ${_candidate})
        set(_lamp_metis_target ${_candidate})
        target_link_libraries(lamp_linalg INTERFACE ${_candidate})
        break()
    endif()
endforeach()
unset(_candidate)

# /openmp:llvm needs its runtime beside each executable on a clean shell.
if(MSVC)
    string(REGEX REPLACE "/Tools/MSVC/.*" "" _lamp_vc_root "${CMAKE_CXX_COMPILER}")
    file(GLOB _lamp_openmp_runtimes
        "${_lamp_vc_root}/Redist/MSVC/*/debug_nonredist/x64/Microsoft.VC143.OpenMP.LLVM/libomp140.x86_64.dll")
    list(SORT _lamp_openmp_runtimes COMPARE NATURAL ORDER DESCENDING)
    if(_lamp_openmp_runtimes)
        list(GET _lamp_openmp_runtimes 0 _lamp_openmp_default)
    endif()
    set(LAMP_OPENMP_RUNTIME "${_lamp_openmp_default}" CACHE FILEPATH "MSVC LLVM OpenMP runtime DLL")
    if(NOT EXISTS "${LAMP_OPENMP_RUNTIME}")
        message(FATAL_ERROR "Set LAMP_OPENMP_RUNTIME to the installed libomp140.x86_64.dll")
    endif()
endif()
function(lamp_copy_linalg_runtime target)
    if(WIN32)
        if(MSVC)
            add_custom_command(TARGET ${target} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${LAMP_OPENMP_RUNTIME}" "$<TARGET_FILE_DIR:${target}>" VERBATIM)
        endif()
        foreach(_dependency ${_lamp_superlu_target} OpenBLAS::OpenBLAS ${_lamp_metis_target})
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
set(_lamp_hermite_sources
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/PolynomialVelocitySlab.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Mesh/NativeLeafCoverage.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/LevelSetHermiteGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteSkeleton.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveHermiteTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Core/Geometry/ActiveQuarticHermiteCompletion.cpp")
list(REMOVE_ITEM LAMP_CORE_SOURCES ${_lamp_hermite_sources})
add_library(lamp_hermite_geometry OBJECT ${_lamp_hermite_sources})
target_include_directories(lamp_hermite_geometry PRIVATE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(lamp_hermite_geometry PRIVATE Eigen3::Eigen)
target_compile_options(lamp_hermite_geometry PRIVATE $<TARGET_PROPERTY:lamp_core,COMPILE_OPTIONS>)
target_compile_definitions(lamp_hermite_geometry PRIVATE $<TARGET_PROPERTY:lamp_core,COMPILE_DEFINITIONS>)
set_target_properties(lamp_hermite_geometry PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION FALSE UNITY_BUILD OFF)
foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    set_target_properties(lamp_hermite_geometry PROPERTIES "INTERPROCEDURAL_OPTIMIZATION_${_cfg}" FALSE)
endforeach()
if(MSVC)
    target_compile_options(lamp_hermite_geometry PRIVATE /fp:strict /GL-)
else()
    target_compile_options(lamp_hermite_geometry PRIVATE
        -fno-fast-math -ffp-contract=off -frounding-math -fno-lto)
endif()
add_library(lamp_core STATIC ${LAMP_CORE_SOURCES} $<TARGET_OBJECTS:lamp_hermite_geometry>)
unset(_lamp_hermite_sources)

target_include_directories(lamp_core PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# The public finite-volume dyadic certificate is consumed by Interface.
# Export its existing header-only dependency with the public Core headers.
target_include_directories(lamp_core PUBLIC
    "${LAMP_MULTIPRECISION_INCLUDE_DIR}"
    "${LAMP_BOOST_CONFIG_INCLUDE_DIR}")
target_link_libraries(lamp_core PUBLIC Eigen3::Eigen)

# -------------------------------
# Solvers layer (top-level FV/EB solvers; depends on Core, no module dependency)
# -------------------------------

file(GLOB_RECURSE LAMP_SOLVERS_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Solvers/*.cpp"
)

add_library(lamp_solvers STATIC ${LAMP_SOLVERS_SOURCES})

target_include_directories(lamp_solvers PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(lamp_solvers PUBLIC lamp_core Eigen3::Eigen lamp_linalg)

# Shared App entry/IO facilities and active named workflows. The main entrypoint
# remains explicit below; retained HG sources stay in place as reference.
file(GLOB_RECURSE LAMP_WORKFLOW_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Apps/*.cpp"
)
list(REMOVE_ITEM LAMP_WORKFLOW_SOURCES "${CMAKE_SOURCE_DIR}/src/Apps/lamp_main.cpp")
# Retained HG route sources remain as reference; production does not compile them.
list(FILTER LAMP_WORKFLOW_SOURCES EXCLUDE REGEX "/(Interface|Electrohydrodynamics)/")
# Exact new VOF-EHD reader only; retained HG EHD sources remain excluded.
list(APPEND LAMP_WORKFLOW_SOURCES
    "${CMAKE_SOURCE_DIR}/src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsInput.cpp")
# Exact passive VOF-EHD serializer; retained HG output remains excluded.
list(APPEND LAMP_WORKFLOW_SOURCES
    "${CMAKE_SOURCE_DIR}/src/Apps/Electrohydrodynamics/IO/VofElectrohydrodynamicsOutput.cpp")
# Exact new VOF-EHD entry; old HG App units stay excluded.
list(APPEND LAMP_WORKFLOW_SOURCES
    "${CMAKE_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsApp.cpp")

# Exact new VOF-EHD joint owner; old HG workflows stay excluded.
list(APPEND LAMP_WORKFLOW_SOURCES
    "${CMAKE_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp")

# The interface state is a real module library. Its independent link check
# below has neither App sources nor the include/Apps search path. Future
# transport services belong to this target, not to a second executable copy.
file(GLOB_RECURSE _lamp_interface_sources CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Interface/*.cpp")
list(FILTER _lamp_interface_sources EXCLUDE REGEX "/Compatibility/")
add_library(lamp_interface STATIC ${_lamp_interface_sources})
target_link_libraries(lamp_interface PRIVATE lamp_compile_options PUBLIC lamp_core PRIVATE lamp_linalg)
unset(_lamp_interface_sources)


# All electric services compile/link without Hydro or App runtime.
file(GLOB_RECURSE LAMP_ELECTRIC_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Electric/*.cpp")
add_library(lamp_electric STATIC ${LAMP_ELECTRIC_SOURCES})
target_link_libraries(lamp_electric PRIVATE lamp_compile_options PUBLIC lamp_solvers lamp_interface)

# Native HYDRO has its own implementation owner, with no legacy Hydro/Interface
# compatibility link. Keep one explicit source list for both ownership and
# exclusion from the frozen Hydro glob below.
set(_lamp_hydro_sources
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroDofLayout.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroIntegralDivergence.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroPressureLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMaxwellLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroInterfaceNormalLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroSinglePhaseSaddle.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroPhysicalDivergence.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroEnergyViscosity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroPhaseGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroRegularizedColor.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroBalancedCsfLoad.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroPhaseViscosity.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMaterialTransport.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMomentumGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroSpatialAdvector.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroRegularAdvectorGeometry.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroCharacteristic.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/Detail/HydroFiniteContraction.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroDualMaterial.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroKepConvection.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroVofAdvector.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroVofStep.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroVofRefinement.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroVofMeshTransfer.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroFrozenMidpoint.cpp"
    "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroTwoPhaseStep.cpp")
add_library(lamp_hydrodynamics STATIC ${_lamp_hydro_sources})
target_link_libraries(lamp_hydrodynamics PRIVATE lamp_compile_options PUBLIC lamp_core lamp_interface lamp_electric)
target_include_directories(lamp_hydrodynamics PRIVATE
    "${LAMP_MULTIPRECISION_INCLUDE_DIR}" "${LAMP_BOOST_CONFIG_INCLUDE_DIR}")
# The derived-scalar fallback decomposes native binary values exactly. Only
# this new implementation source needs explicit strict floating semantics.
if(MSVC)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/Detail/HydroFiniteContraction.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofField.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofChargeTransport.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroKepConvection.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofMaxwellStress.cpp"
        "${CMAKE_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMaxwellLoad.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroVofMeshTransfer.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/Detail/HydroFiniteContraction.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofField.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofChargeTransport.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroKepConvection.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Electric/ElectricVofMaxwellStress.cpp"
        "${CMAKE_SOURCE_DIR}/src/Apps/Electrohydrodynamics/VofElectrohydrodynamicsRun.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroMaxwellLoad.cpp"
        "${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroVofMeshTransfer.cpp"
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()

# BF-CSF visibility reuses the exact Core expansion engine on original
# selected-leaf ray endpoints. Keep this one consumer's arithmetic strict.
if(MSVC)
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroBalancedCsfLoad.cpp"
        PROPERTIES COMPILE_OPTIONS "/fp:strict;/GL-" SKIP_UNITY_BUILD_INCLUSION ON)
else()
    set_source_files_properties("${CMAKE_SOURCE_DIR}/src/Modules/Hydrodynamics/HydroBalancedCsfLoad.cpp"
        PROPERTIES COMPILE_OPTIONS "-fno-fast-math;-ffp-contract=off;-frounding-math;-fno-lto" SKIP_UNITY_BUILD_INCLUSION ON)
endif()

unset(_lamp_hydro_sources)

# Plasma compiles against physical inputs, with no Runtime/App include path.
file(GLOB_RECURSE LAMP_PLASMA_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/src/Modules/Plasma/*.cpp")
add_library(lamp_plasma STATIC ${LAMP_PLASMA_SOURCES})
target_link_libraries(lamp_plasma PRIVATE lamp_compile_options PUBLIC lamp_solvers)

add_executable(LAMP ${LAMP_WORKFLOW_SOURCES}
    "${CMAKE_SOURCE_DIR}/src/Apps/lamp_main.cpp")

target_include_directories(LAMP PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    # Shared entry and named App headers live under include/Apps
    ${CMAKE_SOURCE_DIR}/include/Apps
)

# Keep the executable and its runtime DLLs together. Independent validation
# builds can select their own directory without replacing the shared binary.
set(LAMP_RUNTIME_DIRECTORY "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "LAMP executable and runtime output directory")
set_target_properties(LAMP PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${LAMP_RUNTIME_DIRECTORY}")
foreach(_cfg DEBUG RELEASE RELWITHDEBINFO MINSIZEREL)
    if(_cfg STREQUAL "RELEASE")
        set(_lamp_runtime_dir "${LAMP_RUNTIME_DIRECTORY}")
    else()
        string(TOLOWER "${_cfg}" _lamp_runtime_config)
        set(_lamp_runtime_dir "${LAMP_RUNTIME_DIRECTORY}/${_lamp_runtime_config}")
    endif()
    set_target_properties(LAMP PROPERTIES "RUNTIME_OUTPUT_DIRECTORY_${_cfg}" "${_lamp_runtime_dir}")
endforeach()
unset(_lamp_runtime_dir)
unset(_lamp_runtime_config)

target_link_libraries(LAMP PRIVATE lamp_compile_options lamp_hydrodynamics lamp_interface lamp_electric lamp_plasma)

lamp_copy_linalg_runtime(LAMP)

if(MSVC)
    # Emit a linker map so Windows APPCRASH offsets (e.g. 0x335b0) can be resolved quickly.
    target_link_options(LAMP PRIVATE "/MAP:${CMAKE_BINARY_DIR}/LAMP.map")

    # 代码中使用了 OpenMP reduction(max:...) 等特性；MSVC 传统 /openmp 仅支持较旧规范。
    # 使用 LLVM OpenMP 前端以支持更多 OpenMP 语法（VS 2022 17.x）。
    target_compile_options(lamp_core PRIVATE ${_lamp_openmp_option})
    target_compile_options(lamp_solvers PRIVATE ${_lamp_openmp_option})
endif()
unset(_lamp_openmp_option)

include(CTest)
option(LAMP_BUILD_BENCHMARKS "Build complete research/benchmark drivers under cases/benchmarks" OFF)
if(LAMP_BUILD_BENCHMARKS)
    add_subdirectory(cases/benchmarks)
endif()
if(BUILD_TESTING)
    add_subdirectory(tests)
    # These single-source exact arithmetic contracts are created in tests/.
    # Parent-directory source properties do not apply to that target scope.
    foreach(_lamp_exact_contract lamp_hydro_phase_geometry_contract lamp_geometric_vof_transport_contract lamp_planar_cell_cut_contract lamp_electric_solver_contract)
        if(MSVC)
            target_compile_options(${_lamp_exact_contract} PRIVATE /fp:strict /GL-)
        else()
            target_compile_options(${_lamp_exact_contract} PRIVATE -fno-fast-math -ffp-contract=off -frounding-math -fno-lto)
        endif()
        set_target_properties(${_lamp_exact_contract} PROPERTIES INTERPROCEDURAL_OPTIMIZATION FALSE UNITY_BUILD OFF)
    endforeach()
    unset(_lamp_exact_contract)
endif()
