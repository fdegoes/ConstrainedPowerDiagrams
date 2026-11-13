cmake_minimum_required(VERSION 3.18)

# The /usr/local/bin/cmake on this machine is an x86_64 binary (Rosetta), so
# CMake would otherwise target x86_64. Houdini 20.5 here is arm64.
if(APPLE)
    if(NOT CMAKE_OSX_ARCHITECTURES)
        set(CMAKE_OSX_ARCHITECTURES arm64 CACHE STRING "" FORCE)
    endif()
    # Match the deployment target used by the HDK so the static geogram
    # archive does not trigger "built for newer macOS" link warnings.
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        set(CMAKE_OSX_DEPLOYMENT_TARGET 11.0 CACHE STRING "" FORCE)
    endif()
endif()

project(PowerParticles C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Debug Release RelWithDebInfo MinSizeRel)

# Static geogram gets linked into a shared Houdini plugin.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

option(PP_BUILD_HOUDINI "Build the Houdini SOP plugin" ON)
option(PP_BUILD_TESTS "Build the standalone tests" ON)
# Tuning for the build machine (not portable binaries) and link-time
# optimization of the core library and its executables.
option(PP_NATIVE "Optimize for the build machine's CPU (-mcpu/-march=native)" OFF)
option(PP_LTO "Link-time optimization" OFF)
# Sparse Cholesky of the Newton systems: Apple Accelerate (macOS only),
# SuiteSparse's CHOLMOD (macOS and Linux, needs SUITESPARSE_DIR), or Eigen's
# SimplicialLLT (portable, slowest). Defaults to CHOLMOD, so that macOS and
# Linux builds factorize the same way.
set(PP_SPARSE_SOLVER Cholmod CACHE STRING
    "Sparse Cholesky backend: Accelerate, Cholmod or Eigen")
set_property(CACHE PP_SPARSE_SOLVER PROPERTY STRINGS Accelerate Cholmod Eigen)
if(PP_SPARSE_SOLVER STREQUAL "Accelerate" AND NOT APPLE)
    message(FATAL_ERROR "PP_SPARSE_SOLVER=Accelerate needs macOS")
elseif(NOT PP_SPARSE_SOLVER MATCHES "^(Accelerate|Cholmod|Eigen)$")
    message(FATAL_ERROR "PP_SPARSE_SOLVER must be Accelerate, Cholmod or Eigen (got '${PP_SPARSE_SOLVER}')")
endif()
# SuiteSparse install prefix (include/ and lib/), as in ../geomt.
set(SUITESPARSE_DIR "" CACHE PATH "SuiteSparse install prefix (PP_SPARSE_SOLVER=Cholmod)")
set(PP_HOUDINI_HOME "" CACHE PATH
    "Houdini user dir to install into (default: asked to hython)")

###########
# geogram #
###########
# Built as a static library from the source tree in external/geogram, with
# only the pieces we need. GEOGRAM_SUB_BUILD lets us choose the platform.
set(GEOGRAM_SUB_BUILD ON)
if(APPLE)
    set(VORPALINE_PLATFORM "Darwin-aarch64-clang" CACHE STRING "" FORCE)
elseif(UNIX)
    set(VORPALINE_PLATFORM "Linux64-gcc" CACHE STRING "")
endif()
set(GEOGRAM_LIB_ONLY ON)
set(GEOGRAM_WITH_GRAPHICS OFF)
set(GEOGRAM_WITH_LEGACY_NUMERICS OFF)
set(GEOGRAM_WITH_HLBFGS OFF)
set(GEOGRAM_WITH_TETGEN ON)
set(GEOGRAM_WITH_TRIANGLE ON)
set(GEOGRAM_WITH_LUA OFF)
set(GEOGRAM_WITH_FPG OFF)
set(GEOGRAM_WITH_GARGANTUA OFF)
set(GEOGRAM_WITH_TBB OFF)
set(GEOGRAM_WITH_EXPLORAGRAM OFF)
add_subdirectory(external/geogram EXCLUDE_FROM_ALL)

# geogram compiles its own sources with -Weverything; that noise is not ours
# to fix, so silence warnings for the geogram targets only.
foreach(_geo_target geogram geogram_third_party)
    if(TARGET ${_geo_target})
        target_compile_options(${_geo_target} PRIVATE -w)
    endif()
endforeach()

##################
# exploragram OT #
##################
# Semi-discrete optimal transport (weight solve, air particles) copied from
# BrunoLevy/exploragram@9550cc6, optimal_transport/ only. Newton mode only,
# since geogram is built without HLBFGS.
add_library(exploragram_ot STATIC
    external/exploragram_ot/exploragram/optimal_transport/linear_least_squares.cpp
    external/exploragram_ot/exploragram/optimal_transport/optimal_transport.cpp
    external/exploragram_ot/exploragram/optimal_transport/optimal_transport_2d.cpp
    external/exploragram_ot/exploragram/optimal_transport/optimal_transport_3d.cpp
)
target_include_directories(exploragram_ot PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/exploragram_ot)
target_compile_definitions(exploragram_ot PRIVATE exploragram_EXPORTS)
target_compile_options(exploragram_ot PRIVATE -w)
target_link_libraries(exploragram_ot PUBLIC geogram)

#########
# Eigen #
#########
# Header-only. Defaults to the copy used by ../geomt.
set(EIGEN3_INCLUDE_DIR "$ENV{HOME}/MyLibs/eigen-3.4.0" CACHE PATH "Eigen 3 include dir")
if(NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Core")
    message(FATAL_ERROR "Eigen not found: set EIGEN3_INCLUDE_DIR (currently '${EIGEN3_INCLUDE_DIR}')")
endif()
add_library(pp_eigen INTERFACE)
target_include_directories(pp_eigen SYSTEM INTERFACE "${EIGEN3_INCLUDE_DIR}")

########
# core #
########
# Houdini-free code: depends on geogram (private) and Eigen (public).
add_library(pp_core STATIC
    src/core/cellBackend.h
    src/core/cellBackend.cpp
    src/core/cellMoments.h
    src/core/centroidalDiagram.h
    src/core/centroidalDiagram.cpp
    src/core/constrainedIntegrator.h
    src/core/constrainedIntegrator.cpp
    src/core/domain.h
    src/core/domain.cpp
    src/core/dualSolver.h
    src/core/dualSolver.cpp
    src/core/feasibleVolumes.h
    src/core/feasibleVolumes.cpp
    src/core/geogramContext.h
    src/core/geogramContext.cpp
    src/core/geogramMesh.h
    src/core/geogramMesh.cpp
    src/core/geometry.h
    src/core/jointDual.h
    src/core/jointDual.cpp
    src/core/lineSearch.h
    src/core/lineSearch.cpp
    src/core/softenedIntegrator.h
    src/core/softenedIntegrator.cpp
    src/core/gmIntegrator.h
    src/core/gmIntegrator.cpp
    src/core/powerParticlesIntegrator.h
    src/core/powerParticlesIntegrator.cpp
    src/core/integrator.h
    src/core/integrator.cpp
    src/core/domainBoundary.h
    src/core/domainBoundary.cpp
    src/core/powerDiagram.h
    src/core/powerDiagram.cpp
    src/core/powerCells.h
    src/core/powerCells.cpp
    src/core/shiftedDiagram.h
    src/core/shiftedDiagram.cpp
    src/core/sparseCholesky.h
    src/core/sparseCholesky.cpp
    src/core/tetrahedralize.h
    src/core/tetrahedralize.cpp
    src/core/triangulate2D.h
    src/core/triangulate2D.cpp
    src/core/liftSites.h
    src/core/liftSites.cpp
)
target_include_directories(pp_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(pp_core PUBLIC pp_eigen PRIVATE geogram)
if(PP_NATIVE)
    # Apple silicon: the M1 baseline (common to all M chips; -mcpu=native
    # fails when the build tools run under Rosetta). Other ARM takes -mcpu,
    # x86 -march.
    if(APPLE AND CMAKE_OSX_ARCHITECTURES MATCHES "arm64")
        target_compile_options(pp_core PUBLIC -mcpu=apple-m1)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        target_compile_options(pp_core PUBLIC -mcpu=native)
    else()
        target_compile_options(pp_core PUBLIC -march=native)
    endif()
endif()
if(PP_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT PP_HAS_IPO OUTPUT PP_IPO_ERROR)
    if(PP_HAS_IPO)
        set_property(TARGET pp_core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "PP_LTO: ${PP_IPO_ERROR}")
    endif()
endif()
if(PP_SPARSE_SOLVER STREQUAL "Accelerate")
    target_compile_definitions(pp_core PRIVATE PP_USE_ACCELERATE=1)
    target_link_libraries(pp_core PUBLIC "-framework Accelerate")
elseif(PP_SPARSE_SOLVER STREQUAL "Cholmod")
    find_path(CHOLMOD_INCLUDE_DIR cholmod.h
        HINTS "${SUITESPARSE_DIR}/include" PATH_SUFFIXES suitesparse NO_DEFAULT_PATH)
    find_path(CHOLMOD_INCLUDE_DIR cholmod.h PATH_SUFFIXES suitesparse)
    set(_pp_cholmod_libs)
    foreach(_lib cholmod amd camd colamd ccolamd suitesparseconfig)
        find_library(PP_${_lib}_LIBRARY ${_lib} HINTS "${SUITESPARSE_DIR}/lib" "${SUITESPARSE_DIR}/lib64")
        if(NOT PP_${_lib}_LIBRARY)
            message(FATAL_ERROR "SuiteSparse: lib${_lib} not found (set SUITESPARSE_DIR, currently '${SUITESPARSE_DIR}')")
        endif()
        list(APPEND _pp_cholmod_libs ${PP_${_lib}_LIBRARY})
    endforeach()
    if(NOT CHOLMOD_INCLUDE_DIR)
        message(FATAL_ERROR "SuiteSparse: cholmod.h not found (set SUITESPARSE_DIR, currently '${SUITESPARSE_DIR}')")
    endif()
    message(STATUS "CHOLMOD: ${CHOLMOD_INCLUDE_DIR}, ${PP_cholmod_LIBRARY}")
    target_compile_definitions(pp_core PRIVATE PP_USE_CHOLMOD=1)
    target_include_directories(pp_core SYSTEM PRIVATE "${CHOLMOD_INCLUDE_DIR}")
    target_link_libraries(pp_core PUBLIC ${_pp_cholmod_libs})
    # The shared libraries are found at run time from their install dir.
    get_filename_component(_pp_cholmod_libdir "${PP_cholmod_LIBRARY}" DIRECTORY)
    set(CMAKE_BUILD_RPATH "${CMAKE_BUILD_RPATH};${_pp_cholmod_libdir}")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH};${_pp_cholmod_libdir}")
endif()

#########
# tests #
#########
if(PP_BUILD_TESTS)
    add_executable(testPowerDiagram tests/testPowerDiagram.cpp)
    target_link_libraries(testPowerDiagram PRIVATE pp_core)

    add_executable(testVolumeSolver tests/testVolumeSolver.cpp)
    target_link_libraries(testVolumeSolver PRIVATE pp_core)

    add_executable(testJointDual tests/testJointDual.cpp)
    target_link_libraries(testJointDual PRIVATE pp_core)

    # exploragram's semi-discrete OT is the reference for the volume solve.
    add_executable(testDualSolvers tests/testDualSolvers.cpp)
    target_link_libraries(testDualSolvers PRIVATE pp_core geogram exploragram_ot)

    add_executable(testIntegrators tests/testIntegrators.cpp)
    target_link_libraries(testIntegrators PRIVATE pp_core)
endif()

##############
# benchmarks #
##############
if(PP_BUILD_TESTS)
    add_executable(benchBourne benchmarks/bourne/benchBourne.cpp)
    target_link_libraries(benchBourne PRIVATE pp_core)

    add_executable(checkHessian benchmarks/bourne/checkHessian.cpp)
    target_link_libraries(checkHessian PRIVATE pp_core)

    add_executable(makeSynthetic benchmarks/synthetic/makeSynthetic.cpp)
    target_link_libraries(makeSynthetic PRIVATE pp_core)

    add_executable(runPP2D benchmarks/pp2d/runPP2D.cpp)
    target_link_libraries(runPP2D PRIVATE pp_core)
endif()

###########
# Houdini #
###########
if(PP_BUILD_HOUDINI)
    if(DEFINED ENV{HFS})
        set(PP_HFS "$ENV{HFS}")
    else()
        set(PP_HFS "/Applications/Houdini/Current/Frameworks/Houdini.framework/Versions/Current/Resources")
    endif()
    list(APPEND CMAKE_PREFIX_PATH "${PP_HFS}/toolkit/cmake")
    find_package(Houdini REQUIRED)

    add_library(PowerParticlesHDK SHARED
        src/houdini/domainInput.h
        src/houdini/domainInput.cpp
        src/houdini/feasibleVolumesSOP.h
        src/houdini/feasibleVolumesSOP.cpp
        src/houdini/geometryIO.h
        src/houdini/geometryIO.cpp
        src/houdini/newtonParms.h
        src/houdini/newtonParms.cpp
        src/houdini/plane.h
        src/houdini/plane.cpp
        src/houdini/powerBuilder2DSOP.h
        src/houdini/powerBuilder2DSOP.cpp
        src/houdini/powerDiagramSOP.h
        src/houdini/powerDiagramSOP.cpp
        src/houdini/powerDiagramSolverSOP.h
        src/houdini/powerDiagramSolverSOP.cpp
        src/houdini/tetrahedralizeSOP.h
        src/houdini/tetrahedralizeSOP.cpp
        src/houdini/triangulate2DSOP.h
        src/houdini/triangulate2DSOP.cpp
        src/houdini/plugin.cpp
    )
    target_include_directories(PowerParticlesHDK PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(PowerParticlesHDK PRIVATE pp_core Houdini)

    # Install layout follows ../geomt: the dso goes to <home>/dso and the help
    # card to <home>/help/nodes/sop, where <home> is the Houdini user dir
    # (e.g. ~/Library/Preferences/houdini/20.5). Like
    # houdini_get_default_install_dir() we ask hython for it, but fall back to
    # the standard location when hython cannot run (e.g. no license).
    if(NOT PP_HOUDINI_HOME)
        execute_process(
            COMMAND "${PP_HFS}/bin/hython" -c "print(hou.homeHoudiniDirectory())"
            OUTPUT_VARIABLE _hhome
            RESULT_VARIABLE _hstatus
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        if(NOT _hstatus EQUAL 0 OR NOT IS_DIRECTORY "${_hhome}")
            set(_hver "${Houdini_VERSION_MAJOR}.${Houdini_VERSION_MINOR}")
            if(APPLE)
                set(_hhome "$ENV{HOME}/Library/Preferences/houdini/${_hver}")
            else()
                set(_hhome "$ENV{HOME}/houdini${_hver}")
            endif()
            message(STATUS "hython unavailable, using default Houdini user dir")
        endif()
        set(PP_HOUDINI_HOME "${_hhome}" CACHE PATH
            "Houdini user directory the plugin is installed to" FORCE)
    endif()
    message(STATUS "Houdini install dir: ${PP_HOUDINI_HOME}")

    houdini_configure_target(PowerParticlesHDK INSTDIR "${PP_HOUDINI_HOME}/dso")

    # Help cards (shown by the nodes' "?" button), copied after each build.
    set(_help_dir "${PP_HOUDINI_HOME}/help/nodes/sop")
    add_custom_command(TARGET PowerParticlesHDK POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${_help_dir}"
        COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_CURRENT_SOURCE_DIR}/src/houdini/help/nodes/sop/pp_powerdiagram.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/houdini/help/nodes/sop/pp_triangulate2d.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/houdini/help/nodes/sop/pp_tetrahedralize.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/houdini/help/nodes/sop/pp_powerdiagramsolver.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/houdini/help/nodes/sop/pp_feasiblevolumes.txt"
            "${CMAKE_CURRENT_SOURCE_DIR}/src/houdini/help/nodes/sop/pp_powerbuilder2d.txt"
            "${_help_dir}")
endif()
