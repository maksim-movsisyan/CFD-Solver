cmake_minimum_required(VERSION 3.15)

# ============================================================================
# Default build type (only for single-config generators like Ninja, Makefiles)
# ============================================================================
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type: Debug, Release" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
endif()

# ============================================================================
# ============================================================================
# PROJECT <<CFDSolver>> | Single thread cfd solver for compressible flows
#                       | computations
# ============================================================================
# ============================================================================
project(CFDSolver 
        VERSION 0.1.0
        LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# ============================================================================
# Inter-module optimization (LTO/IPO) for Release builds
# ============================================================================
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)

# ============================================================================
# Output directories
# ============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# ============================================================================
# Library sources
# ============================================================================
file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS "src/*.f90")
list(REMOVE_ITEM LIB_SOURCES 
    "${CMAKE_SOURCE_DIR}/src/A_solver.f90"
    "${CMAKE_SOURCE_DIR}/src/A_mesh_generator.f90"
    "${CMAKE_SOURCE_DIR}/src/A_mg_initialization.f90"
)
add_library(cfd_core STATIC ${LIB_SOURCES})

# ============================================================================
# Compilation options
# ============================================================================
target_compile_options(cfd_core PRIVATE
    $<$<Fortran_COMPILER_ID:GNU>:
        $<$<CONFIG:Release>:-Ofast -march=native -funroll-loops -cpp>
        $<$<CONFIG:Debug>:-g -O0 -fcheck=all -finit-real=nan -Warray-temporaries -fbounds-check -cpp>
    >
)

# ============================================================================
# Turn on inter-module optimization
# ============================================================================
if(ipo_supported)
    set_target_properties(cfd_core PROPERTIES INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

# ============================================================================
# Create executable files
# ============================================================================
function(add_cfd_executable EXE_NAME SOURCE_FILE)
    add_executable(${EXE_NAME} ${SOURCE_FILE})
    target_link_libraries(${EXE_NAME} PRIVATE cfd_core)
endfunction()

add_cfd_executable(solver src/A_solver.f90)			# Main solver executable
add_cfd_executable(mesh_generator src/A_mesh_generator.f90)	# Block-structured mesh generator
add_cfd_executable(mg_initialize src/A_mg_initialization.f90)	# Hybrid initialization prepocessor


message(STATUS "=======================================================")
message(STATUS "CFD Solver Configuration Success!")
message(STATUS "Current Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "=======================================================")
