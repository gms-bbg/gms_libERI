cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0144 NEW)

project(
  libERI
  LANGUAGES Fortran C
  VERSION 0.1.0)

# =============================================================================
# Options
# =============================================================================
set(ERI_GPU_ARCH
    "cc70"
    CACHE STRING "NVIDIA GPU compute capability (e.g. cc70, cc80)")
option(ERI_ENABLE_F "Enable F-function integrals (eric + rys)" OFF)
option(ERI_USE_MPI "MPI support in MDI" ON)
option(ERI_ENABLE_TESTING "Build driver executable" ON)
option(ERI_ENABLE_BENCHMARKS "Build performance benchmarks" OFF)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(_ERI_GPU_DEFAULT OFF)
else()
  set(_ERI_GPU_DEFAULT ON)
endif()
option(ERI_ENABLE_GPU "Enable GPU offloading flags" ${_ERI_GPU_DEFAULT})
option(ERI_ENABLE_MDI "Enable MDI library support" ON)
option(ERI_USE_MPI "MPI support in MDI (only used if ERI_ENABLE_MDI=ON)" ON)
option(ERI_ENABLE_TESTING "Build driver executable" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

# =============================================================================
# Module path for Find*.cmake
# =============================================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# =============================================================================
# Dependencies
# =============================================================================
find_package(
  OpenMP
  COMPONENTS Fortran
  REQUIRED)
if(ERI_ENABLE_MDI)
  find_package(
    MPI
    COMPONENTS Fortran
    REQUIRED)
endif()

# =============================================================================
# Compiler flags (must come before targets)
# =============================================================================
add_subdirectory(cmake)

if(ERI_ENABLE_MDI)
  find_package(MDI REQUIRED)
  set(LIBRARY_TYPE "SHARED")
else()
  set(LIBRARY_TYPE "STATIC")
endif()

# =============================================================================
# Main library target
# =============================================================================
set(main_lib ERI)
add_library(${main_lib} ${LIBRARY_TYPE})

# Populate sources from subdirectories
add_subdirectory(src)

# =============================================================================
# Compile definitions — maps options to preprocessor guards
# =============================================================================

# F-functions
if(ERI_ENABLE_F)
  target_compile_definitions(${main_lib} PUBLIC ENABLE_F)
endif()

# GPU offloading (OpenMP target)
if(ERI_ENABLE_GPU)
  target_compile_definitions(${main_lib} PUBLIC ENABLE_GPU)
endif()

# MDI support
if(ERI_ENABLE_MDI)
  target_compile_definitions(${main_lib} PUBLIC USE_MDI)
endif()

# =============================================================================
# Link dependencies
# =============================================================================
target_link_libraries(${main_lib} PRIVATE OpenMP::OpenMP_Fortran
                                          ${CMAKE_DL_LIBS})
if(ERI_ENABLE_MDI)
  target_link_libraries(${main_lib} PRIVATE MPI::MPI_Fortran mdi::mdi)
endif()

# =============================================================================
# Fortran module output directory
# =============================================================================
set_target_properties(${main_lib} PROPERTIES Fortran_MODULE_DIRECTORY
                                             "${CMAKE_BINARY_DIR}/mod")
target_include_directories(
  ${main_lib} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mod>
                     $<INSTALL_INTERFACE:include/libERI>)

# =============================================================================
# Driver executable
# =============================================================================
if(ERI_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests/fortran)
endif()

# =============================================================================
# Performance benchmarks (manual, not ctest)
# =============================================================================
if(ERI_ENABLE_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

# =============================================================================
# Install
# =============================================================================
include(CMakePackageConfigHelpers)

install(
  TARGETS ${main_lib}
  EXPORT libERITargets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(DIRECTORY "${CMAKE_BINARY_DIR}/mod/" DESTINATION "include/libERI")

install(
  EXPORT libERITargets
  FILE libERITargets.cmake
  NAMESPACE libERI::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/libERI")

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/libERIConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/libERIConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/libERI")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/libERIConfig.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/libERI")

# =============================================================================
# Configuration summary
# =============================================================================
include(cmake/modules/eri_status.cmake)
print_summary()
