# FindMDI.cmake — FetchContent wrapper for the MDI library

set(_lib "mdi")
set(_url "https://github.com/MolSSI-MDI/MDI_Library.git")
set(_rev "v1.4.35")

# MDI is a static lib linked into our shared lib — it must be built with -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Configure MDI build options before FetchContent
set(language
    "Fortran"
    CACHE STRING "" FORCE)
set(libtype
    "STATIC"
    CACHE STRING "" FORCE)
set(mpi
    "${ERI_USE_MPI}"
    CACHE STRING "" FORCE)
set(plugins
    "ON"
    CACHE STRING "" FORCE)
set(python_plugins
    "OFF"
    CACHE STRING "" FORCE)

include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")

my_fetch_package("${_lib}" "${_url}" "${_rev}")

# MDI's own CMakeLists doesn't export the directory containing mdi.mod. The
# Fortran module file lands in ${mdi_BINARY_DIR}/MDI_Library/, so we must expose
# that path to any target that links against mdi.
target_include_directories(
  mdi INTERFACE $<BUILD_INTERFACE:${mdi_BINARY_DIR}/MDI_Library>)

unset(_lib)
unset(_url)
unset(_rev)
