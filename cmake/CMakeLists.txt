# cmake/CMakeLists.txt — Per-compiler flag configuration

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Build type" FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
endif()

# =============================================================================
# Manual Flag Override
# =============================================================================
set(ERI_FC_EXTRA_FLAGS
    ""
    CACHE STRING "Override all Fortran flags (leave empty for auto-detection)")

if(NOT ERI_FC_EXTRA_FLAGS STREQUAL "")
  message(STATUS "Using manual Fortran flags: ${ERI_FC_EXTRA_FLAGS}")
  set(CMAKE_Fortran_FLAGS
      "${ERI_FC_EXTRA_FLAGS}"
      PARENT_SCOPE)
  return()
endif()

# =============================================================================
# GPU flags (composed from ERI options)
# =============================================================================
set(FLAGS_GPU "")

if(ERI_ENABLE_GPU)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(FLAGS_GPU
        "-mp=gpu -target=gpu -acc -stdpar=gpu -gpu=${ERI_GPU_ARCH} -gpu=mem:separate"
    )
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID
                                                       STREQUAL "IntelLLVM")
    set(FLAGS_GPU
        "-fiopenmp -fopenmp-targets=spir64_gen -Xopenmp-target-backend \"-device 12.60.7\""
    )
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    set(FLAGS_GPU
        "-mno-unsafe-fp-atomics -hacc_model=auto_async_none:no_fast_addr:no_deep_copy"
    )
  endif()
endif()

# =============================================================================
# Compiler Flags Configuration (auto-detected)
# =============================================================================

if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  set(FLAGS_STANDARD "-mcmodel=medium -fPIC ${FLAGS_GPU}")
  set(FLAGS_DEBUG "-Mbounds -Mchkptr -Ktrap=fp -g")
  set(FLAGS_RELEASE "-O3 -fast")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(FLAGS_STANDARD
      "-ffree-line-length-none -fbacktrace -fPIC -fallow-argument-mismatch")
  set(FLAGS_DEBUG "-ffpe-trap=invalid,zero,overflow -fcheck=all")
  set(FLAGS_RELEASE "-march=native")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID
                                                     STREQUAL "IntelLLVM")
  set(FLAGS_STANDARD "-fPIC -traceback ${FLAGS_GPU}")
  set(FLAGS_DEBUG "-check all -fpe0")
  set(FLAGS_RELEASE "-O3")

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(FLAGS_STANDARD "-eZ -h omp -rm ${FLAGS_GPU}")
  set(FLAGS_DEBUG "-m 3")
  set(FLAGS_RELEASE "-O3")

else()
  message(WARNING "Unknown Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
  set(FLAGS_STANDARD "")
  set(FLAGS_DEBUG "")
  set(FLAGS_RELEASE "")
endif()

# =============================================================================
# Apply flags to CMake variables (propagate to parent scope)
# =============================================================================
set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} ${FLAGS_STANDARD}"
    PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} ${FLAGS_DEBUG}"
    PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS_RELEASE
    "${CMAKE_Fortran_FLAGS_RELEASE} ${FLAGS_RELEASE}"
    PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO} ${FLAGS_RELEASE}"
    PARENT_SCOPE)
