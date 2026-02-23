function(print_summary)
  message(STATUS "")
  message(STATUS "========== libERI Configuration Summary ==========")
  message(
    STATUS
      "  Fortran compiler   : ${CMAKE_Fortran_COMPILER_ID} (${CMAKE_Fortran_COMPILER})"
  )
  message(
    STATUS "  C compiler         : ${CMAKE_C_COMPILER_ID} (${CMAKE_C_COMPILER})"
  )
  message(STATUS "  Build type         : ${CMAKE_BUILD_TYPE}")
  message(STATUS "  GPU enabled        : ${ERI_ENABLE_GPU}")
  if(ERI_ENABLE_GPU)
    message(STATUS "  GPU arch           : ${ERI_GPU_ARCH}")
  endif()
  message(STATUS "  F-functions        : ${ERI_ENABLE_F}")
  message(STATUS "  MDI                : ${ERI_ENABLE_MDI}")
  if(ERI_ENABLE_MDI)
    message(STATUS "  MDI MPI support    : ${ERI_USE_MPI}")
  endif()
  message(STATUS "  Testing            : ${ERI_ENABLE_TESTING}")
  message(STATUS "  Fortran flags      : ${CMAKE_Fortran_FLAGS}")
  message(STATUS "==================================================")
  message(STATUS "")
endfunction()
