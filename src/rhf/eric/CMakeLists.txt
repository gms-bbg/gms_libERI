# src/rhf/eric/CMakeLists.txt — ERIC integral kernels

target_sources(
  ${main_lib}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/eric_kernels.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int0030_ericgen.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int0031_ericgen.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int0032_ericgen.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int1030_ericgen.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int1031_ericgen.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int1130_ericgen.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/int2030_ericgen.F90)
