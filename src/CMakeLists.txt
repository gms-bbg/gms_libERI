# src/CMakeLists.txt — Core library sources

target_sources(
  ${main_lib}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/liberi_parameters.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/liberi_types.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/liberi_constants.f90
          ${CMAKE_CURRENT_SOURCE_DIR}/liberi_driver.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/liberi_shell_pair.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/liberi_boys.F90)

add_subdirectory(rhf)
add_subdirectory(interface)
