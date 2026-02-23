# tests/fortran/CMakeLists.txt — Standalone Fortran tests

add_executable(test_h2_ccpvdz test_h2_ccpvdz.F90)

target_link_libraries(test_h2_ccpvdz PRIVATE ${main_lib} OpenMP::OpenMP_Fortran)

target_include_directories(test_h2_ccpvdz PRIVATE "${CMAKE_BINARY_DIR}/mod")

add_executable(test_h2o_ccpvdz test_h2o_ccpvdz.F90)

target_link_libraries(test_h2o_ccpvdz PRIVATE ${main_lib}
                                              OpenMP::OpenMP_Fortran)

target_include_directories(test_h2o_ccpvdz PRIVATE "${CMAKE_BINARY_DIR}/mod")

add_test(NAME h2_ccpvdz COMMAND $<TARGET_FILE:test_h2_ccpvdz>)
add_test(NAME h2o_ccpvdz COMMAND $<TARGET_FILE:test_h2o_ccpvdz>)

add_executable(test_h2_sonly test_h2_sonly.F90)

target_link_libraries(test_h2_sonly PRIVATE ${main_lib} OpenMP::OpenMP_Fortran)

target_include_directories(test_h2_sonly PRIVATE "${CMAKE_BINARY_DIR}/mod")

add_test(NAME h2_sonly COMMAND $<TARGET_FILE:test_h2_sonly>)

# Direct interface test
add_executable(test_interface_h2 test_interface_h2.F90)
target_link_libraries(test_interface_h2 PRIVATE ${main_lib} OpenMP::OpenMP_Fortran)
target_include_directories(test_interface_h2 PRIVATE "${CMAKE_BINARY_DIR}/mod")
add_test(NAME interface_h2 COMMAND $<TARGET_FILE:test_interface_h2>)

# H2O/cc-pVTZ test - requires F functions
if(ERI_ENABLE_F)
  add_executable(test_h2o_ccpvtz test_h2o_ccpvtz.F90)
  target_link_libraries(test_h2o_ccpvtz PRIVATE ${main_lib} OpenMP::OpenMP_Fortran)
  target_include_directories(test_h2o_ccpvtz PRIVATE "${CMAKE_BINARY_DIR}/mod")
  add_test(NAME h2o_ccpvtz COMMAND $<TARGET_FILE:test_h2o_ccpvtz>)
  # Copy reference file to build directory
  configure_file(fa_h2o_ccpvtz_reference.dat ${CMAKE_CURRENT_BINARY_DIR}/fa_h2o_ccpvtz_reference.dat COPYONLY)
endif()
