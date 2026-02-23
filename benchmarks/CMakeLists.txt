# benchmarks/CMakeLists.txt — Performance benchmarks (not run by ctest)

add_executable(bench_c2h6 bench_c2h6.F90)
target_link_libraries(bench_c2h6 PRIVATE ${main_lib} OpenMP::OpenMP_Fortran)
if(MATH_LIBRARY)
  target_link_libraries(bench_c2h6 PRIVATE ${MATH_LIBRARY})
endif()
target_include_directories(bench_c2h6 PRIVATE "${CMAKE_BINARY_DIR}/mod")

add_executable(bench_c60 bench_c60.F90)
target_link_libraries(bench_c60 PRIVATE ${main_lib} OpenMP::OpenMP_Fortran)
if(MATH_LIBRARY)
  target_link_libraries(bench_c60 PRIVATE ${MATH_LIBRARY})
endif()
target_include_directories(bench_c60 PRIVATE "${CMAKE_BINARY_DIR}/mod")

# No add_test() — these are manual benchmarks
