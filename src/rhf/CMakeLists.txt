# src/rhf/CMakeLists.txt

add_subdirectory(rot_axis)

if(ERI_ENABLE_F)
  add_subdirectory(eric)
  add_subdirectory(rys)
endif()
