
target_sources(
  ${main_lib}
  PRIVATE 
          ${CMAKE_CURRENT_SOURCE_DIR}/liberi_interface.F90
          ${CMAKE_CURRENT_SOURCE_DIR}/gms_liberi_wrapper.F90
          )

if(ERI_ENABLE_MDI)
  target_sources(${main_lib} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/liberi_mdi_api.F90)
endif()
