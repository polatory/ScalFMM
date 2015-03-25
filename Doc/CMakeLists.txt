# add a target to generate API documentation with Doxygen
OPTION( SCALFMM_BUILD_DOC "Set to ON to build the Doxygen documentation " OFF )
IF(SCALFMM_BUILD_DOC)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(
      doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
#    INSTALL(FILES   ${SCALFMM_BINARY_DIR}/Doc/scalfmm.tag
#      DESTINATION doc/
#      )
#    INSTALL(DIRECTORY   ${SCALFMM_BINARY_DIR}/Doc/html
#      DESTINATION doc/
#      )
  endif(DOXYGEN_FOUND)
endif(SCALFMM_BUILD_DOC)
