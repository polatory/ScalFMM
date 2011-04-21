cmake_minimum_required (VERSION 2.6)

# Project Infos
project(scalfmm CXX)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
IF(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
ENDIF(insource)
#
# Options
OPTION( SCALFMM_USE_CBLAS "Set to ON to build ScaFMM with BLAS" ON )
OPTION( BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( BUILD_UTESTS "Set to ON to build UTests" ON )
#
if(SCALFMM_USE_CBLAS )
OPTION( SCALFMM_USE_MKL "Set to ON to build MKL CBLAS" OFF )
  IF (    SCALFMM_USE_MKL ) 
     SET(CBLAS_LIBRARIES "-L$ENV{MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
  ELSE (    SCALFMM_USE_MKL ) 
      FIND_PACKAGE(BLAS)
      SET(CBLAS_LIBRARIES "-lcblas")
  ENDIF (    SCALFMM_USE_MKL)
ENDIF(SCALFMM_USE_CBLAS)
#
ADD_DEFINITIONS(-O2 -Wall -Wshadow -Wpointer-arith -Wcast-qual)

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/Sources/ScalFMM_config.h.cmake 
		${CMAKE_BINARY_DIR}/Sources/ScalFMM_config.h  )
# Build - lib
add_subdirectory(Sources)

# Build - Tests
MESSAGE( STATUS "BUILD_TESTS = ${BUILD_TESTS}" )
if( ${BUILD_TESTS} STREQUAL "ON" )
	add_subdirectory(Tests)
endif()

# Build - UTests
MESSAGE( STATUS "BUILD_UTESTS = ${BUILD_UTESTS}" )
if( ${BUILD_UTESTS} STREQUAL "ON" )
	add_subdirectory(UTests)
endif()

# Build - doc
add_subdirectory(Doc)
