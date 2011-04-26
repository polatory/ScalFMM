cmake_minimum_required (VERSION 2.6)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

# Options
OPTION( SCALFMM_USE_CBLAS "Set to ON to build ScaFMM with BLAS" ON )
OPTION( SCALFMM_USE_MPI "Set to ON to build ScaFMM with MPI" OFF )
OPTION( BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( BUILD_UTESTS "Set to ON to build UTests" ON )

# MPI option
if( SCALFMM_USE_MPI )
  SET(CMAKE_CXX_COMPILER mpicxx)
endif()

# Project Infos
project(scalfmm CXX)

# Blas option
if( SCALFMM_USE_CBLAS )
  OPTION( SCALFMM_USE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )
  if( SCALFMM_USE_MKL_AS_BLAS )
     SET(CBLAS_LIBRARIES "-L$ENV{MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
  else() 
      FIND_PACKAGE(BLAS)
      SET(CBLAS_LIBRARIES "-lcblas")
  endif()
endif()

# Compile option
ADD_DEFINITIONS(-O2 -Wall -Wshadow -Wpointer-arith -Wcast-qual)


CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake 
		${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h  )

# Build - lib
add_subdirectory(Src)

# Build - Tests
MESSAGE( STATUS "BUILD_TESTS = ${BUILD_TESTS}" )
if( BUILD_TESTS )
	add_subdirectory(Tests)
endif()

# Build - UTests
MESSAGE( STATUS "BUILD_UTESTS = ${BUILD_UTESTS}" )
if( BUILD_UTESTS )
	add_subdirectory(UTests)
endif()

# Build - doc
add_subdirectory(Doc)
