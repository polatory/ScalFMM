cmake_minimum_required (VERSION 2.6)

# Project Infos
project(scalfmm CXX)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

# Options
OPTION( FUSE_CBLAS "Set to ON to build ScaFMM with BLAS" ON )
OPTION( FUSE_MPI "Set to ON to build ScaFMM with MPI" OFF )
OPTION( BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( BUILD_UTESTS "Set to ON to build UTests" ON )

# Blas option
if( FUSE_CBLAS )
  OPTION( FUSE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )
  if( FUSE_MKL_AS_BLAS )
     SET(CBLAS_LIBRARIES "-L$ENV{MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
  else( FUSE_MKL_AS_BLAS ) 
      FIND_PACKAGE(BLAS)
      SET(CBLAS_LIBRARIES "-lcblas")
  endif( FUSE_MKL_AS_BLAS )
endif( FUSE_CBLAS )

# MPI option
if( FUSE_MPI )
  SET(CMAKE_CXX_COMPILER mpicxx)
endif(FUSE_MPI)


# Compile option
ADD_DEFINITIONS(-O2 -Wall -Wshadow -Wpointer-arith -Wcast-qual)

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/Sources/ScalFmmConfig.h.cmake 
		${CMAKE_BINARY_DIR}/Sources/ScalFmmConfig.h  )

# Build - lib
add_subdirectory(Sources)

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
