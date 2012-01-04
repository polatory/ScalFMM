cmake_minimum_required (VERSION 2.6)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

# Options
OPTION( SCALFMM_USE_CBLAS "Set to ON to build ScaFMM with BLAS" OFF )
OPTION( SCALFMM_USE_MPI "Set to ON to build ScaFMM with MPI" ON )
OPTION( SCALFMM_USE_TRACE "Set to ON to print trace or use itac trace" OFF )
OPTION( SCALFMM_BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( SCALFMM_BUILD_UTESTS "Set to ON to build UTests" ON )
OPTION( SCALFMM_BUILD_DEBUG "Set to ON to build in Debug" OFF )
OPTION( SCALFMM_USE_MEM_STATS "Set to ON to profile memory" ON )


# MPI option has to be set before project
if( SCALFMM_USE_MPI )
  SET(CMAKE_CXX_COMPILER mpicxx)
endif()

# Project Infos
project(scalfmm CXX)

# Debug
MESSAGE( STATUS "SCALFMM_BUILD_DEBUG = ${SCALFMM_BUILD_DEBUG}" )
if( SCALFMM_BUILD_DEBUG )
  SET(CMAKE_BUILD_TYPE Debug)
endif()

# Trace
MESSAGE( STATUS "SCALFMM_USE_TRACE = ${SCALFMM_USE_TRACE}" )
if( SCALFMM_USE_TRACE )
  OPTION( SCALFMM_USE_ITAC "Set to ON to use itac trace" OFF )
  if( SCALFMM_USE_ITAC )
    ADD_DEFINITIONS(-I$VT_ROOT/include -trace)
  endif()
endif()

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
ADD_DEFINITIONS(-O2 -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion)

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake 
		${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h  )

# Build - lib
add_subdirectory(Src)

# Build - Tests
MESSAGE( STATUS "SCALFMM_BUILD_TESTS = ${SCALFMM_BUILD_TESTS}" )
if( SCALFMM_BUILD_TESTS )
	add_subdirectory(Tests)
endif()

# Build - UTests
MESSAGE( STATUS "SCALFMM_BUILD_UTESTS = ${SCALFMM_BUILD_UTESTS}" )
if( SCALFMM_BUILD_UTESTS )
	ENABLE_TESTING()
	add_subdirectory(UTests)
endif()

# Use Mem stats
MESSAGE( STATUS "SCALFMM_USE_MEM_STATS = ${SCALFMM_USE_MEM_STATS}" )

# Build - doc
add_subdirectory(Doc)
