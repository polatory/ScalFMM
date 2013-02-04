cmake_minimum_required (VERSION 2.8.3 FATAL_ERROR)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

# Project Infos
project(ScalFMM)
#===========================================================================
# Version Number
#===========================================================================

# ScalFMM version number.  An even minor number corresponds to releases.
SET(ScalFMM_MAJOR_VERSION 1)
SET(ScalFMM_MINOR_VERSION 0)
SET(ScalFMM_BUILD_VERSION 0)
SET(ScalFMM_VERSION  "${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}.${ScalFMM_BUILD_VERSION}"
  )

# MPI option has to be set before project, cannot be changed in the cache!
if( SCALFMM_USE_MPI )
	include(CMakeForceCompiler)
	CMAKE_FORCE_C_COMPILER(mpicc "MPI C Compiler") 
    	CMAKE_FORCE_CXX_COMPILER(mpicxx "MPI C++ Compiler")
else()
	MESSAGE(STATUS "Remove CMake cache and run cmake .. -DSCALFMM_USE_MPI=ON to enable MPI" )
endif(SCALFMM_USE_MPI)
#
# Active language
# -----------------------
ENABLE_LANGUAGE(CXX )
MESSAGE(STATUS " CXX  ${CMAKE_CXX_COMPILER_ID}" )

# Options
OPTION( SCALFMM_USE_MPI "Set to ON to build ScaFMM with MPI" OFF )
OPTION( SCALFMM_USE_BLAS "Set to ON to build ScaFMM with BLAS" ON )
OPTION( SCALFMM_USE_TRACE "Set to ON to print trace or use itac trace" OFF )
OPTION( SCALFMM_BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( SCALFMM_BUILD_UTESTS "Set to ON to build UTests" ON )
OPTION( SCALFMM_BUILD_DEBUG "Set to ON to build in Debug" OFF )
OPTION( SCALFMM_USE_MEM_STATS "Set to ON to profile memory" OFF )
OPTION( SCALFMM_USE_DOUBLE_PRECISION "Set to ON to compile in double precision" ON )
OPTION( SCALFMM_ATTACHE_SOURCE "Set to ON to compile with -g" OFF )
OPTION( SCALFMM_USE_ADDONS "Set to ON to compile add ons" OFF )
#
SET(SCALFMM_LIBRARIES "")

# Test if openmp is here
find_package (OpenMP) 
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Debug
MESSAGE( STATUS "SCALFMM_BUILD_DEBUG = ${SCALFMM_BUILD_DEBUG}" )
if( SCALFMM_BUILD_DEBUG )
  SET(CMAKE_BUILD_TYPE Debug)
  ADD_DEFINITIONS(-O0)
else()
  SET(CMAKE_BUILD_TYPE Release)
  # force -O2 in release
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
  STRING(REPLACE "-O3" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
  # if compiler is intel add -ip
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    ADD_DEFINITIONS(-ip)
  ENDIF()
endif()

# Attach source code to exec
if( SCALFMM_ATTACHE_SOURCE )
  MESSAGE( STATUS "Use -g in compiler flags" )
  ADD_DEFINITIONS(-g)
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
if( SCALFMM_USE_BLAS )
  OPTION( SCALFMM_USE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )
  if( SCALFMM_USE_MKL_AS_BLAS )
     SET(BLAS_LIBRARIES "-L$ENV{MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" CACHE STRING "Set your MKL flags")
  elseif(SCALFMM_USE_EXTERNAL_BLAS) 
     MESSAGE(STATUS "BLAS SET BY EXTERNAL PROGRAM = ${BLAS_LIBRARIES}")
  else()
	INCLUDE(FindBLAS)
	INCLUDE(FindLAPACK)
  endif()
  #SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}")
  SET(SCALFMM_LIBRARIES ${SCALFMM_LIBRARIES} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  #MESSAGE(${SCALFMM_LIBRARIES})
endif()

# Compile option
ADD_DEFINITIONS(-Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion)

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake 
		${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h  )

# Use Mem stats
MESSAGE( STATUS "SCALFMM_USE_MEM_STATS = ${SCALFMM_USE_MEM_STATS}" )

# Add CBLAS
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${CBLAS_LIBRARIES}")

# Build - lib
add_subdirectory(Src)

# Link with scalfmm lib
set(scalfmm_lib scalfmm)

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

# Build - doc
add_subdirectory(Doc)

# Build - Addons
if( SCALFMM_USE_ADDONS )
	file(GLOB all_addons ./Addons/*)

	foreach(addon_name ${all_addons}) 
		MESSAGE(STATUS "Addons += ${addon_name}" )
		add_subdirectory(${addon_name})
	endforeach()
endif()

# Add custom clean command
if("${CMAKE_GENERATOR}" MATCHES "Make")
	ADD_CUSTOM_TARGET(reset
		COMMAND ${CMAKE_MAKE_PROGRAM} clean && rm -r ./CMake* && cmake ..
		WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
		COMMENT "Cleaning Build dir and CMake cache, then do cmake.."
		)
endif()
#
#
#===========================================================================
#    Build and export cmake files 
#===========================================================================
#
MESSAGE(WARNING "ScalFMM_BINARY_DIR: "${Scalfmm_BINARY_DIR})
CONFIGURE_FILE(${ScalFMM_SOURCE_DIR}/ScalFMMConfig.cmake.in
               ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake 
               @ONLY IMMEDIATE
               )
INSTALL(FILES 
  ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake
  DESTINATION lib/
  )
