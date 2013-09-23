 cmake_minimum_required (VERSION 2.8.3 FATAL_ERROR)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

#
# MPI option has to be set before project, cannot be changed in the cache!
if( ScalFMM_USE_MPI )
	include(CMakeForceCompiler)
	CMAKE_FORCE_C_COMPILER(mpicc "MPI C Compiler") 
    	CMAKE_FORCE_CXX_COMPILER(mpicxx "MPI C++ Compiler")
	set(ScalFMM_USE_MPI ON CACHE BOOL "ScalFmm use MPI")
else()
	MESSAGE(STATUS "Remove CMake cache and run cmake .. -DScalFMM_USE_MPI=ON to enable MPI" )
endif(ScalFMM_USE_MPI)

#
# Project Infos
project(ScalFMM)
#===========================================================================
# Version Number
#===========================================================================

# ScalFMM version number.  An even minor number corresponds to releases.
SET(ScalFMM_MAJOR_VERSION 1)
SET(ScalFMM_MINOR_VERSION 1)
SET(ScalFMM_PATCH_VERSION 1)
SET(ScalFMM_VERSION  "${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}.${ScalFMM_PATCH_VERSION}" )

#
# Active language
# -----------------------
ENABLE_LANGUAGE(CXX )
MESSAGE(STATUS " CXX  ${CMAKE_CXX_COMPILER_ID}" )

#
# Options
OPTION( ScalFMM_USE_BLAS "Set to ON to build ScaFMM with BLAS" OFF )
OPTION( ScalFMM_USE_TRACE "Set to ON to print trace or use itac trace" OFF )
OPTION( ScalFMM_BUILD_TESTS "Set to ON to build fonctionnalities Tests" OFF )
OPTION( ScalFMM_BUILD_UTESTS "Set to ON to build UTests" OFF )
OPTION( ScalFMM_BUILD_DEBUG          "Set to ON to build in Debug"              OFF )
OPTION( ScalFMM_USE_MEM_STATS        "Set to ON to profile memory"              OFF )
OPTION( ScalFMM_USE_DOUBLE_PRECISION "Set to ON to compile in double precision" ON  )
OPTION( ScalFMM_ATTACHE_SOURCE       "Set to ON to compile with -g"             OFF )
OPTION( ScalFMM_USE_ADDONS           "Set to ON to compile add ons"             OFF )
OPTION( ScalFMM_USE_SSE              "Set to ON to compile with sse support"    ON  )

# Set scalfmm to default libraries
SET(SCALFMM_LIBRARIES "")

# Test if openmp is here
find_package (OpenMP) 
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS            "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS          "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS   "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Debug
MESSAGE( STATUS "ScalFMM_BUILD_DEBUG          = ${ScalFMM_BUILD_DEBUG}" )
if( ScalfMM_BUILD_DEBUG )
  SET(CMAKE_BUILD_TYPE Debug)
  ADD_DEFINITIONS(-O0)

  # Compile optimization
  IF( APPLE ) 
    SET(SCALFMM_FLAGS_OPTI_DEBUG "-m64  -funroll-loops" CACHE STRING "Set your optimization flags for debug mode.")
  ELSE()
    # Test if not apple and 64bits
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
       SET(SCALFMM_FLAGS_OPTI_DEBUG "-m64 -march=native -funroll-loops" CACHE STRING "Set your optimization flags for debug mode.")
    else()
       SET(SCALFMM_FLAGS_OPTI_DEBUG "-march=native -funroll-loops" CACHE STRING "Set your optimization flags for debug mode.")
    endif()
  ENDIF()
  ADD_DEFINITIONS(${SCALFMM_FLAGS_OPTI_DEBUG})
else()
  SET(CMAKE_BUILD_TYPE Release)
  # force -O2 in release
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
  STRING(REPLACE "-O3" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
  # if compiler is intel add -ip
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    ADD_DEFINITIONS(-ip)
  ENDIF()

  # Compile optimization
  IF( APPLE ) 
    SET(SCALFMM_FLAGS_OPTI_RELEASE "-m64 -Ofast  -funroll-loops" CACHE STRING "Set your optimization flags for release mode.")
  ELSE()
    # Test if not apple and 64bits
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
       SET(SCALFMM_FLAGS_OPTI_RELEASE "-m64 -Ofast -flto -march=native -funroll-loops" CACHE STRING "Set your optimization flags for release mode.")
    else()
       SET(SCALFMM_FLAGS_OPTI_RELEASE "-Ofast -flto -march=native -funroll-loops" CACHE STRING "Set your optimization flags for release mode.")
    endif()
  ENDIF()
  ADD_DEFINITIONS(${SCALFMM_FLAGS_OPTI_RELEASE})
endif()

# Attach source code to exec
MESSAGE(STATUS "ScalFMM_ATTACHE_SOURCE       = ${ScalFMM_ATTACHE_SOURCE}" )
if( ScalFMM_ATTACHE_SOURCE )
  MESSAGE( STATUS "Option -g is used, the code is attached to the binary." )
  ADD_DEFINITIONS(-g)
endif()

# Trace
MESSAGE( STATUS "ScalFMM_USE_TRACE            = ${ScalFMM_USE_TRACE}" )
if( ScalFMM_USE_TRACE )
  OPTION( ScalFMM_USE_ITAC "Set to ON to use itac trace" OFF )
  if( ScaLFMM_USE_ITAC )
    ADD_DEFINITIONS(-I$VT_ROOT/include -trace)
  endif()
endif()

# Blas option
if( ScalFMM_USE_BLAS )
  OPTION( ScalFMM_USE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )
  if( ScalFMM_USE_MKL_AS_BLAS )
     SET(BLAS_LIBRARIES "-L$ENV{MKLROOT}/lib  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core" CACHE STRING "Set your MKL flags")
  elseif(ScalFMM_USE_EXTERNAL_BLAS) 
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
ADD_DEFINITIONS(-Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion -fpic )
# 
# C++ 2011
ADD_DEFINITIONS(-std=c++11)
IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND APPLE)
  ADD_DEFINITIONS(-stdlib=libstdc++)
endif()

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake 
		${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h  )


# Use SSE
MESSAGE( STATUS "ScalFMM_USE_SSE              = ${ScalFMM_USE_SSE}" )
if( ScalFMM_USE_SSE )
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -axAVX")
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -mfpmath=sse")
  ENDIF()
endif()

# Use Mem stats
MESSAGE( STATUS "ScalFMM_USE_MEM_STATS        = ${ScalFMM_USE_MEM_STATS}" )

# Add CBLAS
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${CBLAS_LIBRARIES}")

# Build - lib
add_subdirectory(Src)

# Link with scalfmm lib
set(scalfmm_lib scalfmm)

# Build - Tests
MESSAGE( STATUS "ScalFMM_BUILD_TESTS          = ${ScalFMM_BUILD_TESTS}" )
if( ScalFMM_BUILD_TESTS )
	add_subdirectory(Tests)
endif()

# Build - UTests
MESSAGE( STATUS "ScalFMM_BUILD_UTESTS         = ${ScalFMM_BUILD_UTESTS}" )
if( ScalFMM_BUILD_UTESTS )
	ENABLE_TESTING()
	add_subdirectory(UTests)
endif()

# Build - doc
add_subdirectory(Doc)

# Build - Addons
if( ScalFMM_USE_ADDONS )
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
MESSAGE(STATUS "ScalFMM_BINARY_DIR:          " ${ScalFMM_BINARY_DIR})
CONFIGURE_FILE(${ScalFMM_SOURCE_DIR}/ScalFMMConfig.cmake.in
               ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake 
               @ONLY IMMEDIATE
               )
INSTALL(FILES 
  ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake
  DESTINATION lib/
  )
#
# Add CPACK Module
# ----------------
#
# build a CPack driven installer package
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/Licence.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "${ScalFMM_MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_NINOR "${ScalFMM_MINOR_VERSION}")
set(CPACK_PACKAGE_VERSION_PATCH	"${ScalFMM_PATCH_VERSION}")
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY "ON")
#
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME  "ScalFMM-${ScalFMM_VERSION}") 
SET(CPACK_SOURCE_IGNORE_FILES "\\\\.git;.DS_Store;.*~;/*.aux;/*.idx;/*.log;/*.out;/*.toc;/*.ilg;/*.ind;scalfmm.pro*")
LIST(APPEND CPACK_SOURCE_IGNORE_FILES "${CMAKE_BINARY_DIR}")
#
INCLUDE(CPack)
