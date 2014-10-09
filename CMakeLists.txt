 cmake_minimum_required (VERSION 2.8.3 FATAL_ERROR)
if(CMAKE_MAJOR_VERSION GREATER 1 AND CMAKE_MINOR_VERSION GREATER 7 AND CMAKE_PATCH_VERSION GREATER 11)
  cmake_policy(SET CMP0022 OLD)
endif()
cmake_policy(SET CMP0004 OLD)
# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)
SET(CMAKE_MODULE_PATH  ${CMAKE_SOURCE_DIR}/CMakeModules)

#
# MPI option has to be set before project, cannot be changed in the cache!
if( ScalFMM_USE_MPI )
	include(CMakeForceCompiler)
	CMAKE_FORCE_C_COMPILER(mpicc         "MPI C Compiler") 
    	CMAKE_FORCE_CXX_COMPILER(mpicxx   "MPI C++ Compiler")
	set(ScalFMM_USE_MPI ON CACHE BOOL   "ScalFmm use MPI")
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
SET(ScalFMM_MINOR_VERSION 3)
SET(ScalFMM_PATCH_VERSION 0)
SET(ScalFMM_VERSION  "${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}.${ScalFMM_PATCH_VERSION}" )

#
# Active language
# -----------------------
ENABLE_LANGUAGE(CXX )
#
if( ScalFMM_USE_MPI )
    TRY_COMPILE(COMPILE_INTEL  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_MODULE_PATH}/compileTestIntel.cpp COMPILE_DEFINITIONS  "${CMAKE_CXX_FLAGS}")    
    if (COMPILE_INTEL) 
       SET(CMAKE_CXX_COMPILER_ID  "Intel")
        endif()
        endif()
MESSAGE(STATUS " CXX  ${CMAKE_CXX_COMPILER_ID}" )

#
# Options
OPTION( ScalFMM_USE_BLAS 		"Set to ON to build ScaFMM with BLAS" OFF )
OPTION( ScalFMM_USE_FFT  		"Set to ON to build ScaFMM with FFTW" OFF )
OPTION( ScalFMM_BUILD_TESTS 		"Set to ON to build functionalities Tests" OFF )
OPTION( ScalFMM_BUILD_UTESTS 		"Set to ON to build UTests" OFF )
OPTION( ScalFMM_BUILD_DEBUG          "Set to ON to build in Debug"                 OFF )
OPTION( ScalFMM_USE_LOG              "Set to ON to print output debug information" OFF )
OPTION( ScalFMM_USE_MEM_STATS        "Set to ON to profile memory"                 OFF )
OPTION( ScalFMM_USE_DOUBLE_PRECISION "Set to ON to compile in double precision"    ON  )
OPTION( ScalFMM_ATTACHE_SOURCE       "Set to ON to compile with -g"                OFF )
OPTION( ScalFMM_USE_ADDONS           "Set to ON to compile add ons"                OFF )
OPTION( ScalFMM_USE_SSE              "Set to ON to compile with SSE support"       ON  )
OPTION( ScalFMM_USE_AVX              "Set to ON to compile with AVX support"       OFF  )
OPTION( ScalFMM_USE_ASSERT           "Set to ON to enable safe tests during execution" ON  )
OPTION( ScalFMM_USE_MIC_NATIVE       "Set to ON to compile in native mode for MIC" OFF  )
OPTION( ScalFMM_BUILD_ONLY_LIB        "Set to ON to compile only the lib (examples are not compiled) " OFF  )
#OPTION( ScalFMM_ONLY_DEVEL             "Set to ON to compile Development tools (only scalfmm team)" OFF  )
# Set scalfmm to default libraries
SET(SCALFMM_LIBRARIES "")
SET(ScaLFMM_CXX_FLAGS "")
#
#
# Test if openmp is here
#
find_package (OpenMP) 
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS                  "${CMAKE_C_FLAGS}     ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS              "${CMAKE_CXX_FLAGS}  ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS   "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")	
# 	 SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	else()
		MESSAGE(WARNING "OPENMP NOT FOUND")
endif(OPENMP_FOUND)
MESSAGE(STATUS "ScalFMM_BUILD_DEBUG          = ${ScalFMM_BUILD_DEBUG}" )
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  SET(FLAGS_64bits "-m64")
else()
  SET(FLAGS_64bits "")
endif()
##############################################################################
#                                                       Compile options                                                                        #
##############################################################################
#  -xHost -mfpmath=sse
# -Wall Wnosign-conversion 
#
SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS}  -std=c++11  -fpic  -Wall ") 
#
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # INTEL
    IF (APPLE)
        SET(SSE_FLAGS  "-msse4  -mfpmath=sse")   # -mtune=native -march=native
    ELSE(APPLE)
        SET(AVX_FLAGS  "-mavx")
        SET(SSE_FLAGS  "-axSSE4.2")
    ENDIF(APPLE)
    #-Wshadow -Wpointer-arith -Wcast-qual -Wconversion  -Wall -Wnosign-conversion ")
else()
    # NOT      INTEL
    if(NOT ${ScalFMM_USE_MPI})
        SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS}   -Wzero-as-null-pointer-constant")
	else()
		SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS}   -Wno-literal-suffix")
    endif()
    SET(AVX_FLAGS "-mtune=native -march=native")
IF (APPLE)
    SET(SSE_FLAGS  "-msse3  -mfpmath=sse")   # -mtune=native -march=native
ELSE(APPLE)
    SET(SSE_FLAGS  "-mtune=native -march=native")   # -mtune=native -march=native
 ENDIF(APPLE)
endif()
##############################################################################
#                                                        Debug      flags                                                                        #
##############################################################################
#
if( ScalFMM_BUILD_DEBUG )
#
# Compile Debug flags
#
  SET(CMAKE_BUILD_TYPE Debug)
  ADD_DEFINITIONS(-O0)
  IF( APPLE ) 
    SET(SCALFMM_FLAGS_OPTI_DEBUG "-m64  -funroll-loops" CACHE STRING "Set your optimization flags for debug mode.")
  ELSE(APPLE)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        SET(SCALFMM_FLAGS_OPTI_DEBUG "${FLAGS_64bits} -fp-model strict  -funroll-loops" CACHE STRING "Set your optimization flags for debug mode.")
    else()
        SET(SCALFMM_FLAGS_OPTI_DEBUG "${FLAGS_64bits} -funroll-loops" CACHE STRING "Set your optimization flags for debug mode.")
    endif()
  ENDIF(APPLE)
#  ADD_DEFINITIONS(${SCALFMM_FLAGS_OPTI_DEBUG})
  SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS} ${SCALFMM_FLAGS_OPTI_DEBUG}")
else()
##############################################################################
#                                                            Release flags                                                                        #
##############################################################################
  SET(CMAKE_BUILD_TYPE Release)
  # force -O3 in release
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
  STRING(REPLACE "-O2" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
  # if compiler is intel add -ip
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS} -ip -no-prec-div -no-prec-sqrt")
  ENDIF()
  # Compile optimization
  IF( APPLE ) 
    SET(SCALFMM_FLAGS_OPTI_RELEASE "-m64  -funroll-loops" CACHE STRING "Set your optimization flags for release mode.")
  ELSE(APPLE)
    # Not apple system - Check the compiler flags 
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      SET(SCALFMM_FLAGS_OPTI_RELEASE "${FLAGS_64bits} -fp-model precise -fp-model source -fimf-precision=low  -funroll-loops -ftree-vectorize" CACHE STRING "Set your optimization flags for release mode.")
    else()
      SET(SCALFMM_FLAGS_OPTI_RELEASE "${FLAGS_64bits} -ffast-math -funroll-loops -ftree-vectorize" CACHE STRING "Set your optimization flags for release mode.")
    endif()
  ENDIF(APPLE)
#
  SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS} ${SCALFMM_FLAGS_OPTI_RELEASE}")
endif()
##############################################################################
#                                           Attach source code to exec                                                                   #
##############################################################################
MESSAGE(STATUS "ScalFMM_ATTACHE_SOURCE       = ${ScalFMM_ATTACHE_SOURCE}" )
if( ScalFMM_ATTACHE_SOURCE )
  MESSAGE( STATUS "Option -g is used, the code is attached to the binary." )
#  ADD_DEFINITIONS(-g)
    SET(ScaLFMM_CXX_FLAGS  "${ScaLFMM_CXX_FLAGS} -g")
endif(ScalFMM_ATTACHE_SOURCE)

##############################################################################
#                                               Blas option                                                                                       #  
##############################################################################
# 
if( ScalFMM_USE_BLAS )
  OPTION( ScalFMM_USE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )
  MESSAGE(STATUS "CMAKE_CXX_COMPILER_ID STREQUAL  ${CMAKE_CXX_COMPILER_ID}")
  if( ScalFMM_USE_MKL_AS_BLAS )
        IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
              SET(BLAS_LIBRARIES  "-mkl=sequential" CACHE STRING "Set your MKL flags")
        else()
              SET(BLAS_LIBRARIES
                       "-L$ENV{MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_sequential;-lmkl_core" CACHE STRING "Set your MKL flags")
         ENDIF()
         UNSET(LAPACK_LIBRARIES)
         SET(BLASLAPACK_LIBRARIES ${BLAS_LIBRARIES})
  elseif(ScalFMM_USE_EXTERNAL_BLAS) 
         MESSAGE(STATUS "BLAS SET BY EXTERNAL PROGRAM = ${BLAS_LIBRARIES}")
  else()
	INCLUDE(FindBLAS)
	if (NOT BLAS_FOUND)
	   MESSAGE(FATAL_ERROR : "NO BLAS LIBRARY FOUND")
	endif()
	INCLUDE(FindLAPACK)
	SET(BLASLAPACK_LIBRARIES "${BLAS_LIBRARIES};${LAPACK_LIBRARIES}")
  endif()
  SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};${BLASLAPACK_LIBRARIES}")
  MESSAGE(STATUS "SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")
endif(ScalFMM_USE_BLAS)
#
# FFT option
#
if( ScalFMM_USE_FFT )
  OPTION( ScalFMM_USE_MKL_AS_FFTW "Set to ON to use MKL FFTW" OFF )
  if( ScalFMM_USE_MKL_AS_FFTW )
      if( ScalFMM_USE_MKL_AS_BLAS )
               UNSET(FFT_LIBRARIES CACHE)
       #        UNSET(FFT_INCLUDES CACHE)
                SET(FFT_INCLUDES "$ENV{MKLROOT}/include/fftw" CACHE STRING "Set your MKL flags")
        ELSE()
             SET(FFT_LIBRARIES "-L$ENV{MKLROOT}/lib; -lmkl_intel_lp64; -lmkl_sequential; -lmkl_core; -lpthread; -lm" CACHE STRING "Set your MKL flags")
             SET(FFT_INCLUDES "$ENV{MKLROOT}/include/fftw" CACHE STRING "Set your MKL flags")
            SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${FFT_LIBRARIES}")
        ENDIF()        
  else()
    SET(FFT_LIBRARIES "-lfftw3" CACHE STRING "Use LIBFFTW")
    SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${FFT_LIBRARIES}")
  endif()
  
  SET(SCALFMM_INCLUDES "${SCALFMM_INCLUDES}; ${FFT_INCLUDES}")
  MESSAGE(STATUS "SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")
  MESSAGE(STATUS "SCALFMM_INCLUDES          = ${SCALFMM_INCLUDES}")
  
endif(ScalFMM_USE_FFT)

# Compile option
#ADD_DEFINITIONS(-Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion -fpic )
# 
# C++ 2011
   MESSAGE( STATUS "CMAKE_SYSTEM:   ${CMAKE_SYSTEM}  CMAKE_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID}  CMAKE_CXX_COMPILER ${NAME}")
IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE INTEL_VERSION)
    MESSAGE( STATUS "Intel:   ${INTEL_VERSION}")
  IF(APPLE)
    string(FIND ${CMAKE_CXX_COMPILER}  "icl++" NAME)
    if( ${NAME} GREATER  0)
      if(  ${INTEL_VERSION} EQUAL 15. OR  ${INTEL_VERSION} GREATER 15.0.0) 
        MESSAGE( STATUS "  Intel compiler is icl++ ( version >= 15.0.0)")	
        SET(INTEL_ICL_COMPILER "ON")
      else()
        MESSAGE(FATAL_ERROR  "  Intel compiler should be icl++ ( version >= 15.0.0)")  
      endif()
    else()
      MESSAGE(FATAL_ERROR  "  Intel compiler should be icl++ ( version >= 15.0.0)")
    endif()
  endif()
endif()
##############################################################################
##############################################################################

##################################################################
#                                               Use SSE                                                                     #
##################################################################

MESSAGE( STATUS "ScalFMM_USE_SSE              = ${ScalFMM_USE_SSE}" )
if( ScalFMM_USE_SSE )
	if(NOT EXISTS ${CMAKE_MODULE_PATH}/compileTestSse.cpp)
		message(FATAL_ERROR "The CompileTestSseFile does not exist (${CMAKE_MODULE_PATH}/compileTestSse.cpp)" )
	endif()
	MESSAGE( STATUS "SSE_FLAGS ${SSE_FLAGS}  -- ${CMAKE_CXX_FLAGS}  ")
	TRY_COMPILE(COMPILE_SSE  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_MODULE_PATH}/compileTestSse.cpp COMPILE_DEFINITIONS  "${CMAKE_CXX_FLAGS} ${SSE_FLAGS}" OUTPUT_VARIABLE COMPILE_SSE_OUTPUT)	
	
   	IF(${COMPILE_SSE})		
   	    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${SSE_FLAGS}")
   		     
		TRY_COMPILE(COMPILE_RESULT_VAR  ${CMAKE_CURRENT_BINARY_DIR} 
 						     ${CMAKE_MODULE_PATH}/checkSSEpe.cpp  COMPILE_DEFINITIONS  "${CMAKE_CXX_FLAGS} ${SSE_FLAGS}")
  		IF( NOT ${COMPILE_RESULT_VAR})
    	 		 SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__SSEPE_INTEL_COMPILER")
  		ENDIF()
  #		SET(ScalFMM_USE_AVX OFF)
  	ELSE(${COMPILE_SSE})
  		 MESSAGE(FATAL_ERROR "SSE NOT SUPPORTED ; Set ScalFMM_USE_SSE  to OFF \n Output from test is : ${COMPILE_SSE_OUTPUT}")	
  	ENDIF(${COMPILE_SSE})
endif()
##################################################################
#                                                Use AVX                                                                   #
##################################################################

MESSAGE(STATUS "ScalFMM_USE_AVX               = ${ScalFMM_USE_AVX}")
if(ScalFMM_USE_AVX)
	if(NOT EXISTS ${CMAKE_MODULE_PATH}/compileTestAvx.cpp)
		message(FATAL_ERROR "The CompileTestSseFile does not exist (${CMAKE_MODULE_PATH}/compileTestAvx.cpp)" )
	endif()
	
	TRY_COMPILE(COMPILE_AVX  ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_MODULE_PATH}/compileTestAvx.cpp COMPILE_DEFINITIONS  "${CMAKE_CXX_FLAGS}  ${AVX_FLAGS}" OUTPUT_VARIABLE COMPILE_AVX_OUTPUT)
   	IF(${COMPILE_AVX})
   	      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   ${AVX_FLAGS}")
#  	      SET( ScalFMM_USE_SSE   OFF   FORCE) # ne marche pas
 	 		TRY_COMPILE(COMPILE_RESULT_AVSPE  ${CMAKE_CURRENT_BINARY_DIR} 
 						     ${CMAKE_MODULE_PATH}/checkAVXpe.cpp  COMPILE_DEFINITIONS  "${CMAKE_CXX_FLAGS} ")
	  		IF( NOT ${COMPILE_RESULT_AVSPE})
 		   		  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__AVXPE_INTEL_COMPILER")
  			ENDIF()
 
		MESSAGE(STATUS ${CMAKE_CXX_FLAGS} )
ELSE(${COMPILE_AVX})
  		 MESSAGE(FATAL_ERROR "AVX NOT SUPPORTED ; Set ScalFMM_USE_AVX  to OFF \n Output from test is : ${COMPILE_AVX_OUTPUT} ")
ENDIF(${COMPILE_AVX})
endif(ScalFMM_USE_AVX)
#
# Error if both ScalFMM_USE_AVX AND ScalFMM_USE_SSE are set
#
IF( ScalFMM_USE_AVX AND ScalFMM_USE_SSE)
    MESSAGE(FATAL_ERROR "Check ScalFMM_USE_SSE or ScalFMM_USE_AVX BUT NOT BOTH. ")
ENDIF(ScalFMM_USE_AVX AND ScalFMM_USE_SSE)
##################################################################
#                                           Use  native MIC compilation                                          #
##################################################################
if( ScalFMM_USE_MIC_NATIVE ) 
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic")
else()
 # SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xhost")
ENDIF()
##################################################################
#                                
##################################################################

# Add CBLAS
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${CBLAS_LIBRARIES}")

##################################################################
#                                                                                                                                #                                          
#                                            END SETTING VARIABLES                                              # 
#                                                                                                                                #                                          
##################################################################
#
# Generate  ScalFmmConfig.h 
#
CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake 
        ${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h  )
#
##################################################################
#                                        Build - lib                                                                         #                                          
##################################################################
#
add_subdirectory(Src)
# Link with scalfmm lib
set(scalfmm_lib scalfmm)
#
##################################################################
#                                        Remove compilation  of the drivers                                 #                                          
##################################################################
if( NOT ScalFMM_BUILD_ONLY_LIB ) 
# Build - Examples and drivers
add_subdirectory(Examples)
endif()
##################################################################
#                                        Build tools (Internal use)                                                  #                                          
##################################################################
# Build - UTILs
if( ScalFMM_ONLY_DEVEL ) 
#if( NOT ScalFMM_BUILD_ONLY_LIB ) 
    add_subdirectory(Utils)
endif()

##################################################################
#                                        Build -Tests to see examples                                           #                                          
##################################################################
# Build - Tests
MESSAGE( STATUS "ScalFMM_BUILD_TESTS          = ${ScalFMM_BUILD_TESTS}" )
if( ScalFMM_BUILD_TESTS )
	add_subdirectory(Tests)
endif()

##################################################################
#                                        Build - Unitary and numerical Tests                                 #                                          
##################################################################
# 
MESSAGE( STATUS "ScalFMM_BUILD_UTESTS         = ${ScalFMM_BUILD_UTESTS}" )
if( ScalFMM_BUILD_UTESTS )
	INCLUDE(CTest)
	ENABLE_TESTING()
	add_subdirectory(UTests)
endif()

##################################################################
#                                           Add - doc                                                                     #
##################################################################
#
add_subdirectory(Doc)
##################################################################
#                                           Build - Addons                                                              #
##################################################################
#
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
##################################################################
#                                                                                                                                #                                          
#                                   Build and export cmake files                                                 #  
#                                                                                                                                #                                          
##################################################################

MESSAGE(STATUS "ScalFMM_BINARY_DIR:          " ${ScalFMM_BINARY_DIR})
CONFIGURE_FILE(${ScalFMM_SOURCE_DIR}/ScalFMMConfig.cmake.in
               ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake 
               @ONLY IMMEDIATE
               )
#
INSTALL(FILES 
  ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake
  DESTINATION lib/
  )
  IF( NOT ScalFMM_BUILD_ONLY_LIB) 
    INSTALL(FILES 
        	${ScalFMM_SOURCE_DIR}/Data/test20k.fma
        DESTINATION Data/
  )	
  endif()
  
 #################################################################
#                                                                                                                                #                                          
#                               build a CPack driven installer package                                       # 
#                                                                                                                                #                                          
##################################################################
#
# Use git commit number as   CPACK_PACKAGE_VERSION_PATCH
set(CPACK_PACKAGE_VERSION_PATCH 0)
execute_process(
        COMMAND git rev-list HEAD --count
        OUTPUT_VARIABLE CPACK_PACKAGE_VERSION_PATCH
        RESULT_VARIABLE RET
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
 )
 #STRING_TRIM(PATCH1 PATCH)
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY "ON")
#
set (CPACK_RESOURCE_FILE_LICENSE  "${CMAKE_CURRENT_SOURCE_DIR}/Licence.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "${ScalFMM_MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_MINOR "${ScalFMM_MINOR_VERSION}")
#
SET(PACK_PACKAGE_VERSION "${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME  "ScalFMM-${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}") 
SET(CPACK_SOURCE_IGNORE_FILES "\\\\.git;.DS_Store;.*~;/*.aux;/*.idx;/*.log;/*.out;/*.toc;/*.ilg;/*.ind;scalfmm.pro*;org.eclipse.core.resources.prefs;.cproject;.project")
LIST(APPEND CPACK_SOURCE_IGNORE_FILES  "${CMAKE_BINARY_DIR};${CMAKE_SOURCE_DIR}/Utils/;Notes;Deprecated;/Build*;/noDist/;/Bench;ToRemove;Olivier")
# Uniform;GroupTree;Adaptive;testUnif*;/*Lagrange*")
#LIST(APPEND CPACK_SOURCE_IGNORE_FILES  "Stages;Uniform;Adaptive;testUnif*;*Lagrange*" )
#
INCLUDE(CPack)
#
#
 ##################################################################
#                                                                                                                                 #                                          
#                                               PRINT MESSAGES                                                        # 
#                                                                                                                                 #                                          
 ##################################################################
message(STATUS "CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")
message(STATUS "CPACK_SOURCE_PACKAGE_FILE_NAME = ${CPACK_SOURCE_PACKAGE_FILE_NAME}")
#
 #                    Use Mem stats
MESSAGE( STATUS "ScalFMM_USE_MEM_STATS        = ${ScalFMM_USE_MEM_STATS}" )
#                    Use Log
MESSAGE( STATUS "ScalFMM_USE_LOG              = ${ScalFMM_USE_LOG}" )
#                    Use Assert
MESSAGE( STATUS "ScalFMM_USE_ASSERT           = ${ScalFMM_USE_ASSERT}" )
 #
MESSAGE(STATUS "CMAKE_CXX_FLAGS   =  ${CMAKE_CXX_FLAGS}")
MESSAGE(STATUS "ScaLFMM_CXX_FLAGS   =  ${ScaLFMM_CXX_FLAGS}")
MESSAGE(STATUS "SCALFMM_LIBRARIES   =  ${SCALFMM_LIBRARIES}")
MESSAGE(STATUS "SCALFMM_INCLUDES   =  ${SCALFMM_INCLUDES}")
 ##################################################################
#                                                        END                                                                  # 
 ##################################################################
 
