cmake_minimum_required (VERSION 2.8.3 FATAL_ERROR)

if(CMAKE_MAJOR_VERSION GREATER 1 AND CMAKE_MINOR_VERSION GREATER 7 AND CMAKE_PATCH_VERSION GREATER 11)
  cmake_policy(SET CMP0022 OLD)
endif()

cmake_policy(SET CMP0004 OLD)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

#===========================================================================
# Project Declaration
#===========================================================================
project(ScalFMM C CXX)

# check if compiling into source directories
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
  message(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

set(ScalFMM_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)


#===========================================================================
# Version Number
#===========================================================================

# ScalFMM version number.  An even minor number corresponds to releases.
set(ScalFMM_MAJOR_VERSION 1)
set(ScalFMM_MINOR_VERSION 3)
set(ScalFMM_PATCH_VERSION 0)
set(ScalFMM_VERSION  "${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}.${ScalFMM_PATCH_VERSION}" )

set( MORSE_DISTRIB_DIR "" CACHE PATH "Directory of MORSE distribution")

if (MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_SOURCE_DIR}/CMakeModules/morse/")

  # Add extra cmake module path and initialize morse cmake modules
  # --------------------------------------------------------------
  if(MORSE_DISTRIB_DIR)
      list(APPEND CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/cmake_modules)
      list(APPEND CMAKE_MODULE_PATH "${MORSE_DISTRIB_DIR}/cmake_modules/morse")
      set(MORSE_CMAKE_MODULE_PATH ${MORSE_DISTRIB_DIR}/cmake_modules/morse )
  elseif(EXISTS "${CMAKE_SOURCE_DIR}/CMakeModules/morse/")
      list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules/morse/)
      set(MORSE_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules/morse )
  endif()
  include(MorseInit)

  #
  # Options
  option( ScalFMM_USE_MPI              "Set to ON to build ScaFMM with MPI"          OFF )
  option( ScalFMM_USE_BLAS             "Set to ON to build ScaFMM with BLAS"         OFF )
  option( ScalFMM_USE_FFT              "Set to ON to build ScaFMM with FFTW"         OFF )
  option( ScalFMM_BUILD_TESTS          "Set to ON to build functionalities Tests"    OFF )
  option( ScalFMM_BUILD_UTESTS         "Set to ON to build UTests"                   OFF )
  option( ScalFMM_BUILD_DEBUG          "Set to ON to build in Debug"                 OFF )
  option( ScalFMM_USE_LOG              "Set to ON to print output debug information" OFF )
  option( ScalFMM_USE_MEM_STATS        "Set to ON to profile memory"                 OFF )
  option( ScalFMM_USE_DOUBLE_PRECISION "Set to ON to compile in double precision"    ON  )
  option( ScalFMM_ATTACHE_SOURCE       "Set to ON to compile with -g"                OFF )
  option( ScalFMM_USE_ADDONS           "Set to ON to compile add ons"                OFF )
  option( ScalFMM_USE_SSE              "Set to ON to compile with SSE support"       ON  )
  option( ScalFMM_USE_AVX              "Set to ON to compile with AVX support"       OFF  )
  option( ScalFMM_USE_ASSERT           "Set to ON to enable safe tests during execution" ON  )
  option( ScalFMM_USE_MIC_NATIVE       "Set to ON to compile in native mode for MIC" OFF  )
  option( ScalFMM_BUILD_ONLY_LIB       "Set to ON to compile only the lib (examples are not compiled) " OFF  )
  option( ScalFMM_ONLY_DEVEL           "Set to ON to compile Development tools (only scalfmm team)" ON )
  if( ScalFMM_ONLY_DEVEL )
    option( ScalFMM_USE_STARPU "Set to ON to build ScalFMM with StarPU" OFF )
  endif()
  if( ScalFMM_USE_MPI )
    try_compile(COMPILE_INTEL ${CMAKE_CURRENT_BINARY_DIR}
                ${ScalFMM_CMAKE_MODULE_PATH}/compileTestIntel.cpp
                COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS}")
    if (COMPILE_INTEL)
      set(CMAKE_CXX_COMPILER_ID  "Intel")
    endif()
  endif()
  message(STATUS "CXX  ${CMAKE_CXX_COMPILER_ID}" )
  # Set scalfmm to default libraries
  set(SCALFMM_LIBRARIES "")
  set(ScalFMM_CXX_FLAGS "-std=c++11 -fpic -Wall")
  #
  #
  # Test if openmp is here
  #
  find_package (OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else(OPENMP_FOUND)
    message(WARNING "OPENMP NOT FOUND")
  endif(OPENMP_FOUND)
  message(STATUS "ScalFMM_BUILD_DEBUG = ${ScalFMM_BUILD_DEBUG}" )
  #
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} -m64")
  endif()
  ##############################################################################
  #                           Compile options                                  #
  ##############################################################################
  #  -xHost -mfpmath=sse
  # -Wall Wnosign-conversion
  #
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    # INTEL
    IF (APPLE)
      set(SSE_FLAGS  "-msse4  -mfpmath=sse")   # -mtune=native -march=native
    else(APPLE)
      set(AVX_FLAGS  "-mavx")
      set(SSE_FLAGS  "-axSSE4.2")
    endif(APPLE)
    #-Wshadow -Wpointer-arith -Wcast-qual -Wconversion  -Wall -Wnosign-conversion ")
  else()
    # NOT      INTEL
    if(NOT ScalFMM_USE_MPI)
      include(CheckCCompilerFlag)
      check_c_compiler_flag(-Wzero-as-null-pointer-constant HAS_WZERO_NULL_PTR_FLAG)
      if(HAS_WZERO_NULL_PTR_FLAG)
        set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} -Wzero-as-null-pointer-constant")
      endif()
    else()
      include(CheckCCompilerFlag)
      check_c_compiler_flag(-Wno-literal-suffix HAS_NO_LITERAL_SUFFIX_FLAG)
      if(HAS_NO_LITERAL_SUFFIX_FLAG)
        set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} -Wno-literal-suffix")
      endif()
    endif()
    set(AVX_FLAGS "-mtune=native -march=native")
    IF (APPLE)
      set(SSE_FLAGS  "-msse3  -mfpmath=sse")   # -mtune=native -march=native
    else(APPLE)
      set(SSE_FLAGS  "-mtune=native -march=native")   # -mtune=native -march=native
    endif(APPLE)
  endif()

  ##############################################################################
  #                           FUSE list                                        #
  ##############################################################################
  set(FUSE_LIST "")
  # then do list(APPEND FUSE_LIST "BLAS") to protect from FUSE_BLAS
  list(APPEND FUSE_LIST "MPI")

  ##############################################################################
  #                        Debug      flags                                    #
  ##############################################################################
  #
  if( ScalFMM_BUILD_DEBUG )
    #
    # Compile Debug flags
    #
    set(CMAKE_BUILD_TYPE Debug)
    ADD_DEFINITIONS(-O0)
    if( APPLE )
      set(SCALFMM_FLAGS_OPTI_DEBUG "-m64  -funroll-loops" CACHE STRING
        "Set your optimization flags for debug mode.")
    else(APPLE)
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        set(SCALFMM_FLAGS_OPTI_DEBUG "-fp-model strict  -funroll-loops" CACHE STRING
          "Set your optimization flags for debug mode.")
      else()
        set(SCALFMM_FLAGS_OPTI_DEBUG "-funroll-loops" CACHE STRING
          "Set your optimization flags for debug mode.")
      endif()
    endif(APPLE)
    #  ADD_DEFINITIONS(${SCALFMM_FLAGS_OPTI_DEBUG})
    set(ScalFMM_CXX_FLAGS "${ScalFMM_CXX_FLAGS} ${SCALFMM_FLAGS_OPTI_DEBUG}")
    set(CMAKE_CXX_FLAGS_TYPE "${CMAKE_CXX_FLAGS_DEBUG}")

  else()
    ##############################################################################
    #                          Release flags                                     #
    ##############################################################################
    set(CMAKE_BUILD_TYPE Release)
    # force -O3 in release
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    string(REPLACE "-O2" "-O3" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
    set(CMAKE_CXX_FLAGS_TYPE  "${CMAKE_CXX_FLAGS_RELEASE}")
    # if compiler is intel add -ip
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} -ip")
      # set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} -ip -no-prec-div -no-prec-sqrt")
    endif()
    # Compile optimization
    if(APPLE)
      set(SCALFMM_FLAGS_OPTI_RELEASE "-m64  -funroll-loops" CACHE STRING
        "Set your optimization flags for release mode.")
    else(APPLE)
      # Not apple system - Check the compiler flags
      if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        # set(SCALFMM_FLAGS_OPTI_RELEASE "-fp-model precise -fp-model source -fimf-precision=low  -funroll-loops -ftree-vectorize"
        #   CACHE STRING "Set your optimization flags for release mode.")
        set(SCALFMM_FLAGS_OPTI_RELEASE " " CACHE STRING
          "Set your optimization flags for release mode.")
      else()
        set(SCALFMM_FLAGS_OPTI_RELEASE "-ffast-math -funroll-loops -ftree-vectorize" CACHE STRING
          "Set your optimization flags for release mode.")
      endif()
    endif(APPLE)
    #
    set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} ${SCALFMM_FLAGS_OPTI_RELEASE}")
  endif()
  #
  ##############################################################################
  #                Attach source code to exec                                  #
  ##############################################################################
  message(STATUS "ScalFMM_ATTACHE_SOURCE       = ${ScalFMM_ATTACHE_SOURCE}" )
  if( ScalFMM_ATTACHE_SOURCE )
    message(STATUS "Option -g is used, the code is attached to the binary." )
    set(ScalFMM_CXX_FLAGS  "${ScalFMM_CXX_FLAGS} -g")
  endif(ScalFMM_ATTACHE_SOURCE)
  #
  ##############################################################################
  #                                 find MPI                                   #
  ##############################################################################
  #
  if( ScalFMM_USE_MPI )

    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_CXX_COMPILER=path/to/mpicxx
    # at cmake configure
    if(NOT MPI_C_COMPILER)
        set(MPI_C_COMPILER mpicc)
    endif()
    if(NOT MPI_CXX_COMPILER)
        set(MPI_CXX_COMPILER mpicxx)
    endif()

    find_package(MPI REQUIRED)

    if (MPI_CXX_INCLUDE_PATH)
        include_directories( ${MPI_CXX_INCLUDE_PATH} )
    endif()
    if (MPI_CXX_COMPILE_FLAGS)
        set(ScalFMM_CXX_FLAGS "${ScalFMM_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
    endif()
    if (MPI_CXX_INCLUDE_PATH)
        set(SCALFMM_INCLUDES  "${SCALFMM_INCLUDES}; ${MPI_CXX_INCLUDE_PATH}")
    endif()
    if (MPI_CXX_LINK_FLAGS)
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${MPI_CXX_LINK_FLAGS}")
    endif()
    if (MPI_CXX_LIBRARIES)
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${MPI_CXX_LIBRARIES}")
    endif()

  endif()
  #
  ##############################################################################
  #                             Blas option                                    #
  ##############################################################################
  #
  if( ScalFMM_USE_BLAS )
    message(STATUS "CMAKE_CXX_COMPILER_ID STREQUAL  ${CMAKE_CXX_COMPILER_ID}")

    option( ScalFMM_USE_MKL_AS_BLAS "Set to ON to use MKL CBLAS" OFF )

    if( ScalFMM_USE_MKL_AS_BLAS )

      set(BLA_VENDOR "Intel10_64lp_seq")
      find_package(BLASEXT REQUIRED)

      if(BLAS_FOUND AND BLAS_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
      else()
        message(FATAL_ERROR "A BLAS library is required but has not been found")
      endif()
      unset(LAPACK_LIBRARIES)
      if (BLAS_LIBRARIES)
        set(BLASLAPACK_LIBRARIES ${BLAS_LIBRARIES})
      endif()

    elseif(ScalFMM_USE_EXTERNAL_BLAS)

      message(STATUS "BLAS SET BY EXTERNAL PROGRAM = ${BLAS_LIBRARIES}")

    else()

      find_package(BLASEXT REQUIRED)
      find_package(LAPACKEXT REQUIRED)
      set(BLASLAPACK_LIBRARIES)
      if (LAPACK_LIBRARIES)
        set(BLASLAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
      endif()
      if (BLAS_LIBRARIES)
        list(APPEND BLASLAPACK_LIBRARIES "${BLAS_LIBRARIES}")
      endif()

    endif()

    set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${BLASLAPACK_LIBRARIES}")
    #message(STATUS "SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")

  endif(ScalFMM_USE_BLAS)
  list(APPEND FUSE_LIST "BLAS")
  #
  # FFT option
  #
  if( ScalFMM_USE_FFT )

    option( ScalFMM_USE_MKL_AS_FFTW "Set to ON to use MKL FFTW" OFF )

    message(STATUS "SCALFMM USE FFT Configure:")

    if( ScalFMM_USE_MKL_AS_FFTW )

      message(STATUS "         SCALFMM USE MKL ")

      if( ScalFMM_USE_MKL_AS_BLAS )

        unset(FFT_LIBRARIES CACHE)
        message(STATUS "         SCALFMM USE MKL already defined")
        set(FFT_INCLUDES "$ENV{MKLROOT}/include/fftw" CACHE STRING "Set your MKL flags")

      else(ScalFMM_USE_MKL_AS_BLAS)

        # The package can be used with the following COMPONENTS:
        # MKL, THREADS|OMP and/or SIMPLE|DOUBLE|LONG|QUAD
        # Default is DOUBLE and without THREADS|OMP
        find_package(FFTW REQUIRED COMPONENTS MKL)
        set(FFT_LIBRARIES "-L${FFTW_LIBRARY_DIRS_DEP};" CACHE STRING "Set your MKL flags")
        foreach (fft_lib ${FFTW_LIBRARIES_DEP})
          set(FFT_LIBRARIES "${FFT_LIBRARIES} ${fft_lib};")
        endforeach()
        set(FFT_INCLUDES "${FFTW_INCLUDE_DIRS_DEP}" )
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${FFT_LIBRARIES}")

      endif(ScalFMM_USE_MKL_AS_BLAS)

    else(ScalFMM_USE_MKL_AS_FFTW)

      message(STATUS "         SCALFMM USE FFTW ")
      # The package can be used with the following COMPONENTS:
      # MKL, THREADS|OMP and/or SIMPLE|DOUBLE|LONG|QUAD
      # Default is DOUBLE and without THREADS|OMP
      find_package(FFTW REQUIRED)

      set(FFT_LIBRARIES "-L${FFTW_LIBRARY_DIRS_DEP};" CACHE STRING "Set your MKL flags")
      foreach (fft_lib ${FFTW_LIBRARIES_DEP})
        set(FFT_LIBRARIES "${FFT_LIBRARIES} ${fft_lib};")
      endforeach()
      set(FFT_INCLUDES "${FFTW_INCLUDE_DIRS_DEP}" )
      set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${FFT_LIBRARIES}")

    endif(ScalFMM_USE_MKL_AS_FFTW)

    if (FFT_INCLUDES)
        set(SCALFMM_INCLUDES "${SCALFMM_INCLUDES}; ${FFT_INCLUDES}")
    endif()
    message(STATUS "         SCALFMM_LIBRARIES          = ${SCALFMM_LIBRARIES}")
    message(STATUS "         SCALFMM_INCLUDES          = ${SCALFMM_INCLUDES}")

  endif(ScalFMM_USE_FFT)
  list(APPEND FUSE_LIST "FFT")

  # Compile option
  #ADD_DEFINITIONS(-Wall -Wshadow -Wpointer-arith -Wcast-qual -Wconversion -fpic )
  #
  # C++ 2011
  message( STATUS "CMAKE_SYSTEM:   ${CMAKE_SYSTEM}  CMAKE_CXX_COMPILER_ID ${CMAKE_CXX_COMPILER_ID}  CMAKE_CXX_COMPILER ${NAME}")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE INTEL_VERSION)
    message( STATUS "Intel:   ${INTEL_VERSION}")
    if(APPLE)
      string(FIND ${CMAKE_CXX_COMPILER}  "icl++" NAME)
      if( ${NAME} GREATER  0)
        if(  ${INTEL_VERSION} EQUAL 15. OR  ${INTEL_VERSION} GREATER 15.0.0)
  	message( STATUS "  Intel compiler is icl++ ( version >= 15.0.0)")
  	set(INTEL_ICL_COMPILER "ON")
        else()
  	message(FATAL_ERROR  "  Intel compiler should be icl++ ( version >= 15.0.0)")
        endif()
      else()
        message(FATAL_ERROR  "  Intel compiler should be icl++ ( version >= 15.0.0)")
      endif()
    endif()
  endif()
  #
  ##################################################################
  ##################################################################
  #
  message( STATUS "ScalFMM_USE_STARPU             = ${ScalFMM_USE_STARPU}" )
  if( ScalFMM_USE_STARPU )

    # CUDA could be used with StarPU enabled
    option( ScalFMM_USE_CUDA "Set to ON to use CUDA with StarPU" OFF )
    message( STATUS "ScalFMM_USE_CUDA             = ${ScalFMM_USE_CUDA}" )
    if(ScalFMM_USE_CUDA)
        execute_process(COMMAND nvcc --version ERROR_VARIABLE cuda_error_output OUTPUT_QUIET)
        if(cuda_error_output)
            message( FATAL_ERROR "nvcc is needed with CUDA." )
        endif()
        if(NOT DEFINED CUSTOM_CUDA_FLAGS)
            set( CUSTOM_CUDA_FLAGS "-std=c++11;-arch=sm_20" CACHE
            STRING "Set your CUDA flags, for example : -arch=sm_20;-ptxas-options=-v;-use_fast_math")
        endif()
        # This is needed to remove backslash after space in ADD_CUSTOM_COMMAND
        separate_arguments(CUSTOM_CUDA_FLAGS)
        message( STATUS "CUSTOM_CUDA_FLAGS             = ${CUSTOM_CUDA_FLAGS}" )

        find_package(CUDA REQUIRED)

        if (CUDA_INCLUDE_DIRS)
            include_directories(${CUDA_INCLUDE_DIRS})
        endif()
        if (CUDA_LIBRARIES)
            set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${CUDA_LIBRARIES}")
        endif()
        
        set(CUDA_NEEDED_INCLUDE_DIRS ${CMAKE_BINARY_DIR}/Src)
    endif()

    # Find StarPU with a list of optional components
    set(ScalFMM_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")
    # create list of components in order to make a single call to find_package(starpu...)
    set(STARPU_COMPONENT_LIST "HWLOC")
    if(ScalFMM_USE_MPI)
        list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if(ScalFMM_USE_CUDA)
        list(APPEND STARPU_COMPONENT_LIST "CUDA")
    endif()

    find_package(STARPU ${ScalFMM_STARPU_VERSION} REQUIRED
                 COMPONENTS ${STARPU_COMPONENT_LIST})

    # Append list of libraries and include dirs
    include_directories(${STARPU_INCLUDE_DIRS_DEP})
    foreach (starpu_libdir ${STARPU_LIBRARY_DIRS_DEP})
      if (${starpu_libdir} MATCHES "^ *-L")
          set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${starpu_libdir}")
      else()
          set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; -L${starpu_libdir}")
      endif()
    endforeach()
    foreach (starpu_lib ${STARPU_LIBRARIES_DEP})
      if (EXISTS ${starpu_lib} OR ${starpu_lib} MATCHES "^ *-")
          set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${starpu_lib}")
      else()
          set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; -l${starpu_lib}")
      endif()
    endforeach()

    if (CUDA_LIBRARIES)
        set(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; ${CUDA_LIBRARIES}")
    endif()

    # Message
    message(STATUS "         STARPU_LIBRARIES          = ${STARPU_LIBRARIES}")
    if (STARPU_INCLUDE_DIRS)
        message(STATUS "         STARPU_INCLUDES           = ${STARPU_INCLUDE_DIRS}")
    endif()

    OPTION( ScalFMM_USE_OPENCL      "Set to ON to use OPENCL with StarPU" OFF )
    MESSAGE( STATUS "ScalFMM_USE_OPENCL             = ${ScalFMM_USE_OPENCL}" )
    if(ScalFMM_USE_OPENCL)
        include_directories($ENV{OPENCL_INC})
        SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES}; -L$ENV{OPENCL_LIB}; -lOpenCL")
    endif()
  endif(ScalFMM_USE_STARPU)
  list(APPEND FUSE_LIST "STARPU")
  list(APPEND FUSE_LIST "CUDA")
  list(APPEND FUSE_LIST "OPENCL")

  ##################################################################
  #                         Use SSE                                #
  ##################################################################

  message( STATUS "ScalFMM_USE_SSE              = ${ScalFMM_USE_SSE}" )
  if( ScalFMM_USE_SSE )
    if(NOT EXISTS ${ScalFMM_CMAKE_MODULE_PATH}/compileTestSse.cpp)
      message(FATAL_ERROR "The CompileTestSseFile does not exist (${ScalFMM_CMAKE_MODULE_PATH}/compileTestSse.cpp)" )
    endif()
    message( STATUS "SSE_FLAGS ${SSE_FLAGS}  -- ${CMAKE_CXX_FLAGS}  ")
    try_compile(COMPILE_SSE  ${CMAKE_CURRENT_BINARY_DIR}
                ${ScalFMM_CMAKE_MODULE_PATH}/compileTestSse.cpp
                COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${SSE_FLAGS}"
                OUTPUT_VARIABLE COMPILE_SSE_OUTPUT)

    if(${COMPILE_SSE})
      set(ScalFMM_CXX_FLAGS "${ScalFMM_CXX_FLAGS} ${SSE_FLAGS}")

      try_compile(COMPILE_RESULT_VAR ${CMAKE_CURRENT_BINARY_DIR}
                  ${ScalFMM_CMAKE_MODULE_PATH}/checkSSEpe.cpp
                  COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${SSE_FLAGS}")
      if( NOT ${COMPILE_RESULT_VAR})
        set(__SSEPE_INTEL_COMPILER ON)
      endif()
      #		set(ScalFMM_USE_AVX OFF)
    else(${COMPILE_SSE})
      message(FATAL_ERROR "SSE NOT SUPPORTED ; Set ScalFMM_USE_SSE  to OFF \n Output from test is : ${COMPILE_SSE_OUTPUT}")
    endif(${COMPILE_SSE})
  endif()
  list(APPEND FUSE_LIST "SSE")

  ##################################################################
  #                           Use AVX                              #
  ##################################################################

  message(STATUS "ScalFMM_USE_AVX               = ${ScalFMM_USE_AVX}")
  if(ScalFMM_USE_AVX)
    if(NOT EXISTS ${ScalFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp)
      message(FATAL_ERROR "The CompileTestSseFile does not exist (${ScalFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp)" )
    endif()

    try_compile(COMPILE_AVX ${CMAKE_CURRENT_BINARY_DIR}
                ${ScalFMM_CMAKE_MODULE_PATH}/compileTestAvx.cpp
                COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${AVX_FLAGS}"
                OUTPUT_VARIABLE COMPILE_AVX_OUTPUT)
    if(${COMPILE_AVX})
      set(ScalFMM_CXX_FLAGS "${ScalFMM_CXX_FLAGS}   ${AVX_FLAGS}")
      #             set( ScalFMM_USE_SSE   OFF   FORCE) # ne marche pas
      try_compile(COMPILE_RESULT_AVSPE ${CMAKE_CURRENT_BINARY_DIR}
                  ${ScalFMM_CMAKE_MODULE_PATH}/checkAVXpe.cpp
                  COMPILE_DEFINITIONS "${CMAKE_CXX_FLAGS} ${AVX_FLAGS}")
      if( NOT ${COMPILE_RESULT_AVSPE})
        set(__AVXPE_INTEL_COMPILER ON)
      endif()

      message(STATUS ${CMAKE_CXX_FLAGS} )
    else(${COMPILE_AVX})
      message(FATAL_ERROR "AVX NOT SUPPORTED ; Set ScalFMM_USE_AVX  to OFF \n Output from test is : ${COMPILE_AVX_OUTPUT} ")
    endif(${COMPILE_AVX})
  endif(ScalFMM_USE_AVX)
  list(APPEND FUSE_LIST "AVX")
  #
  # Error if both ScalFMM_USE_AVX AND ScalFMM_USE_SSE are set
  #
  if( ScalFMM_USE_AVX AND ScalFMM_USE_SSE)
    message(FATAL_ERROR "Check ScalFMM_USE_SSE or ScalFMM_USE_AVX BUT NOT BOTH. ")
  endif(ScalFMM_USE_AVX AND ScalFMM_USE_SSE)
  ##################################################################
  #                     Use  native MIC compilation                #
  ##################################################################
  if( ScalFMM_USE_MIC_NATIVE )
    include(CheckCCompilerFlag)
    check_c_compiler_flag(-mmic HAS_MMIC_FLAG)
    if(NOT HAS_MMIC_FLAG)
      message(FATAL_ERROR "MIC NOT SUPPORTED ; Set ScalFMM_USE_MIC_NATIVE to OFF")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmic")
  else()
    # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xhost")
  endif()
  list(APPEND FUSE_LIST "MIC")
  ##################################################################
  #
  ##################################################################
  #
  # Add CBLAS
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${CBLAS_LIBRARIES}")
  # SI CBLAS necessaire utiliser la ligne cidessous
  #set(SCALFMM_LIBRARIES   "${CBLAS_LIBRARIES}; ${SCALFMM_LIBRARIES} ")
  #
  ##################################################################
  #                                                                #
  #                     END SETTING VARIABLES                      #
  #                                                                #
  ##################################################################
  #
  # Generate  ScalFmmConfig.h
  #
  set(ScalFMM_COMPILE_FLAGS "${ScalFMM_CXX_FLAGS} ${CMAKE_CXX_FLAGS_TYPE} ${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${ScalFMM_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS_TYPE} ${CMAKE_CXX_FLAGS}")
  set(ScalFMM_COMPILE_LIBS   "${SCALFMM_LIBRARIES}")
  configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/Src/ScalFmmConfig.h.cmake
                  ${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h  )
  #
  ##################################################################
  #                      Build - lib                               #
  ##################################################################
  #
  add_subdirectory(Src)
  # Link with scalfmm lib
  set(scalfmm_lib scalfmm)
  #
  ##################################################################
  #            Remove compilation  of the drivers                  #
  ##################################################################
  if( NOT ScalFMM_BUILD_ONLY_LIB )
    # Build - Examples and drivers
    add_subdirectory(Examples)
  endif()
  ##################################################################
  #           Build tools (Internal use)                           #
  ##################################################################
  # Build - UTILs
  if( ScalFMM_ONLY_DEVEL )
    #if( NOT ScalFMM_BUILD_ONLY_LIB )
    add_subdirectory(Utils)
  endif()

  ##################################################################
  #                Build -Tests to see examples                    #
  ##################################################################
  # Build - Tests
  message( STATUS "ScalFMM_BUILD_TESTS          = ${ScalFMM_BUILD_TESTS}" )
  if( ScalFMM_BUILD_TESTS )
    add_subdirectory(Tests)
  endif()

  ##################################################################
  #              Build - Unitary and numerical Tests               #
  ##################################################################
  #
  message( STATUS "ScalFMM_BUILD_UTESTS         = ${ScalFMM_BUILD_UTESTS}" )
  if( ScalFMM_BUILD_UTESTS )
    include(CTest)
    enable_testing()
    add_subdirectory(UTests)
  endif()

  ##################################################################
  #                           Add - doc                            #
  ##################################################################
  #
  add_subdirectory(Doc)
  ##################################################################
  #                       Build - Addons                           #
  ##################################################################
  #
  if( ScalFMM_USE_ADDONS )
    file(GLOB all_addons ./Addons/*)
    # Options
    option( ScalFMM_ADDON_CKERNELAPI "Set to ON to build ScaFMM FMM API interface" OFF )

    foreach(addon_name ${all_addons})
      message(STATUS "Addons += ${addon_name}" )
      add_subdirectory(${addon_name})
    endforeach()
  endif()

  # Add custom clean command
  if("${CMAKE_GENERATOR}" MATCHES "Make")
    add_custom_target(reset
      COMMAND ${CMAKE_MAKE_PROGRAM} clean && rm -r ./CMake* && cmake ..
      WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
      COMMENT "Cleaning Build dir and CMake cache, then do cmake.."
      )
  endif()
  #
  ##################################################################
  #                                                                #
  #                  Build and export cmake files                  #
  #                                                                #
  ##################################################################
  #
  message(STATUS "ScalFMM_BINARY_DIR:          " ${ScalFMM_BINARY_DIR})
  configure_file(${ScalFMM_SOURCE_DIR}/ScalFMMConfig.cmake.in
    ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake
    @ONLY IMMEDIATE
    )
  #
  install(FILES
    ${ScalFMM_BINARY_DIR}/ScalFMMConfig.cmake
    DESTINATION lib/
    )
  if( NOT ScalFMM_BUILD_ONLY_LIB)
    install(FILES
      ${ScalFMM_SOURCE_DIR}/Data/test20k.fma
      DESTINATION Data/
      )
  endif()

  ##################################################################
  #                                                                #
  #              build a CPack driven installer package            #
  #                                                                #
  ##################################################################
  #
  # Use git commit number as CPACK_PACKAGE_VERSION_PATCH
  set(CPACK_PACKAGE_VERSION_PATCH 0)
  execute_process(
    COMMAND git rev-list HEAD --count
    OUTPUT_VARIABLE CPACK_PACKAGE_VERSION_PATCH
    RESULT_VARIABLE RET
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  #string_TRIM(PATCH1 PATCH)
  set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY "ON")
  #
  set (CPACK_RESOURCE_FILE_LICENSE  "${CMAKE_CURRENT_SOURCE_DIR}/Licence.txt")
  set(CPACK_PACKAGE_VERSION_MAJOR "${ScalFMM_MAJOR_VERSION}")
  set(CPACK_PACKAGE_VERSION_MINOR "${ScalFMM_MINOR_VERSION}")
  #
  set(PACK_PACKAGE_VERSION "${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_SOURCE_GENERATOR "TGZ")
  set(CPACK_SOURCE_PACKAGE_FILE_NAME  "ScalFMM-${ScalFMM_MAJOR_VERSION}.${ScalFMM_MINOR_VERSION}-${CPACK_PACKAGE_VERSION_PATCH}")
  set(CPACK_SOURCE_IGNORE_FILES "\\\\.git;.DS_Store;.*~;/*.aux;/*.idx;/*.log;/*.out;/*.toc;/*.ilg;/*.ind;scalfmm.pro*;org.eclipse.core.resources.prefs;.cproject;.project")
  list(APPEND CPACK_SOURCE_IGNORE_FILES  "${CMAKE_BINARY_DIR};${CMAKE_SOURCE_DIR}/Utils/;Notes;Deprecated;/Build*;/noDist/;/Bench;ToRemove;Olivier")
  # Uniform;GroupTree;Adaptive;testUnif*;/*Lagrange*")
  #list(APPEND CPACK_SOURCE_IGNORE_FILES  "Stages;Uniform;Adaptive;testUnif*;*Lagrange*" )
  #
  include(CPack)
  #
  #
  ##################################################################
  #                                                                #
  #                      PRINT messages                            #
  #                                                                #
  ##################################################################
  message(STATUS "CPACK_SOURCE_IGNORE_FILES = ${CPACK_SOURCE_IGNORE_FILES}")
  message(STATUS "CPACK_SOURCE_PACKAGE_FILE_NAME = ${CPACK_SOURCE_PACKAGE_FILE_NAME}")
  #
  #                    Use Mem stats
  message(STATUS "ScalFMM_USE_MEM_STATS        = ${ScalFMM_USE_MEM_STATS}" )
  #                    Use Log
  message(STATUS "ScalFMM_USE_LOG              = ${ScalFMM_USE_LOG}" )
  #                    Use Assert
  message(STATUS "ScalFMM_USE_ASSERT           = ${ScalFMM_USE_ASSERT}" )
  #
  #message(STATUS "CMAKE_CXX_FLAGS   =  ${CMAKE_CXX_FLAGS}")
  #set(CMAKE_CXX_FLAGS "${ScalFMM_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

  message(STATUS "CMAKE_CXX_FLAGS     =  ${CMAKE_CXX_FLAGS}")
  message(STATUS "SCALFMM_CXX_FLAGS   =  ${ScalFMM_CXX_FLAGS}")
  message(STATUS "SCALFMM_LIBRARIES   =  ${SCALFMM_LIBRARIES}")
  message(STATUS "SCALFMM_INCLUDES    =  ${SCALFMM_INCLUDES}")
  ##################################################################
  #                               END                              #
  ##################################################################

else(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_SOURCE_DIR}/CMakeModules/morse/")

  message(STATUS "MORSE_DISTRIB_DIR is not set")
  message(STATUS "Please indicate where is located your MORSE distribution directory."
    " This is necessary to find cmake_modules.")
  message(STATUS "If you are a MORSE developer, you normally have access to the svn repository of the MORSE_DISTRIB_DIR: "
    "svn checkout svn+ssh://yourlogin@scm.gforge.inria.fr/svnroot/morse/trunk/morse_distrib")

endif(MORSE_DISTRIB_DIR OR EXISTS "${CMAKE_SOURCE_DIR}/CMakeModules/morse/")
