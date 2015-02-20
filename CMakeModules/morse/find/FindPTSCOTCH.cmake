###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find PTSCOTCH include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PTSCOTCH
#               [REQUIRED]             # Fail with error if ptscotch is not found
#              )
#
#  PTSCOTCH depends on the following libraries:
#   - Threads
#   - MPI
#
# This module finds headers and ptscotch library.
# Results are reported in variables:
#  PTSCOTCH_FOUND            - True if headers and requested libraries were found
#  PTSCOTCH_INCLUDE_DIRS     - ptscotch include directories
#  PTSCOTCH_LIBRARY_DIRS     - Link directories for ptscotch libraries
#  PTSCOTCH_LIBRARIES        - ptscotch component libraries to be linked
#  PTSCOTCH_INCLUDE_DIRS_DEP - ptscotch + dependencies include directories
#  PTSCOTCH_LIBRARY_DIRS_DEP - ptscotch + dependencies link directories
#  PTSCOTCH_LIBRARIES_DEP    - ptscotch libraries + dependencies
#  PTSCOTCH_INTSIZE          - Number of octets occupied by a SCOTCH_Num
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPTSCOTCH=path/to/ptscotch):
#  PTSCOTCH_DIR              - Where to find the base directory of ptscotch
#  PTSCOTCH_INCDIR           - Where to find the header files
#  PTSCOTCH_LIBDIR           - Where to find the library files

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)

if (NOT PTSCOTCH_FOUND)
    set(PTSCOTCH_DIR "" CACHE PATH "Root directory of PTSCOTCH library")
    if (NOT PTSCOTCH_FIND_QUIETLY)
        message(STATUS "A cache variable, namely PTSCOTCH_DIR, has been set to specify the install directory of PTSCOTCH")
    endif()
endif()

# PTSCOTCH may depend on Threads, try to find it
if (NOT THREADS_FOUND)
    if (PTSCOTCH_FIND_REQUIRED)
        find_package(Threads REQUIRED)
    else()
        find_package(Threads)
    endif()
endif()

# PTSCOTCH may depend on MPI, try to find it
if (NOT MPI_FOUND)
    if (PTSCOTCH_FIND_REQUIRED)
        find_package(MPI REQUIRED)
    else()
        find_package(MPI)
    endif()
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
if(WIN32)
    string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
else()
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
endif()
list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)


# Try to find the ptscotch header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(PTSCOTCH_INCDIR)
    set(PTSCOTCH_ptscotch.h_DIRS "PTSCOTCH_ptscotch.h_DIRS-NOTFOUND")
    find_path(PTSCOTCH_ptscotch.h_DIRS
      NAMES ptscotch.h
      HINTS ${PTSCOTCH_INCDIR})
else()
    if(PTSCOTCH_DIR)
        set(PTSCOTCH_ptscotch.h_DIRS "PTSCOTCH_ptscotch.h_DIRS-NOTFOUND")
        find_path(PTSCOTCH_ptscotch.h_DIRS
          NAMES ptscotch.h
          HINTS ${PTSCOTCH_DIR}
          PATH_SUFFIXES include)
    else()
        set(PTSCOTCH_ptscotch.h_DIRS "PTSCOTCH_ptscotch.h_DIRS-NOTFOUND")
        find_path(PTSCOTCH_ptscotch.h_DIRS
          NAMES ptscotch.h
          HINTS ${_inc_env})
    endif()
endif()
mark_as_advanced(PTSCOTCH_ptscotch.h_DIRS)

# If found, add path to cmake variable
# ------------------------------------
if (PTSCOTCH_ptscotch.h_DIRS)
    set(PTSCOTCH_INCLUDE_DIRS "${PTSCOTCH_ptscotch.h_DIRS}")
else ()
    set(PTSCOTCH_INCLUDE_DIRS "PTSCOTCH_INCLUDE_DIRS-NOTFOUND")
    if (NOT PTSCOTCH_FIND_QUIETLY)
        message(STATUS "Looking for ptscotch -- ptscotch.h not found")
    endif()
endif()


# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
else()
    if(APPLE)
        string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
    else()
        string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    endif()
    list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
    list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
endif()
list(REMOVE_DUPLICATES _lib_env)

# Try to find the ptscotch lib in the given paths
# ----------------------------------------------

set(PTSCOTCH_libs_to_find "ptscotch;scotch;scotcherrexit")

# call cmake macro to find the lib path
if(PTSCOTCH_LIBDIR)
    foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
        set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
        find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
            NAMES ${ptscotch_lib}
            HINTS ${PTSCOTCH_LIBDIR})
    endforeach()
else()
    if(PTSCOTCH_DIR)
        foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
            set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
            find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
                NAMES ${ptscotch_lib}
                HINTS ${PTSCOTCH_DIR}
                PATH_SUFFIXES lib lib32 lib64)
        endforeach()
    else()
        foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})
            set(PTSCOTCH_${ptscotch_lib}_LIBRARY "PTSCOTCH_${ptscotch_lib}_LIBRARY-NOTFOUND")
            find_library(PTSCOTCH_${ptscotch_lib}_LIBRARY
                NAMES ${ptscotch_lib}
                HINTS ${_lib_env})
        endforeach()
    endif()
endif()

set(PTSCOTCH_LIBRARIES "")
set(PTSCOTCH_LIBRARY_DIRS "")
# If found, add path to cmake variable
# ------------------------------------
foreach(ptscotch_lib ${PTSCOTCH_libs_to_find})

    if (PTSCOTCH_${ptscotch_lib}_LIBRARY)
        get_filename_component(${ptscotch_lib}_lib_path "${PTSCOTCH_${ptscotch_lib}_LIBRARY}" PATH)
        # set cmake variables
        list(APPEND PTSCOTCH_LIBRARIES "${PTSCOTCH_${ptscotch_lib}_LIBRARY}")
        list(APPEND PTSCOTCH_LIBRARY_DIRS "${${ptscotch_lib}_lib_path}")
    else ()
        list(APPEND PTSCOTCH_LIBRARIES "${PTSCOTCH_${ptscotch_lib}_LIBRARY}")
        if (NOT PTSCOTCH_FIND_QUIETLY)
            message(STATUS "Looking for ptscotch -- lib ${ptscotch_lib} not found")
        endif()
    endif ()

    mark_as_advanced(PTSCOTCH_${ptscotch_lib}_LIBRARY)

endforeach()

# check a function to validate the find
if(PTSCOTCH_LIBRARIES)

    set(REQUIRED_INCDIRS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # PTSCOTCH
    if (PTSCOTCH_INCLUDE_DIRS)
        set(REQUIRED_INCDIRS  "${PTSCOTCH_INCLUDE_DIRS}")
    endif()
    if (PTSCOTCH_LIBRARY_DIRS)
        set(REQUIRED_LIBDIRS "${PTSCOTCH_LIBRARY_DIRS}")
    endif()
    set(REQUIRED_LIBS "${PTSCOTCH_LIBRARIES}")
    # MPI
    if (MPI_FOUND)
        if (MPI_C_INCLUDE_PATH)
            list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
        endif()
        if (MPI_C_LINK_FLAGS)
            list(APPEND REQUIRED_LIBS "${MPI_C_LINK_FLAGS}")
        endif()
        list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
    endif()
    # THREADS
    if(CMAKE_THREAD_LIBS_INIT)
        list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")
    endif()

    # set required libraries for link
    set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
    set(CMAKE_REQUIRED_LIBRARIES)
    foreach(lib_dir ${REQUIRED_LIBDIRS})
        list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
    endforeach()
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")

    # test link
    unset(PTSCOTCH_WORKS CACHE)
    include(CheckFunctionExists)
    check_function_exists(SCOTCH_dgraphInit PTSCOTCH_WORKS)
    mark_as_advanced(PTSCOTCH_WORKS)

    if(PTSCOTCH_WORKS)
        # save link with dependencies
        set(PTSCOTCH_LIBRARIES_DEP "${REQUIRED_LIBS}")
        set(PTSCOTCH_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
        set(PTSCOTCH_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    else()
        if(NOT PTSCOTCH_FIND_QUIETLY)
            message(STATUS "Looking for PTSCOTCH : test of SCOTCH_dgraphInit with PTSCOTCH library fails")
            message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
            message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
            message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
        endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
endif(PTSCOTCH_LIBRARIES)


# Check the size of SCOTCH_Num
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${PTSCOTCH_INCLUDE_DIRS})

include(CheckCSourceRuns)
#stdio.h and stdint.h should be included by scotch.h directly
set(PTSCOTCH_C_TEST_SCOTCH_Num_4 "
#include <stdio.h>
#include <stdint.h>
#include <ptscotch.h>
int main(int argc, char **argv) {
  if (sizeof(SCOTCH_Num) == 4)
    return 0;
  else
    return 1;
}
")

set(PTSCOTCH_C_TEST_SCOTCH_Num_8 "
#include <stdio.h>
#include <stdint.h>
#include <ptscotch.h>
int main(int argc, char **argv) {
  if (sizeof(SCOTCH_Num) == 8)
    return 0;
  else
    return 1;
}
")
check_c_source_runs("${PTSCOTCH_C_TEST_SCOTCH_Num_4}" PTSCOTCH_Num_4)
if(NOT PTSCOTCH_Num_4)
  check_c_source_runs("${PTSCOTCH_C_TEST_SCOTCH_Num_8}" PTSCOTCH_Num_8)
  if(NOT PTSCOTCH_Num_8)
    set(PTSCOTCH_INTSIZE -1)
  else()
    set(PTSCOTCH_INTSIZE 8)
  endif()
else()
  set(PTSCOTCH_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

# check that PTSCOTCH has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTSCOTCH DEFAULT_MSG
                                  PTSCOTCH_LIBRARIES
                                  PTSCOTCH_WORKS)
#
# TODO: Add possibility to check for specific functions in the library
#