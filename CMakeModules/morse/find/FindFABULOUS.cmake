# - Try to find fabulous (C-api)
#
# If this script have difficulties to find fabulous, you can try to help
# it by setting the variable FABULOUS_DIR to the prefix path where fabulous
# was installed
#
# Once done this will define
#  FABULOUS_FOUND - System has fabulous
#  FABULOUS_INCLUDE_DIRS - The fabulous include directories
#  FABULOUS_MODULE_DIRS - The fabulous module directories for Fortran API
#  FABULOUS_LIBRARIES - The libraries needed to use fabulous
#  FABULOUS_DEFINITIONS - Compiler switches required for using fabulous
#
include(FindPackageHandleStandardArgs)

macro(FABULOUS_FIND_LIBRARIES_FROM_PKGCONFIG_RESULTS _prefix _pc_xprefix)
    foreach(_library ${${_pc_xprefix}_LIBRARIES})
        get_filename_component(_library ${_library} NAME_WE)
        unset(_library_path)
        unset(_library_path CACHE)
        find_library(_library_path NAMES ${_library}
            HINTS ${${_pc_xprefix}_LIBDIR} ${${_pc_xprefix}_LIBRARY_DIRS} )
        if (_library_path)
            list(APPEND ${_prefix}_LIBRARIES ${_library_path})
        else()
            message(FATAL_ERROR "Dependency of ${_prefix} '${_library}' NOT FOUND")
        endif()
        unset(_library_path CACHE)
    endforeach()
endmacro()

macro(FABULOUS_CHECK_FUNCTION_EXISTS _prefix _function)
    include(CheckFunctionExists)
    unset(${_prefix}_WORKS)
    unset(${_prefix}_WORKS CACHE)
    set(CMAKE_REQUIRED_LIBRARIES ${${_prefix}_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES ${${_prefix}_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_DEFINITIONS ${${_prefix}_DEFINITIONS})
    check_function_exists(${_function} ${_prefix}_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES "")
    set(CMAKE_REQUIRED_INCLUDES "")
    set(CMAKE_REQUIRED_DEFINITIONS "")
    mark_as_advanced(${_prefix}_WORKS)
endmacro()

find_package(PkgConfig QUIET)

set(ENV_FABULOUS_DIR "$ENV{FABULOUS_DIR}")
set(ENV_FABULOUS_INCDIR "$ENV{FABULOUS_INCDIR}")
set(ENV_FABULOUS_LIBDIR "$ENV{FABULOUS_LIBDIR}")
set(FABULOUS_GIVEN_BY_USER "FALSE")
if ( FABULOUS_DIR OR ENV_FABULOUS_DIR
        OR ( FABULOUS_INCDIR AND FABULOUS_LIBDIR )
        OR ( ENV_FABULOUS_INCDIR AND ENV_FABULOUS_LIBDIR ) )
    set(FABULOUS_GIVEN_BY_USER "TRUE")
endif()

set(FABULOUS_STATIC_FIND_QUIETLY "TRUE")
set(FABULOUS_SHARED_FIND_QUIETLY "TRUE")

if ((NOT FABULOUS_FOUND) AND (NOT FABULOUS_GIVEN_BY_USER) AND PKG_CONFIG_FOUND)
    pkg_check_modules(PC_FABULOUS QUIET fabulous)

    find_path(FABULOUS_STATIC_INCLUDE_DIR NAMES fabulous.h
        HINTS ${PC_FABULOUS_STATIC_INCLUDEDIR} ${PC_FABULOUS_STATIC_INCLUDE_DIRS} )
    find_library(FABULOUS_STATIC_LIBRARY NAMES libfabulous.a
        HINTS ${PC_FABULOUS_STATIC_LIBDIR} ${PC_FABULOUS_STATIC_LIBRARY_DIRS} )

    find_path(FABULOUS_SHARED_INCLUDE_DIR NAMES fabulous.h
        HINTS ${PC_FABULOUS_INCLUDEDIR} ${PC_FABULOUS_INCLUDE_DIRS} )
    find_library(FABULOUS_SHARED_LIBRARY NAMES libfabulous.so
        HINTS ${PC_FABULOUS_LIBDIR} ${PC_FABULOUS_LIBRARY_DIRS} )

    # handle the QUIETLY and REQUIRED arguments and set FABULOUS_FOUND to TRUE
    # if all listed variables are TRUE

    find_package_handle_standard_args(
        FABULOUS_STATIC DEFAULT_MSG
        FABULOUS_STATIC_LIBRARY FABULOUS_STATIC_INCLUDE_DIR)
    mark_as_advanced(FABULOUS_STATIC_INCLUDE_DIR FABULOUS_STATIC_LIBRARY)

    find_package_handle_standard_args(
        FABULOUS_SHARED DEFAULT_MSG
        FABULOUS_SHARED_LIBRARY FABULOUS_SHARED_INCLUDE_DIR)
    mark_as_advanced(FABULOUS_SHARED_INCLUDE_DIR FABULOUS_SHARED_LIBRARY)

    if (FABULOUS_STATIC_FOUND AND NOT FABULOUS_SHARED_FOUND)
        set(FABULOUS_INCLUDE_DIRS ${FABULOUS_STATIC_INCLUDE_DIR} )
        set(FABULOUS_DEFINITIONS ${PC_FABULOUS_STATIC_CFLAGS_OTHER} )
        set(FABULOUS_LIBRARIES "")
        fabulous_find_libraries_from_pkgconfig_results(FABULOUS PC_FABULOUS_STATIC)
    elseif(FABULOUS_SHARED_FOUND)
        set(FABULOUS_INCLUDE_DIRS ${FABULOUS_INCLUDE_DIR} )
        set(FABULOUS_DEFINITIONS ${PC_FABULOUS_CFLAGS_OTHER} )
        set(FABULOUS_LIBRARIES "")
        fabulous_find_libraries_from_pkgconfig_results(FABULOUS PC_FABULOUS)
    endif()
    fabulous_check_function_exists(FABULOUS fabulous_create)
    find_package_handle_standard_args(
        FABULOUS DEFAULT_MSG
        FABULOUS_LIBRARIES FABULOUS_INCLUDE_DIRS FABULOUS_WORKS)
endif()

if ((NOT FABULOUS_FOUND) AND (FABULOUS_GIVEN_BY_USER OR (NOT PKG_CONFIG_FOUND)))

    # Currently the C-api (compiled version) does not depent on chameleon
    # so the library only depends on CBLAS AND LAPACKE

    set(FABULOUS_DEFINITIONS "")
    if ( ( FABULOUS_INCDIR AND FABULOUS_LIBDIR ) OR ( ENV_FABULOUS_INCDIR AND ENV_FABULOUS_LIBDIR ) )

        if ((NOT FABULOUS_LIBDIR) AND (NOT FABULOUS_INCDIR)
                AND (ENV_FABULOUS_INCDIR AND ENV_FABULOUS_LIBDIR) )
            set(FABULOUS_LIBDIR ${ENV_FABULOUS_LIBDIR})
            set(FABULOUS_INCDIR ${ENV_FABULOUS_INCDIR})
        endif()

        find_path(FABULOUS_INCLUDE_DIRS NAMES fabulous.h HINTS ${FABULOUS_INCDIR})
        find_path(FABULOUS_MODULE_DIRS NAMES fabulous_mod.mod HINTS ${FABULOUS_INCDIR}
          PATH_SUFFIXES fabulous)
        find_library(FABULOUS_STATIC_LIBRARY NAMES libfabulous.a HINTS ${FABULOUS_LIBDIR})
        find_library(FABULOUS_SHARED_LIBRARY NAMES libfabulous.so HINTS ${FABULOUS_LIBDIR})
    else()
        if (ENV_FABULOUS_DIR AND NOT FABULOUS_DIR)
            set(FABULOUS_DIR "${ENV_FABULOUS_DIR}" CACHE PATH "Installation prefix where fabulous is installed")
        else()
            set(FABULOUS_DIR "${FABULOUS_DIR}" CACHE PATH "Installation prefix where fabulous is installed")
        endif()

        find_path(FABULOUS_INCLUDE_DIRS NAMES fabulous.h
            HINTS ${FABULOUS_DIR}
            PATH_SUFFIXES include include/fabulous)
        find_path(FABULOUS_MODULE_DIRS NAMES fabulous_mod.mod
            HINTS ${FABULOUS_DIR}
            PATH_SUFFIXES include include/fabulous)
        find_library(FABULOUS_STATIC_LIBRARY NAMES libfabulous.a
            HINTS ${FABULOUS_DIR}
            PATH_SUFFIXES lib lib32 lib64 lib/fabulous lib32/fabulous lib64/fabulous)
        find_library(FABULOUS_SHARED_LIBRARY NAMES libfabulous.so
            HINTS ${FABULOUS_DIR}
            PATH_SUFFIXES lib lib32 lib64 lib/fabulous lib32/fabulous lib64/fabulous)
    endif()

    find_package_handle_standard_args(FABULOUS_STATIC DEFAULT_MSG FABULOUS_STATIC_LIBRARY)
    find_package_handle_standard_args(FABULOUS_SHARED DEFAULT_MSG FABULOUS_SHARED_LIBRARY)
    mark_as_advanced(FABULOUS_STATIC_LIBRARY FABULOUS_SHARED_LIBRARY)

    if (FABULOUS_FIND_REQUIRED)
      find_package(CBLAS REQUIRED)
      find_package(LAPACKE REQUIRED)
    else()
      find_package(CBLAS)
      find_package(LAPACKE)
    endif()

    if (FABULOUS_STATIC_FOUND AND NOT FABULOUS_SHARED_FOUND)
        set(FABULOUS_LIBRARIES ${FABULOUS_STATIC_LIBRARY} ${CBLAS_LIBRARIES} ${LAPACKE_LIBRARIES} "-lstdc++" "-lm")
    elseif(FABULOUS_SHARED_FOUND)
        set(FABULOUS_LIBRARIES ${FABULOUS_SHARED_LIBRARY} ${CBLAS_LIBRARIES} ${LAPACKE_LIBRARIES})
    endif()
    fabulous_check_function_exists(FABULOUS fabulous_create)
    find_package_handle_standard_args(
        FABULOUS DEFAULT_MSG
        FABULOUS_LIBRARIES FABULOUS_INCLUDE_DIRS FABULOUS_WORKS)
endif()
