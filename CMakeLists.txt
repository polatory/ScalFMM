cmake_minimum_required (VERSION 2.6)

# Project Infos
project(scalfmm CXX)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
IF(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
ENDIF(insource)

# Options
OPTION( BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( BUILD_UTESTS "Set to ON to build UTests" ON )

ADD_DEFINITIONS(-O2 -Wall -Wshadow -Wpointer-arith -Wcast-qual)

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
