cmake_minimum_required (VERSION 2.6)

# Project Infos
project(FMB_All CXX)

# Options
OPTION( BUILD_TESTS "Set to ON to build fonctionnalities Tests" ON )
OPTION( BUILD_UTESTS "Set to ON to build UTests" ON )

ADD_DEFINITIONS(-O2 -Wall -W -Wshadow -Wpointer-arith -Wcast-qual)

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
