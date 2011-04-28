cmake_minimum_required(VERSION 2.6)

# Configuration du projet
project(Lib_scalfmm)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

# Searching all cpp file
file(
	GLOB_RECURSE
	source_lib_files
	./*.cpp
)

# Adding cpp files to project
add_library(
	scalfmm
	STATIC
	${source_lib_files}
)

# Adding the entire project dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_BINARY_DIR}/Src 
)

