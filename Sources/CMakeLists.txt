cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(Lib_scalfmm)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

#Configuration de la bibliothèque
file(
	GLOB_RECURSE
	source_lib_files
	./*.cpp
)

add_library(
	scalfmm
	STATIC
	${source_lib_files}
)

INCLUDE_DIRECTORIES(
     ${CMAKE_BINARY_DIR}/Sources 
)

