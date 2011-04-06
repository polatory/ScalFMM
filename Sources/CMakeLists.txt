cmake_minimum_required(VERSION 2.6)

#Configuration du projet
project(FMB++)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

#Configuration de la bibliothèque
file(
	GLOB_RECURSE
	source_lib_files
	./*.cpp
)
add_library(
	FMB
	STATIC
	${source_lib_files}
)
