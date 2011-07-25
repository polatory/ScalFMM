cmake_minimum_required(VERSION 2.6)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

# Configuration du projet
project(Lib_scalfmm)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

# Test if openmp is here
include(FindOpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

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


# Install lib
install( TARGETS scalfmm ARCHIVE DESTINATION lib )

# Install header
SET(my_include_dirs "." "Components" "Containers" "Core" "Extensions" "Files" "Fmb" "Utils")
FOREACH(my_dir ${my_include_dirs})	
	file(GLOB 	 	 	 	
	 	hpp_in_dir
		${my_dir}/*.hpp
	)
	INSTALL( FILES ${hpp_in_dir} DESTINATION include/ScalFmm/${my_dir} )
ENDFOREACH()
INSTALL( FILES "${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h" DESTINATION include/ScalFmm/Utils/${my_dir} )

