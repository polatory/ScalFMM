cmake_minimum_required (VERSION 2.6)

project(UTests_scalfmm CXX)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})

# Test if openmp is here
include(FindOpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Link with fmb lib
set(fmb_lib scalfmm)

# Find all code files
file(	
	GLOB_RECURSE
	source_tests_files
	./*.cpp
)

# Add execs - 1 cpp = 1 exec
foreach(exec ${source_tests_files}) 
	get_filename_component(
		execname ${exec}
		NAME_WE
	)
	add_executable(
		${execname}
		${exec}
	)
	target_link_libraries(
		${execname}
		${fmb_lib}
	)
endforeach(exec)
