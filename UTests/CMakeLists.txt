cmake_minimum_required (VERSION 2.8.5)

# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

project(UTests_scalfmm CXX)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})

# Find all code files
file(	
	GLOB_RECURSE
	source_tests_files
	./*.cpp
)

# Adding the project sources dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_BINARY_DIR}/Src 
)

# Add execs - 1 cpp = 1 exec
foreach(exec ${source_tests_files}) 
	get_filename_component(
		execname ${exec}
		NAME_WE
	)

	set(compile_exec "TRUE")	

	# Test Blas dependency
	file(STRINGS "${exec}" lines_blas REGEX "@FUSE_BLAS")
	if(lines_blas)
		if( NOT SCALFMM_USE_BLAS )
			MESSAGE( STATUS "This needs BLAS = ${exec}" )
			set(compile_exec "FALSE")
		endif()
	endif()

	# Test MPI dependency
	file(STRINGS "${exec}" lines_mpi REGEX "@FUSE_MPI")
	if(lines_mpi)
		if( NOT SCALFMM_USE_MPI )
			MESSAGE( STATUS "This needs MPI = ${exec}" )
			set(compile_exec "FALSE")
		endif()
	endif()

	# Dependency are OK
	if( compile_exec )
		add_executable(
			${execname}
			${exec}
		)
		target_link_libraries(
			${execname}
			${scalfmm_lib}
			${BLAS_LIBRARIES}
			${LAPACK_LIBRARIES}
		)

		add_test(
			${execname}
			${CMAKE_CURRENT_BINARY_DIR}/${execname}
		)
	endif()
endforeach(exec)
