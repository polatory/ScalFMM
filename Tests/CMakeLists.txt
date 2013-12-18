# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

project(Tests_scalfmm CXX)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})
ADD_DEFINITIONS( ${ScaLFMM_CXX_FLAGS})
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
		if( NOT ScalFMM_USE_BLAS )
			MESSAGE( STATUS "This needs BLAS = ${exec}" )
			set(compile_exec "FALSE")
		endif()
	endif()

	# Test FFT dependency
	file(STRINGS "${exec}" lines_fft REGEX "@FUSE_FFT")
	if(lines_fft)
		if( NOT ScalFMM_USE_FFT )
			MESSAGE( STATUS "This needs FFT = ${exec}" )
			set(compile_exec "FALSE")
		endif()
	endif()

	# Test MPI dependency
	file(STRINGS "${exec}" lines_mpi REGEX "@FUSE_MPI")
	if(lines_mpi)
		if( NOT ScalFMM_USE_MPI )
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
		SET_TARGET_PROPERTIES(${execname} PROPERTIES ENABLE_EXPORTS TRUE)
		target_link_libraries(
			${execname}
			${scalfmm_lib}
			${SCALFMM_LIBRARIES}
			)
	endif()
endforeach(exec)
