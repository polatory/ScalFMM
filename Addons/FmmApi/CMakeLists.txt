# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

project(Addons_fmmapi_scalfmm CXX)

# Active language
# -----------------------
ENABLE_LANGUAGE(CXX )
MESSAGE(STATUS " CXX  ${CMAKE_CXX_COMPILER_ID}" )

# Options
OPTION( SCALFMM_ADDON_FMMAPI "Set to ON to build ScaFMM FMM API interface" OFF )

# if ask to build addon
if(SCALFMM_ADDON_FMMAPI)
	# first build lib scalfmmapi
	set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

	# Searching all cpp file
	file( GLOB_RECURSE source_lib_files Src/*.cpp )

	# Adding cpp files to project
	add_library( scalfmmapi STATIC ${source_lib_files} )

	# Add blas library (even if it is set to off)
	target_link_libraries( scalfmmapi scalfmm)

	# Adding the entire project dir as an include dir
	INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/Src )

	# Install lib
	install( TARGETS scalfmmapi ARCHIVE DESTINATION lib )

	# Install header
	SET(my_include_dirs "Src")
	
	file( GLOB hpp_in_dir Src/*.hpp )
	INSTALL( FILES ${hpp_in_dir} DESTINATION include/ScalFmm/FmmApi )

	file( GLOB_RECURSE source_tests_files Tests/*.cpp )
	INCLUDE_DIRECTORIES( ${CMAKE_BINARY_DIR}/Src )

	# Then build test files
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
			add_executable(	${execname} ${exec} )
			# link to scalfmm and scalfmmapi
			target_link_libraries(
				${execname}
				${scalfmm_lib}
				scalfmmapi
				${BLAS_LIBRARIES}
				${LAPACK_LIBRARIES}
				)
		endif()
	endforeach(exec)
endif()

