# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto scalfmm/Build and tapes cmake ../")
endif(insource)

# Configuration du projet
project(Lib_scalfmm)
set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

ADD_DEFINITIONS(${ScaLFMM_CXX_FLAGS})
MESSAGE( ${ScaLFMM_CXX_FLAGS})
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

# Add blas library (even if it is set to off)
target_link_libraries(
	scalfmm
	${scalfmm_lib}
	${SCALFMM_LIBRARIES}
)

# Adding the entire project dir as an include dir
INCLUDE_DIRECTORIES(
     ${CMAKE_SOURCE_DIR}/Src 
     ${CMAKE_BINARY_DIR}/Src 
)


# Install lib
install( TARGETS scalfmm ARCHIVE DESTINATION lib )

# Install header
SET(my_include_dirs "."
  "Arranger"        "Core"           "Utils"
  "Kernels/Chebyshev"       "Components"       "Extensions"  
  "Containers"       "Files"       "Kernels/Spherical"
  "Kernels/P2P" "Kernels/Taylor" "Kernels/ROtation")

FOREACH(my_dir ${my_include_dirs})	
	file(GLOB 	 	 	 	
	 	hpp_in_dir
		${my_dir}/*.hpp ${my_dir}/*.h
	)
	INSTALL( FILES ${hpp_in_dir} DESTINATION include/${my_dir} )
ENDFOREACH()

INSTALL( FILES "${CMAKE_BINARY_DIR}/Src/ScalFmmConfig.h" DESTINATION include/Utils/${my_dir} )

