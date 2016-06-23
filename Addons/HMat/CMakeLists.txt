# check if compiling into source directories
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" insource)
if(insource)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build. Goto ./Build and tapes cmake ../")
endif(insource)

project(ADDONS_HMAT_SCALFMM CXX C)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SCALFMM_CXX_FLAGS}")

# Active language
# -----------------------
ENABLE_LANGUAGE(CXX C)
MESSAGE(STATUS " CXX  ${CMAKE_CXX_COMPILER_ID}" )

# Options
OPTION( SCALFMM_ADDON_HMAT "Set to ON to build ScalFMM HMat interface" OFF )

# if ask to build addon
if(SCALFMM_ADDON_HMAT)
        list(APPEND FUSE_LIST "SCOTCH")
        set(SCOTCH_NOT_FOUND $ENV{SCOTCH_LIB})
        if(NOT SCOTCH_NOT_FOUND) 
            MESSAGE(STATUS " SCOTCH not found " )
        else()
            OPTION( SCALFMM_USE_SCOTCH "Set to ON to use scotch" ON )
            if (SCALFMM_USE_SCOTCH)
               include_directories($ENV{SCOTCH_INC})
               SET(SCALFMM_LIBRARIES "${SCALFMM_LIBRARIES};-L$ENV{SCOTCH_LIB};-lscotch;-lscotcherr")        
            endif()
        endif()

        # first build lib scalfmmhmat
        set(LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

        # Searching all cpp file
        file( GLOB_RECURSE source_lib_files Src/*.cpp )

        # Adding cpp files to project
        add_library( scalfmmhmat ${source_lib_files} )
        
        # Add blas library (even if it is set to off)
        target_link_libraries( scalfmmhmat scalfmm)

        # Adding the entire project dir as an include dir
        INCLUDE_DIRECTORIES(
             ${SCALFMM_BINARY_DIR}/Src
             ${SCALFMM_SOURCE_DIR}/Src
             ${SCALFMM_INCLUDES}
        )

        # Install lib
        install( TARGETS scalfmmhmat ARCHIVE DESTINATION lib )

        # Install headers
        SET(my_include_dirs "Src/Blocks" "Src/Clustering" "Src/Containers" "Src/Utils" "Src/Viewers" "CClusteringLibrary")

        FOREACH(my_dir ${my_include_dirs})  
            file(GLOB               
                hpp_in_dir
                ${my_dir}/*.hpp ${my_dir}/*.h
            )
            INSTALL( FILES ${hpp_in_dir} DESTINATION include/ScalFmm/HMat/${my_dir} )
        ENDFOREACH()

        # Add C Clustering Library
        file( GLOB_RECURSE ccl_lib_files CClusteringLibrary/*.c )
        add_library( cclusteringlib ${ccl_lib_files} )
        INCLUDE_DIRECTORIES(CClusteringLibrary/)
        target_link_libraries( cclusteringlib scalfmm)
        install( TARGETS cclusteringlib ARCHIVE DESTINATION lib )

        # Tests
        file( GLOB_RECURSE source_tests_files Tests/*.cpp )
        INCLUDE_DIRECTORIES( ${SCALFMM_BINARY_DIR}/Src )

        # Then build test files
        SET(hmat_list_execs "")
        foreach(exec ${source_tests_files})
                get_filename_component(
                        execname ${exec}
                        NAME_WE
                )

            set(compile_exec "TRUE")

            foreach(fuse_key ${FUSE_LIST})
                file(STRINGS "${exec}" lines_fuse REGEX "@FUSE_${fuse_key}")
                if(lines_fuse)
                        if( NOT SCALFMM_USE_${fuse_key} )
                                MESSAGE( STATUS "This needs ${fuse_key} = ${exec}" )
                                set(compile_exec "FALSE")
                        endif()
                endif()
            endforeach()

                # Dependency are OK
                if( compile_exec )
                        add_executable(	${execname} ${exec} )
                        # link to scalfmm and scalfmmhmat and cclusteringlib
                        target_link_libraries(
                                ${execname}
                                ${scalfmm_lib}
                                scalfmmhmat
                                cclusteringlib
                                ${SCALFMM_LIBRARIES}
                                )
                         LIST(APPEND hmat_list_execs ${execname})
                endif()
        endforeach(exec)
        
        add_custom_target(hmat DEPENDS ${hmat_list_execs})
endif()
