# CMake convenience macros
#  
# Author: Tomas Hillberg


macro(APPEND_PATHS_UNIQUE LIST)

    list( LENGTH ${LIST} SIZE )
    foreach( elm ${ARGN} )
        if( EXISTS ${elm} )
            get_filename_component( elm2 ${elm} REALPATH )
            if( SIZE EQUAL 0 )               
                set( ${LIST} ${elm2} CACHE INTERNAL "" )
                set( SIZE 1 )
            else()
                list( APPEND ${LIST} ${elm2} )
            endif()
        endif()
    endforeach()

    list( LENGTH ${LIST} SIZE )
    if( SIZE GREATER 0 )
        list( REMOVE_DUPLICATES ${LIST} )
    endif()
endmacro()


macro( APPEND_LIBS_UNIQUE LIST ELEMENT )
    set( prefix "" )
    foreach( elm ${ELEMENT} ${ARGN} )
        if( "${elm}" STREQUAL "" )
        elseif( "${elm}" STREQUAL "optimized" OR "${elm}" STREQUAL "debug" )
            set( prefix "${elm};" )
        else()
            list( FIND "${LIST}" "${elm}" FOUND )
            if( FOUND EQUAL -1 )
                list( APPEND ${LIST} "${prefix}${elm}" )
            endif()
            set( prefix "" )
        endif()
    endforeach()
endmacro()


macro( USE_EXTERNAL ... )
    foreach( ext ${ARGV} )
        file( TO_CMAKE_PATH "${WFWFS_DIR}/cmake/use_${ext}.cmake" EXT_CONFIG )
        if( EXISTS "${EXT_CONFIG}" )
            include( "${EXT_CONFIG}" )
        else()
            file( TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/cmake/use_${ext}.cmake" EXT_CONFIG2 )
            if( EXISTS "${EXT_CONFIG2}" )
                include( "${EXT_CONFIG2}" )
            else()
                message( WARNING "Can't locate configuration file for ${ext} (\"${EXT_CONFIG}\")" )
            endif()
        endif()
    endforeach()
endmacro()


# -----------------------------------------------------------------------------
# A couple of convenience functions 
# -----------------------------------------------------------------------------
macro( WFWFS_SETDIRS )
    if( DEFINED WFWFS_CURRENT_INCLUDES )
        list( REMOVE_DUPLICATES WFWFS_CURRENT_INCLUDES )
    endif()
    if( DEFINED WFWFS_CURRENT_LIBDIRS )
        list( REMOVE_DUPLICATES WFWFS_CURRENT_LIBDIRS )
    endif()
    include_directories( ${WFWFS_CURRENT_INCLUDES} )
    link_directories( ${WFWFS_CURRENT_LIBDIRS} )
endmacro()


macro(WFWFS_ADD_EXECUTABLE target)
    WFWFS_SETDIRS()
    add_executable( ${target} ${ARGN} )
    get_target_property( src ${target} SOURCES )
#    target_link_libraries( ${target} ${WFWFS_CURRENT_LIBRARIES} ${WFWFS_CURRENT_LIBRARIES} ) # WFWFS_CURRENT_LIBRARIES is needed a 2nd time because of circular dependencies
    target_link_libraries( ${target} ${WFWFS_CURRENT_LIBRARIES} )
    set_target_properties( ${target} PROPERTIES DEBUG_POSTFIX "-d" )
endmacro()


macro(WFWFS_ADD_LIBRARY target)
    WFWFS_SETDIRS()
    add_library( ${target} ${ARGN} )
    get_target_property( src ${target} SOURCES )
    # TODO: Eventually the install-dir should really be user configurable.
    set_target_properties( ${target} PROPERTIES DEBUG_POSTFIX "-d" )
    install( TARGETS ${target} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
endmacro()



# -----------------------------------------------------------------------------
# Macro for generic (git) versioning
# The template file will be configured with the following variables defined:
#     REVISION_COMMIT_TIME      (timestamp as string)
#     REVISION_COMMIT_COMMENT   (string)
#     REVISION_BUILD_TIME       (timestamp as string)
#     REVISION_ID               (string, e.g hash-tag)
#     REVISION_DESCRIPTION      (string, id + comment)
# -----------------------------------------------------------------------------
macro( CHECK_REVISION REVNAME REVISION_TEMPLATE_LOCATION REVISION_FILE_LOCATION ) # optional 4th argument = path

    if( IS_DIRECTORY "${ARGN}" ) # default path is current project
        set( GIT_WORKTREE "${ARGN}" )
    else()
        set( GIT_WORKTREE "${PROJECT_SOURCE_DIR}" )
    endif()
    
    find_program( GIT_EXECUTABLE NAMES "git" NO_CMAKE_PATH)
    if( EXISTS "${GIT_EXECUTABLE}" AND EXISTS "${GIT_WORKTREE}/.git" )
        add_custom_target( GIT_CHECK_${REVNAME}
                        COMMAND ${CMAKE_COMMAND}
                            -DREVNAME:STRING="${REVNAME}"
                            -DGIT_DIR:STRING="${GIT_WORKTREE}"
                            -DREVISION_TEMPLATE_LOCATION="${REVISION_TEMPLATE_LOCATION}"
                            -DREVISION_FILE_LOCATION="${REVISION_FILE_LOCATION}"
                            -DBUILD_DIR="${CMAKE_BINARY_DIR}"
                            -P "${WFWFS_DIR}/cmake/revision.cmake"
                            ${CMAKE_BINARY_DIR}
        )
        set_source_files_properties( ${REVISION_FILE_LOCATION} PROPERTIES GENERATED 1 )
    else()
        message( STATUS "Couldn't find git executable and/or .git directory" )
        # Set default values for generating file without git data.
        set( REVISION_COMMIT_TIME "unknown" )
        set( REVISION_BUILD_TIME "unknown" )
        set( REVISION_ID "unknown" )
        set( REVISION_COMMIT_COMMENT "No revision control" )
        set( REVISION_DESCRIPTION "No revision control" )
        set( ${REVNAME}_VERSION_MAJOR 0 )
        set( ${REVNAME}_VERSION_PATCH 0 )
        set( ${REVNAME}_VERSION_COMMIT 0 )
        configure_file( "${REVISION_TEMPLATE_LOCATION}" "${REVISION_FILE_LOCATION}" )
    endif()

endmacro()

