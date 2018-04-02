# project config stuff
#  
# Author: Tomas Hillberg

# simple include guard
if( DEFINED WFWFS_CONFIGURATION_LOADED )
    return()
endif()
set( WFWFS_CONFIGURATION_LOADED True )

# Enable empty ELSE and ENDIF's
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# Load compiler specific configuration
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
    include(${CMAKE_CURRENT_LIST_DIR}/clang.cmake)
elseif( CMAKE_COMPILER_IS_GNUCC )
    include(${CMAKE_CURRENT_LIST_DIR}/gcc.cmake)
else()
    message(FATAL_ERROR "Unsupported compiler")
endif()


include( ${CMAKE_CURRENT_LIST_DIR}/macros.cmake )
    
get_filename_component(WFWFS_DIR ${WFWFS_DIR} REALPATH)

include_directories(${WFWFS_DIR}/src)

message(STATUS "The WFWFS will be installed into: \"${CMAKE_INSTALL_PREFIX}\"")

