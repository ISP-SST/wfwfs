cmake_minimum_required( VERSION  2.8 )

IF( POLICY CMP0048 )  # The ``project()`` command manages VERSION variables. cf. cmake --help-policy CMP0048
    CMAKE_POLICY( SET CMP0048 OLD )   # We manage versions manually
ENDIF()

project( WFWFS )

set( WFWFS_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE PATH "Path to source tree." )
mark_as_advanced( WFWFS_DIR )

include( ${WFWFS_DIR}/cmake/config.cmake )

add_subdirectory( src )


