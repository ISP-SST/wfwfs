USE_EXTERNAL( boost fits gsl pleora threads zlib )

if( NOT DEFINED Boost_FOUND )
    message(STATUS "Boost not found. Try your systems equivalent of \"apt-get install libboost-all-dev\"" )
    return()
endif()

if( DEFINED Fits_FOUND )
    option( WFWFS_WITH_FITS "Build with FITS support" ON )
    if( WFWFS_WITH_FITS )
        message( STATUS "Building with FITS support" )
        add_definitions( -DWFWFS_WITH_FITS )
    endif()
else()
    unset( WFWFS_WITH_FITS CACHE )
    message( WARNING "CFITSIO not found. Try your systems equivalent of \"apt-get install libcfitsio-dev\"" )
    return()
endif()

if( DEFINED Pleora_FOUND )
    option( WFWFS_WITH_PLEORA "Build with Pleora support" ON )
    if( WFWFS_WITH_PLEORA )
        message( STATUS "Building with Pleora support" )
        add_definitions( -DWFWFS_WITH_PLEORA )
    endif()
else()
    unset( WFWFS_WITH_PLEORA CACHE )
    message( WARNING "Pleora not found. Install the Pleora-SDK or give a hint to cmake as: \"-DPleora_DIR=/your/path/to/pleora/\"" )
    return()
endif()

##### revison info ######
set( WFWFS_REVISION_SRC_FILE "${CMAKE_CURRENT_BINARY_DIR}/revision.cpp" )
set_source_files_properties( ${WFWFS_REVISION_SRC_FILE} GENERATED )
CHECK_REVISION( WFWFS "${CMAKE_CURRENT_LIST_DIR}/revision.cpp.template" "${WFWFS_REVISION_SRC_FILE}")

file( GLOB WFWFS_CPP "*.cpp" )
WFWFS_ADD_EXECUTABLE( wfwfs ${WFWFS_CPP} ${WFWFS_REVISION_SRC_FILE} )

if( TARGET GIT_CHECK_WFWFS )
    add_dependencies( wfwfs GIT_CHECK_WFWFS )
endif()

