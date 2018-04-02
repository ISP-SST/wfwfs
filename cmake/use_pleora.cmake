#
# Set input data for FindExternal.cmake
#

set(EXT_NAME "Pleora")
#set(EXT_DEBUG 1)

if( NOT EXISTS "${Pleora_DIR}" ) # check for the GENICAM enviroment to use as hint
    if( IS_DIRECTORY "$ENV{GENICAM_ROOT_V3_0}" )
        set( Pleora_DIR "$ENV{GENICAM_ROOT_V3_0}/../../" )
    elseif( IS_DIRECTORY "$ENV{GENICAM_ROOT}" )
        set( Pleora_DIR "$ENV{GENICAM_ROOT}/../../" )
    endif()
endif()

# add some default locations to the search
set( EXT_HINT ${Pleora_DIR}
              "/opt/pleora/ebus_sdk/*/"
              "/usr/local/pleora/ebus_sdk/*/"
              "/usr/pleora/ebus_sdk/*/"
)

set(EXT_HEADER_FILE PvVersion.h)
set(EXT_VERSION_FILE PvVersion.h)
set(EXT_COMPONENTS PvBase PvDevice PvBuffer PvGenICam PvStream PvSystem)

set(EXT_MAJOR_REGEXP "VERSION_MAJOR")
set(EXT_MINOR_REGEXP "VERSION_MINOR")
set(EXT_PATCH_REGEXP "VERSION_SUB")


# Attempt to locate libs/headers automagically
include("${CMAKE_CURRENT_LIST_DIR}/FindExternal.cmake")


appendPaths()
