# Global configuration for the clang compiler
#
# Author: Tomas Hillberg

message(STATUS "Loading clang-specific configuration. (${CMAKE_CURRENT_LIST_FILE})")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -ferror-limit=1")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lc++")
