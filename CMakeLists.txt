cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17 -m64 -O3")

#set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

#set(DIVISIBLE_INSTALL_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
#set(DIVISIBLE_INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
#set(DIVISIBLE_INSTALL_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)

set(COMBINE_HEADERS_DIR ${PROJECT_SOURCE_DIR}/interface)

#include_directories(${DIVISIBLE_INSTALL_INCLUDE_DIR})
include_directories(${COMBINE_HEADERS_DIR})

#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT REQUIRED COMPONENTS RIO Net)
include_directories($ENV{ROOTSYS}/include/root)

#---Set Link and include directories--------------------------------------------------------------
#include_directories(${ROOT_INCLUDE_DIRS})
#link_directories(${ROOT_LIBRARY_DIR})

##---Set Flags-------------------------------------------------------------------------------------
#add_definitions(${ROOT_DEFINITIONS})
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ROOT_C_FLAGS}")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ROOT_fortran_FLAGS}")

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{ROOTSYS}/etc/cmake)
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#message(STATUS "CMAKE_PREFIX_PATH: " ${CMAKE_PREFIX_PATH})
#message(STATUS "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
#find_package(ROOT REQUIRED COMPONENTS Graf3d Graf Gpad MathCore RIO Hist Tree Net Thread Core)

add_subdirectory(src)
#add_subdirectory(src_test)
