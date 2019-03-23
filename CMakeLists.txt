cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17 -m64 -O2 -L/usr/lib/root -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -lvdt -lRooFit -lRooFitCore -lRooStats -lMinuit -lFoam -lHistFactory -lboost_filesystem -lboost_program_options -lboost_system -rdynamic")

set(COMBINE_HEADERS_DIR ${PROJECT_SOURCE_DIR}/interface)

set(COMBINE_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(COMBINE_INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(COMBINE_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)

include_directories(${COMBINE_HEADERS_DIR})

include_directories($ENV{ROOTSYS}/include/root)

add_subdirectory(src)
add_subdirectory(bin)
