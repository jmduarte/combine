cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++17 -m64 -O2 -L/usr/lib/root -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -lvdt -lRooFit -lRooFitCore -lRooStats -lMinuit -lFoam -lHistFactory -lboost_filesystem -lboost_program_options -lboost_system -rdynamic")


if(NOT COMBINE_HEADERS_DIR)
    set(COMBINE_HEADERS_DIR ${PROJECT_SOURCE_DIR}/include)
endif()

include_directories(${COMBINE_HEADERS_DIR})

include_directories($ENV{ROOTSYS}/include/root)

add_subdirectory(src)
#add_subdirectory(bin)

unset(COMBINE_HEADERS_DIR CACHE)
