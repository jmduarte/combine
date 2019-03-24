cmake_minimum_required(VERSION 3.2)
project(combine CXX)

add_executable(combine combine.cpp)
target_link_libraries(combine PUBLIC Combine)

install(TARGETS combine DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
