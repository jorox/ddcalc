cmake_minimum_required(VERSION 2.8)
project( ddcalc )

include_directories ( ${CMAKE_SOURCE_DIR}/src )
MESSAGE( STATUS "CMAKE_SOURCE_DIR:         " ${CMAKE_SOURCE_DIR} )
MESSAGE( STATUS "CMAKE_INCLUDE_PATH: " ${CMAKE_INCLUDE_PATH} )

add_compile_options(-std=c++11)
set(CMAKE_CXX_STANDARD 11)
file(GLOB SOURCES "src/*.cpp")

add_executable( ddcalc ${SOURCES} )
