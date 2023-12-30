cmake_minimum_required(VERSION 3.16.3)
project(ES)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "-O3")
SET(CMAKE_C_FLAGS "-O3")
file(GLOB source_cpp src/*.cpp)
file(GLOB source_h include/*.h)
file(GLOB source_c src/*.c)

#
# sudo apt install libboost1.74-all-dev
#
find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
find_package(pagmo REQUIRED)

include_directories(include)
include_directories(liblsoda/src)

# Create a target out of the library compilation result
add_custom_target(lsoda_target
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/liblsoda/src
    COMMAND make
    COMMENT "Building liblsoda")

# Create an library target out of the library compilation result
add_library(liblsoda STATIC IMPORTED GLOBAL)
add_dependencies(liblsoda lsoda_target)
set_target_properties(liblsoda
    PROPERTIES
    IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/liblsoda/src/liblsoda.so
    IMPORTED_IMPLIB ${CMAKE_SOURCE_DIR}/liblsoda/src/liblsoda.so
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/liblsoda/src)

add_executable(ES ${source_cpp} ${source_c} src/GRNEDOHelpers.cpp)
target_link_libraries(ES PUBLIC liblsoda Pagmo::pagmo)

add_executable(test
        include/ProblemDescription.h
        include/GRNSeries.h
        src/GRNSeries.cpp
        include/GRNEDOHelpers.h
        src/GRNEDOHelpers.cpp
        src/read.cpp
        src/rk4.c
        src/Individual.cpp
        include/Individual.h
        tests/unit_tests.cpp)
target_link_libraries(test PUBLIC liblsoda Pagmo::pagmo
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_SYSTEM_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

add_executable(test2 tests/test2.cpp)
target_link_libraries(test2 PUBLIC liblsoda)

