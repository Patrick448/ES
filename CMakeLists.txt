cmake_minimum_required(VERSION 3.16.3)
project(ES)

set(CMAKE_CXX_STANDARD 14)

file(GLOB source_cpp src/*.cpp)
file(GLOB source_h include/*.h)
file(GLOB source_c src/*.c)


include_directories(include)
include_directories(liblsoda/src)


# Create a target out of the library compilation result
add_custom_target(lsoda_target
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/liblsoda/src
    COMMAND make
    COMMENT "Building liblsoda")

# Create an library target out of the library compilation result
add_library(liblsoda STATIC IMPORTED GLOBAL)
add_dependencies(liblsoda lsoda_target)
set_target_properties(liblsoda
    PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/liblsoda/src/liblsoda.so
    IMPORTED_IMPLIB ${CMAKE_CURRENT_SOURCE_DIR}/liblsoda/src/liblsoda.so
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/liblsoda/src)

add_executable(ES ${source_cpp} ${source_c})
target_link_libraries(ES PUBLIC liblsoda)

# add_executable(test tests/test.cpp)
# target_link_libraries(test PUBLIC liblsoda)
