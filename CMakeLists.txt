cmake_minimum_required(VERSION 3.16.3)
project(ES)

set(CMAKE_CXX_STANDARD 14)
#SET(CMAKE_CXX_FLAGS "-O3")
#SET(CMAKE_C_FLAGS "-O3")

file(GLOB source_cpp src/*.cpp)
file(GLOB source_h include/*.h)
file(GLOB source_c src/*.c)


include_directories(include)
include_directories(liblsoda/src)

add_custom_target(
   liblsoda
   COMMAND make
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/liblsoda
)

add_executable(ES ${source_cpp} ${source_c})
add_dependencies(ES liblsoda)
target_link_libraries(ES PUBLIC ${CMAKE_SOURCE_DIR}/liblsoda/src/liblsoda.so)

