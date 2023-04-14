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

# function(add_external_library)
#     set(options)
#     set(oneValueArgs TARGET WORKING_DIRECTORY OUTPUT COMMENT)
#     set(multiValueArgs COMMAND INCLUDE_DIRS)
#     cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" ${multiValueArgs}" ${ARGN})

#     # Specify how to build the result of the library
#     add_custom_command(OUTPUT "${ARGS_OUTPUT}"
#         COMMAND ${ARGS_COMMAND}
#         WORKING_DIRECTORY "${ARGS_WORKING_DIRECTORY}"
#         COMMENT "${ARGS_COMMENT}")

#     # Create a target out of the library compilation result
#     add_custom_target(${ARGS_TARGET}_target DEPENDS ${ARGS_OUTPUT})

#     # Create an library target out of the library compilation result
#     add_library(${ARGS_TARGET} STATIC IMPORTED GLOBAL)
#     add_dependencies(${ARGS_TARGET} ${ARGS_TARGET}_target)

#     # Specify where the library is and where to find the headers
#     set_target_properties(${ARGS_TARGET}
#         PROPERTIES
#         IMPORTED_LOCATION "${ARGS_OUTPUT}"
#         INTERFACE_INCLUDE_DIRECTORIES "${ARGS_INCLUDE_DIRS}")
# endfunction()

# # Example
# add_external_library(TARGET YourLib
#     COMMAND           /bin/bash compile_your_lib.sh
#     WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
#     OUTPUT            "output/yourlib.a"
#     INCLUDE_DIRS      "include/a" "include/b"
#     COMMENT           "Building YourLib")

# add_executable(YourExe)
# target_link_libraries(YourExe YourLib)

# add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/liblsoda
#         COMMAND ${ARGS_COMMAND}
#         WORKING_DIRECTORY "${ARGS_WORKING_DIRECTORY}"
#         COMMENT "${ARGS_COMMENT}")

# add_custom_target(
#    liblsoda
#    COMMAND make
#    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/liblsoda
# )

add_executable(ES ${source_cpp} ${source_c})
add_dependencies(ES liblsoda)
target_link_libraries(ES INTERFACE ${CMAKE_SOURCE_DIR}/liblsoda/src/liblsoda.so)

