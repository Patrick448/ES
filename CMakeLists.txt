cmake_minimum_required(VERSION 3.16.3)
project(ES)

set(CMAKE_CXX_STANDARD 14)
#SET(CMAKE_CXX_FLAGS "-O3")
#SET(CMAKE_C_FLAGS "-O3")

include_directories(.)

add_library(liblsoda STATIC IMPORTED)
set_property(TARGET liblsoda PROPERTY
        IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/liblsoda.so")
#[[set_property(TARGET liblsoda PROPERTY
        IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/liblsoda.dll")
set_property(TARGET liblsoda PROPERTY
        IMPORTED_IMPLIB "${CMAKE_SOURCE_DIR}/liblsoda.lib")]]


add_executable(ES
        lsoda.h
        dependencies.h
        ESAlgorithm.cpp
        ESAlgorithm.h
        grn_eval.c
        grn_eval.h
        Individual.cpp
        Individual.h
        main.cpp
        ode.c
        ode.h
        read.cpp
        read.h
        rk4.c
        rk4.h

        )

target_link_libraries(ES PRIVATE liblsoda)
#target_link_libraries(ES "${CMAKE_SOURCE_DIR}/liblsoda.lib")