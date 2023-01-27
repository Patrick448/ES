cmake_minimum_required(VERSION 3.20)
project(ES)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(ES
        main.cpp
        ESAlgorithm.cpp
        ESAlgorithm.h
        Individual.cpp Individual.h)
