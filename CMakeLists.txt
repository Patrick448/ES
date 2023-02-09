cmake_minimum_required(VERSION 3.20)
project(ES)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(ES
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
        rk4.h)
