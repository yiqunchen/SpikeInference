cmake_minimum_required(VERSION 3.14)
project(src)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(src
        funPieceList.cpp
        funPieceList.h
        fpop.cpp
        fpop.h
        selective_inference.cpp
        selective_inference.h
        fpop_inference.h
        fpop_inference.cpp
        utils.h
        utils.cpp
        main.cpp
        test_utils.h)
