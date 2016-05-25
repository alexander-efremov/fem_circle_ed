cmake_minimum_required(VERSION 2.8)

project(fem)

add_definitions(-std=c++0x)

set(HEADER_FILES
    utils.h
    timer.h
    common.h
    consts.h
    )

set(SOURCE_FILES
       solver3.cpp solver4.cpp solver5.cpp consts.cpp)

add_library(fem STATIC ${SOURCE_FILES} ${HEADER_FILES})
