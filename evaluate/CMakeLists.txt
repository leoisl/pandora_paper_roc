cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR pandora1_paper)
project(${PROJECT_NAME_STR} C CXX)

file(GLOB_RECURSE SRC_FILES
        ${PROJECT_SOURCE_DIR}/*.cpp
        ${PROJECT_SOURCE_DIR}/*.h)

add_executable(${PROJECT_NAME} ${SRC_FILES})

