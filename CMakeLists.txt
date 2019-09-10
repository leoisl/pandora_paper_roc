cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR pandora1_paper)
project(${PROJECT_NAME_STR} C CXX)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fPIC)

#add_executable(pangenome_snps evaluate/Pangenome_SNPs.h)
add_subdirectory(tests)
