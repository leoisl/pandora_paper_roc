cmake_minimum_required(VERSION 2.8.7)

set(PROJECT_NAME_STR pandora1_paper_cpp)
project(${PROJECT_NAME_STR} C CXX)

#add_executable(pangenome_snps evaluate/Pangenome_SNPs.h)
add_subdirectory(tests)
