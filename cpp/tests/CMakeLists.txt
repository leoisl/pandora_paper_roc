include(${PROJECT_SOURCE_DIR}/ext/gtest.cmake)

find_package(Boost 1.71.0 COMPONENTS serialization system REQUIRED)

file(GLOB_RECURSE TEST_SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.h
        ${PROJECT_SOURCE_DIR}/src/*.cpp
        ${PROJECT_SOURCE_DIR}/tests/*.cpp)
message(INFO ${TEST_SRC_FILES})
add_executable(pandora1_paper_cpp_test ${TEST_SRC_FILES})

add_dependencies(pandora1_paper_cpp_test gtest)

target_link_libraries(pandora1_paper_cpp_test
        libgmock
        libgtest
        pthread
        ${Boost_LIBRARIES}
)

target_include_directories(pandora1_paper_cpp_test PUBLIC
        ${PROJECT_SOURCE_DIR}/src
        ${Boost_INCLUDE_DIR}
        ${source_dir}/googletest/include #Change here for CLion unit tests inspection
        ${source_dir}/googlemock/include
        )

add_test(NAME pandora1_paper_cpp_test COMMAND pandora1_paper_cpp_test)
