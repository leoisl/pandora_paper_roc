include(${PROJECT_SOURCE_DIR}/ext/gtest.cmake)

find_package(Boost 1.71.0 COMPONENTS serialization system REQUIRED)

file(GLOB_RECURSE TEST_SRC_FILES
        ${PROJECT_SOURCE_DIR}/evaluate/*.h
        ${PROJECT_SOURCE_DIR}/evaluate/*.cpp
        ${PROJECT_SOURCE_DIR}/tests/*.cpp)
add_executable(pandora1_paper_test ${TEST_SRC_FILES})

add_dependencies(pandora1_paper_test gtest)

target_link_libraries(pandora1_paper_test
        libgmock
        libgtest
        pthread
        ${Boost_LIBRARIES}
)

target_include_directories(pandora1_paper_test PUBLIC
        ${PROJECT_SOURCE_DIR}/evaluate
        ${Boost_INCLUDE_DIR}
        ${source_dir}/googletest/include #Change here for CLion unit tests inspection
        ${source_dir}/googlemock/include
        /usr/include/python3.6 #TODO: change this
        )

add_test(NAME pandora1_paper_test COMMAND pandora1_paper_test)
