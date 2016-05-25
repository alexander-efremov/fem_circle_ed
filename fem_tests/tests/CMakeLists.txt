include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runBasicFemTests
        tests.cpp)

target_link_libraries(runBasicFemTests gtest gtest_main)
target_link_libraries(runBasicFemTests fem)
