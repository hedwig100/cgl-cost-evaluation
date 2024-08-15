enable_testing()
include(GoogleTest)
include_directories(${PROJECT_SOURCE_DIR}/lib)

## discrete_log_test
add_executable(discrete_log_test
    discrete_log_test.cpp
)
target_link_libraries(discrete_log_test
    GTest::gtest_main
)
gtest_discover_tests(discrete_log_test)

## elliptic_test
add_executable(elliptic_test
    elliptic_test.cpp
)
target_link_libraries(elliptic_test
    GTest::gtest_main
)
gtest_discover_tests(elliptic_test)

## entangled_basis_test
add_executable(entangled_basis_test
    entangled_basis_test.cpp
)
target_link_libraries(entangled_basis_test
    GTest::gtest_main
)
gtest_discover_tests(entangled_basis_test)

## field_test
add_executable(field_test
    field_test.cpp
)
target_link_libraries(field_test
    GTest::gtest_main
)
gtest_discover_tests(field_test)

## mod_test
add_executable(mod_test
    mod_test.cpp
)
target_link_libraries(mod_test
    GTest::gtest_main
)
gtest_discover_tests(mod_test)

## montgomery_test
add_executable(montgomery_test
    montgomery_test.cpp
)
target_link_libraries(montgomery_test
    GTest::gtest_main
)
gtest_discover_tests(montgomery_test)

## prime_test
add_executable(prime_test
    prime_test.cpp
)
target_link_libraries(prime_test
    GTest::gtest_main
)
gtest_discover_tests(prime_test)

## quadrt_test
add_executable(quadrt_test
    quadrt_test.cpp
)
target_link_libraries(quadrt_test
    GTest::gtest_main
)
gtest_discover_tests(quadrt_test)

## random_field_test
add_executable(random_field_test
    random_field_test.cpp
)
target_link_libraries(random_field_test
    GTest::gtest_main
)
gtest_discover_tests(random_field_test)

## random_test
add_executable(random_test
    random_test.cpp
)
target_link_libraries(random_test
    GTest::gtest_main
)
gtest_discover_tests(random_test)

## sqrt_test
add_executable(sqrt_test
    sqrt_test.cpp
)
target_link_libraries(sqrt_test
    GTest::gtest_main
)
gtest_discover_tests(sqrt_test)

## utils_test
add_executable(utils_test
    utils_test.cpp
)
target_link_libraries(utils_test
    GTest::gtest_main
)
gtest_discover_tests(utils_test)