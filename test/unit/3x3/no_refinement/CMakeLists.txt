set(TEST_NAME "unittest_3x3_unref")
add_executable(${TEST_NAME} test.F90)
target_link_libraries(${TEST_NAME} create_block)
target_link_libraries(${TEST_NAME} common)
target_link_libraries(${TEST_NAME} init)
target_link_libraries(${TEST_NAME} refin)
target_link_libraries(${TEST_NAME} choose)

