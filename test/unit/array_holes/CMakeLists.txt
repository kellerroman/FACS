set(TEST_NAME "array_holes")
add_pfunit_test (${TEST_NAME} ut_array_holes.pf "" "")
# add_test(${TEST_NAME} ./${TEST_NAME})
# add_executable(${TEST_NAME} test.F90)
target_link_libraries(${TEST_NAME} common)

