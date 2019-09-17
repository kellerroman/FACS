add_pfunit_test (refinement_simple ut_refi.pf "" "")
target_link_libraries(refinement_simple unittest_environment)

add_pfunit_test (refinement_complex ut_refic.pf "" "")
target_link_libraries(refinement_complex unittest_environment)
# target_include_directories(my_tests PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

