#--------------------------------------------
#
#   CMAKE CONFIG FILE FOR ADAPT
#
#   Author: Roman Keller
#    Start: 20.01.2019
#     Last: 20.01.2019
#--------------------------------------------

add_library(unittest_environment unittest_environment.F90)
target_link_libraries(unittest_environment create_block)
target_link_libraries(unittest_environment common)
target_link_libraries(unittest_environment init)
target_link_libraries(unittest_environment refin)
target_link_libraries(unittest_environment choose)
target_include_directories (unittest_environment  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

set_target_properties (unittest_environment PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# add_subdirectory(3x3)
add_subdirectory(array_holes)
# add_subdirectory(refine)
