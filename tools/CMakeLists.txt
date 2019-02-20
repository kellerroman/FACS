add_executable(sol2csv sol2csv.F90)
target_link_libraries(sol2csv common)

add_executable(sol2para sol2para.F90)
target_link_libraries(sol2para common)

add_library(create_block create_block.F90)
target_link_libraries(create_block common)
target_include_directories (create_block  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

#target_link_libraries(create_xdmf.exe "hdf5_fortran")

#add_library(module_gridgen module_gridgen.F90)
# this is required so the compiler finds the module when compiling a code using this module
#target_include_directories (sol2csv  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories (sol2para PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
#target_link_libraries(module_gridgen INTERFACE "hdf5_fortran")

#add_library(module_post module_post.F90)
