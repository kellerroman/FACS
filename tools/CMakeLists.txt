add_executable(sol2csv sol2csv.F90)
add_executable(sol2para sol2para.F90)
target_link_libraries(sol2csv common)
target_link_libraries(sol2para common)

set(CMAKE_Fortran_MODULE_DIRECTORY $(CMAKE_BINARY_DIR}/mod_files)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include))
#target_link_libraries(create_xdmf.exe "hdf5_fortran")

#add_library(module_gridgen module_gridgen.F90)
# this is required so the compiler finds the module when compiling a code using this module
#target_include_directories (sol2csv  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories (sol2para PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
#target_link_libraries(module_gridgen INTERFACE "hdf5_fortran")

#add_library(module_post module_post.F90)
