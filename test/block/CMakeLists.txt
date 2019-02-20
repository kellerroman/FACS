add_executable(gridgen_block gridgen.F90)
target_link_libraries(gridgen_block create_block)
target_link_libraries(gridgen_block common)

