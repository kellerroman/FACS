#--------------------------------------------
#
#   CMAKE CONFIG FILE FOR ADAPT
#
#   Author: Roman Keller
#    Start: 20.01.2019
#     Last: 20.01.2019
#--------------------------------------------
cmake_minimum_required(VERSION 3.7)
project (ADAPT)
enable_language (Fortran)

#---------------------
# Set up the compiler
#---------------------
set (CMAKE_Fortran_COMPILER gfortran)

# set the preporcessor directives for compilation
set (PREPROC_DIRECTIVES "-DVERBOSE")
# set Compiler Directives to use with every BUILD_TYPE
 set (GLOBAL_FLAGS "-std=f2008 -ffree-form -Wall -Wextra -Wconversion -Warray-temporaries -Wimplicit-interface -Wsurprising")
# set the compiler flags -ffree-line-length-130
set (CMAKE_Fortran_FLAGS_RELEASE "${PREPROC_DIRECTIVES} ${GLOBAL_FLAGS} -O3")
set (CMAKE_Fortran_FLAGS_DEBUG   "${PREPROC_DIRECTIVES} ${GLOBAL_FLAGS} -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=all -fbounds-check -g -O0")
#set (CMAKE_BUILD_TYPE Debug)
set (CMAKE_BUILD_TYPE Release)

# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND
  "${GIT_EXECUTABLE}" describe --match=NeVeRmAtCh --always --abbrev=40 --dirty --tags
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_SHA1
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the date of the commit
execute_process(COMMAND
  "${GIT_EXECUTABLE}" log -1 --format=%ad --date=local
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_DATE
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the subject of the commit
execute_process(COMMAND
  "${GIT_EXECUTABLE}" log -1 --format=%s
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_COMMIT_SUBJECT
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

#------------------------
# Define the source path
#------------------------
file (GLOB sources src/*.F90)
file (GLOB src_common src/common/*.F90)
file (GLOB src_init src/init/*.F90)
file (GLOB src_ref src/refinement/*.F90)

#--------------------------
# Include HDF5 Library
#--------------------------
set(HDF5_ROOT $ENV{HDF5_HOME})
find_package(HDF5 REQUIRED)

include_directories(${HDF5_INCLUDE_DIRS})

#set(PRE_CONFIGURE_FILE "src/git_version.h.in")
#set(POST_CONFIGURE_FILE "git_version.h")
set(PRE_CONFIGURE_FILE "src/git_version.F90.in")
set(POST_CONFIGURE_FILE "${CMAKE_BINARY_DIR}/git_version.F90")
include(git_watcher.cmake)

#-------------------------
# Define the build & link
#-------------------------

add_library(common ${src_common})
target_include_directories (common  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

add_library(init ${src_init})
target_include_directories (init  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(init common)

add_library(refin ${src_ref})
target_include_directories (refin  PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(refin common)

add_executable (adapt.exe ${POST_CONFIGURE_FILE} ${sources})
#target_link_libraries(heat.exe "${HDF5}" "${HDF5_Fortran_LIBRARIES}")
target_link_libraries(adapt.exe "hdf5" "hdf5_fortran")
target_link_libraries(adapt.exe common)
target_link_libraries(adapt.exe init)
target_link_libraries(adapt.exe refin)

#include_directories(${CMAKE_BINARY_DIR})

add_subdirectory(tools)
add_subdirectory(test)
