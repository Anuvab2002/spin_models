cmake_minimum_required(VERSION 3.16.3)
project("Chaos in spin models" )
#
#add_subdirectory(src)
#add_subdirectory(tests)
#
enable_language(Fortran)
#
#enable_testing()
#
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(GLOB_RECURSE sources  src/*.f90 src/*.h)
#
add_executable(spin ${sources})
#
#file(GLOB_RECURSE lib_sources  src/*.f90)
#list(REMOVE_ITEM lib_sources src/main.f90)
#add_library(ke_lib ${lib_sources})
#
#add_subdirectory(tests)
