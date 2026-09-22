cmake_minimum_required(VERSION 3.16.3)
project("Chaos in spin models" )
#
#add_subdirectory(src)
#add_subdirectory(tests)
#
enable_language(Fortran)
find_package(LAPACK REQUIRED)
#
#enable_testing()
#
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#
file(GLOB_RECURSE lib_sources CONFIGURE_DEPENDS  src/*.f90)
list(REMOVE_ITEM lib_sources src/main.f90)
add_library(spin_lib ${lib_sources})
target_link_libraries(spin_lib PRIVATE ${LAPACK_LIBRARIES})
#
add_executable(spin src/main.f90)
#
target_link_libraries(spin PRIVATE spin_lib)
#
#add_subdirectory(tests)
