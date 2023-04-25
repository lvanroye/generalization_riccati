cmake_minimum_required(VERSION 3.0.0)
project(symspals VERSION 0.1.0 LANGUAGES CXX C Fortran)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp -O0 -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra -Wno-unused-parameter")



add_subdirectory(external)
add_subdirectory(src)

include(CPack)