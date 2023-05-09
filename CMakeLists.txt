cmake_minimum_required(VERSION 3.0.0)
project(symspals VERSION 0.1.0 LANGUAGES CXX C Fortran)
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Profile")
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp -O0 -DPANEL_SIZE=${PANEL_SIZE} -Wall -Wextra -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
set(CMAKE_C_FLAGS_RELEASE "-O3 -march=native")

set(CMAKE_CXX_FLAGS_PROFILE "-O1 -g -pg -fno-omit-frame-pointer -fno-inline-functions")
set(CMAKE_ASM_FLAGS_PROFILE "-O1 -g -pg -fno-omit-frame-pointer -fno-inline-functions")
set(CMAKE_C_FLAGS_PROFILE "-O1 -g -pg -fno-omit-frame-pointer -fno-inline-functions")
set(CMAKE_FORTRAN_FLAGS_PROFILE "-O1 -pg -g -fno-omit-frame-pointer -fno-inline-functions")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(external)
add_subdirectory(src)
