cmake_minimum_required(VERSION 3.0.0)
project(symspals VERSION 0.1.0 LANGUAGES CXX)
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

add_library(symspals SHARED src/symspals/symspals.cpp src/symspals/expressions.hpp src/symspals/sparse_system.hpp)
add_executable(expressions_example src/symspals/expressions_example.cpp)
add_executable(sparse_system_example src/symspals/sparse_system_example.cpp)

include(CPack)

# make the symspals library installable as a cmake package
install(TARGETS symspals
    EXPORT symspalsConfig 
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)
# include the configuration file to the installation
INSTALL(EXPORT symspalsConfig DESTINATION cmake)
