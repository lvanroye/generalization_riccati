add_library(symspals SHARED symspals/symspals.cpp)
target_link_libraries(symspals MUMPS)
add_executable(expressions_example symspals/expressions_example.cpp)
add_executable(sparse_system_example symspals/sparse_system_example.cpp)
target_link_libraries(expressions_example symspals)
target_link_libraries(sparse_system_example symspals)
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