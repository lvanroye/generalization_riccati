add_library(symspals SHARED symspals/symspals.cpp symspals/timing.c)
target_link_directories(symspals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../external/pardiso_shared_lib)
target_link_libraries(symspals MUMPS ma57lib pardiso blas m gomp)
target_include_directories(symspals PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/symspals)
add_executable(expressions_example symspals/expressions_example.cpp)
add_executable(sparse_system_example symspals/sparse_system_example.cpp)
target_link_libraries(expressions_example symspals)
target_link_libraries(sparse_system_example symspals)

add_library(gen_riccati SHARED gen_riccati/gen_riccati.cpp gen_riccati/blasfeo_wrapper.cpp)
target_include_directories(gen_riccati PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/gen_riccati)
target_link_libraries(gen_riccati blasfeo)




add_executable(benchmark_paper benchmark_auxilary/benchmark_paper.cpp)
target_link_libraries(benchmark_paper fatrop symspals blasfeo gen_riccati)

add_executable(benchmark_quadrotor benchmark_auxilary/benchmark_quadrotor.cpp)
target_link_libraries(benchmark_quadrotor fatrop symspals blasfeo gen_riccati)

add_executable(benchmark_random benchmark_auxilary/benchmark_random.cpp)
target_link_libraries(benchmark_random fatrop symspals blasfeo gen_riccati)