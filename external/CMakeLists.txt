add_subdirectory(fatrop)
include_directories("fatrop/fatrop")
set(metis off CACHE BOOL "metis ordering library" FORCE)
set(parallel off CACHE BOOL "parallel" FORCE)
# add_library(ma57lib SHARED ma57_src/ddeps.f ma57_src/ma57d.f ma57_src/fakemetis.f)
add_library(ma57lib STATIC ma57_src/ddeps.f ma57_src/ma57d.f)
add_subdirectory(mumps)
add_subdirectory(metis-4.0.3)
link_directories(pardiso_shared_lib/)
target_link_libraries(ma57lib PRIVATE metis)