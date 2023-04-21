add_subdirectory(fatrop)
include_directories("fatrop/fatrop")
set(metis ON CACHE BOOL "metis ordering library" FORCE)
add_subdirectory(mumps)