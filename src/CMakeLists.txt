add_executable(mumemto pfp_mum.cpp ref_builder.cpp)
target_link_libraries(mumemto common pfp gsacak64 sdsl "-fopenmp")# ${PROFILER_LIB})
target_include_directories(mumemto PUBLIC "../include/")
target_compile_options(mumemto PUBLIC "-std=c++17" "-fopenmp" "-DM64" "-march=native")

#add_executable(viz viz.cpp)
#target_link_libraries(viz common pfp gsacak64 sdsl "-fopenmp" matplot)# ${PROFILER_LIB})
#target_include_directories(viz PUBLIC "../include/")
#target_compile_options(viz PUBLIC "-std=c++17" "-fopenmp" "-DM64" "-march=native")