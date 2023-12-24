add_executable(pfp_mum64 pfp_mum.cpp ref_builder.cpp)
target_link_libraries(pfp_mum64 common pfp gsacak64 sdsl "-fopenmp")
target_include_directories(pfp_mum64 PUBLIC "../include/")
target_compile_options(pfp_mum64 PUBLIC "-std=c++17" "-fopenmp" "-DM64" "-march=native" "-g")




