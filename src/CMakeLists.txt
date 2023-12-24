add_executable(pfp_doc64 pfp_doc.cpp ref_builder.cpp)
target_link_libraries(pfp_doc64 common pfp gsacak64 sdsl "-fopenmp")
target_include_directories(pfp_doc64 PUBLIC "../include/")
target_compile_options(pfp_doc64 PUBLIC "-std=c++17" "-fopenmp" "-DM64" "-march=native" "-g")




