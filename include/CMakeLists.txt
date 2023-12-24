set(COMMON_SOURCES common.hpp)
add_library(common OBJECT ${COMMON_SOURCES})
target_link_libraries(common)

set(PFP_SOURCES dictionary.hpp parse.hpp pfp.hpp)
add_library(pfp OBJECT ${PFP_SOURCES})
target_link_libraries(pfp common sdsl)# divsufsort divsufsort64)