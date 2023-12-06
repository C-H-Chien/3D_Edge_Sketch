
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#> main lib
find_package(Threads REQUIRED)
link_libraries(edge_reconstruction)

add_executable(edge_reconstruction-main main.cpp )

set_target_properties(edge_reconstruction-main PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(edge_reconstruction-main Threads::Threads)

set_property(TARGET edge_reconstruction-main PROPERTY CUDA_SEPARABLE_COMPILATION ON)

