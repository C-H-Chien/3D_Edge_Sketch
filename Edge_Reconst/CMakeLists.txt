#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
set ( control_sources
  
  #> Multiview Geometry Utilities
  util.hpp
  util.cpp
  PairEdgeHypo.hpp
  PairEdgeHypo.cpp
  getReprojectedEdgel.hpp
  getReprojectedEdgel.cpp
  getSupportedEdgels.hpp
  getSupportedEdgels.cpp
  getOrientationList.hpp
  getOrientationList.cpp

  file_reader.hpp
  file_reader.cpp
  edge_mapping.hpp
  edge_mapping.cpp

  EdgeSketch_Core.cpp
  EdgeSketch_Core.hpp
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library( edge_reconstruction ${control_sources} )
target_link_libraries(edge_reconstruction 
        -L/gpfs/data/bkimia/cchien3/yaml-cpp/bin/lib64 yaml-cpp
        -L/usr/lib64 pthread
)

set_target_properties(edge_reconstruction PROPERTIES LINKER_LANGUAGE CXX)
#set_target_properties(edge_reconstruction PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
