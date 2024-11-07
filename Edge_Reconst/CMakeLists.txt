#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
set ( control_sources
  
  #> Multiview Geometry Utilities
  util.hpp
  util.cpp
  PairEdgeHypo.hpp
  PairEdgeHypo.cpp
  getReprojectedEdgel.hpp
  getReprojectedEdgel.cpp
  getQuadrilateral.hpp
  getQuadrilateral.cpp
  getSupportedEdgels.hpp
  getSupportedEdgels.cpp
  getOrientationList.hpp
  getOrientationList.cpp

  file_reader.hpp
  file_reader.cpp
  edge_mapping.hpp
  edge_mapping.cpp
  iteration.hpp
  iteration.cpp

  EdgeSketch_Core.cpp
  EdgeSketch_Core.hpp

  #subpixel_point_set.hpp
  #lemsvpe_CH/vgl_polygon_CH.hpp
  #lemsvpe_CH/vgl_point_2d_CH.hpp
  #lemsvpe_CH/vgl_point_2d_CH.cpp
  #lemsvpe_CH/vgl_polygon_scan_iterator_CH.cpp
  #lemsvpe_CH/vgl_polygon_scan_iterator_CH.hpp
  #lemsvpe_CH/vgl_region_scan_iterator_CH.hpp
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library( edge_reconstruction ${control_sources} )
target_link_libraries(edge_reconstruction 
        #-L/gpfs/runtime/opt/magma-gpu/2.4.0/lib magma
        #-L${CUDAPATH}/lib cublas cudart cusparse
        #-L/gpfs/runtime/opt/openblas/0.3.7/lib openblas
        -L/gpfs/data/bkimia/cchien3/yaml-cpp/bin/lib64 yaml-cpp
        -L/usr/lib64 pthread
)

set_target_properties(edge_reconstruction PROPERTIES LINKER_LANGUAGE CXX)
#set_target_properties(edge_reconstruction PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
