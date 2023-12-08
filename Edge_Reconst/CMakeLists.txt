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

  subpixel_point_set.hpp

  lemsvpe_CH/vgl_polygon_CH.hpp
  lemsvpe_CH/vgl_point_2d_CH.hpp
  lemsvpe_CH/vgl_point_2d_CH.cpp
  lemsvpe_CH/vgl_polygon_scan_iterator_CH.cpp
  lemsvpe_CH/vgl_polygon_scan_iterator_CH.hpp
  lemsvpe_CH/vgl_region_scan_iterator_CH.hpp

  #> GPU
  GPU_Edge_Reconst/gpu_edge_reconst.hpp
  GPU_Edge_Reconst/gpu_kernels.hpp
  GPU_Edge_Reconst/indices.hpp
  GPU_Edge_Reconst/consistency_check.hpp
  GPU_Edge_Reconst/gpu_preprocess_edgels_HYPO1.cu
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_library( edge_reconstruction ${control_sources} )
target_link_libraries(edge_reconstruction 
        PUBLIC OpenMP::OpenMP_CXX
        #-L/gpfs/runtime/opt/magma-gpu/2.4.0/lib magma
        #-L${CUDAPATH}/lib cublas cudart cusparse
        #-L/gpfs/runtime/opt/openblas/0.3.7/lib openblas
        -L/usr/lib64 pthread
)

#set_target_properties(edge_reconstruction PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(edge_reconstruction PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
