#ifndef ITERATION_HPP
#define ITERATION_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <assert.h>
#include <string>
#include <ctime>
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <Eigen/Dense>
#include <algorithm> 
#include <utility>  
#include "../Edge_Reconst/util.hpp"
#include "../Edge_Reconst/PairEdgeHypo.hpp"
#include "../Edge_Reconst/getReprojectedEdgel.hpp"
#include "../Edge_Reconst/getQuadrilateral.hpp"
#include "../Edge_Reconst/getSupportedEdgels.hpp"
#include "../Edge_Reconst/getOrientationList.hpp"
#include "../Edge_Reconst/linearTriangulationUtil.hpp"
#include "../Edge_Reconst/definitions.h"
#include "../Edge_Reconst/lemsvpe_CH/vgl_polygon_CH.hpp"
#include "../Edge_Reconst/lemsvpe_CH/vgl_polygon_scan_iterator_CH.hpp"
#include "../Edge_Reconst/subpixel_point_set.hpp"
#include "../Edge_Reconst/file_reader.hpp"
#include "../Edge_Reconst/edge_mapping.hpp"
#include "../Edge_Reconst/iteration.hpp"

// Forward declarations of any classes or structs used
class EdgeMapping;

Eigen::MatrixXd project3DEdgesToView(
    const Eigen::MatrixXd& edges3D, 
    const Eigen::Matrix3d& R, 
    const Eigen::Vector3d& T, 
    const Eigen::Matrix3d& K,
    const Eigen::Matrix3d& R_hyp01, 
    const Eigen::Vector3d& T_hpy01
);

std::vector<int> findClosestObservedEdges(
    const Eigen::MatrixXd& projectedEdges, 
    const Eigen::MatrixXd& observedEdges, 
    double threshold
);

std::pair<int, int> selectBestViews(
    const std::vector<std::vector<int>>& claimedEdges
);


// CorePipelineOutput struct for storing the output of the core_pipeline function
struct CorePipelineOutput {
    Eigen::MatrixXd Gamma1s; // 3D edge points

    // Mapping structure: frame index -> 2D edge point -> list of corresponding 3D points
    std::map<int, std::unordered_map<Eigen::Vector2d, std::vector<Eigen::Vector3d>, HashEigenVector2d>> frame_to_edge_to_3D_map;
};

#endif  // ITERATION_HPP
