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
#include "../Edge_Reconst/definitions.h"
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

#endif  // ITERATION_HPP
