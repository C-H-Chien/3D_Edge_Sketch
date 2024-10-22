#ifndef ITERATION_HPP
#define ITERATION_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// Forward declarations of any classes or structs used
class EdgeMapping;

// Function Prototypes
Eigen::MatrixXd core_pipeline(
    const Eigen::MatrixXd& Edges_HYPO1, 
    const Eigen::MatrixXd& Edges_HYPO2, 
    const Eigen::Matrix3d& R21, 
    const Eigen::Vector3d& T21, 
    const Eigen::Matrix3d& F, 
    const Eigen::Matrix3d& F12,
    const std::vector<Eigen::Matrix3d>& All_R, 
    const std::vector<Eigen::Vector3d>& All_T,
    const std::vector<Eigen::MatrixXd>& All_Edgels,
    EdgeMapping& edgeMapping, 
    Eigen::MatrixXd& paired_edge,
    const Eigen::Matrix3d& K1,
    const Eigen::Matrix3d& K2,
    const int HYPO1_VIEW_INDX,
    const int HYPO2_VIEW_INDX,
    const double parallelangle,
    const double DIST_THRESH,
    const int MAX_NUM_OF_SUPPORT_VIEWS,
    const int NUM_OF_OPENMP_THREADS
);

Eigen::MatrixXd project3DEdgesToView(
    const Eigen::MatrixXd& edges3D, 
    const Eigen::Matrix3d& R, 
    const Eigen::Vector3d& T, 
    const Eigen::Matrix3d& K
);

std::vector<int> findClosestObservedEdges(
    const Eigen::MatrixXd& projectedEdges, 
    const Eigen::MatrixXd& observedEdges, 
    double threshold
);

std::pair<int, int> selectBestViews(
    const std::vector<std::vector<int>>& claimedEdges, 
    const std::vector<Eigen::Vector3d>& cameraPositions, 
    double baselineThreshold
);

#endif  // ITERATION_HPP
