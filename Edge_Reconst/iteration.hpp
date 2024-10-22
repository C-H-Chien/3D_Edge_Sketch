#ifndef ITERATION_HPP
#define ITERATION_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>

// Forward declarations of any classes or structs used
class EdgeMapping;

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
