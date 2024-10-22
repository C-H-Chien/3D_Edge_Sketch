#include "iteration.hpp"

using namespace MultiviewGeometryUtil;


Eigen::MatrixXd project3DEdgesToView(const Eigen::MatrixXd& edges3D, const Eigen::Matrix3d& R, const Eigen::Vector3d& T, const Eigen::Matrix3d& K) {

    Eigen::MatrixXd edges2D(edges3D.rows(), 2);

    for (int i = 0; i < edges3D.rows(); ++i) {
        Eigen::Vector3d point3D = edges3D.row(i).transpose();
        Eigen::Vector3d point_camera = R * point3D + T;

        // Check if the Z value is zero to avoid division by zero
        if (point_camera(2) == 0) {
            std::cout << "Error: Point " << i << " is located at infinity (Z=0) after camera transformation.\n"<<std::endl;
            continue;  
        }
        
        Eigen::Vector3d point_image = K * point_camera;
        edges2D(i, 0) = point_image(0) / point_image(2);
        edges2D(i, 1) = point_image(1) / point_image(2);
    }
    
    return edges2D;
}




std::vector<int> findClosestObservedEdges(const Eigen::MatrixXd& projectedEdges, const Eigen::MatrixXd& observedEdges, double threshold) {
    std::vector<int> claimedEdges;
    for (int i = 0; i < projectedEdges.rows(); ++i) {
        double minDist = std::numeric_limits<double>::max();
        int closestEdgeIdx = -1;
        for (int j = 0; j < observedEdges.rows(); ++j) {
            double dist = (projectedEdges.row(i) - observedEdges.row(j).head<2>()).norm();
            if (dist < minDist && dist < threshold) {
                minDist = dist;
                closestEdgeIdx = j;
            }
        }
        if (closestEdgeIdx != -1) {
            claimedEdges.push_back(closestEdgeIdx);
        }
    }
    return claimedEdges;
}



std::pair<int, int> selectBestViews(const std::vector<std::vector<int>>& claimedEdges, const std::vector<Eigen::Vector3d>& cameraPositions, double baselineThreshold) {
    int bestView1 = -1, bestView2 = -1;
    int minClaimedEdges = std::numeric_limits<int>::max();

    for (int i = 0; i < claimedEdges.size(); ++i) {
        for (int j = i + 1; j < claimedEdges.size(); ++j) {
            double baseline = (cameraPositions[i] - cameraPositions[j]).norm();
            if (baseline > baselineThreshold) {
                int totalClaimedEdges = claimedEdges[i].size() + claimedEdges[j].size();
                if (totalClaimedEdges < minClaimedEdges) {
                    minClaimedEdges = totalClaimedEdges;
                    bestView1 = i;
                    bestView2 = j;
                }
            }
        }
    }
    return {bestView1, bestView2};
}