#include "iteration.hpp"

using namespace MultiviewGeometryUtil;


Eigen::MatrixXd project3DEdgesToView(const Eigen::MatrixXd& edges3D, const Eigen::Matrix3d& R, const Eigen::Vector3d& T, const Eigen::Matrix3d& K, const Eigen::Matrix3d& R_hyp01, const Eigen::Vector3d& T_hpy01) {

    Eigen::MatrixXd edges2D(edges3D.rows(), 2);

    for (int i = 0; i < edges3D.rows(); ++i) {
        Eigen::Vector3d point3D = edges3D.row(i).transpose();
        Eigen::Vector3d world_point3D = R_hyp01.transpose() * (point3D - T_hpy01);
        Eigen::Vector3d point_camera = R * world_point3D + T;

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



std::pair<int, int> selectBestViews(const std::vector<std::vector<int>>& claimedEdges) {
    std::vector<std::pair<int, int>> frameSupportCounts;

    for (int i = 0; i < claimedEdges.size(); i++) {
        int numSupportedEdges = claimedEdges[i].size();
        frameSupportCounts.push_back({i, numSupportedEdges});
    }

    std::sort(frameSupportCounts.begin(), frameSupportCounts.end(), 
              [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                  return a.second < b.second;
              });

    int bestView1 = frameSupportCounts[0].first;
    int bestView2 = frameSupportCounts[1].first;

    std::cout << "Selected frames with the least supported edges: " << bestView1 << " and " << bestView2 << std::endl;

    return {bestView1, bestView2};
}
