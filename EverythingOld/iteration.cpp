#include "iteration.hpp"

using namespace MultiviewGeometryUtil;


Eigen::MatrixXd project3DEdgesToView(const Eigen::MatrixXd& edges3D, const Eigen::Matrix3d& R, const Eigen::Vector3d& T, const Eigen::Matrix3d& K, const Eigen::Matrix3d& R_hyp01, const Eigen::Vector3d& T_hpy01) {

    Eigen::MatrixXd edges2D(edges3D.rows(), 2);

    for (int i = 0; i < edges3D.rows(); ++i) {
        Eigen::Vector3d point3D = edges3D.row(i).transpose();
        //Eigen::Vector3d world_point3D = R_hyp01.transpose() * (point3D - T_hpy01);
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

int claim_Projected_Edges(const Eigen::MatrixXd& projectedEdges, const Eigen::MatrixXd& observedEdges, double threshold) {
    
    int num_of_claimed_observed_edges = 0;

    //> Loop over all observed edges
    for (int i = 0; i < observedEdges.rows(); ++i) {

        //> Loop over all projected edges
        for (int j = 0; j < projectedEdges.rows(); ++j) {

            //> Calculate the Euclidean distance
            double dist = (projectedEdges.row(j) - observedEdges.row(i).head<2>()).norm();

            //> If the projected edge and the observed edge has Euclidean distance less than the "threshold"
            if (dist < threshold) {
                num_of_claimed_observed_edges++;
                break;
            }
        }
    }

    return num_of_claimed_observed_edges;
}

void select_Next_Best_Hypothesis_Views( 
    const std::vector< int >& claimedEdges, std::vector<Eigen::MatrixXd> All_Edgels, \
    std::pair<int, int> &next_hypothesis_views, double &least_ratio ) 
{
    std::vector<std::pair<int, double>> frameSupportCounts;

    double ratio_claimed_over_unclaimed;
    for (int i = 0; i < claimedEdges.size(); i++) {
        ratio_claimed_over_unclaimed = (double)(claimedEdges[i]) / (double)(All_Edgels[i].rows());
        frameSupportCounts.push_back({i, ratio_claimed_over_unclaimed});
    }

    std::sort(frameSupportCounts.begin(), frameSupportCounts.end(), 
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                  return a.second < b.second;
              });

    int bestView1 = frameSupportCounts[0].first;
    int bestView2 = frameSupportCounts[1].first;
    next_hypothesis_views = std::make_pair(bestView1, bestView2);
    least_ratio = frameSupportCounts[0].second;

    //> log: show the selected hypothesis views
    std::string out_str = "Selected frames with the least supported edges: " + std::to_string(bestView1) + " and " + std::to_string(bestView2);
    LOG_INFOR_MESG(out_str);

    // return {bestView1, bestView2};
}

std::pair<int, int> selectBestViews(const std::vector<std::vector<int>>& claimedEdges, int num_of_claimed_edges) {
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

    //> log: show the selected hypothesis views
    std::string out_str = "Selected frames with the least supported edges: " + std::to_string(bestView1) + " and " + std::to_string(bestView2);
    LOG_INFOR_MESG(out_str);

    return {bestView1, bestView2};
}
