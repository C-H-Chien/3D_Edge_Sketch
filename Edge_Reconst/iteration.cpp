#include "iteration.hpp"
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

using namespace MultiviewGeometryUtil;



Eigen::MatrixXd project3DEdgesToView(const Eigen::MatrixXd& edges3D, const Eigen::Matrix3d& R, const Eigen::Vector3d& T, const Eigen::Matrix3d& K) {
    Eigen::MatrixXd edges2D(edges3D.rows(), 2);
    for (int i = 0; i < edges3D.rows(); ++i) {
        Eigen::Vector3d point3D = edges3D.row(i);
        Eigen::Vector3d point_camera = R * point3D + T;
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