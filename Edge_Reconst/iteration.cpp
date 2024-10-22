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
    const int HYPO2_VIEW_INDX
) {
    std::cout << "Core pipeline start" << std::endl;

    PairEdgeHypothesis::pair_edge_hypothesis       PairHypo;
    GetOrientationList::get_OrientationList        getOre;

    Eigen::MatrixXd OreListdegree    = getOre.getOreList(HYPO1_VIEW_INDX, HYPO2_VIEW_INDX, Edges_HYPO2, All_R, All_T, K1, K2);
    Eigen::MatrixXd OreListBardegree = getOre.getOreListBar(Edges_HYPO1, All_R, All_T, K1, K2, HYPO2_VIEW_INDX, HYPO1_VIEW_INDX);
    
    unsigned nthreads = NUM_OF_OPENMP_THREADS;
    omp_set_num_threads(nthreads);
    double itime = omp_get_wtime();

    Eigen::MatrixXd Gamma1s;

    #pragma omp parallel
    {
        // Local thread-specific variables for supported indices
        std::vector<Eigen::MatrixXd> local_thread_supported_indices;

        // First loop: loop over all edgels from hypothesis view 1
        #pragma omp for schedule(static, nthreads)
        for (int edge_idx = 0; edge_idx < Edges_HYPO1.rows(); edge_idx++) {

            // Edge Boundary Check: not too close to boundary
            if (Edges_HYPO1(edge_idx, 0) < 10 || Edges_HYPO1(edge_idx, 0) > IMGCOLS - 10 || 
                Edges_HYPO1(edge_idx, 1) < 10 || Edges_HYPO1(edge_idx, 1) > IMGROWS - 10) {
                continue;
            }

            // Paired Edge Check: not yet been paired
            if (paired_edge(edge_idx, 0) != -2) {
                continue;
            }

            // Get the current edge from Hypo1
            Eigen::Vector3d pt_edgel_HYPO1;
            pt_edgel_HYPO1 << Edges_HYPO1(edge_idx, 0), Edges_HYPO1(edge_idx, 1), 1;

            // Get angle Thresholds from OreListBar (Degree) 
            double thresh_ore21_1 = OreListBardegree(edge_idx, 0);
            double thresh_ore21_2 = OreListBardegree(edge_idx, 1);

            // Find Edges in Hypo2 Based on Angle Thresholds
            Eigen::MatrixXd HYPO2_idx_raw = PairHypo.getHYPO2_idx_Ore(OreListdegree, thresh_ore21_1, thresh_ore21_2);
            if (HYPO2_idx_raw.rows() == 0) {
                continue;
            }

            // Retrieve Hypo2 Edgels
            Eigen::MatrixXd edgels_HYPO2 = PairHypo.getedgels_HYPO2_Ore(Edges_HYPO2, OreListdegree, thresh_ore21_1, thresh_ore21_2);

            // Correct Edgels in Hypo2 Based on Epipolar Constraints
            Eigen::MatrixXd edgel_HYPO1 = Edges_HYPO1.row(edge_idx);
            Eigen::MatrixXd edgels_HYPO2_corrected = PairHypo.edgelsHYPO2correct(edgels_HYPO2, edgel_HYPO1, F, F12, HYPO2_idx_raw);

            // Process paired edges
            if (HYPO2_idx_raw.rows() == 0 || edgels_HYPO2_corrected.rows() == 0) {
                continue;
            }

            // Store final edge data
            Eigen::MatrixXd Edges_HYPO1_final(edgels_HYPO2_corrected.rows(), 4);
            Edges_HYPO1_final << edgels_HYPO2_corrected.col(0), edgels_HYPO2_corrected.col(1), edgels_HYPO2_corrected.col(2), edgels_HYPO2_corrected.col(3);
            Eigen::MatrixXd Edges_HYPO2_final(edgels_HYPO2_corrected.rows(), 4);
            Edges_HYPO2_final << edgels_HYPO2_corrected.col(4), edgels_HYPO2_corrected.col(5), edgels_HYPO2_corrected.col(6), edgels_HYPO2_corrected.col(7);

            // Stack supported indices for validation views
            std::vector<std::pair<int, Eigen::Vector2d>> validation_support_edges;

            // Loop through validation views to find the supporting edges
            for (int val_idx = 0; val_idx < DATASET_NUM_OF_FRAMES; ++val_idx) {
                if (val_idx == HYPO1_VIEW_INDX || val_idx == HYPO2_VIEW_INDX) {
                    continue;  // Skip hypothesis views
                }

                // Retrieve support index from paired_edge for the current validation view
                int support_idx = paired_edge(edge_idx, val_idx + 2);  // +2 accounts for the first two columns for HYPO1 and HYPO2
                if (support_idx != -2) {
                    // Retrieve the supporting edge from the validation view
                    Eigen::MatrixXd edges_for_val_frame = All_Edgels[val_idx];
                    Eigen::Vector2d supporting_edge = edges_for_val_frame.row(support_idx).head<2>();

                    // Store validation view and the supporting edge
                    validation_support_edges.emplace_back(val_idx, supporting_edge);

                    // Add the supporting edge to the edgeMapping for the 3D edge
                    edgeMapping.add3DToSupportingEdgesMapping(edgel_HYPO1, supporting_edge, val_idx);
                }
            }

            // Calculate 3D edge using triangulation
            std::vector<Eigen::Vector2d> pts;
            pts.push_back(Edges_HYPO1_final.row(0).head<2>());
            pts.push_back(Edges_HYPO2_final.row(0).head<2>());

            std::vector<Eigen::Matrix3d> Rs = { R21 };
            std::vector<Eigen::Vector3d> Ts = { T21 };
            std::vector<double> K1_v = { K1(0, 2), K1(1, 2), K1(0, 0), K1(1, 1) };

            Eigen::Vector3d edge_pt_3D = linearTriangulation(2, pts, Rs, Ts, K1_v);
            if (!edge_pt_3D.hasNaN()) {
                Gamma1s.conservativeResize(Gamma1s.rows() + 1, 3);
                Gamma1s.row(Gamma1s.rows() - 1) = edge_pt_3D;

                // Add the 3D edge and its supporting edges to the mapping
                edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, pts[0], HYPO1_VIEW_INDX);
                edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, pts[1], HYPO2_VIEW_INDX);
            }
        }  // End of first loop

        // Timing for the parallel section
        double ftime = omp_get_wtime();
        double exec_time = ftime - itime;
        std::cout << "It took " << exec_time << " second(s) to finish this round." << std::endl;
    }  // End of OpenMP parallel section

    return Gamma1s;  // Return the Gamma1s matrix with the reconstructed 3D edges
}



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