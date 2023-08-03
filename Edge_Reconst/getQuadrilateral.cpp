#ifndef GETQUADRILATERAL_CPP
#define GETQUADRILATERAL_CPP
// ====================================================================================================
//
// =====================================================================================================
#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>
#include <random>
#include <algorithm>
#include <iomanip>
#include <string.h>
#include <assert.h>
#include <vector>
#include <chrono>

#include "getQuadrilateral.hpp"
#include "util.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace GetQuadrilateral {
    
    get_Quadrilateral::get_Quadrilateral( ) { }
    
    Eigen::MatrixXd get_Quadrilateral::getQuadrilateralPoints(Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd pt_edge_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e2  = {0,1,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[HYPO1_VIEW_INDX];
        Eigen::Vector3d T1  = All_T[HYPO1_VIEW_INDX];
        Eigen::Matrix3d R2  = All_R[HYPO2_VIEW_INDX];
        Eigen::Vector3d T2  = All_T[HYPO2_VIEW_INDX];
        Eigen::Matrix3d R3  = All_R[VALID_INDX];
        Eigen::Vector3d T3  = All_T[VALID_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d R32 = util.getRelativePose_R21(R2, R3);
        Eigen::Vector3d T32 = util.getRelativePose_T21(R2, R3, T2, T3);
        Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
        Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);

        Eigen::Vector3d epipole_met_view1 = {double(e1.transpose() * R21.transpose() *T21) / double(e3.transpose()*R21.transpose()*T21), double(e2.transpose()*R21.transpose()*T21) / double(e3.transpose()*R21.transpose()*T21), 1};
        Eigen::Vector3d epipole_met_view2 = {double(e1.transpose()*T21) / double(e3.transpose()*T21), double(e2.transpose()*T21) / double(e3.transpose()*T21), 1};
        Eigen::Vector3d epipole_pix_view1 = K * epipole_met_view1;
        Eigen::Vector3d epipole_pix_view2 = K * epipole_met_view2;

        double slope_hypo1        = double(pt_edge_HYPO1(1)-epipole_pix_view1(1))/double(pt_edge_HYPO1(0)-epipole_pix_view1(0));
        double intersect_hypo1    = double(pt_edge_HYPO1(1)) - slope_hypo1*double(pt_edge_HYPO1(0));
        Eigen::Vector3d hypo1_pt1 = {double(pt_edge_HYPO1(0)), intersect_hypo1+slope_hypo1*double(pt_edge_HYPO1(0))+DIST_THRESH, 1};
        Eigen::Vector3d hypo1_pt2 = {double(pt_edge_HYPO1(0)), intersect_hypo1+slope_hypo1*double(pt_edge_HYPO1(0))-DIST_THRESH, 1};
        double slope_hypo2        = double(pt_edge_HYPO2(1)-epipole_pix_view2(1))/double(pt_edge_HYPO2(0)-epipole_pix_view2(0));
        double intersect_hypo2    = double(pt_edge_HYPO2(1)) - slope_hypo2*double(pt_edge_HYPO2(0));
        Eigen::Vector3d hypo2_pt1 = {double(pt_edge_HYPO2(0)), intersect_hypo2+slope_hypo2*double(pt_edge_HYPO2(0))+DIST_THRESH, 1};
        Eigen::Vector3d hypo2_pt2 = {double(pt_edge_HYPO2(0)), intersect_hypo2+slope_hypo2*double(pt_edge_HYPO2(0))-DIST_THRESH, 1};

        Eigen::Matrix3d F31   = util.getFundamentalMatrix(K.inverse(), R31, T31);
        Eigen::Matrix3d F32   = util.getFundamentalMatrix(K.inverse(), R32, T32);

        Eigen::Vector3d coeffshypo1_1 = F31*hypo1_pt1;
        Eigen::Vector3d coeffshypo1_2 = F31*hypo1_pt2;
        Eigen::Vector3d coeffshypo2_1 = F32*hypo2_pt1;
        Eigen::Vector3d coeffshypo2_2 = F32*hypo2_pt2;

        Eigen::MatrixXd QuadrilateralPoints;
        QuadrilateralPoints.conservativeResize(4,2);
        QuadrilateralPoints.row(0)<< double(coeffshypo1_1(1)*coeffshypo2_1(2) - coeffshypo2_1(1)*coeffshypo1_1(2))/double(coeffshypo1_1(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_1(1)), double(coeffshypo1_1(2)*coeffshypo2_1(0) - coeffshypo2_1(2)*coeffshypo1_1(0))/double(coeffshypo1_1(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_1(1));
        QuadrilateralPoints.row(1)<< double(coeffshypo1_1(1)*coeffshypo2_2(2) - coeffshypo2_2(1)*coeffshypo1_1(2))/double(coeffshypo1_1(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_1(1)), double(coeffshypo1_1(2)*coeffshypo2_2(0) - coeffshypo2_2(2)*coeffshypo1_1(0))/double(coeffshypo1_1(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_1(1));
        QuadrilateralPoints.row(2)<< double(coeffshypo1_2(1)*coeffshypo2_1(2) - coeffshypo2_1(1)*coeffshypo1_2(2))/double(coeffshypo1_2(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_2(1)), double(coeffshypo1_2(2)*coeffshypo2_1(0) - coeffshypo2_1(2)*coeffshypo1_2(0))/double(coeffshypo1_2(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_2(1));
        QuadrilateralPoints.row(3)<< double(coeffshypo1_2(1)*coeffshypo2_2(2) - coeffshypo2_2(1)*coeffshypo1_2(2))/double(coeffshypo1_2(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_2(1)), double(coeffshypo1_2(2)*coeffshypo2_2(0) - coeffshypo2_2(2)*coeffshypo1_2(0))/double(coeffshypo1_2(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_2(1));
        
        return QuadrilateralPoints;
    }
    
    

}

#endif
