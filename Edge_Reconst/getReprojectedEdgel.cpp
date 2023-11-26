#ifndef GETREPROJECTEDEDGEL_CPP
#define GETREPROJECTEDEDGEL_CPP
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

#include "getReprojectedEdgel.hpp"
#include "util.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace GetReprojectedEdgel {
    
    get_Reprojected_Edgel::get_Reprojected_Edgel( ) { }
    
    Eigen::Vector3d get_Reprojected_Edgel::getTGT_Meters(Eigen::MatrixXd pt_edge, Eigen::Matrix3d K) {
        Eigen::Vector3d pt_edgel;
        pt_edgel << pt_edge(0,0), pt_edge(0,1), 1;
        Eigen::Vector3d Gamma = K.inverse() * pt_edgel;
        Eigen::Vector2d tangent;
        tangent << cos(pt_edge(0,2)), sin(pt_edge(0,2));
        Eigen::Vector3d tgt_to_pixels;
        tgt_to_pixels << (tangent(0)+pt_edgel(0)), (tangent(1)+pt_edgel(1)), 1;
        Eigen::Vector3d tgt_to_meters = K.inverse() * tgt_to_pixels;
        Eigen::Vector3d tgt_meters = tgt_to_meters - Gamma;
        return tgt_meters;
    }

    Eigen::MatrixXd get_Reprojected_Edgel::getGamma3Pos(Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd edgels_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2, Eigen::Matrix3d K3) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::Vector3d pt_edgel_HYPO1;
        pt_edgel_HYPO1 << pt_edge_HYPO1(0,0), pt_edge_HYPO1(0,1), 1;
        Eigen::Vector3d Gamma1 = K1.inverse() * pt_edgel_HYPO1;
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[HYPO1_VIEW_INDX];
        Eigen::Vector3d T1  = All_T[HYPO1_VIEW_INDX];
        Eigen::Matrix3d R2  = All_R[HYPO2_VIEW_INDX];
        Eigen::Vector3d T2  = All_T[HYPO2_VIEW_INDX];
        Eigen::Matrix3d R3  = All_R[VALID_INDX];
        Eigen::Vector3d T3  = All_T[VALID_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
        Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);
        Eigen::MatrixXd edge_pos_gamma3;
        edge_pos_gamma3.conservativeResize(edgels_HYPO2.rows(),2);
        for (int idx_HYPO2 = 0; idx_HYPO2 < edgels_HYPO2.rows(); idx_HYPO2++){
            Eigen::MatrixXd pt_edge_HYPO2 = edgels_HYPO2.row(idx_HYPO2);
            Eigen::Vector3d pt_edgel_HYPO2;
            pt_edgel_HYPO2 << pt_edge_HYPO2(0,0), pt_edge_HYPO2(0,1), 1;
            Eigen::Vector3d Gamma2 = K2.inverse() * pt_edgel_HYPO2;
            double rho1 = (double(e1.transpose() * T21) - double(e3.transpose() * T21) * double(e1.transpose() *Gamma2))/(double(e3.transpose() * R21 * Gamma1)* double(e1.transpose() * Gamma2) - double(e1.transpose() * R21 * Gamma1));
            double rho3 = rho1 * double(e3.transpose() * R31 * Gamma1) + double(e3.transpose()*T31);
            Eigen::Vector3d Gamma3 = 1 / rho3 * (R31 * rho1 * Gamma1 + T31);
            Eigen::Vector3d point3 = K3 * Gamma3;
            edge_pos_gamma3.row(idx_HYPO2) << point3(0), point3(1);
        }
        return edge_pos_gamma3;
    }
    
    Eigen::MatrixXd get_Reprojected_Edgel::getGamma3Tgt(Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd edgels_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::Vector3d pt_edgel_HYPO1;
        pt_edgel_HYPO1 << pt_edge_HYPO1(0,0), pt_edge_HYPO1(0,1), 1;
        Eigen::Vector3d Gamma1 = K1.inverse() * pt_edgel_HYPO1;
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[HYPO1_VIEW_INDX];
        Eigen::Vector3d T1  = All_T[HYPO1_VIEW_INDX];
        Eigen::Matrix3d R2  = All_R[HYPO2_VIEW_INDX];
        Eigen::Vector3d T2  = All_T[HYPO2_VIEW_INDX];
        Eigen::Matrix3d R3  = All_R[VALID_INDX];
        Eigen::Vector3d T3  = All_T[VALID_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
        Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);
        Eigen::Vector3d tgt1_meters = getTGT_Meters(pt_edge_HYPO1, K1);
        Eigen::MatrixXd edge_tgt_gamma3;
        edge_tgt_gamma3.conservativeResize(edgels_HYPO2.rows(),2);
        for (int idx_HYPO2 = 0; idx_HYPO2 < edgels_HYPO2.rows(); idx_HYPO2++){
            Eigen::MatrixXd pt_edge_HYPO2 = edgels_HYPO2.row(idx_HYPO2);
            Eigen::Vector3d pt_edgel_HYPO2;
            pt_edgel_HYPO2 << pt_edge_HYPO2(0,0), pt_edge_HYPO2(0,1), 1;
            Eigen::Vector3d Gamma2 = K2.inverse() * pt_edgel_HYPO2;
            double rho1 = (double(e1.transpose() * T21) - double(e3.transpose() * T21) * double(e1.transpose() *Gamma2))/(double(e3.transpose() * R21 * Gamma1)* double(e1.transpose() * Gamma2) - double(e1.transpose() * R21 * Gamma1));
            double rho3 = rho1 * double(e3.transpose() * R31 * Gamma1) + double(e3.transpose()*T31);
            Eigen::Vector3d Gamma3 = 1 / rho3 * (R31 * rho1 * Gamma1 + T31);
            Eigen::Vector3d tgt2_meters = getTGT_Meters(pt_edge_HYPO2, K2);
            Eigen::Vector3d n1 = tgt1_meters.cross(Gamma1);
            Eigen::Vector3d n2 = R21.transpose() * tgt2_meters.cross(Gamma2);
            Eigen::Vector3d T_v1 = n1.cross(n2) / (n1.cross(n2) ).norm();
            Eigen::Vector3d T_v3 = R31 * T_v1;
            Eigen::Vector3d Tt_v3 = T_v3 - double(e3.transpose()*T_v3)*Gamma3;
            Eigen::Vector3d t_v3  = Tt_v3 / Tt_v3.norm();
            if(IF_ICLNUIM_DATASET == 1){
                t_v3(1) = -1*t_v3(1);
            }
            edge_tgt_gamma3.row(idx_HYPO2) << t_v3(0), t_v3(1);
        }
        return edge_tgt_gamma3;
    }

}

#endif
