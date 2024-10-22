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

#include <stdio.h>
#include <stdlib.h>
//using namespace std;

namespace GetQuadrilateral {
    
    get_Quadrilateral::get_Quadrilateral( ) { }
    
    Eigen::MatrixXd get_Quadrilateral::getQuadrilateralPoints(int hyp01_view_indx, int hyp02_view_indx, Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd pt_edge_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2, Eigen::Matrix3d K3) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e2  = {0,1,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[hyp01_view_indx];
        Eigen::Vector3d T1  = All_T[hyp01_view_indx];
        Eigen::Matrix3d R2  = All_R[hyp02_view_indx];
        Eigen::Vector3d T2  = All_T[hyp02_view_indx];
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
        Eigen::Vector3d epipole_pix_view1 = K1 * epipole_met_view1;
        Eigen::Vector3d epipole_pix_view2 = K2 * epipole_met_view2;

        double slope_hypo1        = double(pt_edge_HYPO1(1)-epipole_pix_view1(1))/double(pt_edge_HYPO1(0)-epipole_pix_view1(0));
        double intersect_hypo1    = double(pt_edge_HYPO1(1)) - slope_hypo1*double(pt_edge_HYPO1(0));
        Eigen::Vector3d hypo1_pt1 = {double(pt_edge_HYPO1(0)), intersect_hypo1+slope_hypo1*double(pt_edge_HYPO1(0))+DIST_THRESH, 1};
        Eigen::Vector3d hypo1_pt2 = {double(pt_edge_HYPO1(0)), intersect_hypo1+slope_hypo1*double(pt_edge_HYPO1(0))-DIST_THRESH, 1};
        double slope_hypo2        = double(pt_edge_HYPO2(1)-epipole_pix_view2(1))/double(pt_edge_HYPO2(0)-epipole_pix_view2(0));
        double intersect_hypo2    = double(pt_edge_HYPO2(1)) - slope_hypo2*double(pt_edge_HYPO2(0));
        Eigen::Vector3d hypo2_pt1 = {double(pt_edge_HYPO2(0)), intersect_hypo2+slope_hypo2*double(pt_edge_HYPO2(0))+DIST_THRESH, 1};
        Eigen::Vector3d hypo2_pt2 = {double(pt_edge_HYPO2(0)), intersect_hypo2+slope_hypo2*double(pt_edge_HYPO2(0))-DIST_THRESH, 1};

        Eigen::Matrix3d F31   = util.getFundamentalMatrix(K3.inverse(), K1.inverse(), R31, T31);
        Eigen::Matrix3d F32   = util.getFundamentalMatrix(K3.inverse(), K2.inverse(), R32, T32);

        Eigen::Vector3d coeffshypo1_1 = F31*hypo1_pt1;
        Eigen::Vector3d coeffshypo1_2 = F31*hypo1_pt2;
        Eigen::Vector3d coeffshypo2_1 = F32*hypo2_pt1;
        Eigen::Vector3d coeffshypo2_2 = F32*hypo2_pt2;

        Eigen::Vector2d inter1((coeffshypo1_1(1)*coeffshypo2_1(2) - coeffshypo2_1(1)*coeffshypo1_1(2))/(coeffshypo1_1(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_1(1)), (coeffshypo1_1(2)*coeffshypo2_1(0) - coeffshypo2_1(2)*coeffshypo1_1(0))/(coeffshypo1_1(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_1(1)));
        Eigen::Vector2d inter2((coeffshypo1_1(1)*coeffshypo2_2(2) - coeffshypo2_2(1)*coeffshypo1_1(2))/(coeffshypo1_1(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_1(1)), (coeffshypo1_1(2)*coeffshypo2_2(0) - coeffshypo2_2(2)*coeffshypo1_1(0))/(coeffshypo1_1(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_1(1)));
        Eigen::Vector2d inter3((coeffshypo1_2(1)*coeffshypo2_1(2) - coeffshypo2_1(1)*coeffshypo1_2(2))/(coeffshypo1_2(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_2(1)), (coeffshypo1_2(2)*coeffshypo2_1(0) - coeffshypo2_1(2)*coeffshypo1_2(0))/(coeffshypo1_2(0)*coeffshypo2_1(1) - coeffshypo2_1(0)*coeffshypo1_2(1)));
        Eigen::Vector2d inter4((coeffshypo1_2(1)*coeffshypo2_2(2) - coeffshypo2_2(1)*coeffshypo1_2(2))/(coeffshypo1_2(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_2(1)), (coeffshypo1_2(2)*coeffshypo2_2(0) - coeffshypo2_2(2)*coeffshypo1_2(0))/(coeffshypo1_2(0)*coeffshypo2_2(1) - coeffshypo2_2(0)*coeffshypo1_2(1)));

        Eigen::Vector2d v1 = inter2 - inter1;
        Eigen::Vector2d v2 = inter3 - inter1;
        Eigen::Vector2d v3 = inter4 - inter1;

        Eigen::Vector3d angle_list(acos(v1.dot(v2)/(v1.norm()*v2.norm())), acos(v1.dot(v3)/(v1.norm()*v3.norm())),acos(v2.dot(v3)/(v2.norm()*v3.norm())));

        Eigen::MatrixXd QuadrilateralPoints;
        QuadrilateralPoints.conservativeResize(4,2);

        if(angle_list(0) == angle_list.maxCoeff()){
            QuadrilateralPoints.row(0)<< inter1(0),inter1(1);
            QuadrilateralPoints.row(1)<< inter2(0),inter2(1);
            QuadrilateralPoints.row(2)<< inter4(0),inter4(1);
            QuadrilateralPoints.row(3)<< inter3(0),inter3(1);
        }else if(angle_list(1) ==angle_list.maxCoeff()){
            QuadrilateralPoints.row(0)<< inter1(0),inter1(1);
            QuadrilateralPoints.row(1)<< inter2(0),inter2(1);
            QuadrilateralPoints.row(2)<< inter3(0),inter3(1);
            QuadrilateralPoints.row(3)<< inter4(0),inter4(1);
        }else{
            QuadrilateralPoints.row(0)<< inter1(0),inter1(1);
            QuadrilateralPoints.row(1)<< inter3(0),inter3(1);
            QuadrilateralPoints.row(2)<< inter2(0),inter2(1);
            QuadrilateralPoints.row(3)<< inter4(0),inter4(1);
        }
        
        return QuadrilateralPoints;
    }
    
    Eigen::MatrixXd get_Quadrilateral::getInliner(int hyp01_view_indx, int hyp02_view_indx, Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd pt_edge_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2, Eigen::Matrix3d K3, Eigen::MatrixXd TO_Edges_VALID) {
        Eigen::MatrixXd QuadrilateralPoints = getQuadrilateralPoints(hyp01_view_indx, hyp02_view_indx, pt_edge_HYPO1, pt_edge_HYPO2, All_R, All_T, VALID_INDX, K1, K2, K3);
        Eigen::Vector2d v1 = {double(QuadrilateralPoints(1,0) - QuadrilateralPoints(0,0)),double(QuadrilateralPoints(1,1) - QuadrilateralPoints(0,1))};
        Eigen::Vector2d v2 = {double(QuadrilateralPoints(2,0) - QuadrilateralPoints(0,0)),double(QuadrilateralPoints(2,1) - QuadrilateralPoints(0,1))};
        Eigen::Vector2d v3 = {double(QuadrilateralPoints(3,0) - QuadrilateralPoints(0,0)),double(QuadrilateralPoints(3,1) - QuadrilateralPoints(0,1))};
        Eigen::Vector3d anglelist;
        anglelist(0) = acos(v1.dot(v2)/(v1.norm()*v2.norm()));
        anglelist(1) = acos(v1.dot(v3)/(v1.norm()*v3.norm()));
        anglelist(2) = acos(v2.dot(v3)/(v2.norm()*v3.norm())); 
        Eigen::Index   maxIndex;
        double maxangle = anglelist.maxCoeff(&maxIndex); 
        
        Eigen::Vector2d inter_idx;
        if(maxIndex == 0){
            inter_idx = {1,2};
        }else if(maxIndex == 1){
            inter_idx = {1,3};
        }else{
            inter_idx = {2,3};
        }
        
        Eigen::Matrix2d V_matrix;
        V_matrix << double(QuadrilateralPoints(inter_idx(0),0) - QuadrilateralPoints(0,0)), double(QuadrilateralPoints(inter_idx(1),0) - QuadrilateralPoints(0,0)), 
                    double(QuadrilateralPoints(inter_idx(0),1) - QuadrilateralPoints(0,1)), double(QuadrilateralPoints(inter_idx(1),1) - QuadrilateralPoints(0,1));
        Eigen::MatrixXd TO_Edges_VALID_XY;
        TO_Edges_VALID_XY.conservativeResize(TO_Edges_VALID.rows(),2);
        TO_Edges_VALID_XY.col(0) = TO_Edges_VALID.col(0);
        TO_Edges_VALID_XY.col(1) = TO_Edges_VALID.col(1);
        Eigen::MatrixXd intersection1;
        intersection1.conservativeResize(TO_Edges_VALID.rows(),2);
        intersection1.col(0) = Eigen::VectorXd::Ones(intersection1.rows())*QuadrilateralPoints(0,0);
        intersection1.col(1) = Eigen::VectorXd::Ones(intersection1.rows())*QuadrilateralPoints(0,1);
        Eigen::MatrixXd Val_ab = V_matrix.inverse() * (TO_Edges_VALID_XY.transpose() - intersection1.transpose());
        
        Eigen::MatrixXd inliner;
        int idx_inliner = 0;
        for (int idx_VALI = 0; idx_VALI < Val_ab.cols(); idx_VALI++){
            if(double(Val_ab(0,idx_VALI)) >= 0 && double(Val_ab(1,idx_VALI)) >= 0 && double(Val_ab(0,idx_VALI)) <= 1 && double(Val_ab(1,idx_VALI)) <= 1){
                inliner.conservativeResize(idx_inliner+1,1);
                inliner.row(idx_inliner) << double(idx_VALI);
                idx_inliner ++;
            }
        }
        return inliner;
    }

}

#endif
