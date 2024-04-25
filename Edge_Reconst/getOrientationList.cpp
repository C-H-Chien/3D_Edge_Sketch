#ifndef GETOREINTATIONLIST_CPP
#define GETOREINTATIONLIST_CPP
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

#include "getOrientationList.hpp"
#include "util.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

#include <stdio.h>
#include <stdlib.h>
// using namespace std;

namespace GetOrientationList {
    
    get_OrientationList::get_OrientationList( ) { }
    
    Eigen::MatrixXd get_OrientationList::getOreListBar(Eigen::MatrixXd Edges_HYPO1, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2, int VALID_INDX, int REFIDX) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::MatrixXd OreListBar_raw;
        OreListBar_raw.conservativeResize(Edges_HYPO1.rows(),2);
        Eigen::MatrixXd OreListBardegree;
        OreListBardegree.conservativeResize(Edges_HYPO1.rows(),2);
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e2  = {0,1,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[REFIDX];
        Eigen::Vector3d T1  = All_T[REFIDX];
        Eigen::Matrix3d R2  = All_R[VALID_INDX];
        Eigen::Vector3d T2  = All_T[VALID_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d F   = util.getFundamentalMatrix(K1.inverse(), K2.inverse(), R21, T21);
        //
        Eigen::Vector3d epipole_met_view1 = {double(e1.transpose() * R21.transpose() *T21) / double(e3.transpose()*R21.transpose()*T21), double(e2.transpose()*R21.transpose()*T21) / double(e3.transpose()*R21.transpose()*T21), 1};
        Eigen::Vector3d epipole_met_view2 = {double(e1.transpose()*T21) / double(e3.transpose()*T21), double(e2.transpose()*T21) / double(e3.transpose()*T21), 1};
        Eigen::Vector3d epipole_pix_view1 = K1 * epipole_met_view1;
        Eigen::Vector3d epipole_pix_view2 = K2 * epipole_met_view2;
        //
        Eigen::MatrixXd DyDx;
        DyDx.conservativeResize(Edges_HYPO1.rows(),2);
        DyDx.col(0) = Edges_HYPO1.col(1) - Eigen::VectorXd::Ones(Edges_HYPO1.rows())*epipole_pix_view1(1);
        DyDx.col(1) = Edges_HYPO1.col(0) - Eigen::VectorXd::Ones(Edges_HYPO1.rows())*epipole_pix_view1(0);
        //
        Eigen::MatrixXd slope_hypo1                     = DyDx.col(0).array() / DyDx.col(1).array();
        Eigen::MatrixXd dist_hypo1pt_epipoleall         = DyDx.col(0).array() * DyDx.col(0).array()+ DyDx.col(1).array() * DyDx.col(1).array();
        Eigen::MatrixXd dist_hypo1pt12_epipoleallsquare = dist_hypo1pt_epipoleall.col(0) - Eigen::VectorXd::Ones(dist_hypo1pt_epipoleall.rows())*delta*delta;
        Eigen::MatrixXd dist_hypo1pt12_epipoleall       = dist_hypo1pt12_epipoleallsquare.array().sqrt();
        //
        Eigen::MatrixXd thetahypo1all        = ((Eigen::VectorXd::Ones(dist_hypo1pt_epipoleall.rows())*delta).array() / (dist_hypo1pt12_epipoleall.col(0).array())).array().asin();
        Eigen::MatrixXd anglehypo1all        = slope_hypo1.col(0).array().atan();
        Eigen::MatrixXd angle_theta_hypo1all;
        angle_theta_hypo1all.conservativeResize(Edges_HYPO1.rows(),2);
        angle_theta_hypo1all.col(0) = anglehypo1all.col(0)+thetahypo1all.col(0); 
        angle_theta_hypo1all.col(1) = anglehypo1all.col(0)-thetahypo1all.col(0);
        //
        Eigen::MatrixXd dxdyhypo1;
        Eigen::MatrixXd dxdyhypo2;
        dxdyhypo1.conservativeResize(Edges_HYPO1.rows(),2);
        dxdyhypo2.conservativeResize(Edges_HYPO1.rows(),2);
        dxdyhypo1.col(0) = (angle_theta_hypo1all.col(0).array().cos()*dist_hypo1pt12_epipoleall.col(0).array()).array().abs();
        dxdyhypo1.col(1) = (angle_theta_hypo1all.col(0).array().sin()*dist_hypo1pt12_epipoleall.col(0).array()).array().abs();
        dxdyhypo2.col(0) = (angle_theta_hypo1all.col(1).array().cos()*dist_hypo1pt12_epipoleall.col(0).array()).array().abs();
        dxdyhypo2.col(1) = (angle_theta_hypo1all.col(1).array().sin()*dist_hypo1pt12_epipoleall.col(0).array()).array().abs();

        //
        for(int idx_dxdy = 0; idx_dxdy < dxdyhypo1.rows(); idx_dxdy++){
            double dxall = DyDx(idx_dxdy,1);
            double dyall = DyDx(idx_dxdy,0);
            double dx1   = dxdyhypo1(idx_dxdy,0);
            double dy1   = dxdyhypo1(idx_dxdy,1);
            double dx2   = dxdyhypo2(idx_dxdy,0);
            double dy2   = dxdyhypo2(idx_dxdy,1);
            //
            if(dxall < 0){
                dxdyhypo1(idx_dxdy,0) = -dx1;
                dxdyhypo2(idx_dxdy,0) = -dx2;
            }
            //
            if(dyall < 0){
                dxdyhypo1(idx_dxdy,1) = -dy1;
                dxdyhypo2(idx_dxdy,1) = -dy2;
            }
        }

        Eigen::MatrixXd hypo1_pt_1all;
        Eigen::MatrixXd hypo1_pt_2all;
        hypo1_pt_1all.conservativeResize(Edges_HYPO1.rows(),3);
        hypo1_pt_2all.conservativeResize(Edges_HYPO1.rows(),3);
        hypo1_pt_1all.col(0) = dxdyhypo1.col(0) + Eigen::VectorXd::Ones(dxdyhypo1.rows())* epipole_pix_view1(0);
        hypo1_pt_1all.col(1) = dxdyhypo1.col(1) + Eigen::VectorXd::Ones(dxdyhypo1.rows())* epipole_pix_view1(1);
        hypo1_pt_1all.col(2) = Eigen::VectorXd::Ones(dxdyhypo1.rows());
        hypo1_pt_2all.col(0) = dxdyhypo2.col(0) + Eigen::VectorXd::Ones(dxdyhypo1.rows())* epipole_pix_view1(0);
        hypo1_pt_2all.col(1) = dxdyhypo2.col(1) + Eigen::VectorXd::Ones(dxdyhypo1.rows())* epipole_pix_view1(1);
        hypo1_pt_2all.col(2) = Eigen::VectorXd::Ones(dxdyhypo1.rows());
        
        Eigen::MatrixXd coeffspt1T = F * hypo1_pt_1all.transpose();
        Eigen::MatrixXd coeffspt1  = coeffspt1T.transpose();
        Eigen::MatrixXd Apixel_1 = coeffspt1.col(0);
        Eigen::MatrixXd Bpixel_1 = coeffspt1.col(1);
        Eigen::MatrixXd Cpixel_1 = coeffspt1.col(2);

        Eigen::MatrixXd coeffspt2T = F * hypo1_pt_2all.transpose();
        Eigen::MatrixXd coeffspt2  = coeffspt2T.transpose();
        Eigen::MatrixXd Apixel_2 = coeffspt2.col(0);
        Eigen::MatrixXd Bpixel_2 = coeffspt2.col(1);
        Eigen::MatrixXd Cpixel_2 = coeffspt2.col(2);
        Eigen::MatrixXd slope_hypo_pt1 = Bpixel_1.col(0).array()/Apixel_1.col(0).array();
        Eigen::MatrixXd ore_list1bar_1 = slope_hypo_pt1.col(0).array().atan()/M_PI*180;

        Eigen::MatrixXd slope_hypo_pt2 = Bpixel_2.col(0).array()/Apixel_2.col(0).array();
        Eigen::MatrixXd ore_list1bar_2 = slope_hypo_pt2.col(0).array().atan()/M_PI*180;
        Eigen::MatrixXd ore_list1bar_all;
        ore_list1bar_all.conservativeResize(ore_list1bar_1.rows(),2);
        for(int idx_ore12 = 0; idx_ore12 < slope_hypo_pt1.rows(); idx_ore12++){
            if(ore_list1bar_1(idx_ore12,0) < 0){
                ore_list1bar_1(idx_ore12,0) = ore_list1bar_1(idx_ore12,0)+180;
            }
            //
            if(ore_list1bar_2(idx_ore12,0) < 0){
                ore_list1bar_2(idx_ore12,0) = ore_list1bar_2(idx_ore12,0)+180;
            }
            if(ore_list1bar_1(idx_ore12,0) < ore_list1bar_2(idx_ore12,0)){
                ore_list1bar_all.row(idx_ore12) << ore_list1bar_1(idx_ore12,0) , ore_list1bar_2(idx_ore12,0);
            }else{
                ore_list1bar_all.row(idx_ore12) << ore_list1bar_2(idx_ore12,0) , ore_list1bar_1(idx_ore12,0);
            }
        }

        // a_hypo2_all = tan(deg2rad(ore_list1bar_all));
        // c_hypo2_all = epipole_pix_view2(2,1) - a_hypo2_all*epipole_pix_view2(1,1);
        Eigen::MatrixXd p1x_hypo_all;
        p1x_hypo_all.conservativeResize(Apixel_1.rows(),4);
        p1x_hypo_all.col(0) = Eigen::VectorXd::Zero(Edges_HYPO1.rows());
        p1x_hypo_all.col(1) = Eigen::VectorXd::Ones(Edges_HYPO1.rows())*imgcols;
        p1x_hypo_all.col(2) = -1*Cpixel_1.col(0).array()/Apixel_1.col(0).array();
        p1x_hypo_all.col(3) = -1*(Cpixel_1 + Bpixel_1*imgrows).array()/Apixel_1.array();

        Eigen::MatrixXd p1y_hypo_all;
        p1y_hypo_all.conservativeResize(Apixel_1.rows(),4);
        p1y_hypo_all.col(0) = -1*Cpixel_1.col(0).array()/Bpixel_1.col(0).array();
        p1y_hypo_all.col(1) = -1*(Cpixel_1 + Apixel_1*imgcols).array()/Bpixel_1.array();
        p1y_hypo_all.col(2) = Eigen::VectorXd::Zero(Edges_HYPO1.rows());
        p1y_hypo_all.col(3) = Eigen::VectorXd::Ones(Edges_HYPO1.rows())*imgrows;

        Eigen::MatrixXd p2x_hypo_all;
        p2x_hypo_all.conservativeResize(Apixel_2.rows(),4);
        p2x_hypo_all.col(0) = Eigen::VectorXd::Zero(Edges_HYPO1.rows());
        p2x_hypo_all.col(1) = Eigen::VectorXd::Ones(Edges_HYPO1.rows())*imgcols;
        p2x_hypo_all.col(2) = -1*Cpixel_2.col(0).array()/Apixel_2.col(0).array();
        p2x_hypo_all.col(3) = -1*(Cpixel_2 + Bpixel_2*imgrows).array()/Apixel_2.array();

        Eigen::MatrixXd p2y_hypo_all;
        p2y_hypo_all.conservativeResize(Apixel_2.rows(),4);
        p2y_hypo_all.col(0) = -1*Cpixel_2.col(0).array()/Bpixel_2.col(0).array();
        p2y_hypo_all.col(1) = -1*(Cpixel_2 + Apixel_2*imgcols).array()/Bpixel_2.array();
        p2y_hypo_all.col(2) = Eigen::VectorXd::Zero(Edges_HYPO1.rows());
        p2y_hypo_all.col(3) = Eigen::VectorXd::Ones(Edges_HYPO1.rows())*imgrows;

        Eigen::MatrixXd p1_all;
        Eigen::MatrixXd p1_dxdy;
        Eigen::MatrixXd p1_xy;
        Eigen::MatrixXd p1_dist;
        Eigen::MatrixXd p1_final_dxdyp;
        Eigen::MatrixXd p1_final_dxdyn;
        p1_all.conservativeResize(p1x_hypo_all.rows(),4);
        p1_dxdy.conservativeResize(p1x_hypo_all.rows(),4);
        p1_xy.conservativeResize(p1x_hypo_all.rows(),2);
        p1_dist.conservativeResize(1,2);
        p1_final_dxdyp.conservativeResize(p1x_hypo_all.rows(),2);
        p1_final_dxdyn.conservativeResize(p1x_hypo_all.rows(),2);

        Eigen::MatrixXd p2_all;
        Eigen::MatrixXd p2_dxdy;
        Eigen::MatrixXd p2_xy;
        Eigen::MatrixXd p2_dist;
        Eigen::MatrixXd p2_final_dxdyp;
        Eigen::MatrixXd p2_final_dxdyn;
        p2_all.conservativeResize(p2x_hypo_all.rows(),4);
        p2_dxdy.conservativeResize(p2x_hypo_all.rows(),4);
        p2_xy.conservativeResize(p2x_hypo_all.rows(),2);
        p2_dist.conservativeResize(1,2);
        p2_final_dxdyp.conservativeResize(p2x_hypo_all.rows(),2);
        p2_final_dxdyn.conservativeResize(p2x_hypo_all.rows(),2);

        Eigen::MatrixXd slope12all;
        slope12all.conservativeResize(4,1);
        for(int idx_pt12 = 0; idx_pt12 < p1_all.rows(); idx_pt12++){
            int column     = 0;
            int columndist = 0;
            if(p1y_hypo_all(idx_pt12,0) >= 0 && p1y_hypo_all(idx_pt12,0) <= imgrows ){
                p1_all(idx_pt12, 0)  = p1x_hypo_all(idx_pt12,0);
                p1_dxdy(idx_pt12, 0) = p1x_hypo_all(idx_pt12,0) - epipole_pix_view2(0);
                p1_all(idx_pt12, 1)  = p1y_hypo_all(idx_pt12,0);
                p1_dxdy(idx_pt12, 1) = p1y_hypo_all(idx_pt12,0) - epipole_pix_view2(1);
                p1_dist(0, 0)        = p1_dxdy(idx_pt12, 0)*p1_dxdy(idx_pt12, 0) + p1_dxdy(idx_pt12, 1)*p1_dxdy(idx_pt12, 1);
                column     += 2;
                columndist += 1;
            }
            if(p1y_hypo_all(idx_pt12,1) >= 0 && p1y_hypo_all(idx_pt12,1) <= imgrows ){
                p1_all(idx_pt12, column)    = p1x_hypo_all(idx_pt12,1);
                p1_dxdy(idx_pt12, column)   = p1x_hypo_all(idx_pt12,1) - epipole_pix_view2(0);
                p1_all(idx_pt12, column+1)  = p1y_hypo_all(idx_pt12,1);
                p1_dxdy(idx_pt12, column+1) = p1y_hypo_all(idx_pt12,1) - epipole_pix_view2(1);
                p1_dist(0, columndist)      = p1_dxdy(idx_pt12, column)*p1_dxdy(idx_pt12, column) + p1_dxdy(idx_pt12, column+1)*p1_dxdy(idx_pt12, column+1);
                column     += 2;
                columndist += 1;
            }
            if(p1x_hypo_all(idx_pt12,2) >= 0 && p1x_hypo_all(idx_pt12,2) <= imgcols ){
                p1_all(idx_pt12, column)    = p1x_hypo_all(idx_pt12,2);
                p1_dxdy(idx_pt12, column)   = p1x_hypo_all(idx_pt12,2) - epipole_pix_view2(0);
                p1_all(idx_pt12, column+1)  = p1y_hypo_all(idx_pt12,2);
                p1_dxdy(idx_pt12, column+1) = p1y_hypo_all(idx_pt12,2) - epipole_pix_view2(1);
                p1_dist(0, columndist)        = p1_dxdy(idx_pt12, column)*p1_dxdy(idx_pt12, column) + p1_dxdy(idx_pt12, column+1)*p1_dxdy(idx_pt12, column+1);
                column += 2;
            }
            if(p1x_hypo_all(idx_pt12,3) >= 0 && p1x_hypo_all(idx_pt12,3) <= imgcols ){
                p1_all(idx_pt12, 2)  = p1x_hypo_all(idx_pt12,3);
                p1_dxdy(idx_pt12, 2) = p1x_hypo_all(idx_pt12,3) - epipole_pix_view2(0);
                p1_all(idx_pt12, 3)  = p1y_hypo_all(idx_pt12,3);
                p1_dxdy(idx_pt12, 3) = p1y_hypo_all(idx_pt12,3) - epipole_pix_view2(1);
                p1_dist(0, 1)        = p1_dxdy(idx_pt12, 2)*p1_dxdy(idx_pt12, 2) + p1_dxdy(idx_pt12, 3)*p1_dxdy(idx_pt12, 3);
            }
            
            if(p1_dist(0, 0) < p1_dist(0, 1)){
                p1_xy.row(idx_pt12) << p1_all(idx_pt12, 0), p1_all(idx_pt12, 1);
            }else{
                p1_xy.row(idx_pt12) << p1_all(idx_pt12, 2), p1_all(idx_pt12, 3);
            }
            
            double delta_x1 = abs(cos(ore_list1bar_all(idx_pt12,0)/180*M_PI)*delta); 
            double slope1_p = 0;
            double slope1_n = 0;
            if(p1_xy(idx_pt12, 0) == 0 || p1_xy(idx_pt12, 0) == imgcols){
                p1_final_dxdyp(idx_pt12, 0) = p1_xy(idx_pt12,0) - epipole_pix_view2(0);
                p1_final_dxdyp(idx_pt12, 1) = p1_xy(idx_pt12,1) - epipole_pix_view2(1) - delta_x1;
                p1_final_dxdyn(idx_pt12, 0) = p1_xy(idx_pt12,0) - epipole_pix_view2(0);
                p1_final_dxdyn(idx_pt12, 1) = p1_xy(idx_pt12,1) - epipole_pix_view2(1) + delta_x1;
                slope1_p                    = p1_final_dxdyp(idx_pt12, 1)/p1_final_dxdyp(idx_pt12, 0);
                slope1_n                    = p1_final_dxdyn(idx_pt12, 1)/p1_final_dxdyn(idx_pt12, 0);
                slope12all(0,0)             = atan(slope1_p)/M_PI*180;
                slope12all(1,0)             = atan(slope1_n)/M_PI*180;
            }else{
                p1_final_dxdyp(idx_pt12, 0) = p1_xy(idx_pt12,0) - epipole_pix_view2(0) - delta_x1;
                p1_final_dxdyp(idx_pt12, 1) = p1_xy(idx_pt12,1) - epipole_pix_view2(1);
                p1_final_dxdyn(idx_pt12, 0) = p1_xy(idx_pt12,0) - epipole_pix_view2(0) + delta_x1;
                p1_final_dxdyn(idx_pt12, 1) = p1_xy(idx_pt12,1) - epipole_pix_view2(1);
                slope1_p                    = p1_final_dxdyp(idx_pt12, 1)/p1_final_dxdyp(idx_pt12, 0);
                slope1_n                    = p1_final_dxdyn(idx_pt12, 1)/p1_final_dxdyn(idx_pt12, 0);
                slope12all(0,0)             = atan(slope1_p)/M_PI*180;
                slope12all(1,0)             = atan(slope1_n)/M_PI*180;
            }
            //
            column     = 0;
            columndist = 0;
            if(p2y_hypo_all(idx_pt12,0) >= 0 && p2y_hypo_all(idx_pt12,0) <= imgrows ){
                p2_all(idx_pt12, 0)  = p2x_hypo_all(idx_pt12,0);
                p2_dxdy(idx_pt12, 0) = p2x_hypo_all(idx_pt12,0) - epipole_pix_view2(0);
                p2_all(idx_pt12, 1)  = p2y_hypo_all(idx_pt12,0);
                p2_dxdy(idx_pt12, 1) = p2y_hypo_all(idx_pt12,0) - epipole_pix_view2(1);
                p2_dist(0, 0)        = p2_dxdy(idx_pt12, 0)*p2_dxdy(idx_pt12, 0) + p2_dxdy(idx_pt12, 1)*p2_dxdy(idx_pt12, 1);
                column     += 2;
                columndist += 1;
            }
            if(p2y_hypo_all(idx_pt12,1) >= 0 && p2y_hypo_all(idx_pt12,1) <= imgrows ){
                p2_all(idx_pt12, column)    = p2x_hypo_all(idx_pt12,1);
                p2_dxdy(idx_pt12, column)   = p2x_hypo_all(idx_pt12,1) - epipole_pix_view2(0);
                p2_all(idx_pt12, column+1)  = p2y_hypo_all(idx_pt12,1);
                p2_dxdy(idx_pt12, column+1) = p2y_hypo_all(idx_pt12,1) - epipole_pix_view2(1);
                p2_dist(0, columndist)      = p2_dxdy(idx_pt12, column)*p2_dxdy(idx_pt12, column) + p2_dxdy(idx_pt12, column+1)*p2_dxdy(idx_pt12, column+1);
                column     += 2;
                columndist += 1;
            }
            if(p2x_hypo_all(idx_pt12,2) >= 0 && p2x_hypo_all(idx_pt12,2) <= imgcols ){
                p2_all(idx_pt12, column)    = p2x_hypo_all(idx_pt12,2);
                p2_dxdy(idx_pt12, column)   = p2x_hypo_all(idx_pt12,2) - epipole_pix_view2(0);
                p2_all(idx_pt12, column+1)  = p2y_hypo_all(idx_pt12,2);
                p2_dxdy(idx_pt12, column+1) = p2y_hypo_all(idx_pt12,2) - epipole_pix_view2(1);
                p2_dist(0, columndist)      = p2_dxdy(idx_pt12, column)*p2_dxdy(idx_pt12, column) + p2_dxdy(idx_pt12, column+1)*p2_dxdy(idx_pt12, column+1);
                column     += 2;
                columndist += 1;
            }
            if(p2x_hypo_all(idx_pt12,3) >= 0 && p2x_hypo_all(idx_pt12,3) <= imgcols ){
                p2_all(idx_pt12, 2)  = p2x_hypo_all(idx_pt12,3);
                p2_dxdy(idx_pt12, 2) = p2x_hypo_all(idx_pt12,3) - epipole_pix_view2(0);
                p2_all(idx_pt12, 3)  = p2y_hypo_all(idx_pt12,3);
                p2_dxdy(idx_pt12, 3) = p2y_hypo_all(idx_pt12,3) - epipole_pix_view2(1);
                p2_dist(0, 1)        = p2_dxdy(idx_pt12, 2)*p2_dxdy(idx_pt12, 2) + p2_dxdy(idx_pt12, 3)*p2_dxdy(idx_pt12, 3);
            }

            if(p2_dist(0, 0) < p2_dist(0, 1)){
                p2_xy.row(idx_pt12) << p2_all(idx_pt12, 0), p2_all(idx_pt12, 1);
            }else{
                p2_xy.row(idx_pt12) << p2_all(idx_pt12, 2), p2_all(idx_pt12, 3);
            }
            
            double delta_x2 = abs(cos(ore_list1bar_all(idx_pt12,1)/180*M_PI)*delta); 
            double slope2_p = 0;
            double slope2_n = 0;
            if(p2_xy(idx_pt12, 0) == 0 || p2_xy(idx_pt12, 0) == imgcols){
                p2_final_dxdyp(idx_pt12, 0) = p2_xy(idx_pt12,0) - epipole_pix_view2(0);
                p2_final_dxdyp(idx_pt12, 1) = p2_xy(idx_pt12,1) - epipole_pix_view2(1) - delta_x2;
                p2_final_dxdyn(idx_pt12, 0) = p2_xy(idx_pt12,0) - epipole_pix_view2(0);
                p2_final_dxdyn(idx_pt12, 1) = p2_xy(idx_pt12,1) - epipole_pix_view2(1) + delta_x2;
                slope2_p                    = p2_final_dxdyp(idx_pt12, 1)/p2_final_dxdyp(idx_pt12, 0);
                slope2_n                    = p2_final_dxdyn(idx_pt12, 1)/p2_final_dxdyn(idx_pt12, 0);
                slope12all(2,0)             = atan(slope2_p)/M_PI*180;
                slope12all(3,0)             = atan(slope2_n)/M_PI*180;
            }else{
                p2_final_dxdyp(idx_pt12, 0) = p2_xy(idx_pt12,0) - epipole_pix_view2(0) - delta_x2;
                p2_final_dxdyp(idx_pt12, 1) = p2_xy(idx_pt12,1) - epipole_pix_view2(1);
                p2_final_dxdyn(idx_pt12, 0) = p2_xy(idx_pt12,0) - epipole_pix_view2(0) + delta_x2;
                p2_final_dxdyn(idx_pt12, 1) = p2_xy(idx_pt12,1) - epipole_pix_view2(1);
                slope2_p                    = p2_final_dxdyp(idx_pt12, 1)/p2_final_dxdyp(idx_pt12, 0);
                slope2_n                    = p2_final_dxdyn(idx_pt12, 1)/p2_final_dxdyn(idx_pt12, 0);
                slope12all(2,0)             = atan(slope2_p)/M_PI*180;
                slope12all(3,0)             = atan(slope2_n)/M_PI*180;
            }
            /*
            auto itmin = std::min_element(std::begin(slope12all), std::end(slope12all));
            auto const posmin = std::distance(std::begin(slope12all), itmin);
            auto itmax = std::max_element(std::begin(slope12all), std::end(slope12all));
            auto const posmax = std::distance(std::begin(slope12all), itmax);
            */
            if(slope12all(0,0) < 0){
                slope12all(0,0) = slope12all(0,0)+180;
            }
            if(slope12all(1,0) < 0){
                slope12all(1,0) = slope12all(1,0)+180;
            }
            if(slope12all(2,0) < 0){
                slope12all(2,0) = slope12all(2,0)+180;
            }
            if(slope12all(3,0) < 0){
                slope12all(3,0) = slope12all(3,0)+180;
            }
            
            /*
            if(VALID_INDX == 39 && REFIDX == HYPO2_VIEW_INDX){
                std::cout << "VALID_INDX: \n" << VALID_INDX <<std::endl;
                std::cout << "p1_xy: \n" << p1_xy.row(idx_pt12) <<std::endl;
                std::cout << "p2_xy: \n" << p2_xy.row(idx_pt12) <<std::endl;
                std::cout << "p1_dxdy: \n" << p1_dxdy.row(idx_pt12) <<std::endl;
                std::cout << "p1_final_dxdyp: \n" << p1_final_dxdyp.row(idx_pt12) <<std::endl;
                std::cout << "p1_final_dxdyn: \n" << p1_final_dxdyn.row(idx_pt12) <<std::endl;
                std::cout << "p2_dxdy: \n" << p2_dxdy.row(idx_pt12) <<std::endl;
                std::cout << "p2_final_dxdyp: \n" << p2_final_dxdyp.row(idx_pt12) <<std::endl;
                std::cout << "p2_final_dxdyn: \n" << p2_final_dxdyn.row(idx_pt12) <<std::endl;
                std::cout << "delta_x1: \n" << delta_x1 <<std::endl;
                std::cout << "delta_x2: \n" << delta_x2 <<std::endl;
                std::cout << "slope12all: \n" << slope12all <<std::endl;
                if (DEBUG == 1) { std::cerr << "\n—=>DEBUG MODE<=—\n"; exit(1); }
            }
            */
            OreListBar_raw(idx_pt12,0) = slope12all.minCoeff();
            OreListBar_raw(idx_pt12,1) = slope12all.maxCoeff();
            
            /*
            if(OreListBar_raw(idx_pt12,0) <= 0){
                OreListBardegree(idx_pt12,0) = OreListBar_raw(idx_pt12,0)+180;
            }else{
                OreListBardegree(idx_pt12,0) = OreListBar_raw(idx_pt12,0);
            }
            if(OreListBar_raw(idx_pt12,1) <= 0){
                OreListBardegree(idx_pt12,1) = OreListBar_raw(idx_pt12,1)+180;
            }else{
                OreListBardegree(idx_pt12,1) = OreListBar_raw(idx_pt12,1);
            }
            */
            
        }
        
        //OreListBar_raw.col(0) = -1*(Eigen::VectorXd::Ones(Edges_HYPO1.rows())*F(0,0)+F(0,1)*slope_hypo1.col(0));
        //OreListBar_raw.col(1) = (Eigen::VectorXd::Ones(Edges_HYPO1.rows())*F(1,0)+F(1,1)*slope_hypo1.col(0));
        //Eigen::MatrixXd OreListBar = OreListBar_raw.col(0).array()/OreListBar_raw.col(1).array();
        //
        //Eigen::MatrixXd OreListBarAtan = OreListBar.col(0).array().atan();
        //
        
        /*
        Eigen::MatrixXd OreListBar_rawp1 = OreListBar_raw.col(0) + Eigen::VectorXd::Ones(OreListBar_raw.rows())*180;
        Eigen::MatrixXd OreListBar_rawp2 = OreListBar_raw.col(1) + Eigen::VectorXd::Ones(OreListBar_raw.rows())*180;
        OreListBardegree.col(0) = (OreListBar_raw.col(0).array() < 0).select(OreListBar_rawp1, OreListBar_raw.col(0));
        OreListBardegree.col(1) = (OreListBar_raw.col(1).array() < 0).select(OreListBar_rawp2, OreListBar_raw.col(1));
        */

        return OreListBar_raw;
    }

    Eigen::MatrixXd get_OrientationList::getOreListBarVali(Eigen::MatrixXd Edges_HYPO1, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2, int VALID_INDX, int REFIDX) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::MatrixXd OreListBar_raw;
        OreListBar_raw.conservativeResize(Edges_HYPO1.rows(),2);
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e2  = {0,1,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[REFIDX];
        Eigen::Vector3d T1  = All_T[REFIDX];
        Eigen::Matrix3d R2  = All_R[VALID_INDX];
        Eigen::Vector3d T2  = All_T[VALID_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d F   = util.getFundamentalMatrix(K2.inverse(), K1.inverse(), R21, T21);

        Eigen::Vector3d epipole_met_view1 = {double(e1.transpose() * R21.transpose() *T21) / double(e3.transpose()*R21.transpose()*T21), double(e2.transpose()*R21.transpose()*T21) / double(e3.transpose()*R21.transpose()*T21), 1};
        // Eigen::Vector3d epipole_met_view2 = {double(e1.transpose()*T21) / double(e3.transpose()*T21), double(e2.transpose()*T21) / double(e3.transpose()*T21), 1};
        Eigen::Vector3d epipole_pix_view1 = K1 * epipole_met_view1;
        // Eigen::Vector3d epipole_pix_view2 = K2 * epipole_met_view2;

        Eigen::MatrixXd DyDx;
        DyDx.conservativeResize(Edges_HYPO1.rows(),2);
        DyDx.col(0) = Edges_HYPO1.col(1) - Eigen::VectorXd::Ones(Edges_HYPO1.rows())*epipole_pix_view1(1);
        DyDx.col(1) = Edges_HYPO1.col(0) - Eigen::VectorXd::Ones(Edges_HYPO1.rows())*epipole_pix_view1(0);
        Eigen::MatrixXd slope_hypo1 = DyDx.col(0).array()/DyDx.col(1).array();

        OreListBar_raw.col(0) = -1*(Eigen::VectorXd::Ones(Edges_HYPO1.rows())*F(0,0)+F(0,1)*slope_hypo1.col(0));
        OreListBar_raw.col(1) = (Eigen::VectorXd::Ones(Edges_HYPO1.rows())*F(1,0)+F(1,1)*slope_hypo1.col(0));
        Eigen::MatrixXd OreListBar = OreListBar_raw.col(0).array()/OreListBar_raw.col(1).array();
        
        Eigen::MatrixXd OreListBarAtan = OreListBar.col(0).array().atan();
        
        Eigen::MatrixXd OreListBardegree = (OreListBarAtan.col(0)*180)/M_PI;
        Eigen::MatrixXd OreListBardegree1 = OreListBardegree.col(0) + Eigen::VectorXd::Ones(OreListBardegree.rows())*180;
        OreListBardegree = (OreListBardegree.array() < 0).select(OreListBardegree1, OreListBardegree);

        return OreListBardegree;
    }

    Eigen::MatrixXd get_OrientationList::getOreList(Eigen::MatrixXd Edges_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::MatrixXd OreListBar_raw;
        OreListBar_raw.conservativeResize(Edges_HYPO2.rows(),2);
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e2  = {0,1,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[HYPO1_VIEW_INDX];
        Eigen::Vector3d T1  = All_T[HYPO1_VIEW_INDX];
        Eigen::Matrix3d R2  = All_R[HYPO2_VIEW_INDX];
        Eigen::Vector3d T2  = All_T[HYPO2_VIEW_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d F   = util.getFundamentalMatrix(K2.inverse(), K1.inverse(), R21, T21);

        //Eigen::Vector3d epipole_met_view1 = {double(e1.transpose() * R21.transpose() *T21) / double(e3.transpose()*R21.transpose()*T21), double(e2.transpose()*R21.transpose()*T21) / double(e3.transpose()*R21.transpose()*T21), 1};
        Eigen::Vector3d epipole_met_view2 = {double(e1.transpose()*T21) / double(e3.transpose()*T21), double(e2.transpose()*T21) / double(e3.transpose()*T21), 1};
        //Eigen::Vector3d epipole_pix_view1 = K1 * epipole_met_view1;
        Eigen::Vector3d epipole_pix_view2 = K2 * epipole_met_view2;

        Eigen::MatrixXd DyDx;
        DyDx.conservativeResize(Edges_HYPO2.rows(),2);
        DyDx.col(0) = Edges_HYPO2.col(1) - Eigen::VectorXd::Ones(Edges_HYPO2.rows())*epipole_pix_view2(1);
        DyDx.col(1) = Edges_HYPO2.col(0) - Eigen::VectorXd::Ones(Edges_HYPO2.rows())*epipole_pix_view2(0);
        Eigen::MatrixXd OreListBar = DyDx.col(0).array()/DyDx.col(1).array();
        
        //OreListBar_raw.col(0) = -1*(Eigen::VectorXd::Ones(Edges_HYPO2.rows())*F(0,0)+F(0,1)*slope_hypo2.col(0));
        //OreListBar_raw.col(1) = (Eigen::VectorXd::Ones(Edges_HYPO2.rows())*F(1,0)+F(1,1)*slope_hypo2.col(0));
        //Eigen::MatrixXd OreListBar = OreListBar_raw.col(0).array()/OreListBar_raw.col(1).array();
        
        Eigen::MatrixXd OreListBarAtan = OreListBar.col(0).array().atan();
        
        Eigen::MatrixXd OreListBardegree = (OreListBarAtan.col(0)*180)/M_PI;
        Eigen::MatrixXd OreListBardegree1 = OreListBardegree.col(0) + Eigen::VectorXd::Ones(OreListBardegree.rows())*180;
        OreListBardegree = (OreListBardegree.array() < 0).select(OreListBardegree1, OreListBardegree);

        return OreListBardegree;
    }

    Eigen::MatrixXd get_OrientationList::getOreListVali(Eigen::MatrixXd Edges_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2, int VALID_INDX, int REFIDX) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::MatrixXd OreListBar_raw;
        OreListBar_raw.conservativeResize(Edges_HYPO2.rows(),2);
        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e2  = {0,1,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Matrix3d R1  = All_R[REFIDX];
        Eigen::Vector3d T1  = All_T[REFIDX];
        Eigen::Matrix3d R2  = All_R[VALID_INDX];
        Eigen::Vector3d T2  = All_T[VALID_INDX];
        Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
        Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
        Eigen::Matrix3d F   = util.getFundamentalMatrix(K2.inverse(), K1.inverse(), R21, T21);

        //Eigen::Vector3d epipole_met_view1 = {double(e1.transpose() * R21.transpose() *T21) / double(e3.transpose()*R21.transpose()*T21), double(e2.transpose()*R21.transpose()*T21) / double(e3.transpose()*R21.transpose()*T21), 1};
        Eigen::Vector3d epipole_met_view2 = {double(e1.transpose()*T21) / double(e3.transpose()*T21), double(e2.transpose()*T21) / double(e3.transpose()*T21), 1};
        //Eigen::Vector3d epipole_pix_view1 = K1 * epipole_met_view1;
        Eigen::Vector3d epipole_pix_view2 = K2 * epipole_met_view2;

        Eigen::MatrixXd DyDx;
        DyDx.conservativeResize(Edges_HYPO2.rows(),2);
        DyDx.col(0) = Edges_HYPO2.col(1) - Eigen::VectorXd::Ones(Edges_HYPO2.rows())*epipole_pix_view2(1);
        DyDx.col(1) = Edges_HYPO2.col(0) - Eigen::VectorXd::Ones(Edges_HYPO2.rows())*epipole_pix_view2(0);
        Eigen::MatrixXd OreListBar = DyDx.col(0).array()/DyDx.col(1).array();
        
        //OreListBar_raw.col(0) = -1*(Eigen::VectorXd::Ones(Edges_HYPO2.rows())*F(0,0)+F(0,1)*slope_hypo2.col(0));
        //OreListBar_raw.col(1) = (Eigen::VectorXd::Ones(Edges_HYPO2.rows())*F(1,0)+F(1,1)*slope_hypo2.col(0));
        //Eigen::MatrixXd OreListBar = OreListBar_raw.col(0).array()/OreListBar_raw.col(1).array();
        
        Eigen::MatrixXd OreListBarAtan = OreListBar.col(0).array().atan();
        
        Eigen::MatrixXd OreListBardegree = (OreListBarAtan.col(0)*180)/M_PI;
        Eigen::MatrixXd OreListBardegree1 = OreListBardegree.col(0) + Eigen::VectorXd::Ones(OreListBardegree.rows())*180;
        OreListBardegree = (OreListBardegree.array() < 0).select(OreListBardegree1, OreListBardegree);

        return OreListBardegree;
    }

}

#endif
