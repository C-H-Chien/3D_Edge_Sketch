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
using namespace std;

namespace GetOrientationList {
    
    get_OrientationList::get_OrientationList( ) { }
    
    Eigen::MatrixXd get_OrientationList::getOreListBar(Eigen::MatrixXd Edges_HYPO1, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2) {
        MultiviewGeometryUtil::multiview_geometry_util util;
        Eigen::MatrixXd OreListBar_raw;
        OreListBar_raw.conservativeResize(Edges_HYPO1.rows(),2);
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
