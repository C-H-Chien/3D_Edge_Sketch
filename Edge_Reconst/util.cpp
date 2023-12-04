#ifndef UTIL_CPP
#define UTIL_CPP
// ====================================================================================================
//
// Modifications
//    Chiang-Heng Chien  23-07-14:   Intiailly Created. Some functions are shifted from my ICCV code.
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
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

#include "util.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace MultiviewGeometryUtil {
    
    multiview_geometry_util::multiview_geometry_util( ) { }

    Eigen::Matrix3d multiview_geometry_util::getSkewSymmetric(Eigen::Vector3d T) {
        Eigen::Matrix3d T_x = (Eigen::Matrix3d() << 0.,  -T(2),   T(1), T(2),  0.,  -T(0), -T(1),  T(0),   0.).finished();
        return T_x;
    }

    Eigen::Matrix3d multiview_geometry_util::getEssentialMatrix( Eigen::Matrix3d R21, Eigen::Vector3d T21 ) {
        //> E21 = (skew_T(T21)*R21);
        Eigen::Matrix3d T21_x = getSkewSymmetric(T21);
        return T21_x * R21;
    }

    Eigen::Matrix3d multiview_geometry_util::getFundamentalMatrix(Eigen::Matrix3d inverse_K1, Eigen::Matrix3d inverse_K2, Eigen::Matrix3d R21, Eigen::Vector3d T21) {
        //> F21 = inv_K'*(skew_T(T21)*R21)*inv_K;
        Eigen::Matrix3d T21_x = getSkewSymmetric(T21);
        return inverse_K2.transpose() * (T21_x * R21) * inverse_K1;
    }

    Eigen::Matrix3d multiview_geometry_util::getRelativePose_R21(Eigen::Matrix3d R1, Eigen::Matrix3d R2) {
        Eigen::Matrix3d R_1; 
        Eigen::Matrix3d R_2;
        R_1 = R1;
        R_2 = R2;
        /*
        if(IF_TLESS_DATASET == 1){
            R_1 = R1.transpose();
            R_2 = R2.transpose();
        }else{
            R_1 = R1;
            R_2 = R2;
        }
        */
        return R_2* R_1.transpose();
    }

    Eigen::Vector3d multiview_geometry_util::getRelativePose_T21(Eigen::Matrix3d R1, Eigen::Matrix3d R2, Eigen::Vector3d T1, Eigen::Vector3d T2) {
        Eigen::Vector3d C1 = -1*R1.transpose() * T1;
        Eigen::Vector3d C2 = -1*R2.transpose() * T2;
        Eigen::Matrix3d R_1; 
        Eigen::Matrix3d R_2;
        R_1 = R1;
        R_2 = R2;
        /*
        if(IF_TLESS_DATASET == 1){ 
            R_1 = R1.transpose();
            R_2 = R2.transpose();
            C1 = -1*C1;
            C2 = -1*C2;
        }else{
            R_1 = R1;
            R_2 = R2;
        }
        */
        return R_2 * (C1 - C2); 
    }
}

#endif
