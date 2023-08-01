#ifndef UTIL_HPP
#define UTIL_HPP
// =============================================================================
//
// Modifications
//    Chiang-Heng Chien  23-07-14:   Intiailly Created for Multiview Geometry
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ==============================================================================
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

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace MultiviewGeometryUtil {
    
    class multiview_geometry_util {

    public:
        multiview_geometry_util();
        
        Eigen::Matrix3d getSkewSymmetric( Eigen::Vector3d T );
        Eigen::Matrix3d getEssentialMatrix( Eigen::Matrix3d R21, Eigen::Vector3d T21 );
        Eigen::Matrix3d getFundamentalMatrix( Eigen::Matrix3d inverse_K, Eigen::Matrix3d R21, Eigen::Vector3d T21 );
        
        Eigen::Matrix3d getRelativePose_R21(Eigen::Matrix3d R1, Eigen::Matrix3d R2);
        Eigen::Vector3d getRelativePose_T21(Eigen::Matrix3d R1, Eigen::Matrix3d R2, Eigen::Vector3d T1, Eigen::Vector3d T2);
        


    private:
        
    };

}


#endif
