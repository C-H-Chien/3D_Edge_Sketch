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
        Eigen::Matrix3d getFundamentalMatrix(Eigen::Matrix3d inverse_K1, Eigen::Matrix3d inverse_K2, Eigen::Matrix3d R21, Eigen::Vector3d T21);
        
        Eigen::Matrix3d getRelativePose_R21(Eigen::Matrix3d R1, Eigen::Matrix3d R2);
        Eigen::Vector3d getRelativePose_T21(Eigen::Matrix3d R1, Eigen::Matrix3d R2, Eigen::Vector3d T1, Eigen::Vector3d T2);
        
        Eigen::Vector3d linearTriangulation(int N, const std::vector<Eigen::Vector2d> pts,  \
                                                   const std::vector<Eigen::Matrix3d> & Rs, \
                                                   const std::vector<Eigen::Vector3d> & Ts, \
                                                   const Eigen::Matrix3d K);

        void getRelativePoses( Eigen::Matrix3d R1, Eigen::Vector3d T1, Eigen::Matrix3d R2, Eigen::Vector3d T2, 
                               Eigen::Matrix3d &R21, Eigen::Vector3d &T21, Eigen::Matrix3d &R12, Eigen::Vector3d &T12 )
        {
            R21 = getRelativePose_R21(R1, R2);
            T21 = getRelativePose_T21(R1, R2, T1, T2);
            R12 = getRelativePose_R21(R2, R1);
            T12 = getRelativePose_T21(R2, R1, T2, T1);  
        }

        Eigen::Vector3d get3DTangentFromTwo2Dtangents( 
            const Eigen::MatrixXd pt_edge_view1, const Eigen::MatrixXd pt_edge_view2,
            const Eigen::Matrix3d K1,  const Eigen::Matrix3d K2,
            const Eigen::Matrix3d R1,  const Eigen::Vector3d T1,
            const Eigen::Matrix3d R2,  const Eigen::Vector3d T2 );

        Eigen::Vector3d transformToWorldCoordinates( const Eigen::Vector3d& point, const Eigen::Matrix3d& R, const Eigen::Vector3d& T) 
        {
            // Apply the inverse transformation to convert the point to world coordinates
            return R.transpose() * (point - T);
        }

        std::vector<double> check_reproj_error(std::vector<Eigen::Vector2d> points_2D, Eigen::Vector3d point_3D, 
                                               std::vector<Eigen::Matrix3d> Rs, std::vector<Eigen::Vector3d> Ts, Eigen::Matrix3d K);
    private:
        
    };

}


#endif
