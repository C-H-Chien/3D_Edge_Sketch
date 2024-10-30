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
        return R_2* R_1.transpose();
    }

    Eigen::Vector3d multiview_geometry_util::getRelativePose_T21(Eigen::Matrix3d R1, Eigen::Matrix3d R2, Eigen::Vector3d T1, Eigen::Vector3d T2) {
        Eigen::Vector3d C1 = -1*R1.transpose() * T1;
        Eigen::Vector3d C2 = -1*R2.transpose() * T2;
        Eigen::Matrix3d R_1; 
        Eigen::Matrix3d R_2;
        R_1 = R1;
        R_2 = R2;
        return R_2 * (C1 - C2); 
    }

    Eigen::Vector3d multiview_geometry_util::linearTriangulation(
        const int N,
        const std::vector<Eigen::Vector2d> pts, 
        const std::vector<Eigen::Matrix3d> & Rs,
        const std::vector<Eigen::Vector3d> & Ts,
        const std::vector<double> & K )
    {
        Eigen::MatrixXd A(2*N, 4);
        Eigen::MatrixXd ATA(4, 4);
        Eigen::Vector4d GAMMA;
        Eigen::Vector3d pt3D; //> Returned variable

        //> Convert points in pixels to points in meters
        std::vector<Eigen::Vector2d> pts_meters;
        Eigen::Matrix3d K_;
        K_<<K[2], 0.,   K[0], 0.,  K[3], K[1], 0.,  0.,   1.;
        // std::cout << "K_: " << K_ <<std::endl;
        for (int p = 0; p < N; p++) {
            Eigen::Vector2d gamma;
            Eigen::Vector3d homo_p{pts[p](0), pts[p](1), 1.0};
            Eigen::Vector3d p_bar = K_.inverse() * homo_p;
            gamma(0) = p_bar(0);
            gamma(1) = p_bar(1);
            pts_meters.push_back(gamma);
        }

        //> We are computing GAMMA under the first camera coordinate,
        //  so R1 is an identity matrix and T1 is all zeros
        A(0,0) = 0.0; A(0,1) = -1.0; A(0,2) = pts_meters[0](1); A(0,3) = 0.0;
        A(1,0) = 1.0; A(1,1) = 0.0; A(1,2) = -pts_meters[0](0); A(1,3) = 0.0;

        int row_cnter = 2;
        for (int p = 0; p < N-1; p++) {
            Eigen::Matrix3d Rp = Rs[p];
            Eigen::Vector3d Tp = Ts[p];
            Eigen::Vector2d mp = pts_meters[p+1];

            // std::cout << "Rp: " << Rp <<std::endl;
            
            double r1 = Rp(0,0), r2 = Rp(0,1), r3 = Rp(0,2), t1 = Tp(0);
            double r4 = Rp(1,0), r5 = Rp(1,1), r6 = Rp(1,2), t2 = Tp(1);
            double r7 = Rp(2,0), r8 = Rp(2,1), r9 = Rp(2,2), t3 = Tp(2);

            A(row_cnter,   0) = mp(1) * r7 - r4;
            A(row_cnter,   1) = mp(1) * r8 - r5; 
            A(row_cnter,   2) = mp(1) * r9 - r6; 
            A(row_cnter,   3) = mp(1) * t3 - t2;
            A(row_cnter+1, 0) = r1 - mp(0) * r7; 
            A(row_cnter+1, 1) = r2 - mp(0) * r8; 
            A(row_cnter+1, 2) = r3 - mp(0) * r9; 
            A(row_cnter+1, 3) = t1 - mp(0) * t3;
            row_cnter += 2;
        }

        //> Solving the homogeneous linear system and divide the first three rows with the last element
        ATA = A.transpose() * A;
        GAMMA = ATA.jacobiSvd(Eigen::ComputeFullV).matrixV().col( ATA.rows() - 1 );
        GAMMA[0] /= GAMMA[3];
        GAMMA[1] /= GAMMA[3];
        GAMMA[2] /= GAMMA[3];
        
        //> Assign GAMMA to the returned point
        pt3D[0] = GAMMA[0];
        pt3D[1] = GAMMA[1];
        pt3D[2] = GAMMA[2];
        
        return pt3D;

        /*
        template<typename matrix_t, typename vector_t>
        void solveNullspaceLU(const matrix_t& A, vector_t& x){
            x = A.fullPivLu().kernel();
            x.normalize();
        }

        template<typename matrix_t, typename vector_t>
        void solveNullspaceQR(const matrix_t& A, vector_t& x){
            auto qr = A.transpose().colPivHouseholderQr();
            matrix_t Q = qr.householderQ();
            x = Q.col(A.rows() - 1);
            x.normalize();
        }

        template<typename matrix_t, typename vector_t>
        void solveNullspaceSVD(const matrix_t& A, vector_t& x){
            x = A.jacobiSvd(Eigen::ComputeFullV).matrixV().col( A.rows() - 1 );
            x.normalize();
        }
        */
    }
}

#endif
