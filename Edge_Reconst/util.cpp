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

    Eigen::Vector3d multiview_geometry_util::get3DTangentFromTwo2Dtangents( 
        const Eigen::MatrixXd pt_edge_view1, const Eigen::MatrixXd pt_edge_view2,
        const Eigen::Matrix3d K1,  const Eigen::Matrix3d K2,
        const Eigen::Matrix3d R1,  const Eigen::Vector3d T1,
        const Eigen::Matrix3d R2,  const Eigen::Vector3d T2 )
    {
        Eigen::Matrix3d R21, R12;
        Eigen::Vector3d T21, T12;

        getRelativePoses(R1, T1, R2, T2, R21, T21, R12, T12);

        Eigen::Vector3d e1  = {1,0,0};
        Eigen::Vector3d e3  = {0,0,1};
        Eigen::Vector3d Gamma1 = K1.inverse() * Eigen::Vector3d(pt_edge_view1(0), pt_edge_view1(1), 1.0);
        Eigen::Vector3d Gamma2 = K2.inverse() * Eigen::Vector3d(pt_edge_view2(0), pt_edge_view2(1), 1.0);

        Eigen::Vector3d tgt1(cos(pt_edge_view1(2)), sin(pt_edge_view1(2)), 0.0);
        Eigen::Vector3d tgt2(cos(pt_edge_view2(2)), sin(pt_edge_view2(2)), 0.0);
        Eigen::Vector3d tgt1_meters = K1.inverse() * tgt1;
        Eigen::Vector3d tgt2_meters = K2.inverse() * tgt2;

        // double rho1 = (double(e1.transpose() * T21) - double(e3.transpose() * T21) * double(e1.transpose() *Gamma2))/(double(e3.transpose() * R21 * Gamma1)* double(e1.transpose() * Gamma2) - double(e1.transpose() * R21 * Gamma1));

        Eigen::Vector3d n1 = tgt1_meters.cross(Gamma1);
        Eigen::Vector3d n2 = R21.transpose() * tgt2_meters.cross(Gamma2);

        //> This 3D tangent is in the first camera coordinate
        Eigen::Vector3d T3D = n1.cross(n2) / (n1.cross(n2) ).norm();

        //> Convert the 3D tangent from the first camera coordinate to the world coordinate
        Eigen::Vector3d tangents_3D_world = R1.transpose() * T3D;

        return tangents_3D_world;
    }

    Eigen::Vector3d multiview_geometry_util::linearTriangulation(
        const int N,
        const std::vector<Eigen::Vector2d> pts, 
        const std::vector<Eigen::Matrix3d> & Rs,
        const std::vector<Eigen::Vector3d> & Ts,
        const Eigen::Matrix3d K )
    {
        Eigen::MatrixXd A(2*N, 4);
        Eigen::MatrixXd ATA(4, 4);
        Eigen::Vector4d GAMMA;
        Eigen::Vector3d pt3D; //> Returned variable

        //> Convert points in pixels to points in meters
        std::vector<Eigen::Vector2d> pts_meters;
        for (int p = 0; p < N; p++) {
            Eigen::Vector2d gamma;
            Eigen::Vector3d homo_p{pts[p](0), pts[p](1), 1.0};
            Eigen::Vector3d p_bar = K.inverse() * homo_p;
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

    std::vector<double> multiview_geometry_util::check_reproj_error(
        std::vector<Eigen::Vector2d> points_2D, Eigen::Vector3d point_3D, 
        std::vector<Eigen::Matrix3d> Rs, std::vector<Eigen::Vector3d> Ts, Eigen::Matrix3d K) 
    {
        std::vector<double> reproj_errs;
        for (int p = 0; p < Rs.size(); p++) {
            Eigen::Matrix3d Rot     = Rs[p];
            Eigen::Vector3d Transl  = Ts[p];
            Eigen::Vector2d obs_pt  = points_2D[p];

            Eigen::Vector3d project_pt = K * (Rot * point_3D + Transl);
            project_pt[0] /= project_pt[2];
            project_pt[1] /= project_pt[2];

            double reproj_err = std::sqrt( (project_pt[0] - obs_pt[0])*(project_pt[0] - obs_pt[0]) + (project_pt[1] - obs_pt[1])*(project_pt[1] - obs_pt[1]) );
            reproj_errs.push_back(reproj_err);
        }
        return reproj_errs;
    }
}

#endif
