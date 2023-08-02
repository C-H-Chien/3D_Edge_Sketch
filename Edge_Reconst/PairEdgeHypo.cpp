#ifndef PAIREDGEHYPO_CPP
#define PAIREDGEHYPO_CPP
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

#include "PairEdgeHypo.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace PairEdgeHypothesis {
    
    pair_edge_hypothesis::pair_edge_hypothesis( ) { }

    Eigen::MatrixXd pair_edge_hypothesis::getAp_Bp(Eigen::MatrixXd Edges_HYPO2, Eigen::Vector3d pt_edgel_HYPO1, Eigen::Matrix3d F ) {
        Eigen::Vector3d coeffs;
        coeffs = F * pt_edgel_HYPO1;
        Eigen::MatrixXd Ap_Bp;
        Ap_Bp.conservativeResize(Edges_HYPO2.rows(),2);
        Ap_Bp.col(0) = coeffs(0) * Edges_HYPO2.col(0);
        Ap_Bp.col(1) = coeffs(1) * Edges_HYPO2.col(1);
        return Ap_Bp;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getAp_Bp_Dist(Eigen::MatrixXd Edges_HYPO2, Eigen::Vector3d pt_edgel_HYPO1, Eigen::Matrix3d F ) {
        Eigen::Vector3d coeffs;
        coeffs = F * pt_edgel_HYPO1;
        Eigen::MatrixXd Ap_Bp;
        Ap_Bp.conservativeResize(Edges_HYPO2.rows(),2);
        Ap_Bp.col(0) = coeffs(0) * Edges_HYPO2.col(0);
        Ap_Bp.col(1) = coeffs(1) * Edges_HYPO2.col(1);
        Eigen::MatrixXd numerOfDist = Ap_Bp.col(0) + Ap_Bp.col(1) + Eigen::VectorXd::Ones(Edges_HYPO2.rows())*coeffs(2);
        Eigen::MatrixXd denomOfDist = Eigen::VectorXd::Ones(Edges_HYPO2.rows())*(coeffs(0)*coeffs(0)+coeffs(1)*coeffs(1));
        denomOfDist = denomOfDist.array().sqrt();
        return numerOfDist.cwiseAbs()/denomOfDist(0);
    }

    Eigen::MatrixXd pair_edge_hypothesis::getHYPO2_idx(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd numerOfDist) {
        int idx_hypopair = 0;
        Eigen::MatrixXd HYPO2_idx;
        for(int idx_HYPO2 = 0; idx_HYPO2 < numerOfDist.rows(); idx_HYPO2++){
            double distance = numerOfDist(idx_HYPO2,0);
            if(distance < DIST_THRESH){
                HYPO2_idx.conservativeResize(idx_hypopair+1,1);
                HYPO2_idx.row(idx_hypopair) << double(idx_HYPO2);
                idx_hypopair++;
            }
        }
        return HYPO2_idx;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getedgels_HYPO2(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd numerOfDist) {
        int idx_hypopair = 0;
        Eigen::MatrixXd edgels_HYPO2;
        for(int idx_HYPO2 = 0; idx_HYPO2 < numerOfDist.rows(); idx_HYPO2++){
            double distance = numerOfDist(idx_HYPO2,0);
            if(distance < DIST_THRESH){
                edgels_HYPO2.conservativeResize(idx_hypopair+1,4);
                edgels_HYPO2.row(idx_hypopair) = Edges_HYPO2.row(idx_HYPO2);
                idx_hypopair++;
            }
        }
        return edgels_HYPO2;
    }

}

#endif
