#ifndef PAIREDGEHYPO_HPP
#define PAIREDGEHYPO_HPP
// =============================================================================
//
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

#include <stdio.h>
#include <stdlib.h>

namespace PairEdgeHypothesis {
    
    class pair_edge_hypothesis {

    public:
        pair_edge_hypothesis();
        
        Eigen::MatrixXd getAp_Bp(Eigen::MatrixXd Edges_HYPO2, Eigen::Vector3d pt_edgel_HYPO1, Eigen::Matrix3d F );
        Eigen::MatrixXd getAp_Bp_Dist(Eigen::MatrixXd Edges_HYPO2, Eigen::Vector3d pt_edgel_HYPO1, Eigen::Matrix3d F );
        Eigen::MatrixXd getHYPO2_idx(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd numerOfDist);
        Eigen::MatrixXd getedgels_HYPO2(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd numerOfDist);
        Eigen::MatrixXd getHYPO2_idx_Ore(Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2);
        Eigen::MatrixXd getedgels_HYPO2_Ore(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2);
        Eigen::MatrixXd getHYPO2_idx_Ore_sted(Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2);
        Eigen::MatrixXd getHYPO2_idx_Ore_fixed(Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2);
        Eigen::MatrixXd getedgels_HYPO2_Ore_fixed(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2);
        Eigen::MatrixXd edgelsHYPO2correct(Eigen::MatrixXd edgels_HYPO2,  Eigen::MatrixXd edgel_HYPO1, Eigen::Matrix3d F21, Eigen::Matrix3d F12, Eigen::MatrixXd HYPO2_idx_raw);

    private:
        
    };

}


#endif
