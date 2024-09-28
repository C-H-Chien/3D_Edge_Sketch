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

#include <stdio.h>
#include <stdlib.h>

//using namespace std;

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
            double dist = numerOfDist(idx_HYPO2,0);
            if(dist < DIST_THRESH){
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
            double dist = numerOfDist(idx_HYPO2,0);
            if(dist < DIST_THRESH){
                edgels_HYPO2.conservativeResize(idx_hypopair+1,4);
                edgels_HYPO2.row(idx_hypopair) = Edges_HYPO2.row(idx_HYPO2);
                idx_hypopair++;
            }
        }
        return edgels_HYPO2;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getHYPO2_idx_Ore(Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2) {
        int idx_hypopair = 0;
        Eigen::MatrixXd HYPO2_idx;
        std::vector<double> Ore_List1Bar(OreListdegree.data(), OreListdegree.data() + OreListdegree.rows());
        auto it = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_1, thresh_ore21_2](double i){return i > thresh_ore21_1 && i <thresh_ore21_2;});
        while (it != std::end(Ore_List1Bar)) {
            HYPO2_idx.conservativeResize(idx_hypopair+1,1);
            HYPO2_idx.row(idx_hypopair) << double(std::distance(std::begin(Ore_List1Bar), it));
            idx_hypopair++;
            it = std::find_if(std::next(it), std::end(Ore_List1Bar), [thresh_ore21_1, thresh_ore21_2](double i){return i > thresh_ore21_1 && i <thresh_ore21_2;});
        }
        return HYPO2_idx;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getedgels_HYPO2_Ore(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2) {
        int idx_hypopair = 0;
        Eigen::MatrixXd edgels_HYPO2;
        std::vector<double> Ore_List1Bar(OreListdegree.data(), OreListdegree.data() + OreListdegree.rows());
        auto it = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_1, thresh_ore21_2](double i){return i > thresh_ore21_1 && i <thresh_ore21_2;});
        while (it != std::end(Ore_List1Bar)) {
            edgels_HYPO2.conservativeResize(idx_hypopair+1,4);
            //edgels_HYPO2.row(idx_hypopair) << double(std::distance(std::begin(Ore_List1Bar), it));
            edgels_HYPO2.row(idx_hypopair) = Edges_HYPO2.row(std::distance(std::begin(Ore_List1Bar), it));
            idx_hypopair++;
            it = std::find_if(std::next(it), std::end(Ore_List1Bar), [thresh_ore21_1, thresh_ore21_2](double i){return i > thresh_ore21_1 && i <thresh_ore21_2;});
        }
        return edgels_HYPO2;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getHYPO2_idx_Ore_sted(Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2) {
        Eigen::MatrixXd HYPO2_idx;
        HYPO2_idx.conservativeResize(2,1);
        std::vector<double> Ore_List1Bar(OreListdegree.data(), OreListdegree.data() + OreListdegree.rows());
        auto itst = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_1](double i){return i > thresh_ore21_1;});
        HYPO2_idx.row(0) << double(std::distance(std::begin(Ore_List1Bar), itst));
        auto ited = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_2](double i){return i > thresh_ore21_2;});
        HYPO2_idx.row(1) << double(std::distance(std::begin(Ore_List1Bar), ited)-1);
        
        return HYPO2_idx;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getHYPO2_idx_Ore_fixed(Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2) {
        Eigen::MatrixXd HYPO2_idx;
        HYPO2_idx.conservativeResize(20,1);
        std::vector<double> Ore_List1Bar(OreListdegree.data(), OreListdegree.data() + OreListdegree.rows());
        auto itst = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_1](double i){return i > thresh_ore21_1;});
        auto ited = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_2](double i){return i > thresh_ore21_2;});
        int idx_start = std::distance(std::begin(Ore_List1Bar), itst);
        int idx_end   = std::distance(std::begin(Ore_List1Bar), ited)-1;
        int midpoint  = std::round((idx_start+idx_end)/2);
        for(int idx = 0; idx < 20; idx++){
            HYPO2_idx.row(idx) << double(midpoint-11+idx);
        }
        return HYPO2_idx;
    }

    Eigen::MatrixXd pair_edge_hypothesis::getedgels_HYPO2_Ore_fixed(Eigen::MatrixXd Edges_HYPO2, Eigen::MatrixXd OreListdegree, double thresh_ore21_1, double thresh_ore21_2) {
        Eigen::MatrixXd edgels_HYPO2;
        edgels_HYPO2.conservativeResize(20,4);
        std::vector<double> Ore_List1Bar(OreListdegree.data(), OreListdegree.data() + OreListdegree.rows());
        auto itst = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_1](double i){return i > thresh_ore21_1;});
        auto ited = std::find_if(std::begin(Ore_List1Bar), std::end(Ore_List1Bar), [thresh_ore21_2](double i){return i > thresh_ore21_2;});
        int idx_start = std::distance(std::begin(Ore_List1Bar), itst);
        int idx_end   = std::distance(std::begin(Ore_List1Bar), ited)-1;
        int midpoint  = std::round((idx_start+idx_end)/2);
        for(int idx = 0; idx < 20; idx++){
            edgels_HYPO2.row(idx) = Edges_HYPO2.row(midpoint-11+idx);
        }        
        return edgels_HYPO2;
    }

    Eigen::MatrixXd pair_edge_hypothesis::edgelsHYPO2correct(Eigen::MatrixXd edgels_HYPO2,  Eigen::MatrixXd edgel_HYPO1, Eigen::Matrix3d F21, Eigen::Matrix3d F12, Eigen::MatrixXd HYPO2_idx_raw){
        Eigen::MatrixXd edgels_HYPO2_corrected;
        Eigen::MatrixXd xy1_H1;
        xy1_H1.conservativeResize(1,3);
        xy1_H1(0,0) = edgel_HYPO1(0,0);
        xy1_H1(0,1) = edgel_HYPO1(0,1);
        xy1_H1(0,2) = 1;
        Eigen::MatrixXd coeffspt1T = F21 * xy1_H1.transpose();
        Eigen::MatrixXd coeffspt1  = coeffspt1T.transpose();
        Eigen::MatrixXd Apixel_1   = coeffspt1.col(0);
        Eigen::MatrixXd Bpixel_1   = coeffspt1.col(1);
        Eigen::MatrixXd Cpixel_1   = coeffspt1.col(2);
        double a1_line  = -Apixel_1(0,0)/Bpixel_1(0,0);
        double b1_line  = -1;
        double c1_line  = -Cpixel_1(0,0)/Bpixel_1(0,0);
        double a_edgeH1    = tan(edgel_HYPO1(0,2));
        double b_edgeH1    = -1;
        double c_edgeH1    = -(a_edgeH1*edgel_HYPO1(0,0)-edgel_HYPO1(0,1));
        // std::cout << "a1_line: " << a1_line <<std::endl;
        // std::cout << "c1_line: " << c1_line <<std::endl;
        // std::cout << "a_edgeH1: " << a_edgeH1 <<std::endl;
        // std::cout << "c_edgeH1: " << c_edgeH1 <<std::endl;
        double idx_correct = 0;
        for(int idx_hypo2 = 0; idx_hypo2 < edgels_HYPO2.rows(); idx_hypo2++){
            double a_edgeH2 = tan(edgels_HYPO2(idx_hypo2,2));
            double b_edgeH2 = -1;
            double c_edgeH2 = -(a_edgeH2*edgels_HYPO2(idx_hypo2,0)-edgels_HYPO2(idx_hypo2,1));
            double x_currH2 = ((b1_line*c_edgeH2-b_edgeH2*c1_line)/(a1_line*b_edgeH2-a_edgeH2*b1_line) + edgels_HYPO2(idx_hypo2,0))/2;
            double y_currH2 = ((c1_line*a_edgeH2-c_edgeH2*a1_line)/(a1_line*b_edgeH2-a_edgeH2*b1_line) + edgels_HYPO2(idx_hypo2,1))/2;
            double dist2    = sqrt((x_currH2 - edgels_HYPO2(idx_hypo2,0))*(x_currH2 - edgels_HYPO2(idx_hypo2,0))+(y_currH2 - edgels_HYPO2(idx_hypo2,1))*(y_currH2 - edgels_HYPO2(idx_hypo2,1)));

            Eigen::MatrixXd xy1_H2;
            xy1_H2.conservativeResize(1,3);
            xy1_H2(0,0) = x_currH2;
            xy1_H2(0,1) = y_currH2;
            xy1_H2(0,2) = 1;
            Eigen::MatrixXd coeffspt2T = F12 * xy1_H2.transpose();
            Eigen::MatrixXd coeffspt2  = coeffspt2T.transpose();
            Eigen::MatrixXd Apixel_2   = coeffspt2.col(0);
            Eigen::MatrixXd Bpixel_2   = coeffspt2.col(1);
            Eigen::MatrixXd Cpixel_2   = coeffspt2.col(2);
            double a2_line  = -Apixel_2(0,0)/Bpixel_2(0,0);
            double b2_line  = -1;
            double c2_line  = -Cpixel_2(0,0)/Bpixel_2(0,0);
            double x_currH1 = (b2_line*c_edgeH1-b_edgeH1*c2_line)/(a2_line*b_edgeH1-a_edgeH1*b2_line);
            double y_currH1 = (c2_line*a_edgeH1-c_edgeH1*a2_line)/(a2_line*b_edgeH1-a_edgeH1*b2_line);
            double dist1    = sqrt((x_currH1 - edgel_HYPO1(0,0))*(x_currH1 - edgel_HYPO1(0,0))+(y_currH1 - edgel_HYPO1(0,1))*(y_currH1 - edgel_HYPO1(0,1)));
            if(dist1 < circleR && dist2 < circleR){
                edgels_HYPO2_corrected.conservativeResize(idx_correct+1,10);
                edgels_HYPO2_corrected.row(idx_correct) << x_currH1, y_currH1, edgel_HYPO1(0,2), edgel_HYPO1(0,3), x_currH2, y_currH2,  edgels_HYPO2(idx_hypo2,2),  edgels_HYPO2(idx_hypo2,3), HYPO2_idx_raw(idx_hypo2), idx_hypo2;
                idx_correct +=1;
            }
            //std::cout << "dist1: " << dist1 << ", dist2: " << dist2 << ", circleR: " << circleR << std::endl;

        }
        return edgels_HYPO2_corrected;
    }


}

#endif
