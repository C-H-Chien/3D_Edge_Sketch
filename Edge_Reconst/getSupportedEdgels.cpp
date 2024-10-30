#ifndef GETSUPPORTEDGELS_CPP
#define GETSUPPORTEDGELS_CPP
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

#include <stdio.h>
#include <stdlib.h>
//using namespace std;

#include "getSupportedEdgels.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace GetSupportedEdgels {
    
    get_SupportedEdgels::get_SupportedEdgels( ) { }
    
    double get_SupportedEdgels::getSupportIdx(Eigen::Vector2d edgels_tgt_reproj, Eigen::MatrixXd Tangents_VALID, Eigen::MatrixXd inliner) {
        double prev_prod = 0;
        double supported_link_indx = -2;
        double ore_threshold       = cos(double(OREN_THRESH)/180*PI);
        // std::cout << "ore_threshold: " << cos(double(OREN_THRESH)/180*PI) <<std::endl;
        for (int idx_inline = 0; idx_inline < inliner.rows(); idx_inline++){
            Eigen::Vector2d target_edges = {Tangents_VALID(inliner(idx_inline),0), Tangents_VALID(inliner(idx_inline),1)};
            double abs_dot_prod = fabs(edgels_tgt_reproj(0)*target_edges(0) + edgels_tgt_reproj(1)*target_edges(1));
            // std::cout << "abs_dot_prod: " << abs_dot_prod <<std::endl;
            if(abs_dot_prod >= ore_threshold && abs_dot_prod > prev_prod){
                //cout << "prev_prod: "<< prev_prod << endl;
                prev_prod = abs_dot_prod;
                supported_link_indx = inliner(idx_inline); 
            }
        }
        return supported_link_indx;
    }

    void get_SupportedEdgels::printAllSupportedIndices(const std::vector<Eigen::MatrixXd> &all_supported_indices) {
        LOG_INFOR_MESG("Printing all supported indices:");
        
        // Loop through all matrices in all_supported_indices
        for (size_t idx = 0; idx < all_supported_indices.size(); ++idx) {
            std::cout << "Supported Indices for frame " << idx << ":\n";
            const Eigen::MatrixXd &matrix = all_supported_indices[idx];
            
            // Loop through each row and column of the matrix
            for (int row = 0; row < matrix.rows(); ++row) {
                for (int col = 0; col < matrix.cols(); ++col) {
                    std::cout << std::fixed << std::setprecision(4) << matrix(row, col) << " ";
                }
                std::cout << std::endl;  // Move to the next line after printing all columns of a row
            }
            std::cout << std::endl;  // Space between different matrices
        }
    }

    void get_SupportedEdgels::printEdge3DToHypothesisAndSupports(
        const std::unordered_map<Eigen::Matrix<double, 3, 1>, 
        std::tuple<Eigen::Vector2d, Eigen::Vector2d, std::vector<std::pair<int, Eigen::Vector2d>>>, 
        EigenMatrixHash>& edge_3D_to_hypothesis_and_supports) 
    {
        for (const auto& [edge_3D, value] : edge_3D_to_hypothesis_and_supports) {
            // Extract 3D edge
            std::cout << "3D Edge: [" << edge_3D(0) << ", " << edge_3D(1) << ", " << edge_3D(2) << "]\n";

            // Extract the hypothesis edges from view 6 and view 8
            const Eigen::Vector2d& hypothesis_edge_view6 = std::get<0>(value);
            const Eigen::Vector2d& hypothesis_edge_view8 = std::get<1>(value);

            std::cout << "    Hypothesis Edge (View 6): [" << hypothesis_edge_view6(0) << ", " << hypothesis_edge_view6(1) << "]\n";
            std::cout << "    Hypothesis Edge (View 8): [" << hypothesis_edge_view8(0) << ", " << hypothesis_edge_view8(1) << "]\n";

            // Extract the supporting edges with validation view numbers
            const std::vector<std::pair<int, Eigen::Vector2d>>& validation_support_edges = std::get<2>(value);

            // Loop through and print each supporting edge and its corresponding validation view number
            for (size_t i = 0; i < validation_support_edges.size(); ++i) {
                const auto& [val_view_num, support_edge] = validation_support_edges[i];
                std::cout << "    Validation View " << val_view_num << ": Supporting Edge = [" 
                        << support_edge(0) << ", " << support_edge(1) << "]\n";
            }
        }
    }

}

#endif
