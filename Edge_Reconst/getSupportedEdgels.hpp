#ifndef GETSUPPORTEDGELS_HPP
#define GETSUPPORTEDGELS_HPP
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

namespace GetSupportedEdgels {
    
    class get_SupportedEdgels {

    public:
        get_SupportedEdgels();
        
        double getSupportIdx(Eigen::Vector2d edgels_tgt_reproj, Eigen::MatrixXd Tangents_VALID, Eigen::MatrixXd inliner);
        void printAllSupportedIndices(const std::vector<Eigen::MatrixXd> &all_supported_indices);
        // void printEdge3DToHypothesisAndSupports(
        //     const std::unordered_map<Eigen::Matrix<double, 3, 1>, 
        //     std::tuple<Eigen::Vector2d, Eigen::Vector2d, std::vector<std::pair<int, Eigen::Vector2d>>>, 
        //     EigenMatrixHash>& edge_3D_to_hypothesis_and_supports
        // );

    private:
        
    };

}


#endif
