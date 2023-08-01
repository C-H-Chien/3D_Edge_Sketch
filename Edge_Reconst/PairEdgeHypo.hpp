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

namespace PairEdgeHypothesis {
    
    class pair_edge_hypothesis {

    public:
        pair_edge_hypothesis();
        
        Eigen::MatrixXd getAp_Bp(Eigen::MatrixXd Edges_HYPO2, Eigen::Vector3d pt_edgel_HYPO1, Eigen::Matrix3d F );
        


    private:
        
    };

}


#endif
