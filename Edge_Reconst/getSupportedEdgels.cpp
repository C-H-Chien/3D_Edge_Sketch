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
        for (int idx_inline = 0; idx_inline < inliner.rows(); idx_inline++){
            Eigen::Vector2d target_edges = {Tangents_VALID(inliner(idx_inline),0), Tangents_VALID(inliner(idx_inline),1)};
            double abs_dot_prod = abs(edgels_tgt_reproj(0)*target_edges(0) + edgels_tgt_reproj(1)*target_edges(1));
            
            if(abs_dot_prod > OREN_THRESH && abs_dot_prod > prev_prod){
                //cout << "prev_prod: "<< prev_prod << endl;
                prev_prod = abs_dot_prod;
                supported_link_indx = inliner(idx_inline); 
            }
        }
        return supported_link_indx;
    }

}

#endif
