#ifndef GETREPROJECTEDEDGEL_CPP
#define GETREPROJECTEDEDGEL_CPP
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

#include "getReprojectedEdgel.hpp"
#include "util.hpp"
#include "definitions.h"

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

namespace GetReprojectedEdgel {
    
    get_Reprojected_Edgel::get_Reprojected_Edgel( ) { }
    
    Eigen::Vector3d get_Reprojected_Edgel::getTGT_Meters(Eigen::MatrixXd pt_edge, Eigen::Matrix3d K) {
        Eigen::Vector3d pt_edgel;
        pt_edgel << pt_edge(0,0), pt_edge(0,1), 1;
        Eigen::Vector3d Gamma = K.inverse() * pt_edgel;
        Eigen::Vector2d tangent;
        tangent << cos(pt_edge(0,2)), sin(pt_edge(0,2));
        Eigen::Vector3d tgt_to_pixels;
        tgt_to_pixels << (tangent(0)+pt_edgel(0)), (tangent(1)+pt_edgel(1)), 1;
        Eigen::Vector3d tgt_to_meters = K.inverse() * tgt_to_pixels;
        Eigen::Vector3d tgt_meters = tgt_to_meters - Gamma;
        return tgt_meters;
    }

}

#endif
