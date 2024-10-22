#ifndef GETQUADRILATERAL_HPP
#define GETQUADRILATERAL_HPP
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

namespace GetQuadrilateral {
    
    class get_Quadrilateral {

    public:
        get_Quadrilateral();
        
        Eigen::MatrixXd getQuadrilateralPoints(int hyp01_view_indx, int hyp02_view_indx, Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd edgels_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2, Eigen::Matrix3d K3);
        Eigen::MatrixXd getInliner(int hyp01_view_indx, int hyp02_view_indx, Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd pt_edge_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2, Eigen::Matrix3d K3, Eigen::MatrixXd TO_Edges_VALID);
        

    private:
        
    };

}


#endif
