#ifndef GETREPROJECTEDEDGEL_HPP
#define GETREPROJECTEDEDGEL_HPP
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
namespace GetReprojectedEdgel {
    
    class get_Reprojected_Edgel {

    public:
        get_Reprojected_Edgel();
        
        Eigen::Vector3d getTGT_Meters(Eigen::MatrixXd pt_edge, Eigen::Matrix3d K);
        Eigen::MatrixXd getGamma3Pos(Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd edgels_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2, Eigen::Matrix3d K3);
        Eigen::MatrixXd getGamma3Tgt(Eigen::MatrixXd pt_edge_HYPO1, Eigen::MatrixXd edgels_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, int VALID_INDX, Eigen::Matrix3d K1, Eigen::Matrix3d K2);


    private:
        
    };

}


#endif
