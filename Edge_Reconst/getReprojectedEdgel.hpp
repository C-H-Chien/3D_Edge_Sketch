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

namespace GetReprojectedEdgel {
    
    class get_Reprojected_Edgel {

    public:
        get_Reprojected_Edgel();
        
        Eigen::Vector3d getTGT_Meters(Eigen::MatrixXd pt_edge, Eigen::Matrix3d K);


    private:
        
    };

}


#endif
