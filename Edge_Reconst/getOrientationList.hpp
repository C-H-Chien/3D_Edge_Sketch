#ifndef GETOREINTATIONLIST_HPP
#define GETOREINTATIONLIST_HPP
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

namespace GetOrientationList {
    
    class get_OrientationList {

    public:
        get_OrientationList( double, int, int );
        
        Eigen::MatrixXd getOreListBar(Eigen::MatrixXd Edges_HYPO1, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2, int VALID_INDX, int REFIDX);
        Eigen::MatrixXd getOreListBarVali(Eigen::MatrixXd Edges_HYPO1, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2, int VALID_INDX, int REFIDX);
        Eigen::MatrixXd getOreList(int hyp01_view_indx, int hyp02_view_indx, Eigen::MatrixXd Edges_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2);
        std::pair<Eigen::MatrixXd, Eigen::Vector2d> getOreListVali(Eigen::MatrixXd Edges_HYPO2, std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, Eigen::Matrix3d K1, Eigen::Matrix3d K2, int VALID_INDX, int REFIDX);
        

    private:
        //> private variables will recieve values from the EdgeSketch_Core class
        int dataset_img_rows;
        int dataset_img_cols;
        double delta;
    };

}


#endif
