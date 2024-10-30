#ifndef EDGESKETCH_CORE_HPP
#define EDGESKETCH_CORE_HPP
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
#include <vector>
#include <stdio.h>
#include <stdlib.h>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

//> YAML file data reader
#include <yaml-cpp/yaml.h>

//> shared class pointers
#include <file_reader.hpp>

    
class EdgeSketch_Core {

public:
    EdgeSketch_Core( YAML::Node );

    double edge_sketch_time;

    void Read_Sketch_Data();
    void Run_3D_Edge_Sketch();
    void Finalize_Edge_Pairs();
    
    bool Skip_this_Edge( const int edge_idx ) {
      //> Edge Boundary Check: not too close to boundary
      if ( Edges_HYPO1(edge_idx,0) < 10 || Edges_HYPO1(edge_idx,0) > Img_Cols-10 || Edges_HYPO1(edge_idx,1) < 10 || Edges_HYPO1(edge_idx,1) > Img_Rows-10)
        return true;
      
      //> Paired Edge Check: not yet been paired
      if ( paired_edge(edge_idx,0) != -2 )
        return true;
      return false;
    }
    
    
private:
    std::shared_ptr<file_reader> Load_Data = nullptr;
    std::shared_ptr<MultiviewGeometryUtil::multiview_geometry_util> util = nullptr;
    YAML::Node Edge_Sketch_Setting_YAML_File;

    //> Input 3D Edge Sketch Settings
    int Num_Of_OMP_Threads;
    int hyp01_view_indx;
    int hyp02_view_indx;
    double Edge_Loc_Pertubation;
    double Orien_Thresh;
    int Max_Num_Of_Support_Views;
    int Edge_Detection_Init_Thresh;
    int Edge_Detection_Final_Thresh;
    double Parallel_Epipolar_Line_Angle_Deg;
    double Reproj_Dist_Thresh;
    int circleR; //> Unknown setting

    //> Input Dataset Settings
    std::string Dataset_Path;
    std::string Dataset_Name;
    std::string Scene_Name;
    int Num_Of_Total_Imgs;
    int Img_Rows;
    int Img_Cols;
    bool Use_Multiple_K;
    double fx;
    double fy;
    double cx;
    double cy;
    std::string Delta_FileName_Str;

    //> Edges and camera intrinsix/extrinsic matrices of all images
    std::vector<Eigen::Matrix3d> All_R;
    std::vector<Eigen::Vector3d> All_T;
    std::vector<Eigen::Matrix3d> All_K;
    std::vector<Eigen::MatrixXd> All_Edgels; 
    Eigen::Matrix3d K;

    //> Edges and camera intrinsix/extrinsic matrices of the two hypothesis views
    Eigen::Matrix3d K_HYPO1;
    Eigen::Matrix3d K_HYPO2;
    Eigen::MatrixXd Edges_HYPO1;
    Eigen::MatrixXd Edges_HYPO2;
    Eigen::Matrix3d Rot_HYPO1;
    Eigen::Matrix3d Rot_HYPO2;
    Eigen::Vector3d Transl_HYPO1;
    Eigen::Vector3d Transl_HYPO2;
    //> Relative poses and fundmental matrices
    Eigen::Matrix3d R21;
    Eigen::Vector3d T21;
    Eigen::Matrix3d F21;
    Eigen::Matrix3d R12;
    Eigen::Vector3d T12;  
    Eigen::Matrix3d F12;

    Eigen::MatrixXd paired_edge;
    Eigen::MatrixXd OreListdegree;
    Eigen::MatrixXd OreListBardegree;
};




#endif
