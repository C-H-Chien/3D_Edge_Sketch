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
#include <ctime>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

//> YAML file data reader
#include <yaml-cpp/yaml.h>

//> shared class pointers
#include "file_reader.hpp"
#include "util.hpp"
#include "PairEdgeHypo.hpp"
#include "getReprojectedEdgel.hpp"
#include "getSupportedEdgels.hpp"
#include "getOrientationList.hpp"
#include "edge_mapping.hpp"
    
class EdgeSketch_Core {

public:
    std::shared_ptr<EdgeMapping> edgeMapping = nullptr;
    std::vector<Eigen::MatrixXd> paired_edge_final_all;

    //> Constructor
    EdgeSketch_Core( YAML::Node );
    void Read_Camera_Data();
    void Read_Edgels_Data();
    void Set_Hypothesis_Views_Camera();
    void Set_Hypothesis_Views_Edgels();
    void Run_3D_Edge_Sketch();
    void Finalize_Edge_Pairs_and_Reconstruct_3D_Edges();
    void Clear_Data();
    void Stack_3D_Edges();
    void Project_3D_Edges_and_Find_Next_Hypothesis_Views();

    //> Destructor
    ~EdgeSketch_Core();
    
    bool Skip_this_Edge( const int edge_idx ) {
      //> Edge Boundary Check: not too close to boundary
      if ( Edges_HYPO1(edge_idx,0) < 10 || Edges_HYPO1(edge_idx,0) > Img_Cols-10 || Edges_HYPO1(edge_idx,1) < 10 || Edges_HYPO1(edge_idx,1) > Img_Rows-10)
        return true;
      
      //> Paired Edge Check: not yet been paired
      if ( paired_edge(edge_idx,0) != -2 )
        return true;
      return false;
    }

    bool is_Epipolar_Wedges_in_Parallel(double thresh_ore31_1, double thresh_ore31_2, double thresh_ore32_1, double thresh_ore32_2, int idx_pair, Eigen::VectorXd &isparallel, Eigen::MatrixXd &supported_indice_current) {
      Eigen::MatrixXd anglediff(4,1);
      anglediff << fabs(thresh_ore31_1 - thresh_ore32_1), fabs(thresh_ore31_1 - thresh_ore32_2), \
                   fabs(thresh_ore31_2 - thresh_ore32_1), fabs(thresh_ore31_2 - thresh_ore32_2);
      if (anglediff.maxCoeff() <= Parallel_Epipolar_Line_Angle_Deg) {
          isparallel.row(idx_pair) << 0;
          supported_indice_current.row(idx_pair) << -2;
          return true;
      }
      else {
        return false;
      }
    }

    Eigen::MatrixXd paired_edge_final;
    double edge_sketch_time;
    int thresh_EDG;
    int hyp01_view_indx;
    int hyp02_view_indx;
    int Edge_Detection_Init_Thresh;
    int Edge_Detection_Final_Thresh;
    int Max_3D_Edge_Sketch_Passes;

    std::vector< Eigen::MatrixXd > all_supported_indices;
    Eigen::MatrixXd Gamma1s;
    Eigen::MatrixXd all_3D_Edges;
    std::vector< int > claimedEdgesList;
    double least_ratio;
    bool enable_aborting_3D_edge_sketch;
    int num_of_nonveridical_edge_pairs;

    //> timer
    double itime, pair_edges_time;
    double finalize_edge_pair_time;
    double find_next_hypothesis_view_time;

    std::vector<int> history_hypothesis_views_index;
    
private:
    //> sharing the classes
    std::shared_ptr<file_reader> Load_Data = nullptr;
    std::shared_ptr<MultiviewGeometryUtil::multiview_geometry_util> util = nullptr;
    std::shared_ptr<PairEdgeHypothesis::pair_edge_hypothesis> PairHypo = nullptr;
    std::shared_ptr<GetReprojectedEdgel::get_Reprojected_Edgel> getReprojEdgel = nullptr;
    std::shared_ptr<GetSupportedEdgels::get_SupportedEdgels> getSupport = nullptr;
    std::shared_ptr<GetOrientationList::get_OrientationList> getOre = nullptr;
    //std::shared_ptr<EdgeMapping> edgeMapping = nullptr;

    Eigen::MatrixXd project3DEdgesToView(const Eigen::MatrixXd& edges3D, const Eigen::Matrix3d& R, const Eigen::Vector3d& T, const Eigen::Matrix3d& K, const Eigen::Matrix3d& R_hyp01, const Eigen::Vector3d& T_hpy01);
    int claim_Projected_Edges(const Eigen::MatrixXd& projectedEdges, const Eigen::MatrixXd& observedEdges, double threshold);
    void select_Next_Best_Hypothesis_Views( 
      const std::vector< int >& claimedEdges, std::vector<Eigen::MatrixXd> All_Edgels,
      std::pair<int, int> &next_hypothesis_views, double &least_ratio, std::vector<int> history_hypothesis_views_index );

    //> YAML file data parser
    YAML::Node Edge_Sketch_Setting_YAML_File;

    //> Input 3D Edge Sketch Settings
    int Num_Of_OMP_Threads;
    double Edge_Loc_Pertubation;
    double Orien_Thresh;
    int Max_Num_Of_Support_Views;
    double Parallel_Epipolar_Line_Angle_Deg;
    double Reproj_Dist_Thresh;
    double Stop_3D_Edge_Sketch_by_Ratio_Of_Claimed_Edges;
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
