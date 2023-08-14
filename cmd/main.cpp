#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <assert.h>
#include <string>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

//> Include functions

#include "../Edge_Reconst/util.hpp"
#include "../Edge_Reconst/PairEdgeHypo.hpp"
#include "../Edge_Reconst/getReprojectedEdgel.hpp"
#include "../Edge_Reconst/getQuadrilateral.hpp"
#include "../Edge_Reconst/getSupportedEdgels.hpp"
#include "../Edge_Reconst/definitions.h"

using namespace std;
using namespace MultiviewGeometryUtil;

// ========================================================================================================================
// main function
//
// Modifications
//    Chiang-Heng Chien  23-07-18    Initially create a blank repository with minor multiview geometry utility functions.
//
//> (c) LEMS, Brown University
//> Yilin Zheng (yilin_zheng@brown.edu)
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =========================================================================================================================


int main(int argc, char **argv) {
  std::fstream Edge_File;
  std::fstream Rmatrix_File;
  std::fstream Tmatrix_File;
  int d = 0;
  int q = 0;
  int file_idx = 1;
  double rd_data;
  std::vector<Eigen::MatrixXd> All_Edgels;
  Eigen::MatrixXd Edgels;
  Eigen::Vector4d row_edge;
  cout << "read file now\n";
  while(file_idx<51){
    std::string Edge_File_Path = "/users/yzhen105/Edge_Based_Reconstruction/datasets/cabinet/Edges/Edge_"+std::to_string(file_idx)+".txt";
    file_idx ++;
    Edge_File.open(Edge_File_Path, std::ios_base::in);
    if (!Edge_File) { 
       std::cerr << "Edge file not existed!\n"; exit(1); 
       }
    else {
      Edgels.resize(1,4);
      while (Edge_File >> rd_data) {
        row_edge(q) = rd_data;
        q++;
        if(q>3){
          Edgels.conservativeResize(d+1,4);
          Edgels.row(d) = row_edge;
          q = 0;
          d++;
        }
      }
      Edge_File.close();
      All_Edgels.push_back(Edgels);
      Edgels = {};
      d = 0;
      q = 0;
    }
  }
  
  cout<< "Edge file loading finished" <<endl;

  std::vector<Eigen::Matrix3d> All_R;
  Eigen::Matrix3d R_matrix;
  Eigen::Vector3d row_R;
  std::string Rmatrix_File_Path = "/users/yzhen105/Edge_Based_Reconstruction/datasets/cabinet/RnT/R_matrix.txt";
  Rmatrix_File.open(Rmatrix_File_Path, std::ios_base::in);
  if (!Rmatrix_File) { 
    std::cerr << "R_matrix file not existed!\n"; exit(1); 
    }
  else {
    while (Rmatrix_File >> rd_data) {
      row_R(q) = rd_data;
      q++;
      if(q>2){
        R_matrix.row(d) = row_R;
        row_R = {};
        q = 0;
        d++;
      }
      if(d>2){
        All_R.push_back(R_matrix);
        R_matrix = {};
        d = 0;
      }
    }
    Rmatrix_File.close();
  }

  cout<< "R matrix loading finished" <<endl;

  std::vector<Eigen::Vector3d> All_T;
  Eigen::Vector3d T_matrix;
  std::string Tmatrix_File_Path = "/users/yzhen105/Edge_Based_Reconstruction/datasets/cabinet/RnT/T_matrix.txt";
  Tmatrix_File.open(Tmatrix_File_Path, std::ios_base::in);
  if (!Tmatrix_File) { 
    std::cerr << "T_matrix file not existed!\n"; exit(1); 
    }
  else {
    while (Tmatrix_File >> rd_data) {
      T_matrix(d) = rd_data;
      d++;
      if(d>2){
        All_T.push_back(T_matrix);
        T_matrix = {};
        d = 0;
      }
    }
    Tmatrix_File.close();
  }
  
  cout<< "T matrix loading finished" <<endl;

  Eigen::Matrix3d K;
  K<< 537.960322000000, 0, 319.183641000000, 0,	539.597659000000,	247.053820000000,0,	0,	1;

  ////////////////////
  // Pipeline start
  ////////////////////

  MultiviewGeometryUtil::multiview_geometry_util util;
  PairEdgeHypothesis::pair_edge_hypothesis       PairHypo;
  GetReprojectedEdgel::get_Reprojected_Edgel     getReprojEdgel;
  GetQuadrilateral::get_Quadrilateral            getQuad;
  GetSupportedEdgels::get_SupportedEdgels        getSupport;
  
  Eigen::MatrixXd Edges_HYPO1 = All_Edgels[HYPO1_VIEW_INDX];
  Eigen::Matrix3d R1          = All_R[HYPO1_VIEW_INDX];
  Eigen::Vector3d T1          = All_T[HYPO1_VIEW_INDX];
  Eigen::MatrixXd Edges_HYPO2 = All_Edgels[HYPO2_VIEW_INDX];
  Eigen::Matrix3d R2          = All_R[HYPO2_VIEW_INDX];
  Eigen::Vector3d T2          = All_T[HYPO2_VIEW_INDX];
  
  Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
  Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
  Eigen::Matrix3d Tx  = util.getSkewSymmetric(T21);
  Eigen::Matrix3d E   = util.getEssentialMatrix(R21, T21);
  Eigen::Matrix3d F   = util.getFundamentalMatrix(K.inverse(), R21, T21);


  Eigen::Vector3d pt_edgel_HYPO1;

  /////////////////////////////////////////////////
  // should be a for loop of edges in hypo 1 here
  /////////////////////////////////////////////////

  Eigen::MatrixXd paired_edge;
  int pair_num = 0;
  for(int edge_idx = 0; edge_idx < Edges_HYPO1.rows(); edge_idx++){
  pt_edgel_HYPO1 << Edges_HYPO1(edge_idx,0), Edges_HYPO1(edge_idx,1), 1;

  Eigen::MatrixXd ApBp = PairHypo.getAp_Bp(Edges_HYPO2, pt_edgel_HYPO1, F);

  Eigen::MatrixXd numerOfDist = PairHypo.getAp_Bp_Dist(Edges_HYPO2, pt_edgel_HYPO1, F);
  Eigen::MatrixXd HYPO2_idx    = PairHypo.getHYPO2_idx(Edges_HYPO2, numerOfDist);
  Eigen::MatrixXd edgels_HYPO2 = PairHypo.getedgels_HYPO2(Edges_HYPO2, numerOfDist);

  ///////////////////////////////////////////////////
  // should be a for loop of validation views here
  ///////////////////////////////////////////////////

  int VALID_idx = 0;
  int stack_idx = 0;
  Eigen::MatrixXd supported_indices;
  supported_indices.conservativeResize(edgels_HYPO2.rows(),48);
  Eigen::MatrixXd supported_indice_current;
  supported_indice_current.conservativeResize(edgels_HYPO2.rows(),1);
  Eigen::MatrixXd supported_indices_stack;
  
  for (int VALID_INDX = 0; VALID_INDX < 50; VALID_INDX++){
    if(VALID_INDX == HYPO1_VIEW_INDX || VALID_INDX == HYPO2_VIEW_INDX){
      continue;
    }
    Eigen::MatrixXd TO_Edges_VALID = All_Edgels[VALID_INDX];
    Eigen::Matrix3d R3             = All_R[VALID_INDX];
    Eigen::Vector3d T3             = All_T[VALID_INDX];
    Eigen::MatrixXd VALI_Orient    = TO_Edges_VALID.col(2);
    Eigen::MatrixXd Tangents_VALID;
    Tangents_VALID.conservativeResize(TO_Edges_VALID.rows(),2);
    Tangents_VALID.col(0)          = (VALI_Orient.array()).cos();
    Tangents_VALID.col(1)          = (VALI_Orient.array()).sin();
    
    Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
    Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);
    
    Eigen::MatrixXd pt_edge = Edges_HYPO1.row(edge_idx);
    Eigen::Vector3d tgt1_meters = getReprojEdgel.getTGT_Meters(pt_edge, K);
    
    Eigen::MatrixXd edge_pos_gamma3 = getReprojEdgel.getGamma3Pos(pt_edge, edgels_HYPO2, All_R, All_T, VALID_INDX, K);
    Eigen::MatrixXd edge_tgt_gamma3 = getReprojEdgel.getGamma3Tgt(pt_edge, edgels_HYPO2, All_R, All_T, VALID_INDX, K);
    
    // Eigen::MatrixXd QuadrilateralPoints = getQuad.getQuadrilateralPoints(pt_edge, edgels_HYPO2.row(20), All_R, All_T, VALID_INDX, K);
    for (int idx_pair = 0; idx_pair < edgels_HYPO2.rows(); idx_pair++){
      Eigen::MatrixXd inliner = getQuad.getInliner(pt_edge, edgels_HYPO2.row(idx_pair), All_R, All_T, VALID_INDX, K, TO_Edges_VALID);
      Eigen::Vector2d edgels_tgt_reproj = {edge_tgt_gamma3(idx_pair,0), edge_tgt_gamma3(idx_pair,1)};
      double supported_link_indx = getSupport.getSupportIdx(edgels_tgt_reproj, Tangents_VALID, inliner);
      supported_indice_current.row(idx_pair) << supported_link_indx;
      if (supported_link_indx != -2){
        supported_indices_stack.conservativeResize(stack_idx+1,2);
        supported_indices_stack.row(stack_idx) << double(idx_pair), double(supported_link_indx);
        stack_idx++;
      }
    }
    supported_indices.col(VALID_idx) << supported_indice_current.col(0);
    VALID_idx++;
  }
  //cout<< VALID_idx << endl;
  //cout << "supported_indices.col(0)" << endl;
  //cout << supported_indices.col(0) << endl;
  //cout << "supported_indices_stack" << endl;
  //cout << supported_indices_stack.block(0,0,50,2) << endl;
  std::vector<double> indices_stack(supported_indices_stack.data(), supported_indices_stack.data() + supported_indices_stack.rows());
  std::vector<double> indices_stack_unique = indices_stack;
  std::sort(indices_stack_unique.begin(), indices_stack_unique.end());
  std::vector<double>::iterator it1;
  it1 = std::unique(indices_stack_unique.begin(), indices_stack_unique.end());
  indices_stack_unique.resize( std::distance(indices_stack_unique.begin(),it1) );
  //cout << "supported_indices_stack" << endl;
  //cout << indices_stack_unique.size() << endl;
  //std::vector<double>::iterator it2;
  Eigen::VectorXd rep_count;
  rep_count.conservativeResize(indices_stack_unique.size(),1);
  for(int unique_idx = 0; unique_idx<indices_stack_unique.size(); unique_idx++){
    rep_count.row(unique_idx) << double(count(indices_stack.begin(), indices_stack.end(), indices_stack_unique[unique_idx]));
  }
  Eigen::VectorXd::Index   maxIndex;
  double max_support = rep_count.maxCoeff(&maxIndex);
  int numofmax = count(rep_count.data(), rep_count.data()+rep_count.size(), max_support);
  //cout << rep_count.row(maxIndex) << endl;
  if( double(max_support) < MAX_NUM_OF_SUPPORT_VIEWS){
    // cout << max_support << endl;
    continue;
  }
  int finalpair;
  if(numofmax == 1){
    finalpair = indices_stack_unique[int(maxIndex)];
    // cout << finalpair << endl;
  }else{
    //TODO: find the final pair among multiple maximum support
    std::vector<double> rep_count_vec(rep_count.data(), rep_count.data() + rep_count.rows());
    std::vector<int> max_index;
    auto start_it = begin(rep_count_vec);
    while (start_it != end(rep_count_vec)) {
      start_it = std::find(start_it, end(rep_count_vec), max_support);
      if (start_it != end(rep_count_vec)) {
        auto const pos = std::distance(begin(rep_count_vec), start_it);
        max_index.push_back(int(pos));
        ++start_it;
      }
    }
    Eigen::Vector3d coeffs;
    coeffs = F * pt_edgel_HYPO1;
    Eigen::MatrixXd Edge_Pts;
    Edge_Pts.conservativeResize(max_index.size(),2);
    for(int maxidx = 0; maxidx<max_index.size(); maxidx++){
      Edge_Pts.row(maxidx) << edgels_HYPO2(indices_stack_unique[max_index[maxidx]], 0),edgels_HYPO2(indices_stack_unique[max_index[maxidx]], 1) ;
    }
    Eigen::VectorXd Ap = coeffs(0)*Edge_Pts.col(0);
    Eigen::VectorXd Bp = coeffs(1)*Edge_Pts.col(1);
    Eigen::VectorXd numDist = Ap + Bp + Eigen::VectorXd::Ones(Ap.rows())*coeffs(2);
    double denomDist = coeffs(0)*coeffs(0) + coeffs(1)*coeffs(1);
    denomDist = sqrt(denomDist);
    Eigen::VectorXd dist = numDist.cwiseAbs()/denomDist;
    //cout << dist << endl;
    Eigen::VectorXd::Index   minIndex;
    double min_dist = dist.minCoeff(&minIndex);
    finalpair = int(indices_stack_unique[max_index[minIndex]]);
    // cout << finalpair << endl;
  }
  // linearTriangulation code already exist
  paired_edge.conservativeResize(pair_num+1,50);
  paired_edge.row(pair_num) << edge_idx, HYPO2_idx(finalpair), supported_indices.row(finalpair);
  pair_num++;
  }
  cout<<paired_edge<<endl;
}
