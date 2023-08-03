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

  MultiviewGeometryUtil::multiview_geometry_util util;
  PairEdgeHypothesis::pair_edge_hypothesis       PairHypo;
  GetReprojectedEdgel::get_Reprojected_Edgel     getReprojEdgel;
  GetQuadrilateral::get_Quadrilateral            getQuad;
  
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
  int edge_idx = 1819;
  pt_edgel_HYPO1 << Edges_HYPO1(edge_idx-1,0), Edges_HYPO1(edge_idx-1,1), 1;

  Eigen::MatrixXd ApBp = PairHypo.getAp_Bp(Edges_HYPO2, pt_edgel_HYPO1, F);
  //cout << "function" << endl;
  //cout << ApBp.block(0,0,6,2) << endl;

  Eigen::MatrixXd numerOfDist = PairHypo.getAp_Bp_Dist(Edges_HYPO2, pt_edgel_HYPO1, F);
  //cout << "distance" << endl;
  //cout << numerOfDist.block(0,0,6,1) << endl;
  //numerOfDist = numerOfDist.cwiseAbs();
  //cout << numerOfDist.block(0,0,6,1) << endl;

  Eigen::MatrixXd HYPO2_idx    = PairHypo.getHYPO2_idx(Edges_HYPO2, numerOfDist);
  Eigen::MatrixXd edgels_HYPO2 = PairHypo.getedgels_HYPO2(Edges_HYPO2, numerOfDist);

  //cout<< "coordinates will be: "<< endl;
  //cout<< HYPO2_idx << endl;
  //cout<< edgels_HYPO2.rows() << endl;
  //cout<< edgels_HYPO2 << endl;
  cout << "meow~ ~o( =∩ω∩= )m" << endl;

  int VALID_INDX = 0;
  Eigen::MatrixXd TO_Edges_VALID = All_Edgels[VALID_INDX];
  Eigen::Matrix3d R3             = All_R[VALID_INDX];
  Eigen::Vector3d T3             = All_T[VALID_INDX];
  Eigen::MatrixXd VALI_Orient    = TO_Edges_VALID.col(2);
  Eigen::MatrixXd Tangents_VALID;
  Tangents_VALID.conservativeResize(TO_Edges_VALID.rows(),2);
  Tangents_VALID.col(0)          = (VALI_Orient.array()).cos();
  Tangents_VALID.col(1)          = (VALI_Orient.array()).sin();
  //cout << Tangents_VALID.block(0,0,10,2) << endl;

  Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
  Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);
  
  Eigen::MatrixXd pt_edge = Edges_HYPO1.row(edge_idx-1);
  /*Eigen::Vector3d pt_edgel;
  pt_edgel << pt_edge(0,0), pt_edge(0,1), 1;
  Eigen::Vector3d Gamma = K.inverse() * pt_edgel;
  Eigen::Vector2d tangent;
  tangent << cos(pt_edge(0,2)), sin(pt_edge(0,2));
  Eigen::Vector3d tgt_to_pixels;
  tgt_to_pixels << (tangent(0)+pt_edgel(0)), (tangent(1)+pt_edgel(1)), 1;
  Eigen::Vector3d tgt_to_meters = K.inverse() * tgt_to_pixels;*/
  Eigen::Vector3d tgt1_meters = getReprojEdgel.getTGT_Meters(pt_edge, K);
  cout<< "tgt1_meters: " << endl;
  cout<< tgt1_meters << endl;

  Eigen::MatrixXd edge_pos_gamma3 = getReprojEdgel.getGamma3Pos(pt_edge, edgels_HYPO2, All_R, All_T, VALID_INDX, K);
  //cout<< "edge_pos_gamma3: " << endl;
  //cout<< edge_pos_gamma3 << endl;

  Eigen::MatrixXd edge_tgt_gamma3 = getReprojEdgel.getGamma3Tgt(pt_edge, edgels_HYPO2, All_R, All_T, VALID_INDX, K);


  Eigen::Vector3d e1 = {1,0,0};
  Eigen::Vector3d e2 = {0,1,0};
  Eigen::Vector3d e3 = {0,0,1};

  Eigen::Vector3d epipole_met_view1 = {double(e1.transpose() * R21.transpose() *T21) / double(e3.transpose()*R21.transpose()*T21), double(e2.transpose()*R21.transpose()*T21) / double(e3.transpose()*R21.transpose()*T21), 1};
  Eigen::Vector3d epipole_met_view2 = {double(e1.transpose()*T21) / double(e3.transpose()*T21), double(e2.transpose()*T21) / double(e3.transpose()*T21), 1};
  Eigen::Vector3d epipole_pix_view1 = K * epipole_met_view1;
  Eigen::Vector3d epipole_pix_view2 = K * epipole_met_view2;

  /*cout<<"epipole_met_view1: "<<endl;
  cout<<epipole_met_view1<<endl;
  cout<<"epipole_met_view2: "<<endl;
  cout<<epipole_met_view2<<endl;
  cout<<"epipole_pix_view1: "<<endl;
  cout<<epipole_pix_view1<<endl;
  cout<<"epipole_pix_view2: "<<endl;
  cout<<epipole_pix_view2<<endl;*/
  
  double slope_hypo1 = double(pt_edge(1)-epipole_pix_view1(1))/double(pt_edge(0)-epipole_pix_view1(0));
  cout << "slope: "<< slope_hypo1 << endl;
  double intersect_hypo1 = double(pt_edge(1)) - slope_hypo1*double(pt_edge(0));
  cout << "interect: "<< intersect_hypo1 << endl;
  Eigen::Vector3d hypo1_pt1       = {double(pt_edge(0)), intersect_hypo1+slope_hypo1*double(pt_edge(0))+DIST_THRESH, 1};
  Eigen::Vector3d hypo1_pt2       = {double(pt_edge(0)), intersect_hypo1+slope_hypo1*double(pt_edge(0))-DIST_THRESH, 1};
  cout << "hypo1_pt1: "<< endl;
  cout << hypo1_pt1 << endl;
  cout << "hypo1_pt2: "<< endl;
  cout << hypo1_pt2 << endl;

  Eigen::MatrixXd QuadrilateralPoints = getQuad.getQuadrilateralPoints(pt_edge, edgels_HYPO2.row(13), All_R, All_T, VALID_INDX, K);
  cout<< "QuadrilateralPoints: " << endl;
  cout<< QuadrilateralPoints << endl;

}
