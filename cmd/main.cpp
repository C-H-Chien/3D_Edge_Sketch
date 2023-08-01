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

#include "/users/yzhen105/Edge_Based_Reconstruction/Edge_Reconst/util.hpp"
#include "/users/yzhen105/Edge_Based_Reconstruction/Edge_Reconst/definitions.h"

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
  //double testEdges[];
  //std::vector< Eigen::Vector4d > Edgels;
  //Eigen::Vector4d test;
  std::vector<vector<vector<double>>> All_Edgels;
  std::vector<vector<double>> Edgels;
  std::vector<double> row_edge;
  cout << "read file now\n";
  while(file_idx<51){
    std::string Edge_File_Path = "/users/yzhen105/Edge_Based_Reconstruction/datasets/cabinet/Edges/Edge_"+std::to_string(file_idx)+".txt";
    //cout << "file idx: " << file_idx << '\n';
    //cout << "file path: " << Edge_File_Path << '\n';
    file_idx ++;
    Edge_File.open(Edge_File_Path, std::ios_base::in);
    if (!Edge_File) { 
       std::cerr << "Edge file not existed!\n"; exit(1); 
       }
    else {
      while (Edge_File >> rd_data) {
        //row_edge(q) = rd_data;
        row_edge.push_back(rd_data);
        q++;
        if(q>3){
          //cout << "row_edge size: " << row_edge.size() << '\n';
          Edgels.push_back(row_edge);
          //Edgels.insert (Edgels.end(),row_edge);
          //cout <<  Edgels[d] << '\n';
          //cout << "Edgels size: " << Edgels.size() << '\n';
          row_edge = {};
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
        //cout << "test size: " << test.size() << '\n';
        R_matrix.row(d) = row_R;
        //Edgels.insert (Edgels.end(),test);
        //cout <<  Edgels[d] << '\n';
        //cout << "Edgels size: " << Edgels.size() << '\n';
        row_R = {};
        q = 0;
        d++;
      }
      if(d>2){
        //cout << "test size: " << test.size() << '\n';
        All_R.push_back(R_matrix);
        //Edgels.insert (Edgels.end(),test);
        //cout <<  Edgels[d] << '\n';
        //cout << "Edgels size: " << Edgels.size() << '\n';
        R_matrix = {};
        d = 0;
      }
    }
    Rmatrix_File.close();
  }
  

  
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
        //cout << "test size: " << test.size() << '\n';
        All_T.push_back(T_matrix);
        //Edgels.insert (Edgels.end(),test);
        //cout <<  Edgels[d] << '\n';
        //cout << "Edgels size: " << Edgels.size() << '\n';
        T_matrix = {};
        d = 0;
      }
    }
    Tmatrix_File.close();
  }
  
  /*cout << All_T[0] <<endl;
  cout <<endl;
  cout << All_T[20] <<endl;
  cout <<endl;
  cout << All_T[49] <<endl;*/

  std::vector<vector<double>> K_matrix = {{537.960322000000, 0, 319.183641000000},{0,	539.597659000000,	247.053820000000},{0,	0,	1}};
  Eigen::Matrix3d K;
  K<< 537.960322000000, 0, 319.183641000000, 0,	539.597659000000,	247.053820000000,0,	0,	1;

  /*cout << K <<endl;
  cout << K.inverse() <<endl;*/

  //Eigen::Vector3d T21_for_test;
  //T21_for_test<< -0.513199164027762,-0.172032399874067,0.176688656780153;
  MultiviewGeometryUtil::multiview_geometry_util util;
  Eigen::Matrix3d R1  = All_R[5];
  Eigen::Matrix3d R2  = All_R[2];
  Eigen::Vector3d T1  = All_T[5];
  Eigen::Vector3d T2  = All_T[2];
  Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
  Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
  Eigen::Matrix3d Tx  = util.getSkewSymmetric(T21);
  Eigen::Matrix3d E   = util.getEssentialMatrix(R21, T21);
  Eigen::Matrix3d F   = util.getFundamentalMatrix(K.inverse(), R21, T21);
  
  cout << "R21 will be:" <<endl;
  cout << R21 <<endl;

  cout << "T21 will be:" <<endl;
  cout << T21 <<endl;

  cout << "Tx will be:" <<endl;
  cout << Tx <<endl;

  cout << "E will be:" <<endl;
  cout << E <<endl;

  cout << "F will be:" <<endl;
  cout << F <<endl;
  
}
