#ifndef FILE_READER_HPP
#define FILE_READER_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "definitions.h"

void readEdgelFiles(std::vector<Eigen::MatrixXd> &All_Edgels, std::fstream &Edge_File, double &rd_data, Eigen::Vector4d &row_edge, int &file_idx, int &d, int &q, int thresh_EDG);

void readHypothesisEdgelFiles(int hyp01_view_indx, int hyp02_view_indx, std::vector<Eigen::MatrixXd> &All_Edgels_H12, std::fstream &Edge_File, double &rd_data, Eigen::Vector4d &row_edge, int &H_idx, int &file_idx, int &d, int &q, int thresh_EDG);

void readRmatrix(std::vector<Eigen::Matrix3d> &All_R, Eigen::Matrix3d &R_matrix, std::fstream &Rmatrix_File, double &rd_data, Eigen::Vector3d &row_R, int &d, int &q);

void readTmatrix(std::vector<Eigen::Vector3d> &All_T, Eigen::Vector3d &T_matrix, std::fstream &Tmatrix_File, double &rd_data, int &d, int &q);

void readK(std::fstream &Kmatrix_File, std::vector<Eigen::Matrix3d> &All_K, Eigen::Matrix3d &K, Eigen::Matrix3d &K_matrix, Eigen::Vector3d &row_K, double &rd_data, int &d, int &q);

#endif // FILE_READER_HPP
