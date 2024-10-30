#ifndef FILE_READER_HPP
#define FILE_READER_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include "definitions.h"

class file_reader{

public:
    file_reader( std::string, std::string, std::string, int );
    void read_All_Edgels( std::vector<Eigen::MatrixXd> &All_Edgels, int thresh_EDG );
    Eigen::MatrixXd read_Edgels_Of_a_File( int file_idx, double thresh_EDG );
    void readRmatrix( std::vector<Eigen::Matrix3d> &All_R );
    void readTmatrix( std::vector<Eigen::Vector3d> &All_T );
    void readK( std::vector<Eigen::Matrix3d> &All_K );

private:
    std::string dataset_name_sequence_path;
    int Dataset_Total_Num_Of_Images;

    //> Paths for the files
    std::string Edge_File_Path_First_Half;
    std::string Rmatrix_File_Path;
    std::string Tmatrix_File_Path;
    std::string Kmatrix_File_Path;
};

#endif // FILE_READER_HPP
