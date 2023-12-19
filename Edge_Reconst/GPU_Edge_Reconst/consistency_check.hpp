#ifndef CONSISTENCY_CHECK_HPP
#define CONSISTENCY_CHECK_HPP

#include <cmath>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

#include "../getReprojectedEdgel.hpp"
#include "../definitions.h"

template< typename T >
bool Check_PreProcess_Edgels_in_HYPO1( Eigen::MatrixXd Edges_HYPO1, Eigen::Matrix3d K1, \
                                       T* host_check_Edgel_H1_in_Meters, T* host_check_Edgel_H1_Normal_Vectors) 
{
    GetReprojectedEdgel::get_Reprojected_Edgel getReprojEdgel;
    Eigen::Vector3d GPU_Computed_gamma;
    Eigen::Vector3d GPU_Computed_normal_vec;
    GPU_Computed_gamma(2) = 1.0;
    int gamma1_valid_counter = 0;
    int normal_vec_valid_counter = 0;
    T avg_gamma1_err = 0.0;
    T avg_normal_vec_err = 0.0;
    T max_gamma1_err = 0.0;
    T max_normal_vec_err = 0.0;
    int max_gamma1_err_indx, max_normal_vec_err_indx;

    for(int edge_idx = 0; edge_idx < Edges_HYPO1.rows(); edge_idx++) {
        Eigen::MatrixXd pt_edge = Edges_HYPO1.row(edge_idx);
        Eigen::Vector3d pt_edgel_HYPO1;
        pt_edgel_HYPO1 << pt_edge(0,0), pt_edge(0,1), 1;
        Eigen::Vector3d Gamma1 = K1.inverse() * pt_edgel_HYPO1;
        Eigen::Vector3d tgt1_meters = getReprojEdgel.getTGT_Meters(pt_edge, K1);
        Eigen::Vector3d n1 = tgt1_meters.cross(Gamma1);

        //> Fetch data from the array capturing data transferred back from the GPU side
        GPU_Computed_gamma(0) = host_check_Edgel_H1_in_Meters[ 2*edge_idx    ];
        GPU_Computed_gamma(1) = host_check_Edgel_H1_in_Meters[ 2*edge_idx + 1];
        GPU_Computed_normal_vec(0) = host_check_Edgel_H1_Normal_Vectors[ 3*edge_idx     ];
        GPU_Computed_normal_vec(1) = host_check_Edgel_H1_Normal_Vectors[ 3*edge_idx + 1 ];
        GPU_Computed_normal_vec(2) = host_check_Edgel_H1_Normal_Vectors[ 3*edge_idx + 2 ];

        //> Compute the error
        T gamma1_err     = (Gamma1 - GPU_Computed_gamma).norm();
        T normal_vec_err = (n1 - GPU_Computed_normal_vec).norm();
        gamma1_valid_counter = ( gamma1_err <= PREPROCESS_CONSISTENCY_CPU_GPU_TOL ) ? (gamma1_valid_counter + 1) : (gamma1_valid_counter);
        normal_vec_valid_counter = ( normal_vec_err <= PREPROCESS_CONSISTENCY_CPU_GPU_TOL ) ? (normal_vec_valid_counter + 1) : (normal_vec_valid_counter);
        avg_gamma1_err += gamma1_err;
        avg_normal_vec_err += normal_vec_err;

        if (gamma1_err > max_gamma1_err) {
            max_gamma1_err = gamma1_err;
            max_gamma1_err_indx = edge_idx;
        }
        if (normal_vec_err > max_normal_vec_err) {
            max_normal_vec_err = normal_vec_err;
            max_normal_vec_err_indx = edge_idx;
        }
    }
    if (gamma1_valid_counter != Edges_HYPO1.rows())     std::cout << "FAILURE TO PASS CONCSISTENCY CHECK: Preprocessing stage for computing gamma1 in GPU!" << std::endl;
    if (normal_vec_valid_counter != Edges_HYPO1.rows()) std::cout << "FAILURE TO PASS CONCSISTENCY CHECK: Preprocessing stage for computing normal vectors in GPU!" << std::endl;

#if DEBUG_CONSISTENCY_CHECK
    std::cout << "Average gamma1 consistency error: " << avg_gamma1_err / (T)Edges_HYPO1.rows() << std::endl;
    std::cout << "Average normal vec consistency error: " << avg_normal_vec_err / (T)Edges_HYPO1.rows() << std::endl;
    std::cout << std::endl;
    std::cout << max_gamma1_err << ", " << max_gamma1_err_indx << std::endl;
    std::cout << max_normal_vec_err << ", " << max_normal_vec_err_indx << std::endl;
#endif

    return ((gamma1_valid_counter == Edges_HYPO1.rows()) && (normal_vec_valid_counter == Edges_HYPO1.rows())) ? (true) : (false);
}

bool Check_Hypothesis_Edgels_Pairs( int* host_Hypothesis_Edgel_Pair_Index ) 
{
    for (int i = 100; i < 110; i++) {
        std::cout << host_Hypothesis_Edgel_Pair_Index[i] << std::endl;
    }
}

void Write_Final_Edge_Pair_Results( int Num_Of_Edgles_HYPO1, int* host_Hypothesis_Edgel_Pair_Index ) 
{
  std::ofstream GPU_Result_File;
  std::string Output_File_Path = OUTPUT_WRITE_FOLDER + "GPU_Final_Result.txt";
  GPU_Result_File.open(Output_File_Path);
  if (!GPU_Result_File) {
    std::cerr << "Unable to open GPU_Final_Result file!\n";
  }
  else {
    for (int i = 0; i < Num_Of_Edgles_HYPO1; i++)
        GPU_Result_File << i << "\t" << host_Hypothesis_Edgel_Pair_Index[i] << "\n";
  }
  GPU_Result_File.close();
}

/*
template< typename T, typename EigenT >
bool Check_Converted_Matrix_Consistency( T *Converted_Matrix, std::vector<EigenT> ThirdParty_Matrix ) {

    int rows = ThirdParty_Matrix[0].rows();
    int cols = ThirdParty_Matrix[0].cols();
    int consistency_count = 0;
    for (int k = 0; k < ThirdParty_Matrix.size(); k++) {
        for (int i = 0; i < rows; i++) {
            for (int  j = 0; j < cols; j++) {
                consistency_count = ( fabs(ThirdParty_Matrix[k](i,j) - Converted_Matrix(i,j,k)) < CONSISTENCY_TOL ) ? \
                                    (consistency_count + 1) : (consistency_count);
            }
        }
    }
    int total_consistency_count = rows * cols * ThirdParty_Matrix.size();
    return ( consistency_count == total_consistency_count ) ? (true) : (false);
}

template< typename T, typename EigenT >
bool Check_Converted_Vector_Consistency( T *Converted_Vector, std::vector<EigenT> ThirdParty_Vector ) {

    int rows = ThirdParty_Vector[0].rows();
    int consistency_count = 0;
    for (int k = 0; k < ThirdParty_Vector.size(); k++) {
        for (int i = 0; i < rows; i++) {
            consistency_count = ( fabs(ThirdParty_Vector[k](i) - Converted_Vector[(i) + (k) * 3]) < CONSISTENCY_TOL ) ? \
                                (consistency_count + 1) : (consistency_count);
        }
    }

    int total_consistency_count = rows * ThirdParty_Vector.size();
    return ( consistency_count == total_consistency_count ) ? (true) : (false);
}*/

#endif    // CONSISTENCY_CHECK_HPP