#ifndef GPU_EDGE_RECONST_HPP
#define GPU_EDGE_RECONST_HPP

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

#include "./consistency_check.hpp"
#include "./gpu_kernels.hpp"

#include "./indices.hpp"
#include "../definitions.h"

template< typename T >
class EdgeReconstGPU {
    int device_id;
    int Num_Of_Edgels_in_HYPO1;
    int Num_Of_Total_Edgels;
    
    T *Rel_R21,     *dev_Rel_R21;
    T *Rel_T21,     *dev_Rel_T21;
    T *Fund21,      *dev_Fund21;
    T *Rel_Rots,    *dev_Rel_Rots;
    T *Rel_Transls, *dev_Rel_Transls;
    T *Funds,       *dev_Funds;
    T *All_Calib,   *dev_All_Calib;

    T *host_check_Edgel_H1_in_Meters;
    T *host_check_Edgel_H1_Normal_Vectors;

    T   *All_Edgels_List,                  *dev_All_Edgels_List;
    T   *Wedge_Angle_Range_H2_to_VALID,    *dev_Wedge_Angle_Range_H2_to_VALID;
    int *Edgel_List_Start_Indx,            *dev_Edgel_List_Start_Indx;
    int *Truncated_Wedge_Start_Indx_HYPO2, *dev_Truncated_Wedge_Start_Indx_HYPO2;
    int *Truncated_Wedge_Start_Indx_VALID, *dev_Truncated_Wedge_Start_Indx_VALID;
    int *Valid_Views_Indices,              *dev_Valid_Views_Indices;

    T *dev_Edgel_H1_in_Meters;
    T *dev_Edgel_H1_Normal_Vectors;

  public:
    //> Setup the timer
    float time_Process_EdgelH1, time_ER;
    cudaEvent_t start, stop;

    //> Constructor
    EdgeReconstGPU(int, int, Eigen::Matrix3d, Eigen::Vector3d, Eigen::Matrix3d, \
                   std::vector<Eigen::Matrix3d>, std::vector<Eigen::Vector3d>, \
                   std::vector<Eigen::Matrix3d>, std::vector<Eigen::Matrix3d>, \
                   std::vector<Eigen::MatrixXd>, std::vector<int>, std::vector<int>, \
                   std::vector<T> );
    //> Destructor
    ~EdgeReconstGPU();

    //> Member functions
    void GPU_PreProcess_Edgels_HYPO1( Eigen::MatrixXd Edges_HYPO1, Eigen::Matrix3d K1, std::vector<Eigen::Matrix3d> All_K );
    void GPU_Edge_Reconstruction_Main();

    //void read_array_from_file(std::string filename, T *rd_data, int first_dim, int second_dim);
    //void write_array_to_file(std::string filename, T *wr_data, int first_dim, int second_dim);
};

template< typename T >
EdgeReconstGPU<T>::EdgeReconstGPU(int device, int Num_Of_Total_Edgels_for_All_Views, \
                                  Eigen::Matrix3d R21, Eigen::Vector3d T21, Eigen::Matrix3d F21, \
                                  std::vector<Eigen::Matrix3d> R31s, std::vector<Eigen::Vector3d> T31s, std::vector<Eigen::Matrix3d> F31s, std::vector<Eigen::Matrix3d> All_K, \
                                  std::vector<Eigen::MatrixXd> All_Edgels, std::vector<int> idx_truncated_start_HYPO2, std::vector<int> idx_truncated_start_VALID, \
                                  std::vector<T> Wedge_Angle_Range_from_HYPO2_to_VALID ):
    device_id(device), Num_Of_Total_Edgels(Num_Of_Total_Edgels_for_All_Views) {
    
    //> Check whether the use of number of thread-blocks is veridical
    Num_Of_Edgels_in_HYPO1 = All_Edgels[HYPO1_VIEW_INDX].rows();
    assert(NUM_OF_THREADBLOCKS == Num_Of_Edgels_in_HYPO1);

    cudaDeviceProp prop;
    cudacheck( cudaGetDeviceProperties(&prop, device_id));
    std::cout << std::endl;
    std::cout << "## GPU Device : " << prop.name << std::endl;
    cudacheck( cudaSetDevice(device_id) );

    //> Create CUDA events
    cudacheck( cudaEventCreate(&start) );
    cudacheck( cudaEventCreate(&stop) );

    //> Time holder
    time_Process_EdgelH1 = 0.0;
    time_ER              = 0.0;

    //> Allocate CPU memory
    Rel_R21     = new T[ 9 ];
    Rel_T21     = new T[ 3 ];
    Fund21      = new T[ 9 ];
    Rel_Rots    = new T[ 9 * R31s.size()  ];
    Funds       = new T[ 9 * R31s.size()  ];
    Rel_Transls = new T[ 3 * T31s.size()  ];
    All_Calib   = new T[ 2 * All_K.size() ];         //> Only {cx, cy} is necessary

    //> Some CPU memory capturing GPU memory data for consistency check and debugging purposes
    host_check_Edgel_H1_in_Meters      = new T[ 2 * Num_Of_Total_Edgels ];
    host_check_Edgel_H1_Normal_Vectors = new T[ 3 * Num_Of_Total_Edgels ];

    Edgel_List_Start_Indx            = new int[ DATASET_NUM_OF_FRAMES ];
    All_Edgels_List                  = new T[ 3 * Num_Of_Total_Edgels ];
    Truncated_Wedge_Start_Indx_HYPO2 = new int[ Num_Of_Edgels_in_HYPO1 ];
    Truncated_Wedge_Start_Indx_VALID = new int[ (DATASET_NUM_OF_FRAMES-2) * Num_Of_Edgels_in_HYPO1 ];
    Wedge_Angle_Range_H2_to_VALID    = new T[ DATASET_NUM_OF_FRAMES-2 ];
    Valid_Views_Indices              = new int[DATASET_NUM_OF_FRAMES-2];

#if DEBUG_GPU
    std::cout << "Finished allocating CPU memory!" << std::endl;
#endif

    //> 1) Relative Poses
    for (int i = 0; i < 3; i++) {
        Rel_T21[i] = T21[i];
        for (int j = 0; j < 3; j++) {
            Rel_R21[(i)*3 + (j)] = R21(i,j);
            Fund21[(i)*3 + (j)]  = F21(i,j);
        }
    }
    for (int k = 0; k < R31s.size(); k++) {
        for (int i = 0; i < 3; i++) {
            Rel_Transls[(i) + (k) * 3] = T31s[k](i);
            for (int j = 0; j < 3; j++) {
                Rel_Rots[(i)*3 + (j) + (k)*9] = R31s[k](i,j);
                Funds[(i)*3 + (j) + (k)*9]    = F31s[k](i,j);
            }
        }
    }
    for (int k = 0; k < All_K.size(); k++) {
        All_Calib[0 + k*(2)] = All_K[k](0,2);   //> cx
        All_Calib[1 + k*(2)] = All_K[k](1,2);   //> cy
    }

/*#if DEBUG_GPU
    bool is_consistent_Rot    = Check_Converted_Matrix_Consistency<T, Eigen::Matrix3d>( All_Rot, All_R    );
    bool is_consistent_Transl = Check_Converted_Vector_Consistency<T, Eigen::Vector3d>( All_Transl, All_T );
    if (is_consistent_Rot)    std::cout << "Rotation Matrix conversion is successful!" << std::endl;
    if (is_consistent_Transl) std::cout << "Translation Vector conversion is successful!" << std::endl;
    for (int k = 0; k < 2; k++) {
        for (int i = 0; i < 4; i++) std::cout << All_Calib[i + k*(4)] << "\t";
        std::cout << std::endl;
    }
#endif*/

    //> 2) Construct a list of number of edgels for each image by stacking the start index for each view in a list
    int cumulative_start_index = 0;
    for (int k = 0; k < DATASET_NUM_OF_FRAMES; k++) {
        Edgel_List_Start_Indx[k] = cumulative_start_index;
        cumulative_start_index += All_Edgels[k].rows();
    }

    //> 3) Stack all edgels in a single list
    for (int k = 0; k < DATASET_NUM_OF_FRAMES; k++) {
        for (int i = 0; i < All_Edgels[k].rows(); i++) {
            for (int j = 0; j < 3; j++) {
                if (k == 0) All_Edgels_List[(i) * 3 + (j)] = All_Edgels[k](i,j);
                else        All_Edgels_List[(i) * 3 + (j) + (Edgel_List_Start_Indx[k])*3] = All_Edgels[k](i,j);
            }
        }
    }

    //> 4) Truncated epipolar wedge indices for each gamma1 on HYPO2 and all the rest of the validation views
    //> 4-1) Truncated epipolar wedge in H2
    for (int i = 0; i < Num_Of_Edgels_in_HYPO1; i++) Truncated_Wedge_Start_Indx_HYPO2[i] = idx_truncated_start_VALID[i];
    //> 4-2) Truncated epipolar wedge in VALID
    for (int k = 0; k < (DATASET_NUM_OF_FRAMES-2); k++) {
        for (int i = 0; i < Num_Of_Edgels_in_HYPO1; i++) {
            Truncated_Wedge_Start_Indx_VALID[(i) + (k)*Num_Of_Edgels_in_HYPO1] = idx_truncated_start_VALID[(i) + (k)*Num_Of_Edgels_in_HYPO1];
        }
    }

    //> 5) Epipolar wedge angle range for all validation views from H2
    for (int i = 0; i < DATASET_NUM_OF_FRAMES-2; i++) Wedge_Angle_Range_H2_to_VALID[i] = Wedge_Angle_Range_from_HYPO2_to_VALID[i];

    //> 6) Get a list of indices for validation views
    int valid_view_count = 0;
    for (int i = 0; i < (DATASET_NUM_OF_FRAMES); i++) {
        if (i == (HYPO1_VIEW_INDX) || i == (HYPO2_VIEW_INDX)) continue;
        else {
            Valid_Views_Indices[valid_view_count] = i;
            valid_view_count++;
        }
    }

    //> Allocate GPU memory
    cudacheck( cudaMalloc((void**)&dev_Rel_R21,     (9) * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_Fund21,      (9) * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_Rel_T21,     (3) * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_Rel_Rots,    (9) * R31s.size() * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_Funds,       (9) * T31s.size() * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_Rel_Transls, (3) * T31s.size() * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_All_Calib,   (2) * All_K.size() * sizeof(T)) );

    cudacheck( cudaMalloc((void**)&dev_All_Edgels_List,                  (3 * Num_Of_Total_Edgels)                       *sizeof(T))   );
    cudacheck( cudaMalloc((void**)&dev_Wedge_Angle_Range_H2_to_VALID,    (DATASET_NUM_OF_FRAMES-2)                       *sizeof(T))   );
    cudacheck( cudaMalloc((void**)&dev_Edgel_List_Start_Indx,            (DATASET_NUM_OF_FRAMES)                         *sizeof(int)) );
    cudacheck( cudaMalloc((void**)&dev_Valid_Views_Indices,              (DATASET_NUM_OF_FRAMES-2)                       *sizeof(int)) );
    cudacheck( cudaMalloc((void**)&dev_Truncated_Wedge_Start_Indx_HYPO2, (Num_Of_Total_Edgels)                           *sizeof(int)) );
    cudacheck( cudaMalloc((void**)&dev_Truncated_Wedge_Start_Indx_VALID, (DATASET_NUM_OF_FRAMES-2)*Num_Of_Edgels_in_HYPO1*sizeof(int)) );

    cudacheck( cudaMalloc((void**)&dev_Edgel_H1_in_Meters,      (2 * Num_Of_Edgels_in_HYPO1)*sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_Edgel_H1_Normal_Vectors, (3 * Num_Of_Edgels_in_HYPO1)*sizeof(T)) );
#if DEBUG_GPU
    std::cout << "Finished allocating GPU memory!" << std::endl;
#endif
}

template< typename T >
void EdgeReconstGPU<T>::GPU_PreProcess_Edgels_HYPO1( 
    Eigen::MatrixXd Edges_HYPO1, Eigen::Matrix3d K1, std::vector<Eigen::Matrix3d> All_K ) 
{
    //> Transfer memory from CPU to GPU, for Edgel in HYPO1 information
    cudacheck( cudaMemcpy(dev_All_Calib,  All_Calib, (2 * All_K.size())*sizeof(T), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_All_Edgels_List,       All_Edgels_List,       (3 * Num_Of_Total_Edgels)*sizeof(T),   cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Edgel_List_Start_Indx, Edgel_List_Start_Indx, (DATASET_NUM_OF_FRAMES)  *sizeof(int), cudaMemcpyHostToDevice) );

    //> Start the CUDA event timer
    cudacheck( cudaEventRecord(start) );

    //> Create and Launch the GPU kernel
    gpu_Process_Edgels_HYPO1( device_id, dev_All_Calib, dev_All_Edgels_List, dev_Edgel_List_Start_Indx, \
                              dev_Edgel_H1_in_Meters, dev_Edgel_H1_Normal_Vectors );

    //> End the CUDA event timer
    cudacheck( cudaEventRecord(stop) );
	cudacheck( cudaEventSynchronize(stop) );
	cudacheck( cudaEventElapsedTime(&time_Process_EdgelH1, start, stop) );
    printf(" ## Cost of Processing Edgels on HYPO1 in GPU = %8.4f (ms)\n", time_Process_EdgelH1 );

#if CHECK_PREPROCESS_CONSISTENCY_CPU_GPU
    //> Copy data from GPU to CPU
    cudacheck( cudaMemcpy(host_check_Edgel_H1_in_Meters,      dev_Edgel_H1_in_Meters,      (2*Num_Of_Edgels_in_HYPO1)*sizeof(T), cudaMemcpyDeviceToHost) );
    cudacheck( cudaMemcpy(host_check_Edgel_H1_Normal_Vectors, dev_Edgel_H1_Normal_Vectors, (3*Num_Of_Edgels_in_HYPO1)*sizeof(T), cudaMemcpyDeviceToHost) );
    bool pass_consistency_check = Check_PreProcess_Edgels_in_HYPO1( Edges_HYPO1, K1, host_check_Edgel_H1_in_Meters, host_check_Edgel_H1_Normal_Vectors);
    if (pass_consistency_check) std::cout << " (GPU Preprocessing stage has passed the computation consistency check!) " << std::endl;
#endif
}
/*
template< typename T >
void EdgeReconstGPU<T>::GPU_Edge_Reconstruction_Main() {
    //> Transfer memory from CPU to GPU
    cudacheck( cudaMemcpy(dev_Rel_R21,     Rel_R21,     (9)*sizeof(T), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Fund21,      Fund21,      (9)*sizeof(T), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Rel_T21,     Rel_T21,     (3)*sizeof(T), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Rel_Rots,    Rel_Rots,    (9*(DATASET_NUM_OF_FRAMES-2))*sizeof(T), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Funds,       Funds,       (9*(DATASET_NUM_OF_FRAMES-2))*sizeof(T), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Rel_Transls, Rel_Transls, (3*(DATASET_NUM_OF_FRAMES-2))*sizeof(T), cudaMemcpyHostToDevice) );

    cudacheck( cudaMemcpy(dev_Wedge_Angle_Range_H2_to_VALID,    Wedge_Angle_Range_H2_to_VALID,    (DATASET_NUM_OF_FRAMES-2)*sizeof(T),   cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Valid_Views_Indices,              Valid_Views_Indices,              (DATASET_NUM_OF_FRAMES-2)*sizeof(int), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Truncated_Wedge_Start_Indx_HYPO2, Truncated_Wedge_Start_Indx_HYPO2, (Num_Of_Total_Edgels)    *sizeof(int), cudaMemcpyHostToDevice) );
    cudacheck( cudaMemcpy(dev_Truncated_Wedge_Start_Indx_VALID, Truncated_Wedge_Start_Indx_VALID, (DATASET_NUM_OF_FRAMES-2)*Num_Of_Edgels_in_HYPO1*sizeof(int), cudaMemcpyHostToDevice) );
    
    //> Start the CUDA event timer
    cudacheck( cudaEventRecord(start) );

    //> Create and Launch the GPU kernel
    gpu_Edge_Reconstruction_template<T>(device_id, dev_Rel_R21, dev_Rel_T21, dev_Fund21, dev_All_Calib, \
                                        dev_Rel_Rots, dev_Rel_Transls, dev_Funds, \
                                        dev_All_Edgels_List, dev_Edgel_List_Start_Indx, \
                                        dev_Truncated_Wedge_Start_Indx_HYPO2, dev_Truncated_Wedge_Start_Indx_VALID, \
                                        dev_Valid_Views_Indices, dev_Wedge_Angle_Range_H2_to_VALID \
                                        );

    //> End the CUDA event timer
    cudacheck( cudaEventRecord(stop) );
	cudacheck( cudaEventSynchronize(stop) );
	cudacheck( cudaEventElapsedTime(&time_ER, start, stop) );
    printf(" ## GPU Edge Reconstruction time = %8.4f (ms)\n", time_ER );
}*/

// ===================================== Destructor =======================================
template< typename T >
EdgeReconstGPU<T>::~EdgeReconstGPU () {

    //> Free CPU Memory
    delete [] Rel_R21;
    delete [] Rel_T21;
    delete [] Rel_Rots;
    delete [] Rel_Transls;
    delete [] All_Calib;

    delete [] Fund21;
    delete [] Funds;

    delete [] Edgel_List_Start_Indx;
    delete [] All_Edgels_List;
    delete [] Valid_Views_Indices;
    delete [] Truncated_Wedge_Start_Indx_HYPO2;
    delete [] Truncated_Wedge_Start_Indx_VALID;
    delete [] Wedge_Angle_Range_H2_to_VALID;

    delete [] host_check_Edgel_H1_in_Meters;
    delete [] host_check_Edgel_H1_Normal_Vectors;

    //> Free GPU Memory
    cudacheck( cudaFree(dev_Rel_R21) );
    cudacheck( cudaFree(dev_Rel_T21) );
    cudacheck( cudaFree(dev_Rel_Rots) );
    cudacheck( cudaFree(dev_Rel_Transls) );
    cudacheck( cudaFree(dev_All_Calib) );

    cudacheck( cudaFree(dev_Fund21) );
    cudacheck( cudaFree(dev_Funds) );

    cudacheck( cudaFree(dev_Edgel_List_Start_Indx) );
    cudacheck( cudaFree(dev_All_Edgels_List) );
    cudacheck( cudaFree(dev_Valid_Views_Indices) );
    cudacheck( cudaFree(dev_Truncated_Wedge_Start_Indx_HYPO2) );
    cudacheck( cudaFree(dev_Truncated_Wedge_Start_Indx_VALID) );
    cudacheck( cudaFree(dev_Wedge_Angle_Range_H2_to_VALID) );

    cudacheck( cudaFree(dev_Edgel_H1_in_Meters) );
    cudacheck( cudaFree(dev_Edgel_H1_Normal_Vectors) );

    //> Destroy CUDA Events
    cudacheck( cudaEventDestroy(start) );
    cudacheck( cudaEventDestroy(stop) );
}

#endif    // GPU_EDGE_RECONST_HPP
