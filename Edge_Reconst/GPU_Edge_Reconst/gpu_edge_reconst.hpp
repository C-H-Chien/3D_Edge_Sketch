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

#include "./indices.hpp"
#include "../definitions.h"

//> CUDA error check
#define cudacheck( a )  do { \
                            cudaError_t e = a; \
                            if(e != cudaSuccess) { \
                                printf("\033[1;31m"); \
                                printf("Error in %s:%d %s\n", __func__, __LINE__, cudaGetErrorString(e)); \
                                printf("\033[0m"); \
                            }\
                        } while(0)

template< typename T >
class EdgeReconstGPU {
    int device_id;
    int Num_Of_Edgels_H1;
    
    T *All_Rot,    *dev_All_Rot;
    T *All_Transl, *dev_All_Transl;
    T *All_Calib,  *dev_All_Calib;

    T *Edgels_H1,      *dev_Edgels_H1;
    T *OreList_Deg,    *dev_OreList_Deg;
    T *OreListBar_Deg, *dev_OreListBar_Deg;

  public:
    //> Setup the timer
    float time_ER;
    cudaEvent_t start, stop;

    EdgeReconstGPU(int, std::vector<Eigen::Matrix3d>, std::vector<Eigen::Vector3d>, std::vector<Eigen::Matrix3d>, \
                   Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd);
    ~EdgeReconstGPU();

    //> Member functions
    void GPU_Edge_Reconstruction_Main();

    //void read_array_from_file(std::string filename, T *rd_data, int first_dim, int second_dim);
    //void write_array_to_file(std::string filename, T *wr_data, int first_dim, int second_dim);
};

template< typename T >
EdgeReconstGPU<T>::EdgeReconstGPU(int device, \
                                  std::vector<Eigen::Matrix3d> All_R, std::vector<Eigen::Vector3d> All_T, std::vector<Eigen::Matrix3d> All_K, \
                                  Eigen::MatrixXd Edges_HYPO1, Eigen::MatrixXd OreListBardegree, Eigen::MatrixXd OreListdegree):
    device_id(device) {
    
    cudaDeviceProp prop;
    cudacheck( cudaGetDeviceProperties(&prop, device_id));
    std::cout << "## GPU Device : " << prop.name << std::endl;
    cudacheck( cudaSetDevice(device_id) );

    //> Create CUDA events
    cudacheck( cudaEventCreate(&start) );
    cudacheck( cudaEventCreate(&stop) );

    //> Time holder
    time_ER = 0.0;

    //> Allocate CPU memory
    All_Rot        = new T[ 3 * 3 * All_R.size() ];
    All_Transl     = new T[ 3 * All_T.size()     ];
    All_Calib      = new T[ 3 * 3 * All_K.size() ];

    Edgels_H1      = new T[ Edges_HYPO1.rows() * 3  ];
    OreList_Deg    = new T[ OreListdegree.rows()    ];
    OreListBar_Deg = new T[ OreListBardegree.rows() ];

    //> Need to assign and organize values from Eigen arrays to standard C++ arrays
    //> 1) Camera Intrinsic and Extrinsic Matrices
    for (int k = 0; k < All_R.size(); k++) {
        for (int i = 0; i < 3; i++) {
            All_Transl[(i) + (k) * 3] = (T)All_T[k](i);
            for (int  j = 0; j < 3; j++) {
                All_Rot(i,j,k)   = (T)All_R[k](i,j);
                All_Calib(i,j,k) = (T)All_K[k](i,j);
            }
        }
    }

#if DEBUG_GPU
    bool is_consistent_Rot    = Check_Converted_Matrix_Consistency<T, Eigen::Matrix3d>( All_Rot, All_R    );
    bool is_consistent_Calib  = Check_Converted_Matrix_Consistency<T, Eigen::Matrix3d>( All_Calib, All_K  );
    bool is_consistent_Transl = Check_Converted_Vector_Consistency<T, Eigen::Vector3d>( All_Transl, All_T );
    if (is_consistent_Rot)    std::cout << "Rotation Matrix conversion is successful!" << std::endl;
    if (is_consistent_Calib)  std::cout << "Calibration Matrix conversion is successful!" << std::endl;
    if (is_consistent_Transl) std::cout << "Translation Vector conversion is successful!" << std::endl;
#endif

    //> 2) Edgels_H1, OreList_Deg, and OreListBar_Deg
    int edge_count_in_scope = 0;
    for(int edge_idx = 0; edge_idx < Edges_HYPO1.rows(); edge_idx++) {
        if(Edges_HYPO1(edge_idx,0) < 10 || Edges_HYPO1(edge_idx,0) > imgcols-10 || Edges_HYPO1(edge_idx,1) < 10 || Edges_HYPO1(edge_idx,1) > imgrows-10){
            continue;
        }
        else {
            Edgels_H1(edge_count_in_scope, 0) = Edges_HYPO1(edge_idx, 0);       //> edgel x
            Edgels_H1(edge_count_in_scope, 1) = Edges_HYPO1(edge_idx, 1);       //> edgel y
            Edgels_H1(edge_count_in_scope, 2) = Edges_HYPO1(edge_idx, 2);       //> edgel \theta

            OreList_Deg[ edge_count_in_scope ]    = OreListdegree(edge_idx, 0);
            OreListBar_Deg[ edge_count_in_scope ] = OreListBardegree(edge_idx, 0);
            edge_count_in_scope++;
        }
    }

    Num_Of_Edgels_H1 = edge_count_in_scope;
    std::cout << "Number of Edgels in HYPO1 (Originally from TOED):         " << Edges_HYPO1.rows() << std::endl;
    std::cout << "Number of Edgels in HYPO1 (within the valid image scope): " << Num_Of_Edgels_H1 << std::endl;

    //> Allocate GPU memory
    cudacheck( cudaMalloc((void**)&dev_All_Rot,    (3 * 3 * All_R.size()) * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_All_Transl, (3 * All_T.size())     * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_All_Calib,  (3 * 3 * All_K.size()) * sizeof(T)) );

    cudacheck( cudaMalloc((void**)&dev_Edgels_H1,       (Num_Of_Edgels_H1 * 3) * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_OreList_Deg,     (Num_Of_Edgels_H1)     * sizeof(T)) );
    cudacheck( cudaMalloc((void**)&dev_OreListBar_Deg,  (Num_Of_Edgels_H1)     * sizeof(T)) );
}

template< typename T >
void EdgeReconstGPU<T>::GPU_Edge_Reconstruction_Main() {
    //> Transfer memory from CPU to GPU
    
}

// ===================================== Destructor =======================================
template< typename T >
EdgeReconstGPU<T>::~EdgeReconstGPU () {

    //> Free CPU Memory
    delete [] All_Rot;
    delete [] All_Transl;
    delete [] All_Calib;

    delete [] Edgels_H1;
    delete [] OreList_Deg;
    delete [] OreListBar_Deg;

    //> Free GPU Memory
    cudacheck( cudaFree(dev_All_Rot) );
    cudacheck( cudaFree(dev_All_Transl) );
    cudacheck( cudaFree(dev_All_Calib) );

    cudacheck( cudaFree(dev_Edgels_H1) );
    cudacheck( cudaFree(dev_OreList_Deg) );
    cudacheck( cudaFree(dev_OreListBar_Deg) );

    //> Destroy CUDA Events
    cudacheck( cudaEventDestroy(start) );
    cudacheck( cudaEventDestroy(stop) );
    
}

#endif    // GPU_EDGE_RECONST_HPP