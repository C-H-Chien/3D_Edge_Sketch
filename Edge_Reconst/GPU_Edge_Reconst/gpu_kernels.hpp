#ifndef GPU_KERNELS_HPP
#define GPU_KERNELS_HPP

#include <stdio.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

//#ifdef __cplusplus
//extern "C" {
//#endif

//> Single Precision
void gpu_Process_Edgels_HYPO1(
        int device_id, float* dev_All_Calib, 
        float* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx,
        float* dev_Edgel_H1_in_Meters, float* dev_Edgel_H1_Normal_Vectors 
);

//> Double Precision
void gpu_Process_Edgels_HYPO1(
        int device_id, double* dev_All_Calib, 
        double* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx,
        double* dev_Edgel_H1_in_Meters, double* dev_Edgel_H1_Normal_Vectors 
);

//> Single Precision
void gpu_Edge_Reconstruction(
        int device_id,
        float* dev_Rel_R21, float* dev_Rel_T21, float* dev_All_Calib,
        float* dev_Rel_Rots, float* dev_Rel_Transls, 
        float* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx, 
        int *dev_Truncated_Wedge_Start_Indx_HYPO2, 
        int *dev_Truncated_Wedge_Start_Indx_VALID,
        int* dev_Valid_Views_Indices, float* dev_Wedge_Angle_Range_H2_to_VALID,
        float *dev_Edgel_H1_in_Meters, float* dev_Edgel_H1_Normal_Vectors,
        int* dev_Hypothesis_Edgel_Pair_Index
);

//> Double Precision
void gpu_Edge_Reconstruction(
        int device_id,
        double* dev_Rel_R21, double* dev_Rel_T21, double* dev_All_Calib,
        double* dev_Rel_Rots, double* dev_Rel_Transls, 
        double* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx, 
        int *dev_Truncated_Wedge_Start_Indx_HYPO2, 
        int *dev_Truncated_Wedge_Start_Indx_VALID,
        int* dev_Valid_Views_Indices, double* dev_Wedge_Angle_Range_H2_to_VALID,
        double *dev_Edgel_H1_in_Meters, double* dev_Edgel_H1_Normal_Vectors,
        int* dev_Hypothesis_Edgel_Pair_Index
);

//#ifdef __cplusplus
//}
//#endif

#endif // GPU_KERNELS_HPP
