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

//#ifdef __cplusplus
//}
//#endif

#endif // GPU_KERNELS_HPP
