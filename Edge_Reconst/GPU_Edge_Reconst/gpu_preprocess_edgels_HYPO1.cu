#ifndef GPU_PREPROCESS_EDGLES_HYPO1_HPP
#define GPU_PREPROCESS_EDGLES_HYPO1_HPP

#include <iostream>
#include <math.h>
#include "./indices.hpp"
#include "../definitions.h"
#include "./gpu_kernels.hpp"

template< typename T, int Num_Of_Threads_Per_Block >
__global__
void
gpu_Process_Edgels_HYPO1_Kernel( T* dev_All_Calib, T* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx,
                                 T* dev_Edgel_H1_in_Meters, T* dev_Edgel_H1_Normal_Vectors ) 
{
    const int tid = threadIdx.x;
    const int bid  = blockIdx.x;

    //> Shared memory ptrs
    //extern __shared__ double sdata[];

    int offset_indx_HYPO1;
    T r_gamma1_pixels[3] = {0.0};
    T r_gamma1_meters[2] = {0.0};
    T r_tangent[2] = {0.0};
    T r_N_vec[3] = {0.0};
        
    //> Read calibration matrix, {cx, cy} of the hypothesis views
    //> TODO: maybe use GPU constant memory to broadcast {cx, cy} to all threads
    T r_CX = dev_All_Calib[ HYPO1_VIEW_INDX*(2)     ];
    T r_CY = dev_All_Calib[ HYPO1_VIEW_INDX*(2) + 1 ];

    //> Read gamma1 on HYPO1 image from the global memory, according to the thread-block ID
    //  Store in register file so that every thread holds this data.
    offset_indx_HYPO1 = dev_Edgel_List_Start_Indx[HYPO1_VIEW_INDX] * 3;

    //> Get Edgel on HYPO1 from global memory to register file
    r_gamma1_pixels[0] = dev_All_Edgels_List[ offset_indx_HYPO1 + (bid)*(Num_Of_Threads_Per_Block)*3 + (tid)*3     ];
    r_gamma1_pixels[1] = dev_All_Edgels_List[ offset_indx_HYPO1 + (bid)*(Num_Of_Threads_Per_Block)*3 + (tid)*3 + 1 ];
    r_gamma1_pixels[2] = dev_All_Edgels_List[ offset_indx_HYPO1 + (bid)*(Num_Of_Threads_Per_Block)*3 + (tid)*3 + 2 ];

    //> Convert from pixels to meters
    r_gamma1_meters[0] = (r_gamma1_pixels[0] - r_CX) / (CALIB_FX);
    r_gamma1_meters[1] = (r_gamma1_pixels[1] - r_CY) / (CALIB_FY);

    //> Compute the tangent in pixels and convert it to meters
    r_tangent[0] = r_gamma1_pixels[0] + cos(r_gamma1_pixels[2]);
    r_tangent[1] = r_gamma1_pixels[1] + sin(r_gamma1_pixels[2]);

    r_tangent[0] = (r_tangent[0] - r_CX) / (CALIB_FX);
    r_tangent[1] = (r_tangent[1] - r_CY) / (CALIB_FY);

    //> tangent in meters
    r_tangent[0] -= r_gamma1_meters[0];
    r_tangent[1] -= r_gamma1_meters[1];

    //> Compute the normal vector, i.e., the cross product of tangent and gamma1, all in meters (t1 x gamma1)
    r_N_vec[0] = (r_tangent[1]);
    r_N_vec[1] = -(r_tangent[0]);
    r_N_vec[2] = r_tangent[0]*r_gamma1_meters[1] - r_tangent[1]*r_gamma1_meters[0];

    //> Write back to the global memory
    //> (i) gamma1 in meters
    dev_Edgel_H1_in_Meters[ (bid)*(Num_Of_Threads_Per_Block)*2 + (tid)*2    ] = r_gamma1_meters[0];
    dev_Edgel_H1_in_Meters[ (bid)*(Num_Of_Threads_Per_Block)*2 + (tid)*2 + 1] = r_gamma1_meters[1];
    //> (ii) Normal vectors
    dev_Edgel_H1_Normal_Vectors[ (bid)*(Num_Of_Threads_Per_Block)*3 + (tid)*3    ] = r_N_vec[0];
    dev_Edgel_H1_Normal_Vectors[ (bid)*(Num_Of_Threads_Per_Block)*3 + (tid)*3 + 1] = r_N_vec[1];
    dev_Edgel_H1_Normal_Vectors[ (bid)*(Num_Of_Threads_Per_Block)*3 + (tid)*3 + 2] = r_N_vec[2];
}

//> Template for either SP32 or DP64
template< typename T >
void gpu_Process_Edgels_HYPO1_template(
        int device_id,
        T* dev_All_Calib, T* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx,
        T* dev_Edgel_H1_in_Meters, T* dev_Edgel_H1_Normal_Vectors )
{
    const int Num_Of_Threads_Per_Block = NUM_OF_WARPS_PER_BLOCK_PREPROCESS*WARP_SIZE;

    //> Kernel parameters
    dim3 grid(NUM_OF_THREADBLOCKS_PREPROCESS, 1, 1);
    dim3 threads(Num_Of_Threads_Per_Block, 1, 1);

    int shmem = 0;
    //shmem += sizeof(T) * (1);   //> cx
    //shmem += sizeof(T) * (1);   //> cy

    //> Get max. dynamic shared memory on the GPU
    int nthreads_max, shmem_max = 0;
    cudacheck( cudaDeviceGetAttribute(&nthreads_max, cudaDevAttrMaxThreadsPerBlock, device_id) );
#if CUDA_VERSION >= 9000
    cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, device_id) );
    if (shmem <= shmem_max) {
        cudacheck( cudaFuncSetAttribute(gpu_Process_Edgels_HYPO1_Kernel<T, Num_Of_Threads_Per_Block>, \
                                        cudaFuncAttributeMaxDynamicSharedMemorySize, shmem) );
    }
#else
    cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlock, device_id) );
#endif

    if ( shmem > shmem_max ) printf("error: kernel %s requires too many threads or too much shared memory\n", __func__);

    //> GPU Kernel Argument
    void *kernel_args[] = {&dev_All_Calib, \
                           &dev_All_Edgels_List,    &dev_Edgel_List_Start_Indx, \
                           &dev_Edgel_H1_in_Meters, &dev_Edgel_H1_Normal_Vectors};

    //> Launch the GPU kernel
    cudacheck( cudaLaunchKernel((void*)gpu_Process_Edgels_HYPO1_Kernel<T, Num_Of_Threads_Per_Block>, grid, threads, kernel_args, shmem, NULL) );
}

//> Single Precision (SP32)
void gpu_Process_Edgels_HYPO1(
        int device_id, float* dev_All_Calib, 
        float* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx,
        float* dev_Edgel_H1_in_Meters, float* dev_Edgel_H1_Normal_Vectors )
{
    gpu_Process_Edgels_HYPO1_template<float>(
        device_id, dev_All_Calib, 
        dev_All_Edgels_List, dev_Edgel_List_Start_Indx,
        dev_Edgel_H1_in_Meters, dev_Edgel_H1_Normal_Vectors );
}

//> Double Precision (DP64)
void gpu_Process_Edgels_HYPO1(
        int device_id, double* dev_All_Calib, 
        double* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx,
        double* dev_Edgel_H1_in_Meters, double* dev_Edgel_H1_Normal_Vectors )
{
    gpu_Process_Edgels_HYPO1_template<double>(
        device_id, dev_All_Calib, 
        dev_All_Edgels_List, dev_Edgel_List_Start_Indx,
        dev_Edgel_H1_in_Meters, dev_Edgel_H1_Normal_Vectors );
}

#endif    // GPU_PREPROCESS_EDGLES_HYPO1_HPP
