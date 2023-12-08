#ifndef GPU_EDGE_RECONST_MAIN_HPP
#define GPU_EDGE_RECONST_MAIN_HPP

#include <iostream>
#include "../definitions.h"
#include "./dev_functions.cuh"

template< typename T >
__global__
void
gpu_Edge_Reconstruction_Kernel(
        T* dev_Rel_R21, T* dev_Rel_T21, T* dev_All_Calib,
        T* dev_Rel_Rots, T* dev_Rel_Transls,
        T* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx, 
        int *dev_Truncated_Wedge_Start_Indx_HYPO2, 
        int *dev_Truncated_Wedge_Start_Indx_VALID,
        int* dev_Valid_Views_Indices, T* dev_Wedge_Angle_Range_H2_to_VALID,
        T* dev_Edgel_H1_in_Meters, T* dev_Edgel_H1_Normal_Vectors,
        int* dev_Hypothesis_Edgel_Pair_Index)
{
    const int tid = threadIdx.x;
    const int bid  = blockIdx.x;

    //> Shared memory ptrs
    extern __shared__ double sdata[];
    T* s_gamma_HYPO2        = (T*)sdata;                                       //> gamma2 on the HYPO2
    T* s_gamma_VALID        = s_gamma_HYPO2 + (TRUNCATED_WEDGE_RANGE*3);       //> gamma3 on the validation view
    T* s_tangent_VALID      = s_gamma_VALID + (TRUNCATED_WEDGE_RANGE*3);       //> tangents of gamma3 on VALID
    T* s_R21                = s_tangent_VALID + (TRUNCATED_WEDGE_RANGE*2);     //> Relative rotation matrix HYPO2 w.r.t. HYPO1
    T* s_R31                = s_R21 + (9);                                     //> Relative rotation matrix VALID w.r.t. HYPO1
    T* s_T21                = s_R31 + (9);                                     //> Relative translation vector HYPO2 w.r.t. HYPO1
    T* s_T31                = s_T21 + (3);                                     //> Relative translation vector VALID w.r.t. HYPO1
    T* s_gamma1_pixel       = s_T31 + (3);
    T* s_gamma1_meter       = s_gamma1_pixel + (3);
    T* s_gamma1_normal_vec  = s_gamma1_meter + (2);
    T* s_Calib_HYPO1        = s_gamma1_normal_vec + (3);
    T* s_Calib_HYPO2        = s_gamma1_normal_vec + (2);
    T* s_Calib_VALID        = s_gamma1_normal_vec + (2);
    int* s_Num_Of_Supports  = (int*)(s_Calib_VALID + 2);
    int* s_Num_Of_Supports_ = s_Num_Of_Supports + TRUNCATED_WEDGE_RANGE;
    
    int vi = 0, gi = 0, valid_view_index;
    int truncated_wedge_start_idx_HYPO2, truncated_wedge_start_indx_VALID;
    int offset_indx_HYPO1, offset_indx_HYPO2, offset_indx_VALID;
    T r_gamma2_meter[2] = {0.0};

    //>>>>>>>>>>>>>>>>>>>> START OF FETCHING DATA >>>>>>>>>>>>>>>>>>>
    //> Read relative poses (HYPO2 w.r.t. HYPO1)
    if (tid < 9) s_R21[ tid ] = dev_Rel_R21[ tid ];
    if (tid < 3) s_T21[ tid ] = dev_Rel_T21[ tid ];
        
    //> Read calibration matrix, {cx, cy} of the hypothesis views
    s_Calib_HYPO1[0] = dev_All_Calib[ HYPO1_VIEW_INDX*(2)     ];
    s_Calib_HYPO1[1] = dev_All_Calib[ HYPO1_VIEW_INDX*(2) + 1 ];
    s_Calib_HYPO2[0] = dev_All_Calib[ HYPO2_VIEW_INDX*(2)     ];
    s_Calib_HYPO2[1] = dev_All_Calib[ HYPO2_VIEW_INDX*(2) + 1 ];

    //> Read gamma1 on HYPO1 image from the global memory, according to the thread-block ID
    //  Store in register file so that every thread holds this data.
    offset_indx_HYPO1 = dev_Edgel_List_Start_Indx[HYPO1_VIEW_INDX] * 3;
    if (tid < 3) { 
        s_gamma1_pixel[tid] = dev_All_Edgels_List[ offset_indx_HYPO1 + (bid)*3 + (tid) ];
        s_gamma1_normal_vec[tid] = dev_Edgel_H1_Normal_Vectors[ (bid)*3 + (tid) ];
    }
    if (tid < 2) s_gamma1_meter[tid] = dev_Edgel_H1_in_Meters[ bid*2 + tid ];

    //> Get all gamma2 on the truncated epipolar wedge from global memory to shared memory in parallel
    offset_indx_HYPO2 = dev_Edgel_List_Start_Indx[HYPO2_VIEW_INDX] * 3;
    truncated_wedge_start_idx_HYPO2 = dev_Truncated_Wedge_Start_Indx_HYPO2[bid];
    s_gamma_HYPO2[tid*3    ] = dev_All_Edgels_List[ offset_indx_HYPO2 + truncated_wedge_start_idx_HYPO2 + tid*3     ];
    s_gamma_HYPO2[tid*3 + 1] = dev_All_Edgels_List[ offset_indx_HYPO2 + truncated_wedge_start_idx_HYPO2 + tid*3 + 1 ];
    s_gamma_HYPO2[tid*3 + 2] = dev_All_Edgels_List[ offset_indx_HYPO2 + truncated_wedge_start_idx_HYPO2 + tid*3 + 2 ];
    //>>>>>>>>>>>>>>>>>>>> END OF FETCHING DATA >>>>>>>>>>>>>>>>>>>

    //> Compute gamma2 in meter
    r_gamma2_meter[0] = (s_gamma_HYPO2[tid*3    ] - s_Calib_HYPO2[0]) / (CALIB_FX);
    r_gamma2_meter[1] = (s_gamma_HYPO2[tid*3 + 1] - s_Calib_HYPO2[1]) / (CALIB_FY);

    T r_mapped_gamma3_pixel[2] = {0.0};
    T r_mapped_tangent3[3] = {0.0};
    T dist = 0.0;
    T tangent_dot_prod = 0.0;
    bool is_Supported_by_VALID = 0;
    s_Num_Of_Supports[tid] = 0;
    //s_Edgel_Indx_HYPO2[tid] = truncated_wedge_start_idx_HYPO2 + tid;

    //> Loop over all validation views
    #pragma unroll
    for (vi = 0; vi < (DATASET_NUM_OF_FRAMES-2); vi++) {
        valid_view_index = dev_Valid_Views_Indices[vi];
        truncated_wedge_start_indx_VALID = dev_Truncated_Wedge_Start_Indx_VALID[ (bid) + vi*NUM_OF_THREADBLOCKS ];

        //> Get the relative pose of VALID w.r.t. HYPO1 and {cx, cy}
        if (tid < 9) s_R31[ tid ] = dev_Rel_Rots[ vi*(9) + tid ];
        if (tid < 3) s_T31[ tid ] = dev_Rel_Transls[ vi*(3) + tid ];
        
        s_Calib_VALID[0] = dev_All_Calib[ valid_view_index*(2)     ];
        s_Calib_VALID[1] = dev_All_Calib[ valid_view_index*(2) + 1 ];

        //> Get all gamma3 on the truncated epipolar wedge from global memory to shared memory in parallel
        offset_indx_VALID = dev_Edgel_List_Start_Indx[valid_view_index] * 3;
        s_gamma_VALID[tid*3    ] = dev_All_Edgels_List[ offset_indx_VALID + truncated_wedge_start_indx_VALID + tid*3     ];
        s_gamma_VALID[tid*3 + 1] = dev_All_Edgels_List[ offset_indx_VALID + truncated_wedge_start_indx_VALID + tid*3 + 1 ];
        s_gamma_VALID[tid*3 + 2] = dev_All_Edgels_List[ offset_indx_VALID + truncated_wedge_start_indx_VALID + tid*3 + 2 ];
        s_tangent_VALID[tid*2    ] = cos(s_gamma_VALID[tid*3 + 2]);
        s_tangent_VALID[tid*2 + 1] = sin(s_gamma_VALID[tid*3 + 2]);

        //> Map gamma1 and gamma2 to gamma3 for both position (r_mapped_gamma3_pixel) and tangent (r_mapped_tangent3)
        Map_Two_Correspondences_To_Third_View<T>( tid, s_R21, s_T21, s_R31, s_T31, \
                                                  s_gamma1_meter, s_gamma1_normal_vec, \
                                                  r_gamma2_meter, s_gamma_HYPO2, \
                                                  r_mapped_gamma3_pixel, r_mapped_tangent3, \
                                                  s_Calib_HYPO2, s_Calib_VALID );

        //> Compare (i) the position {s_gamma_VALID, r_mapped_gamma3_pixel}, and 
        //          (ii) the dot product of tangents {s_tangent_VALID, r_mapped_tangent3}
        //          for all s_gamma_VALID
        #pragma unroll
        for (gi = 0; gi < (TRUNCATED_WEDGE_RANGE); gi++) {
            dist = sqrt( (s_gamma_VALID[gi*3]   - r_mapped_gamma3_pixel[0])*(s_gamma_VALID[gi*3  ] - r_mapped_gamma3_pixel[0]) \
                       + (s_gamma_VALID[gi*3+1] - r_mapped_gamma3_pixel[1])*(s_gamma_VALID[gi*3+1] - r_mapped_gamma3_pixel[1]) );

            tangent_dot_prod = fabs(s_tangent_VALID[gi*2]*r_mapped_tangent3[0] + s_tangent_VALID[gi*2 + 1]*r_mapped_tangent3[1] );
            //is_Supported_by_VALID = ( dist <= GAMMA3_POS_TOL_PIXEL && tangent_dot_prod >= GAMMA3_TANGENT_DOT_PRODUCT_TOL ) ? (true) : (false);
            is_Supported_by_VALID = ( dist <= GAMMA3_POS_TOL_PIXEL ) ? (true) : (false);

            if ( bid == 310 & tid == 10 ) {
                printf("%f, %f\n", dist, tangent_dot_prod);
            }
        }

        s_Num_Of_Supports[ tid ] = (is_Supported_by_VALID) ? (s_Num_Of_Supports[ tid ] + 1) : (s_Num_Of_Supports[ tid ]);
    }

    s_Num_Of_Supports_[ tid ] = s_Num_Of_Supports[ tid ];
    magma_max_reduce< int, TRUNCATED_WEDGE_RANGE >( tid, s_Num_Of_Supports );

    //> Write the index of the paired edgel in HYPO2 to the global memory
    if ( s_Num_Of_Supports[0] >= MAX_NUM_OF_SUPPORT_VIEWS ) {
        if ( s_Num_Of_Supports[0] == s_Num_Of_Supports_[ tid ] )
            dev_Hypothesis_Edgel_Pair_Index[ bid ] = truncated_wedge_start_idx_HYPO2 + tid;
    }
    else {
        dev_Hypothesis_Edgel_Pair_Index[ bid ] = -2;
    }
}

//extern "C"
//> Template for either SP32 or DP64
template< typename T >
void gpu_Edge_Reconstruction_template(
        int device_id,
        T* dev_Rel_R21, T* dev_Rel_T21, T* dev_All_Calib,
        T* dev_Rel_Rots, T* dev_Rel_Transls, 
        T* dev_All_Edgels_List, int* dev_Edgel_List_Start_Indx, 
        int *dev_Truncated_Wedge_Start_Indx_HYPO2, 
        int *dev_Truncated_Wedge_Start_Indx_VALID,
        int* dev_Valid_Views_Indices, T* dev_Wedge_Angle_Range_H2_to_VALID,
        T *dev_Edgel_H1_in_Meters, T* dev_Edgel_H1_Normal_Vectors,
        int* dev_Hypothesis_Edgel_Pair_Index )
{
    //> Kernel parameters
    dim3 grid(NUM_OF_THREADBLOCKS, 1, 1);
    //dim3 threads(WARP_SIZE, 1, 1);
    dim3 threads(TRUNCATED_WEDGE_RANGE, 1, 1);

    int shmem = 0;
    shmem += sizeof(T) * (TRUNCATED_WEDGE_RANGE * 3);   //> s_gamma_HYPO2
    shmem += sizeof(T) * (TRUNCATED_WEDGE_RANGE * 3);   //> s_gamma_VALID
    shmem += sizeof(T) * (TRUNCATED_WEDGE_RANGE * 2);   //> s_tangent_VALID
    shmem += sizeof(T) * (9);                           //> R21
    shmem += sizeof(T) * (9);                           //> R31
    shmem += sizeof(T) * (3);                           //> T21
    shmem += sizeof(T) * (3);                           //> T31
    shmem += sizeof(T) * (3);                           //> s_gamma1_pixel
    shmem += sizeof(T) * (2);                           //> s_gamma1_meter
    shmem += sizeof(T) * (3);                           //> s_gamma1_normal_vec
    shmem += sizeof(T) * (2);                           //> s_Calib_HYPO1
    shmem += sizeof(T) * (2);                           //> s_Calib_HYPO2
    shmem += sizeof(T) * (2);                           //> s_Calib_VALID
    shmem += sizeof(int) * (TRUNCATED_WEDGE_RANGE);     //> s_Num_Of_Supports
    shmem += sizeof(int) * (TRUNCATED_WEDGE_RANGE);     //> s_Num_Of_Supports_

    //> Get max. dynamic shared memory on the GPU
    int nthreads_max, shmem_max = 0;
    cudacheck( cudaDeviceGetAttribute(&nthreads_max, cudaDevAttrMaxThreadsPerBlock, device_id) );
#if CUDA_VERSION >= 9000
    cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlockOptin, device_id) );
    if (shmem <= shmem_max) {
        cudacheck( cudaFuncSetAttribute(gpu_Edge_Reconstruction_Kernel<T>, \
                                        cudaFuncAttributeMaxDynamicSharedMemorySize, shmem) );
    }
#else
    cudacheck( cudaDeviceGetAttribute (&shmem_max, cudaDevAttrMaxSharedMemoryPerBlock, device_id) );
#endif

    if ( shmem > shmem_max ) printf("error: kernel %s requires too many threads or too much shared memory\n", __func__);

    //> GPU Kernel Argument
    void *kernel_args[] = {&dev_Rel_R21, &dev_Rel_T21, &dev_All_Calib, \
                           &dev_Rel_Rots, &dev_Rel_Transls, \
                           &dev_All_Edgels_List, &dev_Edgel_List_Start_Indx, \
                           &dev_Truncated_Wedge_Start_Indx_HYPO2, &dev_Truncated_Wedge_Start_Indx_VALID, \
                           &dev_Valid_Views_Indices, &dev_Wedge_Angle_Range_H2_to_VALID, \
                           &dev_Edgel_H1_in_Meters, &dev_Edgel_H1_Normal_Vectors, \
                           &dev_Hypothesis_Edgel_Pair_Index };

    //> Launch the GPU kernel
    cudacheck( cudaLaunchKernel((void*)gpu_Edge_Reconstruction_Kernel<T>, grid, threads, kernel_args, shmem, NULL) );
}

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
        int* dev_Hypothesis_Edgel_Pair_Index )
{
    gpu_Edge_Reconstruction_template<float>(
        device_id,
        dev_Rel_R21, dev_Rel_T21, dev_All_Calib,
        dev_Rel_Rots, dev_Rel_Transls, 
        dev_All_Edgels_List, dev_Edgel_List_Start_Indx, 
        dev_Truncated_Wedge_Start_Indx_HYPO2, 
        dev_Truncated_Wedge_Start_Indx_VALID,
        dev_Valid_Views_Indices, dev_Wedge_Angle_Range_H2_to_VALID,
        dev_Edgel_H1_in_Meters, dev_Edgel_H1_Normal_Vectors,
        dev_Hypothesis_Edgel_Pair_Index );
}

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
        int* dev_Hypothesis_Edgel_Pair_Index )
{
    gpu_Edge_Reconstruction_template<double>(
        device_id,
        dev_Rel_R21, dev_Rel_T21, dev_All_Calib,
        dev_Rel_Rots, dev_Rel_Transls, 
        dev_All_Edgels_List, dev_Edgel_List_Start_Indx, 
        dev_Truncated_Wedge_Start_Indx_HYPO2, 
        dev_Truncated_Wedge_Start_Indx_VALID,
        dev_Valid_Views_Indices, dev_Wedge_Angle_Range_H2_to_VALID,
        dev_Edgel_H1_in_Meters, dev_Edgel_H1_Normal_Vectors,
        dev_Hypothesis_Edgel_Pair_Index );
}


#endif    // GPU_EDGE_RECONST_MAIN_HPP