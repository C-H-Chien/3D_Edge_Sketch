#ifndef DEV_FUNCTIONS_CUH
#define DEV_FUNCTIONS_CUH
// ============================================================================
// Device functions
//
// Change Logs
//    Chien  23-12-07:   Initially created. 
//    Chien  23-12-17:   Fix data dependency issue in mapping a third edge from 
//                       two correspondence edge. 
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// ============================================================================
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <math.h>
#include "../definitions.h"

//> included
#include <cuda_runtime.h>

template< typename T, int n >
__device__ void
magma_max_reduce( int i, T* x )
{
//> Credit: MAGMA library; Brown University LEMS lab in collaboration with MAGMA library developers. 
//  Ahmad Abdelfattah and Chiang-Heng Chien, Nov. 2021
    __syncthreads();
    if ( n > 1024 ) { if ( i < 1024 && i + 1024 < n ) { x[i] = max( x[i], x[i+1024] ); }  __syncthreads(); }
    if ( n >  512 ) { if ( i <  512 && i +  512 < n ) { x[i] = max( x[i], x[i+ 512] ); }  __syncthreads(); }
    if ( n >  256 ) { if ( i <  256 && i +  256 < n ) { x[i] = max( x[i], x[i+ 256] ); }  __syncthreads(); }
    if ( n >  128 ) { if ( i <  128 && i +  128 < n ) { x[i] = max( x[i], x[i+ 128] ); }  __syncthreads(); }
    if ( n >   64 ) { if ( i <   64 && i +   64 < n ) { x[i] = max( x[i], x[i+  64] ); }  __syncthreads(); }
    if ( n >   32 ) { if ( i <   32 && i +   32 < n ) { x[i] = max( x[i], x[i+  32] ); }  __syncthreads(); }
    // probably don't need __syncthreads for < 16 threads
    // because of implicit warp level synchronization.
    if ( n >   16 ) { if ( i <   16 && i +   16 < n ) { x[i] = max( x[i], x[i+  16] ); }  __syncthreads(); }
    if ( n >    8 ) { if ( i <    8 && i +    8 < n ) { x[i] = max( x[i], x[i+   8] ); }  __syncthreads(); }
    if ( n >    4 ) { if ( i <    4 && i +    4 < n ) { x[i] = max( x[i], x[i+   4] ); }  __syncthreads(); }
    if ( n >    2 ) { if ( i <    2 && i +    2 < n ) { x[i] = max( x[i], x[i+   2] ); }  __syncthreads(); }
    if ( n >    1 ) { if ( i <    1 && i +    1 < n ) { x[i] = max( x[i], x[i+   1] ); }  __syncthreads(); }
}

template< typename T >
__device__ __inline__ void
Map_Two_Correspondences_To_Third_View( const int tx, T* s_R21, T* s_T21, T* s_R31, T* s_T31, \
                                       T* s_gamma1_meter, T* s_N1_vec, \
                                       T r_gamma2_meter[2], T* s_gamma2_pixel, \
                                       T r_mapped_gamma3_pixel[2], T r_mapped_tangent3[3], \
                                       T* s_Calib_HYPO2, T* s_Calib_VALID )
{
  T n_e1, n_e2, n_e3;
  n_e1 = (s_T21[0] - s_T21[2]*r_gamma2_meter[0]) * (s_R31(0,0)*s_gamma1_meter[0] + s_R31(0,1)*s_gamma1_meter[1] + s_R31(0,2));
  n_e1 += ((s_R21(2,0)*s_gamma1_meter[0] + s_R21(2,1)*s_gamma1_meter[1] + s_R21(2,2))*r_gamma2_meter[0] \
          - (s_R21(0,0)*s_gamma1_meter[0] + s_R21(0,1)*s_gamma1_meter[1] + s_R21(0,2)))*s_T31[0];

  n_e2 = (s_T21[0] - s_T21[2]*r_gamma2_meter[0]) * (s_R31(1,0)*s_gamma1_meter[0] + s_R31(1,1)*s_gamma1_meter[1] + s_R31(1,2));
  n_e2 += ((s_R21(2,0)*s_gamma1_meter[0] + s_R21(2,1)*s_gamma1_meter[1] + s_R21(2,2))*r_gamma2_meter[0] \
          - (s_R21(0,0)*s_gamma1_meter[0] + s_R21(0,1)*s_gamma1_meter[1] + s_R21(0,2)))*s_T31[1];

  n_e3 = (s_T21[0] - s_T21[2]*r_gamma2_meter[0]) * (s_R31(2,0)*s_gamma1_meter[0] + s_R31(2,1)*s_gamma1_meter[1] + s_R31(2,2));
  n_e3 += ((s_R21(2,0)*s_gamma1_meter[0] + s_R21(2,1)*s_gamma1_meter[1] + s_R21(2,2))*r_gamma2_meter[0] \
          - (s_R21(0,0)*s_gamma1_meter[0] + s_R21(0,1)*s_gamma1_meter[1] + s_R21(0,2)))*s_T31[2];

  //> Reuse the memory of r_mapped_gamma3_pixel to compute the mapped gamma3 in meters
  r_mapped_gamma3_pixel[0] = n_e1 / n_e3; //> in meters!
  r_mapped_gamma3_pixel[1] = n_e2 / n_e3; //> in meters!
  
  //> Use the memory from r_mapped_tangent3 to compute the tangent in pixels for gamma2 and convert it to meters
  r_mapped_tangent3[0] = s_gamma2_pixel[tx*3]     + cos(s_gamma2_pixel[tx*3 + 2]);  //> This is x of (p2 + t2) in pixels
  r_mapped_tangent3[1] = s_gamma2_pixel[tx*3 + 1] + sin(s_gamma2_pixel[tx*3 + 2]);  //> This is y of (p2 + t2) in pixels

  r_mapped_tangent3[0] = (r_mapped_tangent3[0] - s_Calib_HYPO2[0]) / (CALIB_FX);  //> This is x of (p2 + t2) in meters
  r_mapped_tangent3[1] = (r_mapped_tangent3[1] - s_Calib_HYPO2[1]) / (CALIB_FY);  //> This is y of (p2 + t2) in meters

  //> tangent in meters
  r_mapped_tangent3[0] -= r_gamma2_meter[0];  //> This is x of t2 in meters
  r_mapped_tangent3[1] -= r_gamma2_meter[1];  //> This is y of t2 in meters

  //> Compute the normal vector, i.e., the cross product of tangent and gamma2, all in meters (t2 x gamma2)
  //  Here, we reuse the memory from n_e1, n_e2, and n_e3 as the normal vector n2 which is the cross product of r_mapped_tangent3 and r_gamma2_meter
  n_e1 = (r_mapped_tangent3[1]);
  n_e2 = -(r_mapped_tangent3[0]);
  n_e3 = r_mapped_tangent3[0]*r_gamma2_meter[1] - r_mapped_tangent3[1]*r_gamma2_meter[0];
  //> Now, (r_mapped_tangent3[0], r_mapped_tangent3[1], r_mapped_tangent3[2]) is the "normal vector" (n2) of the plane spanned by gamma2 and tangent2
  r_mapped_tangent3[0] = s_R21(0,0)*n_e1 + s_R21(1,0)*n_e2 + s_R21(2,0)*n_e3;
  r_mapped_tangent3[1] = s_R21(0,1)*n_e1 + s_R21(1,1)*n_e2 + s_R21(2,1)*n_e3;
  r_mapped_tangent3[2] = s_R21(0,2)*n_e1 + s_R21(1,2)*n_e2 + s_R21(2,2)*n_e3;

  //> Now, reuse the memory from r_mapped_tangent3 to compute the cross product of s_N1_vec and r_mapped_tangent3
  //  so that (n_e1, n_e2, n_e3) is cross(n1, n2)
  n_e1 = s_N1_vec[1]*r_mapped_tangent3[2] - s_N1_vec[2]*r_mapped_tangent3[1]; //> This is x of T_v1
  n_e2 = s_N1_vec[2]*r_mapped_tangent3[0] - s_N1_vec[0]*r_mapped_tangent3[2]; //> This is y of T_v1
  n_e3 = s_N1_vec[0]*r_mapped_tangent3[1] - s_N1_vec[1]*r_mapped_tangent3[0]; //> This is z of T_v1

  //> We reuse the memory of from n_e1, n_e2, and n_e3 again for computing Tangent_v3
  r_mapped_tangent3[0] = s_R31(0,0)*n_e1 + s_R31(0,1)*n_e2 + s_R31(0,2)*n_e3; //> This is x of T_v3
  r_mapped_tangent3[1] = s_R31(1,0)*n_e1 + s_R31(1,1)*n_e2 + s_R31(1,2)*n_e3; //> This is y of T_v3
  r_mapped_tangent3[2] = s_R31(2,0)*n_e1 + s_R31(2,1)*n_e2 + s_R31(2,2)*n_e3; //> This is z of T_v3

  //> Here we reuse the memory r_mapped_tangent3 for computing the projected tangent on the VALID view
  r_mapped_tangent3[0] -= r_mapped_tangent3[2]*r_mapped_gamma3_pixel[0]; //> x of projected t3
  r_mapped_tangent3[1] -= r_mapped_tangent3[2]*r_mapped_gamma3_pixel[1]; //> y of projected t3

  //> Reuse the memory from r_mapped_tangent3 to compute the norm
  r_mapped_tangent3[2] = sqrt(r_mapped_tangent3[0]*r_mapped_tangent3[0] + r_mapped_tangent3[1]*r_mapped_tangent3[1]); //> norm of projected t3

  //> Finally, this is the normalized mapped t3!
  r_mapped_tangent3[0] /= r_mapped_tangent3[2];
  r_mapped_tangent3[1] /= r_mapped_tangent3[2];

  //> Mapped gamma3 in pixels!!
  r_mapped_gamma3_pixel[0] = r_mapped_gamma3_pixel[0] * (CALIB_FX) + s_Calib_VALID[0];
  r_mapped_gamma3_pixel[1] = r_mapped_gamma3_pixel[1] * (CALIB_FY) + s_Calib_VALID[1];
}

#endif