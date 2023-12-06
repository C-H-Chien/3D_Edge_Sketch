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

//#include "gpu_kernels.hpp"
#include "./indices.hpp"
#include "../definitions.h"

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
                //if (fabs(ThirdParty_Matrix[k](i,j) - Converted_Matrix(i,j,k)) >= CONSISTENCY_TOL) {
                //    std::cout << ThirdParty_Matrix[k](i,j) << ", " << Converted_Matrix(i,j,k) << ", " << fabs(ThirdParty_Matrix[k](i,j) - Converted_Matrix(i,j,k)) << std::endl;
                //}
            }
        }
    }
    
    std::cout << ThirdParty_Matrix.size() << std::endl;
    std::cout << "(" << rows << ", " << cols << ")" << std::endl;
    std::cout << consistency_count << std::endl;
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
}

#endif    // CONSISTENCY_CHECK_HPP