#ifndef DEFINITION_H
#define DEFINITION_H

//> Put thresholds here as macros 
#define REPO_DIR                   std::string("/gpfs/data/bkimia/cchien3/Edge_Based_Reconstruction/")
#define OUTPUT_WRITE_FOLDER        std::string("/gpfs/data/bkimia/cchien3/Edge_Based_Reconstruction/outputs_write_files/")

//> CH: Change this number to the number of CPU-cores you want your code to be processed by
#define NUM_OF_OPENMP_THREADS      (8)

#define DATASET_NUM_OF_FRAMES      (50)

#define HYPO1_VIEW_INDX            (5)
#define HYPO2_VIEW_INDX            (15)

#define DIST_THRESH                (2)
#define OREN_THRESH                (0.99)

#define MAX_NUM_OF_SUPPORT_VIEWS   (4)

#define imgcols                    (400)
#define imgrows                    (400)

#define PERCENT_EPIPOLE            (0.025)

#define IF_ICLNUIM_DATASET         (0)
#define IF_MULTIPLE_K              (1)

#define PI                         (3.1415926)
#define DEBUG                      (0)

//> GPU Implementations
//> PART I: PREPROCESSING ON EDGELS IN HYPO1
#define NUM_OF_THREADBLOCKS_PREPROCESS          (100)
#define NUM_OF_WARPS_PER_BLOCK_PREPROCESS       (1)
#define WARP_SIZE                               (32)      //> Constant. Must not change!
#define CHECK_PREPROCESS_CONSISTENCY_CPU_GPU    (true)
#define PREPROCESS_CONSISTENCY_CPU_GPU_TOL      (1e-5)

//> PART II: PAIRING EDGLES FROM HYPO1 AND HYPO2
#define NUM_OF_THREADBLOCKS        (3200)     //> Must be identical to the number of edgels in HYPO1
#define NUM_OF_THREADS_PER_BLOCK   (32)
#define GAMMA_INDEX_RANGE          (10)
#define TRUNCATED_WEDGE_RANGE      (2*GAMMA_INDEX_RANGE)

#define CALIB_FX                   (1075.65091572)
#define CALIB_FY                   (1073.90347929)

#define DEBUG_GPU                  (true)
#define CONSISTENCY_TOL            (1e-4)

//> CUDA error check
#define cudacheck( a )  do { \
                            cudaError_t e = a; \
                            if(e != cudaSuccess) { \
                                printf("\033[1;31m"); \
                                printf("Error in %s:%d %s\n", __func__, __LINE__, cudaGetErrorString(e)); \
                                printf("\033[0m"); \
                            }\
                        } while(0)


#endif
