//> Macros

//> 3D Edge Sketch Settings
#define NUM_OF_OPENMP_THREADS         (32)                 //> Number of CPU cores to run 3D Edge Sketch in parallel
//#define HYPO1_VIEW_INDX             (6)                  //> Index of the 1st hypothesis view
//#define HYPO2_VIEW_INDX             (8)                  //> Index of the 2nd hypothesis view
#define delta                         (0.3)                //> edge location perturbation
#define deltastr                      std::string("03")    //> string of delta value (used as part of the output file name)
#define OREN_THRESH                   (15)                 //> \Delta \theta: orientation threshold
#define MAX_NUM_OF_SUPPORT_VIEWS      (4)                  //> N: minimal number of validation views supporting a hypothesis edge pair
#define THRESHEDG                     (1)                  //> third-order edge threshold for 1st round of multi-thresholding
#define THRESEDGFORALL                (1)                  //> third-order edge threshold for last round of multi-thresholding
#define PARALLEL_EPIPOLAR_LINE_ANGLE  (15)                 //> decide if the picked hypothesis view should be abandoned by the angle between two epipolar lines
#define CIRCLER                       (55)                 //> threshold determining whether the distance between edges in hypothesis view 1 and 2 is small enough for them to be considered a match
#define NOT1STROUND                   (0)                  //> not sure what it is
#define DIST_THRESH                   (2)                  //> distance threshold between observed edge and reprojected edge on the validation view (probably not in use)

//> Dataset specific
#define DATASET_PATH                  std::string("/gpfs/data/bkimia/Datasets/")
#define DATASET_NAME                  std::string("ABC-NEF")   //> ABC-NEF/
#define SCENE_NAME                    std::string("00000006")
#define DATASET_NUM_OF_FRAMES         (50)
#define IMGCOLS                       (800)
#define IMGROWS                       (800)
#define IF_MULTIPLE_K                 (0)

//> Write to the files
#define WRITE_3D_EDGES                (false)

//> Print out in terminal
#define SHOW_EDGE_SKETCH_SETTINGS     (false)

//> Constant values (no change)
#define PI                            (3.1415926)

//> Debugging purpose
#define DEBUG                      (0)
#define DEBUG_READ_FILES           (false)
#define DEBUG_PAIRED_EDGES         (true)
#define SHOW_DATA_LOADING_INFO     (false)
#define SHOW_OMP_NUM_OF_THREADS    (true)

//> Some useful macros
#define LOG_INFOR_MESG(info_msg)        printf("\033[1;32m[INFO] %s\033[0m\n", std::string(info_msg).c_str() );
#define LOG_FILE_ERROR(err_msg)         printf("\033[1;31m[ERROR] File %s not found!\033[0m\n", std::string(err_msg).c_str() );
#define LOG_ERROR(err_msg)              printf("\033[1;31m[ERROR] %s\033[0m", std::string(err_msg).c_str() );
#define LOG_DATA_LOAD_ERROR(err_msg)    printf("\033[1;31m[DATA LOAD ERROR] %s not loaded successfully!\033[0m\n", std::string(err_msg).c_str() );

