//> Macros

//> 3D Edge Sketch Settings
#define NUM_OF_OPENMP_THREADS      (32)                 //> Number of CPU cores to run 3D Edge Sketch in parallel
#define HYPO1_VIEW_INDX            (6)                  //> Index of the 1st hypothesis view
#define HYPO2_VIEW_INDX            (8)                  //> Index of the 2nd hypothesis view
#define delta                      (0.3)                //> edge location perturbation
#define deltastr                   std::string("03")    //> string of delta value (used as part of the output file name)
#define OREN_THRESH                (15)                 //> \Delta \theta: orientation threshold
#define MAX_NUM_OF_SUPPORT_VIEWS   (4)                  //> N: minimal number of validation views supporting a hypothesis edge pair
#define threshEDG                  (32)                  //> third-order edge threshold for 1st round of multi-thresholding
#define threshEDGforall            (1)                  //> third-order edge threshold for last round of multi-thresholding
#define parallelangle              (15)                 //> decide if the picked hypothesis view should be abandoned by the angle between two epipolar lines
#define circleR                    (55)                 //> not sure what it is
#define NOT1STROUND                (0)                  //> not sure what it is
#define DIST_THRESH                (2)                  //> distance threshold between observed edge and reprojected edge on the validation view (probably not in use)

//> Dataset specific
#define DATASET_NAME               std::string("ABC-NEF")   //> ABC-NEF
#define SCENE_NAME                 std::string("00000006")
#define DATASET_NUM_OF_FRAMES      (50)
#define imgcols                    (800)
#define imgrows                    (800)
#define IF_MULTIPLE_K              (0)

//> Constant values (no change)
#define PI                         (3.1415926)

//> Debugging purpose
#define DEBUG                      (0)
#define DEBUG_READ_FILES           (true)

//> Some useful macros
#define LOG_INFOR_MESG(info_msg)        printf("\033[1;32m[INFO] %s\033[0m\n", std::string(info_msg).c_str() );
#define LOG_FILE_ERROR(err_msg)         printf("\033[1;31m[ERROR] File %s not found!\033[0m\n", std::string(err_msg).c_str() );
#define LOG_ERROR(err_msg)              printf("\033[1;31m[ERROR] %s\033[0m\n", std::string(err_msg).c_str() );
#define LOG_DATA_LOAD_ERROR(err_msg)    printf("\033[1;31m[DATA LOAD ERROR] %s not loaded successfully!\033[0m\n", std::string(err_msg).c_str() );

