//> Macros

//> Write to the files
#define WRITE_3D_EDGES                (true)

//> Print out in terminal
#define SHOW_EDGE_SKETCH_SETTINGS     (false)

//> Debugging purpose
#define DEBUG                      (0)
#define DEBUG_READ_FILES           (false)
#define DEBUG_PAIRED_EDGES         (false)
#define SHOW_DATA_LOADING_INFO     (false)
#define SHOW_OMP_NUM_OF_THREADS    (true)

//> Constant values (no change)
#define PI                            (3.1415926)

//> Some useful macros
#define LOG_INFOR_MESG(info_msg)        printf("\033[1;32m[INFO] %s\033[0m\n", std::string(info_msg).c_str() );
#define LOG_FILE_ERROR(err_msg)         printf("\033[1;31m[ERROR] File %s not found!\033[0m\n", std::string(err_msg).c_str() );
#define LOG_ERROR(err_msg)              printf("\033[1;31m[ERROR] %s\033[0m", std::string(err_msg).c_str() );
#define LOG_DATA_LOAD_ERROR(err_msg)    printf("\033[1;31m[DATA LOAD ERROR] %s not loaded successfully!\033[0m\n", std::string(err_msg).c_str() );

