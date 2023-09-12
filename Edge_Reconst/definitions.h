#ifndef DEFINITION_H
#define DEFINITION_H

//> Put thresholds here as macros 
#define REPO_DIR                   std::string("/users/yzhen105/Edge_Based_Reconstruction/")
#define OUTPUT_WRITE_FOLDER        std::string("/users/yzhen105/Edge_Based_Reconstruction/outputs_write_files/")

//> CH: Change this number to the number of CPU-cores you want your code to be processed by
#define NUM_OF_OPENMP_THREADS      (4)

#define DEBUG                      (1)

#define DIST_THRESH                (2)

#define OREN_THRESH                (0.9995)

#define HYPO1_VIEW_INDX            (5)

#define HYPO2_VIEW_INDX            (2)

#define MAX_NUM_OF_SUPPORT_VIEWS   (4)

#define imgcols                    (640)

#define imgrows                    (480)

#define DATASET_NUM_OF_FRAMES      (50)

#endif
