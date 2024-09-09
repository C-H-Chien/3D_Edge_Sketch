#ifndef DEFINITION_H
#define DEFINITION_H

//> Put thresholds here as macros 

#define DATASET_NAME               std::string("ABC-NEF")
#define SCENE_NAME                 std::string("00000162")

//> CH: Change this number to the number of CPU-cores you want your code to be processed by
#define NUM_OF_OPENMP_THREADS      (32)

#define DEBUG                      (0)

#define DATASET_NUM_OF_FRAMES      (50)

#define HYPO1_VIEW_INDX            (11)
#define HYPO2_VIEW_INDX            (27) 
#define HYPO1_VIEW_INDX_NEXT       (19)
#define HYPO2_VIEW_INDX_NEXT       (28) 

#define delta                      (0.3)

#define deltastr                   std::string("03")

#define DIST_THRESH                (2)

#define OREN_THRESH                (15)

#define MAX_NUM_OF_SUPPORT_VIEWS   (4)

#define imgcols                    (800)

#define imgrows                    (800)

#define IF_ICLNUIM_DATASET         (0)

#define IF_MULTIPLE_K              (0)

#define PI                         (3.1415926)

#define parallelangle              (15)

#define circleR                    (55)

#define threshEDG                  (32)

#define threshEDGforall            (1)

#define NOT1STROUND                (0)

#endif
