#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <assert.h>
#include <string>
#include <ctime>
//> Inluce OpenMP here
#include <omp.h>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

//> YAML file data reader 
#include <yaml-cpp/yaml.h>

//> Include functions
#include "../Edge_Reconst/util.hpp"
#include "../Edge_Reconst/PairEdgeHypo.hpp"
#include "../Edge_Reconst/getReprojectedEdgel.hpp"
#include "../Edge_Reconst/getSupportedEdgels.hpp"
#include "../Edge_Reconst/getOrientationList.hpp"
#include "../Edge_Reconst/definitions.h"
#include "../Edge_Reconst/EdgeSketch_Core.hpp"
#include "../Edge_Reconst/file_reader.hpp"
#include "../Edge_Reconst/edge_mapping.hpp"

using namespace MultiviewGeometryUtil;

// ========================================================================================================================
// main function
//
// Modifications
//    Chien  24-07-06    Yilin finalizes the implementation of 3D edge sketch for reconstructing 3D edges from only one pair 
//                       of hypothesis view.
//    Chien  24-10-24    Qiwu continues the implementation enabling 3D edge sketch from multiple pairs of hypothesis images
//                       selected based on the projecting the 3D edges to 2D images.
//    Chien  24-11-07    Finish organizing the code so that everything is clean.
//
//> (c) LEMS, Brown University
//> Yilin Zheng (yilin_zheng@brown.edu)
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
//> Qiwu Zhang (qiwu_zhang@brown.edu)
// =========================================================================================================================

int main(int argc, char **argv) {

#if DELETE_ALL_FILES_UNDER_OUTPUTS
  std::string command = "rm -f ../../" + OUTPUT_FOLDER_NAME + "/*";
  if (system(command.c_str()) == 0) {
      std::cout << "Files under ../../outputs/ are removed successfully." << std::endl;
  } else {
      std::cerr << "Error removing files under ../../outputs/" << std::endl;
  }
#endif

  //> YAML file path
  std::string Edge_Sketch_Settings_Path = "../../3D_Edge_Sketch_Settings.yaml";

  YAML::Node Edge_Sketch_Settings_Map;
  try {
		Edge_Sketch_Settings_Map = YAML::LoadFile(Edge_Sketch_Settings_Path);
#if SHOW_EDGE_SKETCH_SETTINGS
		std::cout << std::endl << Edge_Sketch_Settings_Map << std::endl << std::endl;
#endif
	}
	catch (const std::exception& e) {
		std::cerr << "Exception: " << e.what() << std::endl;
    return 0;
	}

  //> Constructor
  EdgeSketch_Core MWV_Edge_Rec( Edge_Sketch_Settings_Map );

  //> Read camera intrinsic and extrinsic matrices
  MWV_Edge_Rec.Read_Camera_Data();

  //> Iteratively picking hypothesis view pairs to reconstruct 3D edges
  int edge_sketch_pass_count = 0;
  while ( edge_sketch_pass_count < MWV_Edge_Rec.Max_3D_Edge_Sketch_Passes ) {

    //> log: show the selected hypothesis views
    std::string out_str = "Selected views for hypotheses are " + std::to_string(MWV_Edge_Rec.hyp01_view_indx) \
                          + " and " + std::to_string(MWV_Edge_Rec.hyp02_view_indx);
    LOG_INFOR_MESG(out_str);

    //> Setup some constant data used throughout the 3D edge sketch
    MWV_Edge_Rec.Set_Hypothesis_Views_Camera();

    //> Multiple edge thresholding
    for (int edge_thresh = MWV_Edge_Rec.Edge_Detection_Init_Thresh; edge_thresh >= MWV_Edge_Rec.Edge_Detection_Final_Thresh; edge_thresh/=2) {
      std::cout << "       - Edge Threshold = " << edge_thresh << std::endl;
      MWV_Edge_Rec.thresh_EDG = edge_thresh;

      //> Load edges with specific third-order edge threshold
      MWV_Edge_Rec.Read_Edgels_Data();
      MWV_Edge_Rec.Set_Hypothesis_Views_Edgels();

      //> Hypothesis-Validation process
      MWV_Edge_Rec.Run_3D_Edge_Sketch();
    }
    
    //> Finalize hypothesis edge pairs for a two-view triangulation
    MWV_Edge_Rec.Finalize_Edge_Pairs_and_Reconstruct_3D_Edges();
    
    //> Stack all 3D edges located in the world coordinate
    MWV_Edge_Rec.Stack_3D_Edges();

    //> Find the next hypothesis view pairs, if any
    MWV_Edge_Rec.Project_3D_Edges_and_Find_Next_Hypothesis_Views();
    std::cout << "       - Average ratio of claimed 2D edges: " << MWV_Edge_Rec.avg_ratio << std::endl;

    MWV_Edge_Rec.Clear_Data();

    edge_sketch_pass_count++;

    if (MWV_Edge_Rec.enable_aborting_3D_edge_sketch)
      break;
  }



  double total_time = MWV_Edge_Rec.pair_edges_time + MWV_Edge_Rec.finalize_edge_pair_time + MWV_Edge_Rec.find_next_hypothesis_view_time;
  std::string out_time_str = "Total computation time: " + std::to_string(total_time) + " (s)";
  LOG_TIMEIMGS(out_time_str);
  std::cout << "     - Time for Pairing up Edges:              " << MWV_Edge_Rec.pair_edges_time << " (s)" << std::endl;
  std::cout << "     - Time for Finalizing Edge Pairs:         " << MWV_Edge_Rec.finalize_edge_pair_time << " (s)" << std::endl;
  std::cout << "     - Time for Finding Next Hypothesis Views: " << MWV_Edge_Rec.find_next_hypothesis_view_time << " (s)" << std::endl;

  LOG_INFOR_MESG("3D Edge Sketch is Finished!");
  return 0;
}
