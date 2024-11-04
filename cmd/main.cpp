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

#include <set>


//> Include functions
#include "../Edge_Reconst/util.hpp"
#include "../Edge_Reconst/PairEdgeHypo.hpp"
#include "../Edge_Reconst/getReprojectedEdgel.hpp"
#include "../Edge_Reconst/getQuadrilateral.hpp"
#include "../Edge_Reconst/getSupportedEdgels.hpp"
#include "../Edge_Reconst/getOrientationList.hpp"
#include "../Edge_Reconst/linearTriangulationUtil.hpp"
#include "../Edge_Reconst/definitions.h"

//> Added by CH: Efficient Bucketing Method
#include "../Edge_Reconst/lemsvpe_CH/vgl_polygon_CH.hpp"
#include "../Edge_Reconst/lemsvpe_CH/vgl_polygon_scan_iterator_CH.hpp"
#include "../Edge_Reconst/subpixel_point_set.hpp"

//Qiwu
#include "../Edge_Reconst/file_reader.hpp"
#include "../Edge_Reconst/edge_mapping.hpp"
#include "../Edge_Reconst/iteration.hpp"

// using namespace std;
using namespace MultiviewGeometryUtil;

// ========================================================================================================================
// main function
//
// Modifications
//    Chien  24-07-06    Yilin finalizes the implementation of 3D edge sketch for reconstructing 3D edges from only one pair 
//                       of hypothesis view.
//    Chien  24-10-24    Qiwu continues the implementation enabling 3D edge sketch from multiple pairs of hypothesis images
//                       selected based on the projecting the 3D edges to 2D images.
//
//> (c) LEMS, Brown University
//> Yilin Zheng (yilin_zheng@brown.edu)
//> Qiwu Zhang (qiwu_zhang@brown.edu)
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =========================================================================================================================

Eigen::Vector3d transformToWorldCoordinates(const Eigen::Vector3d& point, 
                                            const Eigen::Matrix3d& R, 
                                            const Eigen::Vector3d& T) {
    // Apply the inverse transformation to convert the point to world coordinates
    return R.transpose() * (point - T);
}




void printAllSupportedIndices(const std::vector<Eigen::MatrixXd> &all_supported_indices) {
    std::cout << "Printing all supported indices:\n";
    
    // Loop through all matrices in all_supported_indices
    for (size_t idx = 0; idx < all_supported_indices.size(); ++idx) {
        std::cout << "Supported Indices for frame " << idx << ":\n";
        const Eigen::MatrixXd &matrix = all_supported_indices[idx];
        
        // Loop through each row and column of the matrix
        for (int row = 0; row < matrix.rows(); ++row) {
            for (int col = 0; col < matrix.cols(); ++col) {
                std::cout << std::fixed << std::setprecision(4) << matrix(row, col) << " ";
            }
            std::cout << std::endl;  // Move to the next line after printing all columns of a row
        }
        std::cout << std::endl;  // Space between different matrices
    }
}


void writeMatrixToFile(const std::string& filename, const Eigen::MatrixXd& matrix) {
    std::ofstream file(filename, std::ios::app); // Open in append mode
    if (file.is_open()) {
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                file << std::fixed << std::setprecision(4) << matrix(i, j) << " ";
            }
            file << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file " << filename << "\n";
    }
}


CorePipelineOutput core_pipeline(
  file_reader& Data_Reader,
  const int hyp01_view_indx, 
  const int hyp02_view_indx,
  const std::vector<Eigen::Matrix3d>& All_R, 
  const std::vector<Eigen::Vector3d>& All_T,
  const std::vector<Eigen::Matrix3d>& All_K,
  const Eigen::Matrix3d K,
  double rd_data,
  int d,
  int q)
{
    Eigen::Vector2d target_point_hypo1(520.5543, 426.4240);
    Eigen::Vector2d target_point_hypo2(519.6100, 407.3082);

  std::cout<< "pipeline start" <<std::endl;
  
  clock_t tstart, tend;
  double itime, ftime, exec_time, totaltime=0;
  int thresh_EDG = THRESHEDG;

  std::vector< Eigen::MatrixXd > all_supported_indices;
  EdgeMapping edgeMapping;
  //> Multi-thresholding!!!!!!
  while(thresh_EDG >= THRESEDGFORALL) {
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< READ FILES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    std::fstream Edge_File;

    std::vector<Eigen::MatrixXd> All_Edgels; 
    std::vector<Eigen::MatrixXd> All_Edgels_H12;  

    Eigen::Vector4d row_edge;
    int file_idx = 0;
    int H_idx = 0;
    
    // Read edgel files
    Data_Reader.readEdgelFiles(All_Edgels, Edge_File, rd_data, row_edge, file_idx, d, q, thresh_EDG);  
    Data_Reader.readHypothesisEdgelFiles(hyp01_view_indx, hyp02_view_indx, All_Edgels_H12, Edge_File, rd_data, row_edge, H_idx, file_idx, d, q, thresh_EDG);
  
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< PREPROCESSING >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    clock_t tstart_pre, tend_pre;
    tstart_pre = clock();
    MultiviewGeometryUtil::multiview_geometry_util util;
    PairEdgeHypothesis::pair_edge_hypothesis       PairHypo;
    GetReprojectedEdgel::get_Reprojected_Edgel     getReprojEdgel;
    GetQuadrilateral::get_Quadrilateral            getQuad;
    GetSupportedEdgels::get_SupportedEdgels        getSupport;
    GetOrientationList::get_OrientationList        getOre;
    
    // Assign variables required for Hypo1 and Hypo2
    Eigen::MatrixXd Edges_HYPO1 = All_Edgels_H12[0];
    Eigen::Matrix3d R1          = All_R[hyp01_view_indx];
    Eigen::Vector3d T1          = All_T[hyp01_view_indx];
    Eigen::MatrixXd Edges_HYPO2 = All_Edgels_H12[1];
    Eigen::Matrix3d R2          = All_R[hyp02_view_indx];
    Eigen::Vector3d T2          = All_T[hyp02_view_indx];
    Eigen::MatrixXd Edges_HYPO1_all = All_Edgels[hyp01_view_indx];
    Eigen::MatrixXd Edges_HYPO2_all = All_Edgels[hyp02_view_indx];
    // deal with multiple K scenario
    Eigen::Matrix3d K1;
    Eigen::Matrix3d K2;
    if(IF_MULTIPLE_K == 1){
      K1 = All_K[hyp01_view_indx];
      K2 = All_K[hyp02_view_indx];
    }else{
      K1 = K;
      K2 = K;
    }
    // Relative pose calculation
    Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
    Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
    Eigen::Matrix3d Tx  = util.getSkewSymmetric(T21);
    Eigen::Matrix3d E   = util.getEssentialMatrix(R21, T21);
    Eigen::Matrix3d F21 = util.getFundamentalMatrix(K1.inverse(), K2.inverse(), R21, T21);
    Eigen::Matrix3d R12 = util.getRelativePose_R21(R2, R1);
    Eigen::Vector3d T12 = util.getRelativePose_T21(R2, R1, T2, T1);  
    Eigen::Matrix3d F12 = util.getFundamentalMatrix(K2.inverse(), K1.inverse(), R12, T12);  
    // Initializations for paired edges between Hypo1 and Hypo 2
    Eigen::MatrixXd paired_edge = Eigen::MatrixXd::Constant(Edges_HYPO1_all.rows(),50,-2);
    // Compute epipolar wedges angles rance between Hypo1 and Hypo2 and find the angle range 1, defining a valid range for matching edges in Hypo2
    Eigen::MatrixXd OreListdegree    = getOre.getOreList(hyp01_view_indx, hyp02_view_indx, Edges_HYPO2, All_R, All_T, K1, K2);
    Eigen::MatrixXd OreListBardegree = getOre.getOreListBar(Edges_HYPO1, All_R, All_T, K1, K2, hyp02_view_indx, hyp01_view_indx);
    
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CORE PIPELINE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
  
    std::cout<< "Threshold is : " << thresh_EDG << std::endl;
    itime = omp_get_wtime();

  #pragma omp parallel
  {
    //> Local array stacking all supported indices
    std::vector<Eigen::MatrixXd> local_thread_supported_indices;
    int edge_idx;

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< First loop: loop over all edgels from hypothesis view 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    //<<<<<<<<<<< Identify pairs of edge, correct the positions of the edges from Hypo2, and store the paired edges >>>>>>>>>>>>>>>>//
    #pragma omp for schedule(static, NUM_OF_OPENMP_THREADS)
    for (edge_idx = 0; edge_idx < Edges_HYPO1.rows() ; edge_idx++) {

      //Edge Boundary Check: not too close to boundary
      if(Edges_HYPO1(edge_idx,0) < 10 || Edges_HYPO1(edge_idx,0) > IMGCOLS-10 || Edges_HYPO1(edge_idx,1) < 10 || Edges_HYPO1(edge_idx,1) > IMGROWS-10){
        continue;
      }
      //Paired Edge Check: not yet been paired
      if(paired_edge(edge_idx,0) != -2){
        continue;
      }
      //Get the current edge from Hypo1
      Eigen::Vector3d pt_edgel_HYPO1;
      pt_edgel_HYPO1 << Edges_HYPO1(edge_idx,0), Edges_HYPO1(edge_idx,1), 1;

      //Get angle Thresholds from OreListBar (Degree) 
      double thresh_ore21_1 = OreListBardegree(edge_idx,0);
      double thresh_ore21_2 = OreListBardegree(edge_idx,1);
      
      //Find Edges in Hypo2 Based on Angle Thresholds
      Eigen::MatrixXd HYPO2_idx_raw = PairHypo.getHYPO2_idx_Ore(OreListdegree, thresh_ore21_1, thresh_ore21_2);
      if (HYPO2_idx_raw.rows() == 0){
        continue;
      }

      //Retrieve Hypo2 Edgels
      Eigen::MatrixXd edgels_HYPO2 = PairHypo.getedgels_HYPO2_Ore(Edges_HYPO2, OreListdegree, thresh_ore21_1, thresh_ore21_2);
      //Correct Edgels in Hypo2 Based on Epipolar Constraints
      Eigen::MatrixXd edgel_HYPO1  = Edges_HYPO1.row(edge_idx);
      Eigen::MatrixXd edgels_HYPO2_corrected = PairHypo.edgelsHYPO2correct(edgels_HYPO2, edgel_HYPO1, F21, F12, HYPO2_idx_raw);
      //Organize the Final Edge Data
      Eigen::MatrixXd Edges_HYPO1_final(edgels_HYPO2_corrected.rows(),4);
      Edges_HYPO1_final << edgels_HYPO2_corrected.col(0), edgels_HYPO2_corrected.col(1), edgels_HYPO2_corrected.col(2), edgels_HYPO2_corrected.col(3);
      Eigen::MatrixXd Edges_HYPO2_final(edgels_HYPO2_corrected.rows(),4);
      Edges_HYPO2_final << edgels_HYPO2_corrected.col(4), edgels_HYPO2_corrected.col(5), edgels_HYPO2_corrected.col(6), edgels_HYPO2_corrected.col(7);

      //Store the Hypo2 Indices
      Eigen::MatrixXd HYPO2_idx(edgels_HYPO2_corrected.rows(),1); 
      HYPO2_idx << edgels_HYPO2_corrected.col(8);
      if (HYPO2_idx.rows() == 0){
        continue;
      }

      if (HYPO2_idx.rows() > 0) {
        {
            // Append matrices to file
            writeMatrixToFile("/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/48_edge.txt", Edges_HYPO1_final);
            writeMatrixToFile("/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/43_edge.txt", Edges_HYPO2_final);
        }
      }

    



      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Second loop:loop over all validation views >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
      // Initializations for all validation views
      int VALID_idx = 0;
      int stack_idx = 0;
      Eigen::MatrixXd supported_indices;
      supported_indices.conservativeResize(edgels_HYPO2.rows(),DATASET_NUM_OF_FRAMES-2);
      Eigen::MatrixXd supported_indice_current;
      supported_indice_current.conservativeResize(edgels_HYPO2.rows(),1);
      Eigen::MatrixXd supported_indices_stack;
      
      bool isempty_link = true;

      for (int VALID_INDX = 0; VALID_INDX < DATASET_NUM_OF_FRAMES; VALID_INDX++) {
        //> Skip the two hypothesis views
        if(VALID_INDX == hyp01_view_indx || VALID_INDX == hyp02_view_indx)
          continue;

        // Get camera pose and other info for current validation view
        Eigen::MatrixXd TO_Edges_VALID = All_Edgels[VALID_INDX];
        Eigen::Matrix3d R3             = All_R[VALID_INDX];
        Eigen::Vector3d T3             = All_T[VALID_INDX];
        Eigen::MatrixXd VALI_Orient    = TO_Edges_VALID.col(2);
        Eigen::MatrixXd Tangents_VALID;
        Tangents_VALID.conservativeResize(TO_Edges_VALID.rows(),2);
        Tangents_VALID.col(0)          = (VALI_Orient.array()).cos();
        Tangents_VALID.col(1)          = (VALI_Orient.array()).sin();
        Eigen::Matrix3d K3 = (IF_MULTIPLE_K == 1) ? All_K[VALID_INDX] : K;

        // Calculate relative pose
        Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
        Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);
        
        // Calculate the angle range for epipolar lines (Hypo1 --> Vali)
        Eigen::MatrixXd pt_edge            = Edges_HYPO1_final;
        Eigen::MatrixXd edge_tgt_gamma3    = getReprojEdgel.getGamma3Tgt(hyp01_view_indx, hyp02_view_indx, pt_edge, Edges_HYPO2_final, All_R, All_T, VALID_INDX, K1, K2);
        // Compute Orientation Lists for Hypo1 and Validation View (Vali)
        Eigen::MatrixXd OreListBardegree31 = getOre.getOreListBar(pt_edge, All_R, All_T, K1, K3, VALID_INDX, hyp01_view_indx);
        Eigen::MatrixXd OreListdegree31    = getOre.getOreListVali(TO_Edges_VALID, All_R, All_T, K1, K3, VALID_INDX, hyp01_view_indx);

        // Find all the edges fall inside epipolar wedge on validation view (Hypo1 --> Vali)
        Eigen::MatrixXd OreListBardegree32 = getOre.getOreListBar(Edges_HYPO2_final, All_R, All_T, K2, K3, VALID_INDX, hyp02_view_indx);
        Eigen::MatrixXd OreListdegree32    = getOre.getOreListVali(TO_Edges_VALID, All_R, All_T, K2, K3, VALID_INDX, hyp02_view_indx);
        // Identify Parallel Edges
        Eigen::VectorXd isparallel         = Eigen::VectorXd::Ones(Edges_HYPO2_final.rows());

        //std::cout<<"total number of hypothesis edge pairs is: "<<Edges_HYPO2_final.rows()<<std::endl;
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Third loop: loop over each edge from Hypo2 <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>//
        bool debug_print = false;
        std::set<int> written_frames;

        for (int idx_pair = 0; idx_pair < Edges_HYPO2_final.rows(); idx_pair++) {
          double epsilon = 1e-4; 
          if (std::abs(Edges_HYPO1_final(idx_pair, 0) - target_point_hypo1(0)) < epsilon &&
              std::abs(Edges_HYPO1_final(idx_pair, 1) - target_point_hypo1(1)) < epsilon) {
                debug_print = true;
          }
                
          double thresh_ore31_1 = OreListBardegree31(idx_pair,0);
          double thresh_ore31_2 = OreListBardegree31(idx_pair,1);
          double thresh_ore32_1 = OreListBardegree32(idx_pair,0);
          double thresh_ore32_2 = OreListBardegree32(idx_pair,1);
          
          Eigen::MatrixXd vali_idx31 = PairHypo.getHYPO2_idx_Ore(OreListdegree31, thresh_ore31_1, thresh_ore31_2);
          Eigen::MatrixXd edgels_31  = PairHypo.getedgels_HYPO2_Ore(TO_Edges_VALID, OreListdegree31, thresh_ore31_1, thresh_ore31_2);

          // Find all the edges fall inside epipolar wedge on validation view (Hypo2 --> Vali)
          Eigen::MatrixXd vali_idx32 = PairHypo.getHYPO2_idx_Ore(OreListdegree32, thresh_ore32_1, thresh_ore32_2);
          Eigen::MatrixXd edgels_32  = PairHypo.getedgels_HYPO2_Ore(TO_Edges_VALID, OreListdegree32, thresh_ore32_1, thresh_ore32_2);



          if (debug_print && written_frames.find(VALID_INDX) == written_frames.end()) {
              // Add frame to the set to mark it as written
              written_frames.insert(VALID_INDX);

              // Write threshold values to the file
              std::ofstream thresh_file("/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/threshold_ranges.txt", std::ios::app);
              if (thresh_file.is_open()) {
                  thresh_file << "frame #: " << VALID_INDX << "\n";
                  thresh_file << "Thresh_ore21_1: " << thresh_ore21_1 << ", Thresh_ore21_2: " << thresh_ore21_2 << "\n\n";
                  thresh_file.close();
              } else {
                  std::cerr << "Error opening file to save threshold values.\n";
              }

              // Save edgels_31 to file
              std::ofstream edgels31_file("/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_31.txt", std::ios::app);
              if (edgels31_file.is_open()) {
                  for (int i = 0; i < edgels_31.rows(); ++i) {
                      edgels31_file << edgels_31(i, 0) << " " << edgels_31(i, 1) << "\n";
                  }
                  edgels31_file << "\n";  // Separate entries by an empty line
                  edgels31_file.close();
              } else {
                  std::cerr << "Error opening file to save edgels_31 values.\n";
              }

              // Save edgels_32 to file
              std::ofstream edgels32_file("/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/edgels_32.txt", std::ios::app);
              if (edgels32_file.is_open()) {
                  for (int i = 0; i < edgels_32.rows(); ++i) {
                      edgels32_file << edgels_32(i, 0) << " " << edgels_32(i, 1) << "\n";
                  }
                  edgels32_file << "\n";  // Separate entries by an empty line
                  edgels32_file.close();
              } else {
                  std::cerr << "Error opening file to save edgels_32 values.\n";
              }
          }


 
          // Check if the two wedges could be considered as parallel to each other
          Eigen::MatrixXd anglediff(4,1);
          anglediff << fabs(thresh_ore31_1 - thresh_ore32_1), 
                      fabs(thresh_ore31_1 - thresh_ore32_2),
                      fabs(thresh_ore31_2 - thresh_ore32_1),
                      fabs(thresh_ore31_2 - thresh_ore32_2);
          if(anglediff.maxCoeff() <= PARALLEL_EPIPOLAR_LINE_ANGLE) {
            isparallel.row(idx_pair) << 0;
            supported_indice_current.row(idx_pair) << -2;
            continue;
          }

          // Find all the edges fall inside the two epipolar wedges intersection on validation view 
          // (Hypo1 --> Vali) && (Hypo2 --> Vali)
          // This vector will store the indices that are common between vali_idx31 and vali_idx32
          std::vector<double> v_intersection;
          // Convert Eigen Matrices to Standard Vectors
          std::vector<double> v1(vali_idx31.data(), vali_idx31.data() + vali_idx31.rows());
          std::vector<double> v2(vali_idx32.data(), vali_idx32.data() + vali_idx32.rows());
          // Find Intersection
          set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v_intersection));
          // Convert Intersection to Eigen Vector
          Eigen::VectorXd idxVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v_intersection.data(), v_intersection.size());
          // Create Eigen Matrix Representation of Inliers
          Eigen::MatrixXd inliner(idxVector);

          // Print inliner before calling getSupportIdx
          if (debug_print){
            std::cout << "Inliner indices before getSupportIdx:\n" << inliner << std::endl;
          }
          
          // Calculate orientation of gamma 3
          Eigen::Vector2d edgels_tgt_reproj = {edge_tgt_gamma3(idx_pair,0), edge_tgt_gamma3(idx_pair,1)};
          // Get support from validation view for this gamma 3
          double supported_link_indx = getSupport.getSupportIdx(edgels_tgt_reproj, Tangents_VALID, inliner, debug_print);

          if (supported_link_indx >= TO_Edges_VALID.rows()) {
            LOG_ERROR("Something fishy here!\n");
            std::cout << "(VALID_INDX, supported_link_indx, TO_Edges_VALID.rows())) = (" << VALID_INDX << ", " << supported_link_indx << ", " << TO_Edges_VALID.rows() << ")" << std::endl;
          }

          // Get the supporting edge idx from this validation view
          supported_indice_current.row(idx_pair) << supported_link_indx;

          if (supported_link_indx != -2) {
            supported_indices_stack.conservativeResize(stack_idx+1,2);
            supported_indices_stack.row(stack_idx) << double(idx_pair), double(supported_link_indx);
            isempty_link = false;
            stack_idx++;
          }

          ///////////////////////debug /////////////////////////////////////////////
          if (std::abs(Edges_HYPO1_final(idx_pair, 0) - target_point_hypo1(0)) < epsilon &&
              std::abs(Edges_HYPO1_final(idx_pair, 1) - target_point_hypo1(1)) < epsilon) {
              std::cout << "!!!!Appending!!!!!" << std::endl;
              // File to save the supporting edge points
              std::ofstream support_edge_file("/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/supporting_edge_points.txt", std::ios::app);
              if (!support_edge_file.is_open()) {
                  std::cout << "Failed to open the supporting edge file for writing.\n"<<std::endl;
              }
              int support_idx = supported_indice_current(idx_pair);
              std::cout<<"validation view "<<VALID_INDX<< ",   support idx is: "<<support_idx<<std::endl;
              if (support_idx != -2 && support_idx < TO_Edges_VALID.rows()) {
                  Eigen::Vector2d supporting_edge = TO_Edges_VALID.row(support_idx).head<2>();
                  support_edge_file << "Validation View " << VALID_INDX << ": "
                                    << supporting_edge(0) << " " << supporting_edge(1) << "\n";
              }

              // Close the file after processing all validation views for this idx_pair
              support_edge_file.close();
          }
          debug_print = false;
        ///////////////////////debug /////////////////////////////////////////////
        }
        supported_indices.col(VALID_idx) << supported_indice_current.col(0);

        VALID_idx++;
      } 
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  End of second loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

      //> Now, for each local thread, stack the supported indices
      local_thread_supported_indices.push_back(supported_indices);


      //Check for Empty Supported Indices
      if (isempty_link) {
        continue;
      }
      
      //Create a Stack of Supported Indices
      std::vector<double> indices_stack(supported_indices_stack.data(), supported_indices_stack.data() + supported_indices_stack.rows());
      std::vector<double> indices_stack_unique = indices_stack;
      std::sort(indices_stack_unique.begin(), indices_stack_unique.end());
      std::vector<double>::iterator it1;
      it1 = std::unique(indices_stack_unique.begin(), indices_stack_unique.end());
      indices_stack_unique.resize( std::distance(indices_stack_unique.begin(),it1) );

      Eigen::VectorXd rep_count;
      rep_count.conservativeResize(indices_stack_unique.size(),1);

      //Count the Occurrence of Each Unique Index
      for(int unique_idx = 0; unique_idx<indices_stack_unique.size(); unique_idx++){
        rep_count.row(unique_idx) << double(count(indices_stack.begin(), indices_stack.end(), indices_stack_unique[unique_idx]));
      }

      //Find Maximum Support
      Eigen::VectorXd::Index   maxIndex;
      double max_support = rep_count.maxCoeff(&maxIndex);
      int numofmax = std::count(rep_count.data(), rep_count.data()+rep_count.size(), max_support);
      
      //Threshold Check
      if( max_support < MAX_NUM_OF_SUPPORT_VIEWS){
        continue;
      }
      int finalpair = -2;
      if(numofmax > 1){
        std::vector<double> rep_count_vec(rep_count.data(), rep_count.data() + rep_count.rows());
        //std::cout<< "here"<<std::endl;
        std::vector<int> max_index;
        auto start_it = begin(rep_count_vec);
        while (start_it != end(rep_count_vec)) {
          start_it = std::find(start_it, end(rep_count_vec), max_support);
          if (start_it != end(rep_count_vec)) {
            auto const pos = std::distance(begin(rep_count_vec), start_it);
            max_index.push_back(int(pos));
            ++start_it;
          }
        }

        //Select the Final Paired Edge
        Eigen::Vector3d coeffs;
        coeffs = F21 * pt_edgel_HYPO1;
        Eigen::MatrixXd Edge_Pts;
        Edge_Pts.conservativeResize(max_index.size(),2);
        for(int maxidx = 0; maxidx<max_index.size(); maxidx++){
          Edge_Pts.row(maxidx) << Edges_HYPO2_final(indices_stack_unique[max_index[maxidx]], 0),Edges_HYPO2_final(indices_stack_unique[max_index[maxidx]], 1) ;
        }
        Eigen::VectorXd Ap = coeffs(0)*Edge_Pts.col(0);
        Eigen::VectorXd Bp = coeffs(1)*Edge_Pts.col(1);
        Eigen::VectorXd numDist = Ap + Bp + Eigen::VectorXd::Ones(Ap.rows())*coeffs(2);
        double denomDist = coeffs(0)*coeffs(0) + coeffs(1)*coeffs(1);
        denomDist = sqrt(denomDist);
        Eigen::VectorXd dist = numDist.cwiseAbs()/denomDist;
        Eigen::VectorXd::Index   minIndex;
        double min_dist = dist.minCoeff(&minIndex);
        if(min_dist > DIST_THRESH){
          continue;
        }
        finalpair = int(indices_stack_unique[max_index[minIndex]]);
      }
      else {
        finalpair = int(indices_stack_unique[int(maxIndex)]);
      }

      //> paired_edge is a row vector continaing (hypo1 edge index), (paired hypo2 edge index), (paired validation edge indices)
      paired_edge.row(edge_idx) << edge_idx, HYPO2_idx(finalpair), supported_indices.row(finalpair);
      //std::cout << "paired_edge.row(edge_idx): \n" << paired_edge.row(edge_idx) <<std::endl;
    }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  End of first loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

    //> A critical session to stack all local supported indices
    #pragma omp critical
    all_supported_indices.insert(all_supported_indices.end(), local_thread_supported_indices.begin(), local_thread_supported_indices.end());
  } //> End of pragma omp parallel
  ftime = omp_get_wtime();
  exec_time = ftime - itime;
  totaltime += exec_time;
    
    //> Finalize Paired Edges
    tstart = clock();
    int pair_num = 0;
    std::vector<int> valid_pair_index;
    for(int pair_idx = 0; pair_idx < paired_edge.rows(); pair_idx++) {
      if(paired_edge(pair_idx,0) != -2 && paired_edge(pair_idx,0) != -3) {
        valid_pair_index.push_back(pair_idx);
        pair_num++;
      }
    }
    Eigen::MatrixXd paired_edge_final = Eigen::MatrixXd::Constant(pair_num,50,-2);
    for(int i = 0; i < pair_num; i++){
      paired_edge_final.row(i) << paired_edge.row(valid_pair_index[i]);
    }
    std::string info_str = "Number of pairs is: " + std::to_string(pair_num);
    LOG_INFOR_MESG(info_str);

#if DEBUG_PAIRED_EDGES
    std::ofstream debug_file_paired_edges;
    std::string Output_File_Path_Paired_Edges = "../../outputs/paired_edge_final.txt";
    debug_file_paired_edges.open(Output_File_Path_Paired_Edges);
    debug_file_paired_edges << paired_edge_final;
    debug_file_paired_edges.close();
#endif
    
    Eigen::MatrixXd Gamma1s;
    // gamma3 calculation, next view index pre
    if(thresh_EDG == 1) {
      std::vector<Eigen::Matrix3d> Rs;
      Rs.push_back(R21);
      std::vector<Eigen::Vector3d> Ts;
      Ts.push_back(T21);
      std::vector<double> K1_v;
      K1_v.push_back(K1(0,2));
      K1_v.push_back(K1(1,2));
      K1_v.push_back(K1(0,0));
      K1_v.push_back(K1(1,1));
      LOG_INFOR_MESG("Finalizing edge pairs");
      std::cout << paired_edge_final.rows() << std::endl;
      Gamma1s.conservativeResize(paired_edge_final.rows(),3);
      for (int pair_idx = 0; pair_idx < paired_edge_final.rows(); pair_idx++) {
        // std::cout << pair_idx << ", ";
        // std::cout << paired_edge_final(pair_idx,0) << ", " << paired_edge_final(pair_idx,1) << ", ";
        Eigen::MatrixXd edgel_HYPO1   = Edges_HYPO1_all.row(int(paired_edge_final(pair_idx,0)));  //> edge index in hypo 1
        Eigen::MatrixXd edgel_HYPO2   = Edges_HYPO2_all.row(int(paired_edge_final(pair_idx,1)));  //> edge index in hypo 2
        // Eigen::MatrixXd HYPO2_idx_raw = Edges_HYPO1.row(int(paired_edge_final(pair_idx,0)));      //> WTF is this? Edges_HYPO1 = Edges_HYPO1_all
        Eigen::MatrixXd HYPO2_idx_raw = Edges_HYPO2.row(int(paired_edge_final(pair_idx,1)));

        // std::cout << "P0, ";
        Eigen::MatrixXd edgels_HYPO2_corrected = PairHypo.edgelsHYPO2correct(edgel_HYPO2, edgel_HYPO1, F21, F12, HYPO2_idx_raw);

        // std::cout << "P1, ";

        if (HYPO2_idx_raw.rows() == 0 || edgels_HYPO2_corrected.rows() == 0) {
          std::cout << "No valid matches found for edge " << pair_idx << " at threshold " << thresh_EDG << std::endl;
          continue;
        }

        Eigen::MatrixXd Edges_HYPO1_final(edgels_HYPO2_corrected.rows(),4);
        Edges_HYPO1_final << edgels_HYPO2_corrected.col(0), edgels_HYPO2_corrected.col(1), edgels_HYPO2_corrected.col(2), edgels_HYPO2_corrected.col(3);
        Eigen::MatrixXd Edges_HYPO2_final(edgels_HYPO2_corrected.rows(),4);
        Edges_HYPO2_final << edgels_HYPO2_corrected.col(4), edgels_HYPO2_corrected.col(5), edgels_HYPO2_corrected.col(6), edgels_HYPO2_corrected.col(7);

        // std::cout << "P2, ";
        
        Eigen::Vector2d pt_H1 = Edges_HYPO1_final.row(0);
        Eigen::Vector2d pt_H2 = Edges_HYPO2_final.row(0);
        std::vector<Eigen::Vector2d> pts;
        pts.push_back(pt_H1);
        pts.push_back(pt_H2);

        //> The resultant edge_pt_3D is 3D edges under the first hypothesis view coordinate
        Eigen::Vector3d edge_pt_3D = linearTriangulation(2, pts, Rs, Ts, K1_v);

        // std::cout << "P3, ";

        if (edge_pt_3D.hasNaN()) {
          LOG_ERROR("NaN values detected in edge_pt_3D for pair_idx: ");
          Gamma1s.row(pair_idx)<< 0, 0, 0;  //> TBD
          std::cerr << pair_idx << std::endl;
          continue;
        }

        Gamma1s.row(pair_idx) << edge_pt_3D(0), edge_pt_3D(1), edge_pt_3D(2);
        edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, pt_H1, hyp01_view_indx);
        edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, pt_H2, hyp02_view_indx);

        // std::cout << "P4, ";

        ///////////////////////////////// Add support from validation views /////////////////////////////////
        std::vector<std::pair<int, Eigen::Vector2d>> validation_support_edges;

        // Loop through validation views to find the supporting edges
        int val_indx_in_paired_edge_array = 2; // +2 accounts for the first two columns for HYPO1 and HYPO2
        for (int val_idx = 0; val_idx < DATASET_NUM_OF_FRAMES; ++val_idx) {
            if (val_idx == hyp01_view_indx || val_idx == hyp02_view_indx) {
              continue;  // Skip hypothesis views
            }

            // Retrieve support index from paired_edge for the current validation view
            int support_idx = paired_edge_final(pair_idx, val_indx_in_paired_edge_array);  
            if (support_idx != -2) {
                // Retrieve the supporting edge from the validation view
                Eigen::MatrixXd edges_for_val_frame = All_Edgels[val_idx];

                if (edges_for_val_frame.rows() <= support_idx) {
                  LOG_ERROR("Something buggy here!\n");
                  std::cout << "(pair_idx, val_idx, edges_for_val_frame.rows(), support_idx) = (" << pair_idx << ", " << val_idx << ", " << edges_for_val_frame.rows() << ", " << support_idx << ")" << std::endl;
                  Eigen::Vector2d supporting_edge = edges_for_val_frame.row(support_idx).head<2>();
                }

                Eigen::Vector2d supporting_edge = edges_for_val_frame.row(support_idx).head<2>();

                // Store validation view and the supporting edge
                validation_support_edges.emplace_back(val_idx, supporting_edge);

                // Add the supporting edge to the edgeMapping for the 3D edge
                edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, supporting_edge, val_idx);
            }
            val_indx_in_paired_edge_array++;
        }

        // std::cout << "P5, " << paired_edge_final.rows() << std::endl;
        
        ///////////////////////////////// Add support from validation views /////////////////////////////////
      }
      std::cout << "Finish for-loop" << std::endl;
    }

    tend = clock() - tstart; 
    //std::cout << "It took "<< double(tend)/double(CLOCKS_PER_SEC) <<" second(s) to generate the final edge pair for output file."<<std::endl;
    std::cout << "Number of pairs found till this round: " << paired_edge_final.rows()<<std::endl;

    thresh_EDG /= 2; 
    std::vector<Eigen::MatrixXd> All_Edgels_H12_1;  
    if (thresh_EDG >= 1) {
      while(file_idx < 3) {
        if(file_idx == 1){
          H_idx = hyp01_view_indx;
        }else{
          H_idx = hyp02_view_indx;
        }
       std::string Edge_File_PathH12 = "../../datasets/" + DATASET_NAME + "/" + SCENE_NAME + "/Edges/Edge_" + std::to_string(H_idx)+"_t"+std::to_string(thresh_EDG)+".txt"; 
#if DEBUG_READ_FILES
        std::cout << Edge_File_PathH12 << std::endl;
#endif
        file_idx ++;
        Eigen::MatrixXd Edgels; //> Declare locally, ensuring the memory addresses are different for different frames
        Edge_File.open(Edge_File_PathH12, std::ios_base::in);
        if (!Edge_File) { 
          LOG_FILE_ERROR("Edge file not existed!\n"); exit(1); 
        }
        else {
          Edgels.resize(1,4);
          while (Edge_File >> rd_data) {
            row_edge(q) = rd_data;
            q++;
            if(q>3){
              Edgels.conservativeResize(d+1,4);
              Edgels.row(d) = row_edge;
              q = 0;
              d++;
            }
          }
          Edge_File.close();
          All_Edgels_H12_1.push_back(Edgels);
          //Edgels = {};
          d = 0;
          q = 0;
        }
      }
      file_idx = 1;
      Edges_HYPO1 = All_Edgels_H12_1[0];
      Edges_HYPO2 = All_Edgels_H12_1[1];

      for (int pair_idx = 0; pair_idx < paired_edge.rows(); pair_idx++) {
        if (paired_edge(pair_idx,0) != -2 && paired_edge(pair_idx,0) != -3 ) {
          Edges_HYPO2.row(int(paired_edge(pair_idx,1))) << 0,0,0,0;
          pair_num++;
        }
      }
    }
    else {
      Eigen::Matrix3d R_ref = All_R[hyp01_view_indx];
      Eigen::Vector3d T_ref = All_T[hyp01_view_indx];

      Eigen::MatrixXd Gamma1s_world(Gamma1s.rows(), 3);
      for (int i = 0; i < Gamma1s.rows(); ++i) {
          Eigen::Vector3d point_camera = Gamma1s.row(i).transpose();
          Eigen::Vector3d point_world = transformToWorldCoordinates(point_camera, R_ref, T_ref);
          Gamma1s_world.row(i) = point_world.transpose();
      }

      std::cout<< "pipeline finished" <<std::endl;
      std::cout << "It took "<< totaltime <<" second(s) to finish the whole pipeline."<<std::endl;
      std::ofstream myfile2;
      std::string Output_File_Path2 = "../../outputs/Gamma1s_" + DATASET_NAME + "_" + SCENE_NAME + "_" + std::to_string(hyp01_view_indx)+"n"+std::to_string(hyp02_view_indx)+"_t32to"+std::to_string(thresh_EDG) + "_delta" + deltastr +"_theta" + std::to_string(OREN_THRESH) + "_N" + std::to_string(MAX_NUM_OF_SUPPORT_VIEWS) + ".txt";
      std::cout << Output_File_Path2 << std::endl;
      myfile2.open (Output_File_Path2);
      myfile2 << Gamma1s_world;
      myfile2.close();
      Gamma1s = Gamma1s_world;
    }


    CorePipelineOutput output;
    output.Gamma1s = Gamma1s;

    // Directly add to the frame-to-edge mapping while processing edges
    for (const auto& [edge_3D, supporting_edges] : edgeMapping.edge_3D_to_supporting_edges) {
      for (const auto& [supporting_edge, image_number] : supporting_edges) {
        // Directly update the frame-to-edge mapping
        edgeMapping.add3DToFrameMapping(edge_3D, supporting_edge, image_number);
      }
    }

    // Set the frame-to-edge mapping in the output structure
    output.frame_to_edge_to_3D_map = edgeMapping.frame_to_edge_to_3D_map;

    std::cout << "pipeline finished" << std::endl;
    std::cout << "It took " << totaltime << " second(s) to finish the whole pipeline." << std::endl;

    return output;
  } //> while-loop
}



int main(int argc, char **argv) {

    int hyp01_view_indx  = 47;
    int hyp02_view_indx  = 42;

    //> initiate data reader class object
    file_reader Data_Reader; 

    std::vector<Eigen::Matrix3d> All_R;
    std::vector<Eigen::Vector3d> All_T;
    std::vector<Eigen::Matrix3d> All_K;
    std::vector<Eigen::MatrixXd> All_Edgels; 

    std::fstream Rmatrix_File, Tmatrix_File, Kmatrix_File;
    std::fstream Edge_File;
    
    Eigen::Matrix3d R_matrix, K_matrix;
    Eigen::Vector3d row_R, row_K, T_matrix;
    Eigen::Matrix3d K;
    Eigen::Vector4d row_edge;
    Eigen::MatrixXd all_Edges_3D;

    // Read all required matrices (rotation, translation, and camera matrices)
    double rd_data;
    int d = 0, q = 0;
    int file_idx = 0;

    // Cumulative storage for all frames and edges to 3D mappings across iterations
    std::map<int, std::unordered_map<Eigen::Vector2d, std::vector<Eigen::Vector3d>, HashEigenVector2d>> cumulative_frame_to_edge_to_3D_map;

    ///////////////////// find worst frames according to every 3d edges not only the previous one ////////////////////
    Data_Reader.readRmatrix(All_R, R_matrix, Rmatrix_File, rd_data, row_R, d, q);
    Data_Reader.readTmatrix(All_T, T_matrix, Tmatrix_File, rd_data, d, q);
    Data_Reader.readK(Kmatrix_File, All_K, K, K_matrix, row_K, rd_data, d, q);

    unsigned nthreads = NUM_OF_OPENMP_THREADS;
    omp_set_num_threads(nthreads);
    int ID = omp_get_thread_num();
    std::cout << "Using " << nthreads << " threads for OpenMP parallelization." << std::endl;
    std::cout << "nthreads: " << nthreads << "." << std::endl;

    bool sufficient_reconstruct = false;
    int iteration = 0;

    //while (!sufficient_reconstruct) {
    while (iteration < 1) {

      std::cout << "Iteration " << iteration << ": Selected views for hypotheses are " << hyp01_view_indx << " and " << hyp02_view_indx << std::endl;

      CorePipelineOutput result = core_pipeline(Data_Reader, hyp01_view_indx, hyp02_view_indx, All_R, All_T, All_K, K, rd_data, d, q);

      Eigen::MatrixXd Edges_3D = result.Gamma1s;

      // Merge the current iteration's frame-to-edge-to-3D map into the cumulative map
      for (const auto& [frame, edge_to_3D_map] : result.frame_to_edge_to_3D_map) {
          for (const auto& [edge_2D, corresponding_3D_edges] : edge_to_3D_map) {
              // Add to the cumulative map for each frame and each 2D edge
              auto& cumulative_edges = cumulative_frame_to_edge_to_3D_map[frame][edge_2D];
              cumulative_edges.insert(cumulative_edges.end(), corresponding_3D_edges.begin(), corresponding_3D_edges.end());
          }
      }

      if (all_Edges_3D.rows() == 0) {
        all_Edges_3D = Edges_3D;
      } 
      else {
        all_Edges_3D.conservativeResize(all_Edges_3D.rows() + Edges_3D.rows(), 3);
        all_Edges_3D.bottomRows(Edges_3D.rows()) = Edges_3D;
      }

      // Project the 3D edges to each view 
      std::vector<Eigen::MatrixXd> projectedEdgesList;
      for (int i = 0; i < DATASET_NUM_OF_FRAMES; i++) {
          // Project 3D edges to view i
          Eigen::MatrixXd projectedEdges = project3DEdgesToView(all_Edges_3D, All_R[i], All_T[i], K, All_R[hyp01_view_indx], All_T[hyp01_view_indx]);
          projectedEdgesList.push_back(projectedEdges);
      }

      Data_Reader.readEdgelFiles(All_Edgels, Edge_File, rd_data, row_edge, file_idx, d, q, 1);  

      std::vector<std::vector<int>> claimedEdgesList;
      double threshold = 1.0;  


      int num_frames_with_sufficient_edges = 0;
      for (int i = 0; i < projectedEdgesList.size(); ++i) {
          std::vector<int> claimedEdges = findClosestObservedEdges(projectedEdgesList[i], All_Edgels[i], threshold);
          claimedEdgesList.push_back(claimedEdges);
          
          if (All_Edgels[i].rows() > 0 && (claimedEdges.size() >= All_Edgels[i].rows())) {
              num_frames_with_sufficient_edges++;
          }
      }

      sufficient_reconstruct = (num_frames_with_sufficient_edges == DATASET_NUM_OF_FRAMES);


      // Use the selectBestViews function to determine the two frames with the least supported edges
      std::pair<int, int> bestViews = selectBestViews(claimedEdgesList);
      
      // Assign the best views to the hypothesis indices
      hyp01_view_indx = bestViews.first;
      hyp02_view_indx = bestViews.second;

      if(iteration==5){
        
        std::cout<<"printing to file"<<std::endl;
        std::string outputFilePath = "/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/projected_vs_observed_edges.txt";
        std::ofstream edgeFile(outputFilePath);
        if (!edgeFile.is_open()) {
            std::cerr << "Error opening file to write projected and observed edges." << std::endl;
            return -1;
        }

        for (int frameIdx = 0; frameIdx < DATASET_NUM_OF_FRAMES; ++frameIdx) {
            edgeFile << "Frame " << frameIdx + 1 << "\n";
            edgeFile << "Projected Edges:\n";
            for (int i = 0; i < projectedEdgesList[frameIdx].rows(); i++) {
                edgeFile << projectedEdgesList[frameIdx](i, 0) << " " << projectedEdgesList[frameIdx](i, 1) << "\n";
            }
            edgeFile << "Observed Edges:\n";
            for (int i = 0; i < All_Edgels[frameIdx].rows(); i++) {
                edgeFile << All_Edgels[frameIdx](i, 0) << " " << All_Edgels[frameIdx](i, 1) << "\n";
            }
        }
        edgeFile.close();
      }

      All_Edgels.clear();
      projectedEdgesList.clear();
      claimedEdgesList.clear();
      d = 0;
      q = 0;
      file_idx = 0;
      iteration++;
    }

    std::string allEdges3DFilePath = "/gpfs/data/bkimia/zqiwu/3D/3D_Edge_Sketch/outputs/all_edges_3d.txt";
    std::ofstream allEdges3DFile(allEdges3DFilePath);
    if (!allEdges3DFile.is_open()) {
        std::cerr << "Error opening file to write all edges 3D data." << std::endl;
        return -1;
    }

    // Write the content of all_Edges_3D as a nested list
    allEdges3DFile << "[\n"; 
    for (int i = 0; i < all_Edges_3D.rows(); i++) {
        allEdges3DFile << "  ["; 
        allEdges3DFile << all_Edges_3D(i, 0) << ", "
                       << all_Edges_3D(i, 1) << ", "
                       << all_Edges_3D(i, 2);
        allEdges3DFile << "]"; 
        if (i != all_Edges_3D.rows() - 1) {
            allEdges3DFile << ",\n"; 
        } else {
            allEdges3DFile << "\n"; 
        }
    }
    allEdges3DFile << "]\n"; 

    allEdges3DFile.close();

    // for (const auto& [frame, edge_to_3D_map] : cumulative_frame_to_edge_to_3D_map) {
    //     std::cout << "Frame " << frame << ":\n";
    //     for (const auto& [edge_2D, corresponding_3D_edges] : edge_to_3D_map) {
    //         std::cout << "  Edge (" << edge_2D(0) << ", " << edge_2D(1) << ") -> 3D Edges: ";
    //         for (const auto& edge_3D : corresponding_3D_edges) {
    //             std::cout << "(" << edge_3D(0) << ", " << edge_3D(1) << ", " << edge_3D(2) << "), ";
    //         }
    //         std::cout << std::endl;
    //     }
    // }


    LOG_INFOR_MESG("3D Edge Sketch is Finished!");

    return 0;
} 
