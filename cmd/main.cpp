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
//    Chiang-Heng Chien  23-07-18    Initially create a blank repository with minor multiview geometry utility functions.
//    Chiang-Heng Chien  23-08-19    Add edited bucketing method from VXL. If this repo is to be integrated under LEMSVPE 
//                                   scheme, simply use vgl library from VXL. This is very easy to be handled.
//    Chiang-Heng Chien  23-08-27    Add bucket building and fetching edgels from bucket coordinate code for efficient edgel
//                                   accessing from buckets inside a given quadrilateral.
//
//> (c) LEMS, Brown University
//> Yilin Zheng (yilin_zheng@brown.edu)
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =========================================================================================================================
void printEdge3DToHypothesisAndSupports(
    const std::unordered_map<Eigen::Matrix<double, 3, 1>, 
                             std::tuple<Eigen::Vector2d, Eigen::Vector2d, std::vector<std::pair<int, Eigen::Vector2d>>>, 
                             EigenMatrixHash>& edge_3D_to_hypothesis_and_supports) 
{
    for (const auto& [edge_3D, value] : edge_3D_to_hypothesis_and_supports) {
        // Extract 3D edge
        std::cout << "3D Edge: [" << edge_3D(0) << ", " << edge_3D(1) << ", " << edge_3D(2) << "]\n";

        // Extract the hypothesis edges from view 6 and view 8
        const Eigen::Vector2d& hypothesis_edge_view6 = std::get<0>(value);
        const Eigen::Vector2d& hypothesis_edge_view8 = std::get<1>(value);

        std::cout << "    Hypothesis Edge (View 6): [" << hypothesis_edge_view6(0) << ", " << hypothesis_edge_view6(1) << "]\n";
        std::cout << "    Hypothesis Edge (View 8): [" << hypothesis_edge_view8(0) << ", " << hypothesis_edge_view8(1) << "]\n";

        // Extract the supporting edges with validation view numbers
        const std::vector<std::pair<int, Eigen::Vector2d>>& validation_support_edges = std::get<2>(value);

        // Loop through and print each supporting edge and its corresponding validation view number
        for (size_t i = 0; i < validation_support_edges.size(); ++i) {
            const auto& [val_view_num, support_edge] = validation_support_edges[i];
            std::cout << "    Validation View " << val_view_num << ": Supporting Edge = [" 
                      << support_edge(0) << ", " << support_edge(1) << "]\n";
        }
    }
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


void getInteriorBuckets(
  const vgl_polygon_CH<double> &p, bool boundary_in, 
  std::vector<Eigen::Vector2d> &InteriorBucketCoordinates
)
{
  // iterate
  vgl_polygon_scan_iterator_CH<double> it(p, boundary_in); 

  //std::cout << "Interior points:\n";
  for (it.reset(); it.next(); ) {
    int y = it.scany();
    for (int x = it.startx(); x <= it.endx(); ++x) {
      //std::cout << "(" << x << "," << y << ") ";
      Eigen::Vector2d Bucket_Coordinate;
      Bucket_Coordinate << x, y;
      InteriorBucketCoordinates.push_back( Bucket_Coordinate );
    }
  }
  //std::cout << std::endl;
}

void getEdgelsFromInteriorQuadrilateral( 
  const subpixel_point_set &sp_pts, 
  const std::vector<Eigen::Vector2d> &InteriorBucketCoordinates,
  std::vector< unsigned > &Edgel_Indices 
)
{
  //> Traverse all buckets inside the quadrilateral
  //std::cout << "Number of interior bucket coordinates: " << InteriorBucketCoordinates.size()<<std::endl;
  //std::cout<<"bucket coordinates(starts from 0) are shown below: "<<std::endl;
  for (int bi = 0; bi < InteriorBucketCoordinates.size(); bi++) {
    
    unsigned const i_col = InteriorBucketCoordinates[bi](0);  //> x
    unsigned const i_row = InteriorBucketCoordinates[bi](1);  //> y

    //std::cout<< "coordinate " << bi << ": "<< i_col<< ", "<< i_row <<std::endl;
    //std::cout<< i_col << ", "<< i_row << ";" <<std::endl;


    //> Ignore if bucket coordinate exceeds image boundary
    if (i_row >= sp_pts.nrows() || i_col >= sp_pts.ncols()) continue;

    //> Traverse all edgels inside the bucket
    for (unsigned k = 0; k < sp_pts.cells()[i_row][i_col].size(); ++k) {
      unsigned const p2_idx = sp_pts.cells()[i_row][i_col][k];
      //std::cout<< "inlier edge index(starts from 0): " << p2_idx <<std::endl;
      Edgel_Indices.push_back(p2_idx);
    }
    //std::cout << "sp_pts.cells()[i_row][i_col].size(): " << sp_pts.cells()[i_row][i_col].size() <<std::endl;
  }
  //std::cout << "number of edges in this quadrilateral found by bucketing: "<< Edgel_Indices.size() <<std::endl;
}




Eigen::MatrixXd core_pipeline(
  const int hyp01_view_indx, 
  const int hyp02_view_indx,
  const std::vector<Eigen::Matrix3d>& All_R, 
  const std::vector<Eigen::Vector3d>& All_T,
  const std::vector<Eigen::Matrix3d>& All_K,
  const Eigen::Matrix3d K,
  double rd_data,
  int d,
  int q){


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
    readEdgelFiles(All_Edgels, Edge_File, rd_data, row_edge, file_idx, d, q, thresh_EDG);  
    readHypothesisEdgelFiles(hyp01_view_indx, hyp02_view_indx, All_Edgels_H12, Edge_File, rd_data, row_edge, H_idx, file_idx, d, q, thresh_EDG);
  
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
    Eigen::MatrixXd Edges_pairs = All_Edgels[hyp01_view_indx];
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
    Eigen::Matrix3d F   = util.getFundamentalMatrix(K1.inverse(), K2.inverse(), R21, T21);
    Eigen::Matrix3d R12 = util.getRelativePose_R21(R2, R1);
    Eigen::Vector3d T12 = util.getRelativePose_T21(R2, R1, T2, T1);  
    Eigen::Matrix3d F12 = util.getFundamentalMatrix(K2.inverse(), K1.inverse(), R12, T12);  
    // Initializations for paired edges between Hypo1 and Hypo 2
    Eigen::MatrixXd paired_edge = Eigen::MatrixXd::Constant(Edges_pairs.rows(),50,-2);
    // Compute epipolar wedges angles rance between Hypo1 and Hypo2 and find the angle range 1, defining a valid range for matching edges in Hypo2
    Eigen::MatrixXd OreListdegree    = getOre.getOreList(hyp01_view_indx, hyp02_view_indx, Edges_HYPO2, All_R, All_T, K1, K2);
    Eigen::MatrixXd OreListBardegree = getOre.getOreListBar(Edges_HYPO1, All_R, All_T, K1, K2, hyp02_view_indx, hyp01_view_indx);
    
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< CORE PIPELINE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
  
    std::cout<< "Threshold is : " << thresh_EDG << std::endl;
   
    unsigned nthreads = NUM_OF_OPENMP_THREADS;
    omp_set_num_threads(nthreads);
    int ID = omp_get_thread_num();
    itime = omp_get_wtime();
    std::cout << "Using " << nthreads << " threads for OpenMP parallelization." << std::endl;
    std::cout << "nthreads: " << nthreads << "." << std::endl;

  #pragma omp parallel
  {
    //> Local array stacking all supported indices
    std::vector<Eigen::MatrixXd> local_thread_supported_indices;

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< First loop: loop over all edgels from hypothesis view 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
    //<<<<<<<<<<< Identify pairs of edge, correct the positions of the edges from Hypo2, and store the paired edges >>>>>>>>>>>>>>>>//
    #pragma omp for schedule(static, nthreads)
    for (int edge_idx = 0; edge_idx < Edges_HYPO1.rows() ; edge_idx++) {

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
      Eigen::MatrixXd edgels_HYPO2_corrected = PairHypo.edgelsHYPO2correct(edgels_HYPO2, edgel_HYPO1, F, F12, HYPO2_idx_raw);
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
        Eigen::MatrixXd OreListBardegree31 = getOre.getOreListBar(pt_edge, All_R, All_T, K1, K3, VALID_INDX, hyp01_view_indx);
        Eigen::MatrixXd OreListdegree31    = getOre.getOreListVali(TO_Edges_VALID, All_R, All_T, K1, K3, VALID_INDX, hyp01_view_indx);

        // Find all the edges fall inside epipolar wedge on validation view (Hypo1 --> Vali)
        Eigen::MatrixXd OreListBardegree32 = getOre.getOreListBar(Edges_HYPO2_final, All_R, All_T, K2, K3, VALID_INDX, hyp02_view_indx);
        Eigen::MatrixXd OreListdegree32    = getOre.getOreListVali(TO_Edges_VALID, All_R, All_T, K2, K3, VALID_INDX, hyp02_view_indx);
        Eigen::VectorXd isparallel         = Eigen::VectorXd::Ones(Edges_HYPO2_final.rows());

        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Third loop: loop over each edge from Hypo2 <<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>//
        for (int idx_pair = 0; idx_pair < Edges_HYPO2_final.rows(); idx_pair++) {
          double thresh_ore31_1 = OreListBardegree31(idx_pair,0);
          double thresh_ore31_2 = OreListBardegree31(idx_pair,1);
          double thresh_ore32_1 = OreListBardegree32(idx_pair,0);
          double thresh_ore32_2 = OreListBardegree32(idx_pair,1);
          
          Eigen::MatrixXd vali_idx31 = PairHypo.getHYPO2_idx_Ore(OreListdegree31, thresh_ore31_1, thresh_ore31_2);
          Eigen::MatrixXd edgels_31  = PairHypo.getedgels_HYPO2_Ore(TO_Edges_VALID, OreListdegree31, thresh_ore31_1, thresh_ore31_2);

          // Find all the edges fall inside epipolar wedge on validation view (Hypo2 --> Vali)
          Eigen::MatrixXd vali_idx32 = PairHypo.getHYPO2_idx_Ore(OreListdegree32, thresh_ore32_1, thresh_ore32_2);
          Eigen::MatrixXd edgels_32  = PairHypo.getedgels_HYPO2_Ore(TO_Edges_VALID, OreListdegree32, thresh_ore32_1, thresh_ore32_2);
          
          // Check if the two wedges could be considered as parallel to each other
          Eigen::MatrixXd anglediff(4,1);
          anglediff << fabs(thresh_ore31_1 - thresh_ore32_1), 
                      fabs(thresh_ore31_1 - thresh_ore32_2),
                      fabs(thresh_ore31_2 - thresh_ore32_1),
                      fabs(thresh_ore31_2 - thresh_ore32_2);
          if(anglediff.maxCoeff() <= parallelangle){
            isparallel.row(idx_pair) << 0;
          }

          // Find all the edges fall inside the two epipolar wedges intersection on validation view 
          // (Hypo1 --> Vali) && (Hypo2 --> Vali)
          std::vector<double> v_intersection;
          std::vector<double> v1(vali_idx31.data(), vali_idx31.data() + vali_idx31.rows());
          std::vector<double> v2(vali_idx32.data(), vali_idx32.data() + vali_idx32.rows());
          set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v_intersection));
          Eigen::VectorXd idxVector = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v_intersection.data(), v_intersection.size());
          Eigen::MatrixXd inliner(idxVector);
          
          // Calculate orientation of gamma 3
          Eigen::Vector2d edgels_tgt_reproj = {edge_tgt_gamma3(idx_pair,0), edge_tgt_gamma3(idx_pair,1)};
          // Get support from validation view for this gamma 3
          double supported_link_indx = getSupport.getSupportIdx(edgels_tgt_reproj, Tangents_VALID, inliner);

          // Get the supporting edge idx from this validation view (if isnotparallel)
          if (isparallel(idx_pair,0) != 0){
            supported_indice_current.row(idx_pair) << supported_link_indx;
          }else{
            supported_indice_current.row(idx_pair) << -2;
          }
          if (supported_link_indx != -2 && isparallel(idx_pair,0) != 0){
            supported_indices_stack.conservativeResize(stack_idx+1,2);
            supported_indices_stack.row(stack_idx) << double(idx_pair), double(supported_link_indx);
            isempty_link = false;
            stack_idx++;
          }
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
      if( double(max_support) < MAX_NUM_OF_SUPPORT_VIEWS){
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
        coeffs = F * pt_edgel_HYPO1;
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
        //std::cout << dist <<std::endl;
        Eigen::VectorXd::Index   minIndex;
        double min_dist = dist.minCoeff(&minIndex);
        // std::cout << "min_dist: "<< min_dist <<std::endl;
        if(min_dist > DIST_THRESH){
          continue;
        }
        finalpair = int(indices_stack_unique[max_index[minIndex]]);
      }
      else {
        finalpair = int(indices_stack_unique[int(maxIndex)]);
        // std::cout << "single: " << finalpair <<std::endl;
      }
      // linearTriangulation code already exist
      paired_edge.row(edge_idx) << edge_idx, HYPO2_idx(finalpair), supported_indices.row(finalpair);
      //std::cout << "paired_edge.row(edge_idx): \n" << paired_edge.row(edge_idx) <<std::endl;
    }
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  End of first loop >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//

    //> A critical session to stack all local supported indices
    #pragma omp critical
    all_supported_indices.insert(all_supported_indices.end(), local_thread_supported_indices.begin(), local_thread_supported_indices.end());

    //OpenMP Parallelization Time Reporting
    // #if defined(_OPENMP)
      ftime = omp_get_wtime();
      exec_time = ftime - itime;
      totaltime += exec_time;
      std::cout << "It took "<< exec_time <<" second(s) to finish this round."<<std::endl;
      std::cout << "End of using OpenMP parallelization." << std::endl;
    // #endif
  } //> End of pragma omp parallel
    
    //> CH: Make pair_edge locally, and merge them to a global variable once the for loop is finished.
    // .....
    //Finalize Paired Edges
    tstart = clock();
    int pair_num = 0;
    Eigen::MatrixXd paired_edge_final;
    for(int pair_idx = 0; pair_idx < paired_edge.rows(); pair_idx++){
      if(paired_edge(pair_idx,0) != -2 && paired_edge(pair_idx,0) != -3){
        paired_edge_final.conservativeResize(pair_num+1,50);
        paired_edge_final.row(pair_num) << paired_edge.row(pair_idx);
        pair_num++;
      }
    }
    //std::cout << "number of pairs is: "<< pair_num << std::endl;
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
        Eigen::MatrixXd edgel_HYPO1   = Edges_HYPO1_all.row(int(paired_edge_final(pair_idx,0)));
        Eigen::MatrixXd edgels_HYPO2  = Edges_HYPO2_all.row(int(paired_edge_final(pair_idx,1)));
        Eigen::MatrixXd HYPO2_idx_raw = Edges_HYPO1.row(int(paired_edge_final(pair_idx,0)));


        Eigen::MatrixXd edgels_HYPO2_corrected = PairHypo.edgelsHYPO2correct(edgels_HYPO2, edgel_HYPO1, F, F12, HYPO2_idx_raw);

        if (HYPO2_idx_raw.rows() == 0 || edgels_HYPO2_corrected.rows() == 0) {
          std::cout << "No valid matches found for edge " << pair_idx << " at threshold " << thresh_EDG << std::endl;
          continue;
        }

        Eigen::MatrixXd Edges_HYPO1_final(edgels_HYPO2_corrected.rows(),4);
        Edges_HYPO1_final << edgels_HYPO2_corrected.col(0), edgels_HYPO2_corrected.col(1), edgels_HYPO2_corrected.col(2), edgels_HYPO2_corrected.col(3);
        Eigen::MatrixXd Edges_HYPO2_final(edgels_HYPO2_corrected.rows(),4);
        Edges_HYPO2_final << edgels_HYPO2_corrected.col(4), edgels_HYPO2_corrected.col(5), edgels_HYPO2_corrected.col(6), edgels_HYPO2_corrected.col(7);
        
        Eigen::Vector2d pt_H1 = Edges_HYPO1_final.row(0);
        Eigen::Vector2d pt_H2 = Edges_HYPO2_final.row(0);
        std::vector<Eigen::Vector2d> pts;
        pts.push_back(pt_H1);
        pts.push_back(pt_H2);

        Eigen::Vector3d edge_pt_3D = linearTriangulation(2, pts, Rs,Ts,K1_v);

        if (edge_pt_3D.hasNaN()) {
          std::cerr << "Error: NaN values detected in edge_pt_3D for pair_idx: " << pair_idx << std::endl;
          continue;
        }

        Gamma1s.row(pair_idx)<< edge_pt_3D(0),edge_pt_3D(1),edge_pt_3D(2);
        edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, pt_H1, hyp01_view_indx);
        edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, pt_H2, hyp02_view_indx);

        ///////////////////////////////// Add support from validation views /////////////////////////////////
        std::vector<std::pair<int, Eigen::Vector2d>> validation_support_edges;

        // Loop through validation views to find the supporting edges
        for (int val_idx = 0; val_idx < DATASET_NUM_OF_FRAMES; ++val_idx) {
            if (val_idx == hyp01_view_indx || val_idx == hyp02_view_indx) {
                continue;  // Skip hypothesis views
            }

            // Retrieve support index from paired_edge for the current validation view
            int support_idx = paired_edge_final(pair_idx, val_idx + 2);  // +2 accounts for the first two columns for HYPO1 and HYPO2
            if (support_idx != -2) {
                // Retrieve the supporting edge from the validation view
                Eigen::MatrixXd edges_for_val_frame = All_Edgels[val_idx];
                Eigen::Vector2d supporting_edge = edges_for_val_frame.row(support_idx).head<2>();

                // Store validation view and the supporting edge
                validation_support_edges.emplace_back(val_idx, supporting_edge);

                // Add the supporting edge to the edgeMapping for the 3D edge
                edgeMapping.add3DToSupportingEdgesMapping(edge_pt_3D, supporting_edge, val_idx);
            }
        }
        ///////////////////////////////// Add support from validation views /////////////////////////////////
      }
    }

    tend = clock() - tstart; 
    std::cout << "It took "<< double(tend)/double(CLOCKS_PER_SEC) <<" second(s) to generate the final edge pair for output file."<<std::endl;
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
      std::cout<< "pipeline finished" <<std::endl;
      std::cout << "It took "<< totaltime <<" second(s) to finish the whole pipeline."<<std::endl;
      std::ofstream myfile2;
      std::string Output_File_Path2 = "../../outputs/Gamma1s_" + DATASET_NAME + "_" + SCENE_NAME + "_" + std::to_string(hyp01_view_indx)+"n"+std::to_string(hyp02_view_indx)+"_t32to"+std::to_string(thresh_EDG) + "_delta" + deltastr +"_theta" + std::to_string(OREN_THRESH) + "_N" + std::to_string(MAX_NUM_OF_SUPPORT_VIEWS) + ".txt";
      std::cout << Output_File_Path2 << std::endl;
      myfile2.open (Output_File_Path2);
      myfile2 << Gamma1s;
      myfile2.close();
    }

    // After all operations in the main function, add the following to print the structure

    // std::cout << "Printing all 3D edges and their corresponding supporting edges:\n";

    // for (const auto& [edge_3D, supporting_edges] : edgeMapping.edge_3D_to_supporting_edges) {
    //     // Print the 3D edge
    //     std::cout << "3D Edge: [" << edge_3D(0) << ", " << edge_3D(1) << ", " << edge_3D(2) << "]\n";

    //     // Print all supporting 2D edges and their image numbers
    //     for (const auto& [supporting_edge, image_number] : supporting_edges) {
    //         std::cout << "    Supporting Edge: [" << supporting_edge(0) << ", " << supporting_edge(1) 
    //                   << "] from Image " << image_number << "\n";
    //     }
    // }
    return Gamma1s;
  }
}



int main(int argc, char **argv) {

    int hyp01_view_indx  = 6;
    int hyp02_view_indx  = 8;

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

    // Read all required matrices (rotation, translation, and camera matrices)
    double rd_data;
    int d = 0, q = 0;
    int file_idx = 0;

    readRmatrix(All_R, R_matrix, Rmatrix_File, rd_data, row_R, d, q);
    readTmatrix(All_T, T_matrix, Tmatrix_File, rd_data, d, q);
    readK(Kmatrix_File, All_K, K, K_matrix, row_K, rd_data, d, q);

    for (int iteration = 0; iteration < 1; iteration++) {
        // Call the core pipeline function to reconstruct 3D edges
        Eigen::MatrixXd Edges_3D = core_pipeline(hyp01_view_indx, hyp02_view_indx, All_R, All_T, All_K, K, rd_data, d, q);

        // Project the 3D edges to each view 
        std::vector<Eigen::MatrixXd> projectedEdgesList;
        for (int i = 0; i < DATASET_NUM_OF_FRAMES; i++) {
            // Project 3D edges to view i
            Eigen::MatrixXd projectedEdges = project3DEdgesToView(Edges_3D, All_R[i], All_T[i], K, All_R[hyp01_view_indx], All_T[hyp01_view_indx]);
            projectedEdgesList.push_back(projectedEdges);
        }

        readEdgelFiles(All_Edgels, Edge_File, rd_data, row_edge, file_idx, d, q, 1);  

        std::vector<std::vector<int>> claimedEdgesList;
        double threshold = 2.0;  

        // Find the claimed edges for each frame
        for (int i = 0; i < projectedEdgesList.size(); ++i) {
            std::vector<int> claimedEdges = findClosestObservedEdges(projectedEdgesList[i], All_Edgels[i], threshold);
            claimedEdgesList.push_back(claimedEdges);
        }

        // Use the selectBestViews function to determine the two frames with the least supported edges
        std::pair<int, int> bestViews = selectBestViews(claimedEdgesList);
        
        // Assign the best views to the hypothesis indices
        hyp01_view_indx = bestViews.first;
        hyp02_view_indx = bestViews.second;

        std::cout << "Iteration " << iteration+1 << ": Selected views for hypotheses are " << hyp01_view_indx << " and " << hyp02_view_indx << std::endl;
    }

    return 0;
}
