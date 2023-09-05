#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <assert.h>
#include <string>
#include <ctime>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

//> Include functions

#include "../Edge_Reconst/util.hpp"
#include "../Edge_Reconst/PairEdgeHypo.hpp"
#include "../Edge_Reconst/getReprojectedEdgel.hpp"
#include "../Edge_Reconst/getQuadrilateral.hpp"
#include "../Edge_Reconst/getSupportedEdgels.hpp"
#include "../Edge_Reconst/definitions.h"

//> Added by CH: Efficient Bucketing Method
#include "../Edge_Reconst/lemsvpe_CH/vgl_polygon_CH.hpp"
#include "../Edge_Reconst/lemsvpe_CH/vgl_polygon_scan_iterator_CH.hpp"
#include "../Edge_Reconst/subpixel_point_set.hpp"

using namespace std;
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
  cout << "Number of interior bucket coordinates: " << InteriorBucketCoordinates.size()<<endl;
  cout<<"bucket coordinates(starts from 0) are shown below: "<< endl;
  for (int bi = 0; bi < InteriorBucketCoordinates.size(); bi++) {
    
    unsigned const i_col = InteriorBucketCoordinates[bi](0);  //> x
    unsigned const i_row = InteriorBucketCoordinates[bi](1);  //> y

    //cout<< "coordinate " << bi << ": "<< i_col<< ", "<< i_row << endl;
    cout<< i_col << ", "<< i_row << ";" << endl;


    //> Ignore if bucket coordinate exceeds image boundary
    if (i_row >= sp_pts.nrows() || i_col >= sp_pts.ncols()) continue;

    //> Traverse all edgels inside the bucket
    for (unsigned k = 0; k < sp_pts.cells()[i_row][i_col].size(); ++k) {
      unsigned const p2_idx = sp_pts.cells()[i_row][i_col][k];
      //cout<< "inlier edge index(starts from 0): " << p2_idx << endl;
      Edgel_Indices.push_back(p2_idx);
    }
    //cout << "sp_pts.cells()[i_row][i_col].size(): " << sp_pts.cells()[i_row][i_col].size() << endl;
  }
  
  cout << "number of edges in this quadrilateral found by bucketing: "<< Edgel_Indices.size() << endl;
}

int main(int argc, char **argv) {
  std::fstream Edge_File;
  std::fstream Rmatrix_File;
  std::fstream Tmatrix_File;
  int d = 0;
  int q = 0;
  int file_idx = 1;
  double rd_data;
  std::vector<Eigen::MatrixXd> All_Edgels;  
  //Eigen::MatrixXd Edgels;
  Eigen::Vector4d row_edge;

  //> All_Bucketed_Imgs stores all "bucketed" images
  //> (A "bucketed image" means that edgels are inserted to the buckets of that image)
  std::vector< subpixel_point_set > All_Bucketed_Imgs;
  cout << "read edges file now\n";
  while(file_idx < DATASET_NUM_OF_FRAMES+1) {
    std::string Edge_File_Path = REPO_DIR + "datasets/cabinet/Edges/Edge_"+std::to_string(file_idx)+".txt";
    file_idx ++;
    Eigen::MatrixXd Edgels; //> Declare locally, ensuring the memory addresses are different for different frames
    Edge_File.open(Edge_File_Path, std::ios_base::in);
    if (!Edge_File) { 
       std::cerr << "Edge file not existed!\n"; exit(1); 
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
      All_Edgels.push_back(Edgels);

      //> Building bucket grid for the current frame
      subpixel_point_set bucketed_img( Edgels );
      bucketed_img.build_bucketing_grid( imgrows, imgcols );
      All_Bucketed_Imgs.push_back( bucketed_img );

      //Edgels = {};
      d = 0;
      q = 0;
    }
  }
  
  cout<< "Edge file loading finished" <<endl;

  std::vector<Eigen::Matrix3d> All_R;
  Eigen::Matrix3d R_matrix;
  Eigen::Vector3d row_R;
  std::string Rmatrix_File_Path = REPO_DIR + "datasets/cabinet/RnT/R_matrix.txt";
  Rmatrix_File.open(Rmatrix_File_Path, std::ios_base::in);
  if (!Rmatrix_File) { 
    std::cerr << "R_matrix file not existed!\n"; exit(1); 
    }
  else {
    while (Rmatrix_File >> rd_data) {
      row_R(q) = rd_data;
      q++;
      if(q>2){
        R_matrix.row(d) = row_R;
        row_R = {};
        q = 0;
        d++;
      }
      if(d>2){
        All_R.push_back(R_matrix);
        R_matrix = {};
        d = 0;
      }
    }
    Rmatrix_File.close();
  }

  cout<< "R matrix loading finished" <<endl;

  std::vector<Eigen::Vector3d> All_T;
  Eigen::Vector3d T_matrix;
  std::string Tmatrix_File_Path = REPO_DIR + "datasets/cabinet/RnT/T_matrix.txt";
  Tmatrix_File.open(Tmatrix_File_Path, std::ios_base::in);
  if (!Tmatrix_File) { 
    std::cerr << "T_matrix file not existed!\n"; exit(1); 
    }
  else {
    while (Tmatrix_File >> rd_data) {
      T_matrix(d) = rd_data;
      d++;
      if(d>2){
        All_T.push_back(T_matrix);
        T_matrix = {};
        d = 0;
      }
    }
    Tmatrix_File.close();
  }
  
  cout<< "T matrix loading finished" <<endl;

  Eigen::Matrix3d K;
  K<< 537.960322000000, 0, 319.183641000000, 0,	539.597659000000,	247.053820000000,0,	0,	1;

  ////////////////////
  // Pipeline start
  ////////////////////

  MultiviewGeometryUtil::multiview_geometry_util util;
  PairEdgeHypothesis::pair_edge_hypothesis       PairHypo;
  GetReprojectedEdgel::get_Reprojected_Edgel     getReprojEdgel;
  GetQuadrilateral::get_Quadrilateral            getQuad;
  GetSupportedEdgels::get_SupportedEdgels        getSupport;
  
  Eigen::MatrixXd Edges_HYPO1 = All_Edgels[HYPO1_VIEW_INDX];
  Eigen::Matrix3d R1          = All_R[HYPO1_VIEW_INDX];
  Eigen::Vector3d T1          = All_T[HYPO1_VIEW_INDX];
  Eigen::MatrixXd Edges_HYPO2 = All_Edgels[HYPO2_VIEW_INDX];
  Eigen::Matrix3d R2          = All_R[HYPO2_VIEW_INDX];
  Eigen::Vector3d T2          = All_T[HYPO2_VIEW_INDX];
  
  Eigen::Matrix3d R21 = util.getRelativePose_R21(R1, R2);
  Eigen::Vector3d T21 = util.getRelativePose_T21(R1, R2, T1, T2);
  Eigen::Matrix3d Tx  = util.getSkewSymmetric(T21);
  Eigen::Matrix3d E   = util.getEssentialMatrix(R21, T21);
  Eigen::Matrix3d F   = util.getFundamentalMatrix(K.inverse(), R21, T21);


  Eigen::Vector3d pt_edgel_HYPO1;

  /////////////////////////////////////////////////
  // should be a for loop of edges in hypo 1 here
  /////////////////////////////////////////////////

  Eigen::MatrixXd paired_edge;
  int pair_num = 0;
  int mod1 = 0;
  int mod2 = 0;
  int mod3 = 1;
  cout<< "pipeline start" <<endl;
  //clock_t tstart, tstart1, tend;
  clock_t tstart, tend;
  tstart = clock();
  bool should_break = false;
  //for(int edge_idx = 0; edge_idx < Edges_HYPO1.rows(); edge_idx++){
  for(int edge_idx = 2001; edge_idx < 2002; edge_idx++){
      //cout<<edge_idx<<endl;
    /*mod1 = (edge_idx+1)%10;
    if( mod1 == 0){
      cout << mod3 << ". "<< flush;
      mod3 ++;
    }
    mod2 = (edge_idx+1)%500;
    if(mod2 == 0){
      mod3 = 1;
      cout<< " "<< edge_idx <<endl;
    }*/

    if(Edges_HYPO1(edge_idx,0) < 10 || Edges_HYPO1(edge_idx,0) > imgcols-10 || Edges_HYPO1(edge_idx,1) < 10 || Edges_HYPO1(edge_idx,1) > imgrows-10){
      continue;
    }
    pt_edgel_HYPO1 << Edges_HYPO1(edge_idx,0), Edges_HYPO1(edge_idx,1), 1;

    Eigen::MatrixXd ApBp = PairHypo.getAp_Bp(Edges_HYPO2, pt_edgel_HYPO1, F);

    Eigen::MatrixXd numerOfDist = PairHypo.getAp_Bp_Dist(Edges_HYPO2, pt_edgel_HYPO1, F);
    Eigen::MatrixXd HYPO2_idx    = PairHypo.getHYPO2_idx(Edges_HYPO2, numerOfDist);
    Eigen::MatrixXd edgels_HYPO2 = PairHypo.getedgels_HYPO2(Edges_HYPO2, numerOfDist);

    ///////////////////////////////////////////////////
    // should be a for loop of validation views here
    ///////////////////////////////////////////////////

    int VALID_idx = 0;
    int stack_idx = 0;
    Eigen::MatrixXd supported_indices;
    supported_indices.conservativeResize(edgels_HYPO2.rows(),48);
    Eigen::MatrixXd supported_indice_current;
    supported_indice_current.conservativeResize(edgels_HYPO2.rows(),1);
    Eigen::MatrixXd supported_indices_stack;
    
    // cout<< "run here 1" << endl;
    bool isempty = true;
    //tstart = clock();
    for (int VALID_INDX = 38; VALID_INDX < DATASET_NUM_OF_FRAMES; VALID_INDX++){
      if(VALID_INDX == HYPO1_VIEW_INDX || VALID_INDX == HYPO2_VIEW_INDX){
        continue;
      }
      
      Eigen::MatrixXd TO_Edges_VALID = All_Edgels[VALID_INDX];
      Eigen::Matrix3d R3             = All_R[VALID_INDX];
      Eigen::Vector3d T3             = All_T[VALID_INDX];
      Eigen::MatrixXd VALI_Orient    = TO_Edges_VALID.col(2);
      Eigen::MatrixXd Tangents_VALID;
      Tangents_VALID.conservativeResize(TO_Edges_VALID.rows(),2);
      Tangents_VALID.col(0)          = (VALI_Orient.array()).cos();
      Tangents_VALID.col(1)          = (VALI_Orient.array()).sin();
      
      Eigen::Matrix3d R31 = util.getRelativePose_R21(R1, R3);
      Eigen::Vector3d T31 = util.getRelativePose_T21(R1, R3, T1, T3);
      
      Eigen::MatrixXd pt_edge = Edges_HYPO1.row(edge_idx);
      Eigen::Vector3d tgt1_meters = getReprojEdgel.getTGT_Meters(pt_edge, K);
      
      Eigen::MatrixXd edge_pos_gamma3 = getReprojEdgel.getGamma3Pos(pt_edge, edgels_HYPO2, All_R, All_T, VALID_INDX, K);
      Eigen::MatrixXd edge_tgt_gamma3 = getReprojEdgel.getGamma3Tgt(pt_edge, edgels_HYPO2, All_R, All_T, VALID_INDX, K);
      
      // Eigen::MatrixXd QuadrilateralPoints = getQuad.getQuadrilateralPoints(pt_edge, edgels_HYPO2.row(20), All_R, All_T, VALID_INDX, K);

      
      //>>>>>>>>>>>>>> START OF FETCHING EDGEL IDS FROM A QUADRILATERAL >>>>>>>>>>>>>>

      Eigen::MatrixXd QuadrilateralPoints = getQuad.getQuadrilateralPoints(pt_edge, edgels_HYPO2.row(47), All_R, All_T, VALID_INDX, K);
      cout<< "QuadrilateralPoints: " << endl;
      cout<< QuadrilateralPoints << endl;

      //std::cout << "================================================" << std::endl;
      //std::cout << round(QuadrilateralPoints(0,0)) << "\t" << round(QuadrilateralPoints(0,1)) << std::endl;
      //std::cout << round(QuadrilateralPoints(1,0)) << "\t" << round(QuadrilateralPoints(1,1)) << std::endl;
      //std::cout << round(QuadrilateralPoints(2,0)) << "\t" << round(QuadrilateralPoints(2,1)) << std::endl;
      //std::cout << round(QuadrilateralPoints(3,0)) << "\t" << round(QuadrilateralPoints(3,1)) << std::endl;

      vgl_polygon_CH<double> Quadrilateral_Corner_Pts_vgl(1);
      Quadrilateral_Corner_Pts_vgl.push_back_(round(QuadrilateralPoints(0,0)), round(QuadrilateralPoints(0,1)));
      Quadrilateral_Corner_Pts_vgl.push_back_(round(QuadrilateralPoints(1,0)), round(QuadrilateralPoints(1,1)));
      Quadrilateral_Corner_Pts_vgl.push_back_(round(QuadrilateralPoints(2,0)), round(QuadrilateralPoints(2,1)));
      Quadrilateral_Corner_Pts_vgl.push_back_(round(QuadrilateralPoints(3,0)), round(QuadrilateralPoints(3,1)));

      //> Get bucket coordinates inside the quadrilateral
      std::vector<Eigen::Vector2d> InteriorBucketCoordinates;
      getInteriorBuckets(Quadrilateral_Corner_Pts_vgl, true, InteriorBucketCoordinates);

      //> Fetch edgel IDs of the corresponding buckets
      std::vector< unsigned > Edgel_Indices;
      getEdgelsFromInteriorQuadrilateral( All_Bucketed_Imgs[VALID_INDX], InteriorBucketCoordinates, Edgel_Indices);

      //> Print out a list of edgel IDs and their subpixel coordinates on the validation view locating inside the quadrilateral, if there are any
      //> Edgel ID: Edgel_Indices[ei]
      //> Edgel coordinate (x, y) = (TO_Edges_VALID(Edgel_Indices[ei], 0), TO_Edges_VALID(Edgel_Indices[ei], 1))
      if (Edgel_Indices.size() > 0) {
        cout << "inlier edges are shown below: " << endl;
        cout << "index: (x,y), index starts from 0" << endl;
        for (int ei = 0; ei < Edgel_Indices.size(); ei++) {
          cout << Edgel_Indices[ei]+1 << ";"<<endl;
          //std::cout << Edgel_Indices[ei] << ": (";
          //std::cout << TO_Edges_VALID(Edgel_Indices[ei], 0) << ", " << TO_Edges_VALID(Edgel_Indices[ei], 1) << ")" << std::endl;
        }
        std::cout << std::endl;
        should_break = true;
      }

      //std::cout << "================================================" << std::endl;
      //>>>>>>>>>>>>>> END OF FETCHING EDGEL IDS FROM A QUADRILATERAL >>>>>>>>>>>>>>
      
      
      //tstart1 = clock();
      //for (int idx_pair = 0; idx_pair < edgels_HYPO2.rows(); idx_pair++){
        int idx_pair = 47;
        Eigen::MatrixXd inliner = getQuad.getInliner(pt_edge, edgels_HYPO2.row(idx_pair), All_R, All_T, VALID_INDX, K, TO_Edges_VALID);
        Eigen::Vector2d edgels_tgt_reproj = {edge_tgt_gamma3(idx_pair,0), edge_tgt_gamma3(idx_pair,1)};
        double supported_link_indx = getSupport.getSupportIdx(edgels_tgt_reproj, Tangents_VALID, inliner);
        supported_indice_current.row(idx_pair) << supported_link_indx;
        if (supported_link_indx != -2){
          supported_indices_stack.conservativeResize(stack_idx+1,2);
          supported_indices_stack.row(stack_idx) << double(idx_pair), double(supported_link_indx);
          isempty = false;
          stack_idx++;
        }
        cout << "Number of inlier found using old method: " << inliner.size() << endl;

        for(int oldi = 0; oldi < inliner.size(); oldi++){
          std::cout << inliner(oldi) << ": (";
          std::cout << TO_Edges_VALID(inliner(oldi), 0) << ", " << TO_Edges_VALID(inliner(oldi), 1) << ")" << std::endl;
        }
        //cout << inliner << endl;
        should_break = true;
      if (should_break) break;
      //}
      //tend = clock() - tstart1; 
      //cout << "It took "<< double(tend)/double(CLOCKS_PER_SEC) <<" second(s) to get support from one validation view."<< endl;
      supported_indices.col(VALID_idx) << supported_indice_current.col(0);
      VALID_idx++;
      
    }

    //if (should_break) break; //> CH: Feel free to remove it!

    //tend = clock() - tstart; 
    //cout << "It took "<< double(tend)/double(CLOCKS_PER_SEC) <<" second(s) to get support from validation views."<< endl;
    if(isempty){
      continue;
    }
    // cout<< "run here 2" << endl;
    //cout<< VALID_idx << endl;
    //cout << "supported_indices.col(0)" << endl;
    //cout << supported_indices.col(0) << endl;
    //cout << "supported_indices_stack" << endl;
    //cout << supported_indices_stack.block(0,0,50,2) << endl;
    std::vector<double> indices_stack(supported_indices_stack.data(), supported_indices_stack.data() + supported_indices_stack.rows());
    std::vector<double> indices_stack_unique = indices_stack;
    std::sort(indices_stack_unique.begin(), indices_stack_unique.end());
    std::vector<double>::iterator it1;
    it1 = std::unique(indices_stack_unique.begin(), indices_stack_unique.end());
    indices_stack_unique.resize( std::distance(indices_stack_unique.begin(),it1) );
    //cout << "supported_indices_stack" << endl;
    //cout << indices_stack_unique.size() << endl;
    //std::vector<double>::iterator it2;
    Eigen::VectorXd rep_count;
    rep_count.conservativeResize(indices_stack_unique.size(),1);
    // cout<< "run here 3" << endl;
    for(int unique_idx = 0; unique_idx<indices_stack_unique.size(); unique_idx++){
      rep_count.row(unique_idx) << double(count(indices_stack.begin(), indices_stack.end(), indices_stack_unique[unique_idx]));
    }
    // cout<< "run here 4" << endl;
    // cout<< rep_count << endl;
    Eigen::VectorXd::Index   maxIndex;
    double max_support = rep_count.maxCoeff(&maxIndex);
    int numofmax = count(rep_count.data(), rep_count.data()+rep_count.size(), max_support);
    //cout << rep_count.row(maxIndex) << endl;
    // cout<< "run here 5" << endl;
    if( double(max_support) < MAX_NUM_OF_SUPPORT_VIEWS){
      // cout << max_support << endl;
      continue;
    }
    int finalpair = -2;
    if(numofmax > 1){
      std::vector<double> rep_count_vec(rep_count.data(), rep_count.data() + rep_count.rows());
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
      Eigen::Vector3d coeffs;
      coeffs = F * pt_edgel_HYPO1;
      Eigen::MatrixXd Edge_Pts;
      Edge_Pts.conservativeResize(max_index.size(),2);
      for(int maxidx = 0; maxidx<max_index.size(); maxidx++){
        Edge_Pts.row(maxidx) << edgels_HYPO2(indices_stack_unique[max_index[maxidx]], 0),edgels_HYPO2(indices_stack_unique[max_index[maxidx]], 1) ;
      }
      Eigen::VectorXd Ap = coeffs(0)*Edge_Pts.col(0);
      Eigen::VectorXd Bp = coeffs(1)*Edge_Pts.col(1);
      Eigen::VectorXd numDist = Ap + Bp + Eigen::VectorXd::Ones(Ap.rows())*coeffs(2);
      double denomDist = coeffs(0)*coeffs(0) + coeffs(1)*coeffs(1);
      denomDist = sqrt(denomDist);
      Eigen::VectorXd dist = numDist.cwiseAbs()/denomDist;
      //cout << dist << endl;
      Eigen::VectorXd::Index   minIndex;
      double min_dist = dist.minCoeff(&minIndex);
      if(min_dist > DIST_THRESH){
        continue;
      }
      finalpair = int(indices_stack_unique[max_index[minIndex]]);
      // cout << finalpair << endl;
    }else{
      finalpair = int(indices_stack_unique[int(maxIndex)]);
      // cout << finalpair << endl;
    }
    // linearTriangulation code already exist
    paired_edge.conservativeResize(pair_num+1,50);
    paired_edge.row(pair_num) << edge_idx, HYPO2_idx(finalpair), supported_indices.row(finalpair);
    pair_num++;
  }
  cout<< "pipeline finished" <<endl;
  tend = clock() - tstart; 
  cout << "It took "<< double(tend)/double(CLOCKS_PER_SEC) <<" second(s) to find a pair."<< endl;
  cout << "Number of pairs found: " << paired_edge.rows()<<endl;
  ofstream myfile1;
  std::string Output_File_Path = OUTPUT_WRITE_FOLDER + "pairededge6n3_quadsize2.txt";
  myfile1.open (Output_File_Path);
  myfile1 << paired_edge;
  myfile1.close();
}
