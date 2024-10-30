#ifndef EDGESKETCH_CORE_CPP
#define EDGESKETCH_CORE_CPP
// =============================================================================================================================
//
// ChangLogs
//    
//
//> (c) LEMS, Brown University
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =============================================================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <vector>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

//> YAML file data reader
#include <yaml-cpp/yaml.h>

#include "definitions.hpp"

//> Constructor
EdgeSketch_Core::EdgeSketch_Core(YAML::Node Edge_Sketch_Setting_File)
    : Edge_Sketch_Setting_YAML_File(Edge_Sketch_Setting_File)
{
    //> Parse data from the YAML file
    //> (1) 3D Edge Sketch Settings
    Num_Of_OMP_Threads                  = Edge_Sketch_Setting_YAML_File["Num_Of_OMP_Threads"].as<int>();
    hyp01_view_indx                     = Edge_Sketch_Setting_YAML_File["Init_Hypo1_View_Index"].as<int>();
    hyp02_view_indx                     = Edge_Sketch_Setting_YAML_File["Init_Hypo2_View_Index"].as<int>();
    Edge_Loc_Pertubation                = Edge_Sketch_Setting_YAML_File["delta"].as<double>();
    Orien_Thresh                        = Edge_Sketch_Setting_YAML_File["delta_theta"].as<double>();
    Max_Num_Of_Support_Views            = Edge_Sketch_Setting_YAML_File["Max_Num_Of_Support_Views"].as<int>();
    Edge_Detection_Init_Thresh          = Edge_Sketch_Setting_YAML_File["Multi_Thresh_Init_Thresh"].as<int>();
    Edge_Detection_Final_Thresh         = Edge_Sketch_Setting_YAML_File["Multi_Thresh_Final_Thresh"].as<int>();
    Parallel_Epipolar_Line_Angle_Deg    = Edge_Sketch_Setting_YAML_File["Parallel_Epipolar_Line_Angle"].as<double>();
    Reproj_Dist_Thresh                  = Edge_Sketch_Setting_YAML_File["Reproj_Dist_Thresh"].as<double>();
    circleR                             = Edge_Sketch_Setting_YAML_File["circleR"].as<double>(); //> Unknown setting
    //> (2) Dataset Settings
    Dataset_Path                        = Edge_Sketch_Setting_YAML_File["Dataset_Path"].as<std::string>();
    Dataset_Name                        = Edge_Sketch_Setting_YAML_File["Dataset_Name"].as<std::string>();
    Scene_Name                          = Edge_Sketch_Setting_YAML_File["Scene_Name"].as<std::string>();
    Num_Of_Total_Imgs                   = Edge_Sketch_Setting_YAML_File["Total_Num_Of_Images"].as<int>();
    Img_Rows                            = Edge_Sketch_Setting_YAML_File["Image_Rows"].as<int>();
    Img_Cols                            = Edge_Sketch_Setting_YAML_File["Image_Cols"].as<int>();
    Use_Multiple_K                      = Edge_Sketch_Setting_YAML_File["Use_Multiple_K"].as<bool>();
    fx                                  = Edge_Sketch_Setting_YAML_File["fx"].as<double>();
    fy                                  = Edge_Sketch_Setting_YAML_File["fy"].as<double>();
    cx                                  = Edge_Sketch_Setting_YAML_File["cx"].as<double>();
    cy                                  = Edge_Sketch_Setting_YAML_File["cy"].as<double>();
    //> (3) Other Settings
    Delta_FileName_Str                  = Edge_Sketch_Setting_YAML_File["deltastr"].as<std::string>();

    //> Initialization
    edge_sketch_time = 0.0;

    //> Class objects
    util = std::shared_ptr<MultiviewGeometryUtil::multiview_geometry_util>(new multiview_geometry_util());

    //> Set up OpenMP threads
    omp_set_num_threads(Num_Of_OMP_Threads);
    // int ID = omp_get_thread_num();
#if SHOW_OMP_NUM_OF_THREADS
    std::cout << "Using " << Num_Of_OMP_Threads << " threads for OpenMP parallelization." << std::endl;
#endif
}

void EdgeSketch_Core::Read_Sketch_Data() {
    Load_Data = std::shared_ptr<file_reader>(new file_reader(Dataset_Path, Dataset_Name, Scene_Name));
    
    //> Read absolute camera rotation matrices (all under world coordinate)
    Load_Data->Data_Reader.readRmatrix( All_R );

    //> Read absolute camera translation vectors (all under world coordinate)
    Load_Data->Data_Reader.readTmatrix( All_T );

    //> Read camera intrinsic matrix
    if (Use_Multiple_K)
        Load_Data->Data_Reader.readK( All_K );
    else
        K << fx, 0,	cx, 0, fy, cy, 0, 0, 1;

    //> Read edgels detected at a specific threshold 
    Load_Data->read_All_Edgels( All_Edgels, Edge_Detection_Init_Thresh );
}

void EdgeSketch_Core::Set_Hypothesis_Views() {
    Edges_HYPO1     = All_Edgels[hyp01_view_indx];
    Edges_HYPO2     = All_Edgels[hyp02_view_indx];
    Rot_HYPO1       = All_R[hyp01_view_indx];
    Rot_HYPO2       = All_R[hyp02_view_indx];
    Transl_HYPO1    = All_T[hyp01_view_indx];
    Transl_HYPO2    = All_T[hyp02_view_indx];

    if (Use_Multiple_K) {
        K_HYPO1 = All_K[hyp01_view_indx];
        K_HYPO2 = All_K[hyp02_view_indx];
    }
    else {
        K_HYPO1 = K;
        K_HYPO2 = K;
    }
    util->getRelativePoses( Rot_HYPO1, Transl_HYPO1, Rot_HYPO2, Transl_HYPO2, R21, T21, R12, T12 );
    F21 = util->getFundamentalMatrix(K_HYPO1.inverse(), K_HYPO2.inverse(), R21, T21); 
    F12 = util->getFundamentalMatrix(K_HYPO2.inverse(), K_HYPO1.inverse(), R12, T12);

    //> Initialize a list of paired edges between HYPO1 and HYPO2
    paired_edge         = Eigen::MatrixXd::Constant(Edges_HYPO1.rows(), Num_Of_Total_Imgs, -2);

    //> Compute epipolar wedge angles between HYPO1 and HYPO2 and valid angle range in HYPO1 for fast indexing from edges of HYPO2
    OreListdegree       = getOre.getOreList(hyp01_view_indx, hyp02_view_indx, Edges_HYPO2, All_R, All_T, K1, K2);
    OreListBardegree    = getOre.getOreListBar(Edges_HYPO1, All_R, All_T, K1, K2, hyp02_view_indx, hyp01_view_indx);
}


EdgeSketch_Core::Run_3D_Edge_Sketch() {

    #pragma omp parallel
    {
        //> Local array stacking all supported indices
        std::vector<Eigen::MatrixXd> local_thread_supported_indices;
        int edge_idx;

        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< First loop: loop over all edgels from hypothesis view 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
        //<<<<<<<<<<< Identify pairs of edge, correct the positions of the edges from Hypo2, and store the paired edges >>>>>>>>>>>>>>>>//
        #pragma omp for schedule(static, NUM_OF_OPENMP_THREADS)
        for (edge_idx = 0; edge_idx < Edges_HYPO1.rows() ; edge_idx++) {

            if ( Skip_this_Edge( edge_idx ) ) continue;
            
            //> TODO: Summarize the following piece of code into get HYPO1 edgel and the corresponding HYPO2 edgels
            
            //> Get the current edge from HYPO1
            Eigen::Vector3d pt_edgel_HYPO1;
            pt_edgel_HYPO1 << Edges_HYPO1(edge_idx,0), Edges_HYPO1(edge_idx,1), 1;

            //Get angle Thresholds from OreListBar (in Degree) 
            double thresh_ore21_1 = OreListBardegree(edge_idx, 0);
            double thresh_ore21_2 = OreListBardegree(edge_idx, 1);

            //> Find the corresponding edgel in HYPO2 based on the epipolar angle
            Eigen::MatrixXd HYPO2_idx_raw = PairHypo.getHYPO2_idx_Ore(OreListdegree, thresh_ore21_1, thresh_ore21_2);
            if (HYPO2_idx_raw.rows() == 0) continue;
            //> Retrieve Hypo2 Edgels
            Eigen::MatrixXd edgels_HYPO2 = PairHypo.getedgels_HYPO2_Ore(Edges_HYPO2, OreListdegree, thresh_ore21_1, thresh_ore21_2);
            //> Correct Edgels in Hypo2 Based on Epipolar Constraints
            Eigen::MatrixXd edgels_HYPO2_corrected = PairHypo.edgelsHYPO2correct(edgels_HYPO2, Edges_HYPO1.row(edge_idx), F21, F12, HYPO2_idx_raw);
            //> Organize the Final Edge Data
            Eigen::MatrixXd Edges_HYPO1_final(edgels_HYPO2_corrected.rows(), 4);
            Edges_HYPO1_final << edgels_HYPO2_corrected.col(0), edgels_HYPO2_corrected.col(1), edgels_HYPO2_corrected.col(2), edgels_HYPO2_corrected.col(3);
            Eigen::MatrixXd Edges_HYPO2_final(edgels_HYPO2_corrected.rows(), 4);
            Edges_HYPO2_final << edgels_HYPO2_corrected.col(4), edgels_HYPO2_corrected.col(5), edgels_HYPO2_corrected.col(6), edgels_HYPO2_corrected.col(7);

            //> Store the Hypo2 Indices
            Eigen::MatrixXd HYPO2_idx(edgels_HYPO2_corrected.rows(), 1); 
            HYPO2_idx << edgels_HYPO2_corrected.col(8);
            if (HYPO2_idx.rows() == 0) continue;

            int VALID_idx = 0;
            int stack_idx = 0;
            Eigen::MatrixXd supported_indices;
            supported_indices.conservativeResize(edgels_HYPO2.rows(), Num_Of_Total_Imgs-2);
            Eigen::MatrixXd supported_indice_current;
            supported_indice_current.conservativeResize(edgels_HYPO2.rows(),1);
            Eigen::MatrixXd supported_indices_stack;

            bool isempty_link = true;

            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Second loop:loop over all validation views >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
            for (int VALID_INDX = 0; VALID_INDX < DATASET_NUM_OF_FRAMES; VALID_INDX++) {
                //> Skip the two hypothesis views
                if(VALID_INDX == hyp01_view_indx || VALID_INDX == hyp02_view_indx) continue;

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
                if(anglediff.maxCoeff() <= Parallel_Epipolar_Line_Angle_Deg) {
                    isparallel.row(idx_pair) << 0;
                    supported_indice_current.row(idx_pair) << -2;
                    continue;
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

                if (supported_link_indx >= TO_Edges_VALID.rows()) {
                    LOG_ERROR("Something fishy here!\n");
                    std::cout << "(VALID_INDX, supported_link_indx, TO_Edges_VALID.rows())) = (" << VALID_INDX << ", " << supported_link_indx << ", " << TO_Edges_VALID.rows() << ")" << std::endl;
                }

                // Get the supporting edge idx from this validation view
                supported_indice_current.row(idx_pair) << supported_link_indx;
                // supported_indice_current.row(idx_pair) << (isparallel(idx_pair,0) != 0) ? supported_link_indx : (-2);
                // if (isparallel(idx_pair,0) != 0){
                //   supported_indice_current.row(idx_pair) << supported_link_indx;
                // }
                // else{
                //   supported_indice_current.row(idx_pair) << -2;
                // }
                if (supported_link_indx != -2) {
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

            //> Check for Empty Supported Indices
            if (isempty_link) continue;

            //Create a Stack of Supported Indices
            std::vector<double> indices_stack(supported_indices_stack.data(), supported_indices_stack.data() + supported_indices_stack.rows());
            std::vector<double> indices_stack_unique = indices_stack;
            std::sort(indices_stack_unique.begin(), indices_stack_unique.end());
            std::vector<double>::iterator it1;
            it1 = std::unique(indices_stack_unique.begin(), indices_stack_unique.end());
            indices_stack_unique.resize( std::distance(indices_stack_unique.begin(),it1) );

            //> Count the Occurrence of Each Unique Index
            Eigen::VectorXd rep_count;
            rep_count.conservativeResize(indices_stack_unique.size(),1);
            for (int unique_idx = 0; unique_idx<indices_stack_unique.size(); unique_idx++) {
                rep_count.row(unique_idx) << double(count(indices_stack.begin(), indices_stack.end(), indices_stack_unique[unique_idx]));
            }

            //> Find the maximal number of supports from validation views and check if this number is over the threshold
            Eigen::VectorXd::Index maxIndex;
            double max_support = rep_count.maxCoeff(&maxIndex);
            if( max_support < Max_Num_Of_Support_Views )
                continue;
            
            int finalpair = -2;
            int numofmax = std::count(rep_count.data(), rep_count.data()+rep_count.size(), max_support);
            if (numofmax > 1) {
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
                    Edge_Pts.row(maxidx) << Edges_HYPO2_final(indices_stack_unique[max_index[maxidx]], 0), \
                                            Edges_HYPO2_final(indices_stack_unique[max_index[maxidx]], 1);
                }
                Eigen::VectorXd Ap = coeffs(0)*Edge_Pts.col(0);
                Eigen::VectorXd Bp = coeffs(1)*Edge_Pts.col(1);
                Eigen::VectorXd numDist = Ap + Bp + Eigen::VectorXd::Ones(Ap.rows())*coeffs(2);
                double denomDist = coeffs(0)*coeffs(0) + coeffs(1)*coeffs(1);
                denomDist = sqrt(denomDist);
                Eigen::VectorXd dist = numDist.cwiseAbs()/denomDist;
                Eigen::VectorXd::Index   minIndex;
                double min_dist = dist.minCoeff(&minIndex);
                //> ignore if the distance from the reprojected edge to the epipolar line is greater than some threshold
                if (min_dist > Reproj_Dist_Thresh) continue;

                finalpair = int(indices_stack_unique[max_index[minIndex]]);
            }
            else {
                finalpair = int(indices_stack_unique[int(maxIndex)]);
            }

            //> paired_edge is a row vector continaing (hypo1 edge index), (paired hypo2 edge index), (paired validation edge indices)
            paired_edge.row(edge_idx) << edge_idx, HYPO2_idx(finalpair), supported_indices.row(finalpair);
        }
        //> A critical session to stack all local supported indices
        #pragma omp critical
        all_supported_indices.insert(all_supported_indices.end(), local_thread_supported_indices.begin(), local_thread_supported_indices.end());
    }
}

void EdgeSketch_Core::Finalize_Edge_Pairs() {

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
        Eigen::MatrixXd edgel_HYPO1   = Edges_HYPO1.row(int(paired_edge_final(pair_idx,0)));  //> edge index in hypo 1
        Eigen::MatrixXd edgel_HYPO2   = Edges_HYPO2.row(int(paired_edge_final(pair_idx,1)));  //> edge index in hypo 2
        // Eigen::MatrixXd HYPO2_idx_raw = Edges_HYPO1.row(int(paired_edge_final(pair_idx,0)));      //> WTF is this? Edges_HYPO1 = Edges_HYPO1
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
        Eigen::Vector3d edge_pt_3D = util.linearTriangulation(2, pts, Rs, Ts, K1_v);

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

}

EdgeSketch_Core::~EdgeSketch_Core() { }

#endif
