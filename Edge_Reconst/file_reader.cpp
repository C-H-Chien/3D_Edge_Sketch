#include "file_reader.hpp"

// Function to read edgel files
void readEdgelFiles(std::vector<Eigen::MatrixXd> &All_Edgels, std::fstream &Edge_File, double &rd_data, Eigen::Vector4d &row_edge, int &file_idx, int &d, int &q, int thresh_EDG) {
    std::cout << "Reading edgel files ...\n";
    
    //> Looping over all frames and read corresponding third-order edgels
    while(file_idx < DATASET_NUM_OF_FRAMES) {
        std::string Edge_File_Path = "../../datasets/" + DATASET_NAME + "/" + SCENE_NAME + "/Edges/Edge_" + std::to_string(file_idx) + "_t" + std::to_string(thresh_EDG) + ".txt";
        #if DEBUG_READ_FILES
            std::cout << Edge_File_Path << std::endl;
        #endif
        file_idx++;
        Eigen::MatrixXd Edgels; //> Declare locally, ensuring the memory addresses are different for different frames
        Edge_File.open(Edge_File_Path, std::ios_base::in);
        if (!Edge_File) {
            LOG_FILE_ERROR("Edge file not existed!"); exit(1);
        }else {
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
            d = 0;
            q = 0;
        }
    }
    file_idx = 1;
    std::cout << "Edge file loading finished" << std::endl;
}




void readHypothesisEdgelFiles(int hyp01_view_indx, int hyp02_view_indx, std::vector<Eigen::MatrixXd> &All_Edgels_H12, std::fstream &Edge_File, double &rd_data, Eigen::Vector4d &row_edge, int &H_idx, int &file_idx, int &d, int &q, int thresh_EDG){
    while(file_idx < 3) {
        if(file_idx == 1){
            H_idx = hyp01_view_indx;
        }else{
            H_idx = hyp02_view_indx;
        }
        std::string Edge_File_PathH12 = "../../datasets/" + DATASET_NAME + "/" + SCENE_NAME + "/Edges/Edge_"+std::to_string(H_idx)+"_t" + std::to_string(thresh_EDG) + ".txt"; 
        #if DEBUG_READ_FILES
            std::cout << Edge_File_PathH12 << std::endl;
        #endif
        file_idx ++;
        Eigen::MatrixXd Edgels; //> Declare locally, ensuring the memory addresses are different for different frames
        Edge_File.open(Edge_File_PathH12, std::ios_base::in);
        if (!Edge_File) { 
            LOG_FILE_ERROR("Edge file not existed!"); exit(1); 
        }else{
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
            All_Edgels_H12.push_back(Edgels);
            d = 0;
            q = 0;
        }
    }
    file_idx = 1;
    std::cout<< "HYPO1 and HYPO2 Edge file loading finished" <<std::endl;
}


void readRmatrix(std::vector<Eigen::Matrix3d> &All_R, Eigen::Matrix3d &R_matrix, std::fstream &Rmatrix_File, double &rd_data, Eigen::Vector3d &row_R, int &d, int &q){
  
  std::string Rmatrix_File_Path = "../../datasets/" + DATASET_NAME + "/" + SCENE_NAME + "/RnT/R_matrix.txt";
  Rmatrix_File.open(Rmatrix_File_Path, std::ios_base::in);
  if (!Rmatrix_File) { 
    LOG_FILE_ERROR("R_matrix file not existed!"); exit(1); 
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
  std::cout<< "R matrix loading finished" <<std::endl;
}

void readTmatrix(std::vector<Eigen::Vector3d> &All_T, Eigen::Vector3d &T_matrix, std::fstream &Tmatrix_File, double &rd_data, int &d, int &q){

    std::string Tmatrix_File_Path = "../../datasets/" + DATASET_NAME + "/" + SCENE_NAME + "/RnT/T_matrix.txt";
    Tmatrix_File.open(Tmatrix_File_Path, std::ios_base::in);
    if (!Tmatrix_File) { 
        LOG_FILE_ERROR("T_matrix file not existed!"); exit(1); 
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
  
  std::cout<< "T matrix loading finished" <<std::endl;
}


void readK(std::fstream &Kmatrix_File, std::vector<Eigen::Matrix3d> &All_K, Eigen::Matrix3d &K, Eigen::Matrix3d &K_matrix, Eigen::Vector3d &row_K, double &rd_data, int &d, int &q){
    if(IF_MULTIPLE_K == 1){
        std::string Kmatrix_File_Path = "../../datasets/" + DATASET_NAME + "/" + SCENE_NAME + "/RnT/K_matrix.txt";
        Kmatrix_File.open(Kmatrix_File_Path, std::ios_base::in);
    if (!Kmatrix_File) { 
      LOG_FILE_ERROR("K_matrix file not existed!"); exit(1);
    }else {
      while (Kmatrix_File >> rd_data) {
        row_K(q) = rd_data;
        q++;
        if(q>2){
          K_matrix.row(d) = row_K;
          row_K = {};
          q = 0;
          d++;
        }
        if(d>2){
          All_K.push_back(K_matrix);
          K_matrix = {};
          d = 0;
        }
      }
      Kmatrix_File.close();
    }
  }else {
    if (DATASET_NAME == "ABC-NEF") 
      K << 1111.11136542426,	0,	399.500000000000, 0,	1111.11136542426,	399.500000000000, 0,	0,	1;
    else if (DATASET_NAME == "Replica")
      K << 600,	0,	599.500000000000, 0,	600,	339.500000000000, 0,	0,	1;
  }
  std::cout<< "K matrix loading finished" <<std::endl;
}
