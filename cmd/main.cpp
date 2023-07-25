#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <assert.h>
#include <string>

//> Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>
using namespace std;

// ========================================================================================================================
// main function
//
// Modifications
//    Chiang-Heng Chien  23-07-18    Initially create a blank repository with minor multiview geometry utility functions.
//
//> (c) LEMS, Brown University
//> Yilin Zheng (yilin_zheng@brown.edu)
//> Chiang-Heng Chien (chiang-heng_chien@brown.edu)
// =========================================================================================================================


int main(int argc, char **argv) {
  std::string File_Path = "/users/yzhen105/Edge_Based_Reconstruction/Edge_Reconst/testfile.txt";
  std::fstream data_File;
  int d = 0;
  int q = 0;
  double rd_data;
  //double testEdges[];
  //std::vector< Eigen::Vector4d > Edgels;
  //Eigen::Vector4d test;
  std::vector<vector<double>> Edgels;
  std::vector<double> test;
  cout << "test size: "   << test.size()   << '\n';
  cout << "Edgels size: " << Edgels.size() << '\n';
  cout << "read file now\n";
  data_File.open(File_Path, std::ios_base::in);
  if (!data_File) { 
      std::cerr << "file not existed!\n"; exit(1); 
      }
  else {
      while (data_File >> rd_data) {
          //test(q) = rd_data;
          test.push_back(rd_data);
          q++;
          if(q>3){
            //cout << "test size: " << test.size() << '\n';
            Edgels.push_back(test);
            //Edgels.insert (Edgels.end(),test);
            //cout <<  Edgels[d] << '\n';
            //cout << "Edgels size: " << Edgels.size() << '\n';
            test = {};
            q = 0;
            d++;
          }


      }
      
      
  }
  cout << "print Edgels now\n";
  
  cout << Edgels[0][0] <<endl;
  cout << Edgels[0][1] <<endl;

  cout << Edgels[101][0] <<endl;
  cout << Edgels[101][1] <<endl;
}
