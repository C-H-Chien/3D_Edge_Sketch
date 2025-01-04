#ifndef EDGE_MAPPING_H
#define EDGE_MAPPING_H

#include <unordered_map>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include "util.hpp"

// Hash function for Eigen::Vector3d to use in unordered_map
struct HashEigenVector3d {
    std::size_t operator()(const Eigen::Vector3d& vec) const {
        std::size_t seed = 0;
        for (int i = 0; i < vec.size(); ++i) {
            seed ^= std::hash<double>()(vec[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

// Hash function for Eigen::Vector2d to use in unordered_map
struct HashEigenVector2d {
    std::size_t operator()(const Eigen::Vector2d& vec) const {
        std::size_t h1 = std::hash<double>()(vec(0));
        std::size_t h2 = std::hash<double>()(vec(1));
        return h1 ^ (h2 << 1); // Combine hashes
    }
};

// EdgeMapping class for storing 3D-to-2D edge mappings
class EdgeMapping {

public:
    //> Constructor
    EdgeMapping() {
        util = std::shared_ptr<MultiviewGeometryUtil::multiview_geometry_util>(new MultiviewGeometryUtil::multiview_geometry_util());
    }

    std::unordered_map<Eigen::Vector3d, std::vector<std::pair<Eigen::Vector3d, int>>, HashEigenVector3d> edge_3D_to_supporting_edges;
    std::map<int, std::unordered_map<Eigen::Vector2d, std::vector<Eigen::Vector3d>, HashEigenVector2d>> frame_to_edge_to_3D_map;

    void add3DToSupportingEdgesMapping(const Eigen::Vector3d &edge_3D, const Eigen::Vector3d &supporting_edge, int image_number) 
    {
        edge_3D_to_supporting_edges[edge_3D].emplace_back(supporting_edge, image_number);
    }

    void add3DToFrameMapping(const Eigen::Vector3d& edge_3D, const Eigen::Vector2d& supporting_edge, int frame) 
    {
        frame_to_edge_to_3D_map[frame][supporting_edge].push_back(edge_3D);
    }

    void clear_edge_linking_data() {
        edge_3D_to_supporting_edges.clear();
        frame_to_edge_to_3D_map.clear();
    }

    void write_3D_edge_and_2D_edge_linkings( Eigen::Matrix3d R_ref, Eigen::Vector3d T_ref, int edge_sketch_pass_count )
    {
        std::ofstream file_edge_linking;
        std::string Output_File_Edge_Linking_Path = "../../" + OUTPUT_FOLDER_NAME + "/3D_2D_edge_linkings_pass" + std::to_string(edge_sketch_pass_count) + ".txt";
        file_edge_linking.open(Output_File_Edge_Linking_Path);
        if (!file_edge_linking) {
            LOG_FILE_ERROR(Output_File_Edge_Linking_Path); 
        }
        else {
            //> Loop over all 3D edge points
            for (const auto& [edge_3D, supporting_edges] : edge_3D_to_supporting_edges) {
                //> Transform the 3D edges from the first hypothesis view coordinate (Gamma1s) to the world coordinate (Gamma1s_world)
                Eigen::Vector3d edge_3D_world = util->transformToWorldCoordinates(edge_3D, R_ref, T_ref);
                file_edge_linking << "-1\t-1\t-1\t-1\n";
                for (int ew = 0; ew < 3; ew++) {
                    file_edge_linking << edge_3D_world(ew);
                    file_edge_linking << "\t";
                }
                file_edge_linking << "1.0\n";
                
                //> Print the 3D edge
                // std::cout << "3D Edge: [" << edge_3D_world(0) << ", " << edge_3D_world(1) << ", " << edge_3D_world(2) << "]\n";

                // Print all supporting 2D edges and their image numbers
                for (const auto& [supporting_edge, image_number] : supporting_edges) {
                    for (int ei = 0; ei < 3; ei++) {
                        file_edge_linking << supporting_edge(ei);
                        file_edge_linking << "\t";
                    }
                    file_edge_linking << image_number;
                    file_edge_linking << "\n";
                    // std::cout << "    Supporting Edge: [" << supporting_edge(0) << ", " << supporting_edge(1) << "] from Image " << image_number << "\n";
                }
            }
            file_edge_linking.close();
        }
    }

private:
    std::shared_ptr<MultiviewGeometryUtil::multiview_geometry_util> util = nullptr;
};

#endif  // EDGE_MAPPING_H
