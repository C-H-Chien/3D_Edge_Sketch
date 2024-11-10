#ifndef EDGE_MAPPING_H
#define EDGE_MAPPING_H

#include <unordered_map>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>

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
    std::unordered_map<Eigen::Vector3d, std::vector<std::pair<Eigen::Vector2d, int>>, HashEigenVector3d> edge_3D_to_supporting_edges;
    std::map<int, std::unordered_map<Eigen::Vector2d, std::vector<Eigen::Vector3d>, HashEigenVector2d>> frame_to_edge_to_3D_map;

    void add3DToSupportingEdgesMapping(const Eigen::Vector3d &edge_3D, const Eigen::Vector2d &supporting_edge, int image_number);
    void add3DToFrameMapping(const Eigen::Vector3d& edge_3D, const Eigen::Vector2d& supporting_edge, int frame);


};

#endif  // EDGE_MAPPING_H
