#ifndef EDGE_MAPPING_H
#define EDGE_MAPPING_H

#include <unordered_map>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>

struct HashEigenVector3d {
    std::size_t operator()(const Eigen::Vector3d& vec) const {
        std::size_t seed = 0;
        for (int i = 0; i < vec.size(); ++i) {
            seed ^= std::hash<double>()(vec[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};


class EdgeMapping {
public:

    std::unordered_map<Eigen::Vector3d, std::vector<std::pair<Eigen::Vector2d, int>>, HashEigenVector3d> edge_3D_to_supporting_edges;

    void add3DToSupportingEdgesMapping(const Eigen::Vector3d &edge_3D, const Eigen::Vector2d &supporting_edge, int image_number);
};

#endif  // EDGE_MAPPING_H
