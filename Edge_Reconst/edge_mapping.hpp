#ifndef EDGE_MAPPING_H
#define EDGE_MAPPING_H

#include <unordered_map>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>


class EdgeMapping {
public:
    // Map from hypothesis edge index to corresponding 3D edge
    std::unordered_map<int, Eigen::Vector3d> hypothesis_to_3D_edge;

    // Member functions
    void add3DToSupportingEdgesMapping(const Eigen::Vector3d &edge_3D, const Eigen::Vector2d &supporting_edge, int image_number);
};

#endif  // EDGE_MAPPING_H
