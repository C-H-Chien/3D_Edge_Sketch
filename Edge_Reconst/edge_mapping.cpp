#include "edge_mapping.hpp"

void EdgeMapping::add3DToSupportingEdgesMapping(const Eigen::Vector3d &edge_3D, const Eigen::Vector2d &supporting_edge, int image_number) {
    edge_3D_to_supporting_edges[edge_3D].emplace_back(supporting_edge, image_number);
}

void EdgeMapping::add3DToFrameMapping(const Eigen::Vector3d& edge_3D, const Eigen::Vector2d& supporting_edge, int frame) {
    frame_to_edge_to_3D_map[frame][supporting_edge].push_back(edge_3D);
}
