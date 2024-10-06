#include "edge_mapping.hpp"

// Adds a 3D edge to the map along with its supporting 2D edge and image number
void EdgeMapping::add3DToSupportingEdgesMapping(const Eigen::Vector3d &edge_3D, const Eigen::Vector2d &supporting_edge, int image_number) {
    edge_3D_to_supporting_edges[edge_3D].emplace_back(supporting_edge, image_number);
}

// Maps a hypothesis edge index to a corresponding 3D edge
void EdgeMapping::addHypothesisTo3DMapping(int hypothesis_edge_idx, const Eigen::Vector3d& edge_3D) {
    hypothesis_to_3D_edge[hypothesis_edge_idx] = edge_3D;
}
