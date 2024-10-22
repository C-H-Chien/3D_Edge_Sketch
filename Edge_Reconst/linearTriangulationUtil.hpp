// linearTriangulationUtil.hpp
#ifndef LINEAR_TRIANGULATION_UTIL_HPP
#define LINEAR_TRIANGULATION_UTIL_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>

// Function declaration for linear triangulation
Eigen::Vector3d linearTriangulation(int N,
                                    const std::vector<Eigen::Vector2d> pts, 
                                    const std::vector<Eigen::Matrix3d> & Rs,
                                    const std::vector<Eigen::Vector3d> & Ts,
                                    const std::vector<double> & K);

#endif // LINEAR_TRIANGULATION_UTIL_HPP
