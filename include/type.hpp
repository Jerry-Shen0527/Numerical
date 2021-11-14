#pragma once

#include <iostream>
#include <Eigen/Eigen>
using Float = double;

using Point2f = Eigen::Matrix<Float, 2, 1, 0>;
using Vector = Eigen::Matrix<Float, -1, 1, 0>;
using Matrix = Eigen::Matrix<Float, -1, -1, 0>;

const Float Pi = 3.14159265358979323846264;