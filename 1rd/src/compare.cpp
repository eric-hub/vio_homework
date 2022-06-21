#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>
#include "sophus/se3.hpp"

int main(int argc, char **argv) {
    // 沿Z轴转90度的旋转矩阵
    Eigen::Matrix3d R = Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d(0, 0, 1)).toRotationMatrix();
    Eigen::Vector3d w(0.01, 0.02, 0.03);

    Eigen::Matrix3d R_exp = R * Sophus::SO3d::exp(w).matrix();

    Eigen::Quaterniond q(R);

    Eigen::Quaterniond w_q(1, 0.5 * w(0), 0.5 * w(1), 0.5 * w(2));
    w_q.normalize(); //归一化

    Eigen::Quaterniond q_q = q * w_q;

    Eigen::Matrix3d r_q_q = q_q.toRotationMatrix();

    std::cout.precision(3);
    std::cout << R_exp << std::endl;
    std::cout << r_q_q << std::endl;

    return 0;
}
