#include "sophus/se3.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <iostream>

int main(int argc, char **argv) {
  Eigen::Vector3d w(0.01, 0.02, 0.03);

  Eigen::Quaternion<double> qq(1, 0.5 * w(0), 0.5 * w(1), 0.5 * w(2));
  std::cout << "delta_q:" << std::endl << qq.coeffs() << std::endl;

  Eigen::Quaterniond delta_qq(qq.normalized());

  std::cout << "norm_detla_q:" << std::endl << delta_qq.coeffs() << std::endl;

  Sophus::SO3d delta_SO3(delta_qq);
  Eigen::Matrix3d R =
      Eigen::AngleAxisd(M_PI_4, Eigen::Vector3d(0, 0, 1)).toRotationMatrix();

  std::cout << "R:" << std::endl << R << std::endl;
  Sophus::SO3d SO3_R(R);

  Eigen::Quaterniond q(R);
  std::cout << "q:" << std::endl << q.coeffs() << std::endl;

  Sophus::SO3d update_SO3 = SO3_R * delta_SO3;
  std::cout << "qupdate_R:" << std::endl << update_SO3.matrix() << std::endl;

  Eigen::Quaterniond update_q = q * delta_qq;
  Eigen::Matrix3d update_qR(update_q);
  std::cout << "qupdate_qR:" << std::endl << update_qR << std::endl;

  return 0;
}