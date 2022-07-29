//
// Created by hyj on 18-11-11.
//
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

struct Pose {
    Pose(Eigen::Matrix3d R, Eigen::Vector3d t)
        : Rwc(R), qwc(R), twc(t){};
    Eigen::Matrix3d Rwc;
    Eigen::Quaterniond qwc;
    Eigen::Vector3d twc;

    Eigen::Vector2d uv; // 这帧图像观测到的特征坐标
};
int main() {
    int poseNums = 10;
    double radius = 8;
    double fx = 1.;
    double fy = 1.;
    std::vector<Pose> camera_pose;
    for (int n = 0; n < poseNums; ++n) {
        double theta = n * 2 * M_PI / (poseNums * 4); // 1/4 圆弧
        // 绕 z轴 旋转
        Eigen::Matrix3d R;
        R = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
        Eigen::Vector3d t = Eigen::Vector3d(radius * cos(theta) - radius, radius * sin(theta), 1 * sin(2 * theta));
        camera_pose.push_back(Pose(R, t));
    }

    // 随机数生成 1 个 三维特征点
    std::default_random_engine generator;
    std::uniform_real_distribution<double> xy_rand(-4, 4.0);
    std::uniform_real_distribution<double> z_rand(8., 10.);
    double tx = xy_rand(generator);
    double ty = xy_rand(generator);
    double tz = z_rand(generator);

    Eigen::Vector3d Pw(tx, ty, tz);
    // 这个特征从第三帧相机开始被观测，i=3
    int start_frame_id = 0;
    int end_frame_id = poseNums;
    std::normal_distribution<double> noise(0., 4. / 1000.);
    for (int i = start_frame_id; i < end_frame_id; ++i) {
        Eigen::Matrix3d Rcw = camera_pose[i].Rwc.transpose();
        Eigen::Vector3d Pc = Rcw * (Pw - camera_pose[i].twc);

        double x = Pc.x();
        double y = Pc.y();
        double z = Pc.z();

        camera_pose[i].uv = Eigen::Vector2d(x / z + noise(generator), y / z + noise(generator));
    }

    /// TODO::homework; 请完成三角化估计深度的代码
    // 遍历所有的观测数据，并三角化
    Eigen::Vector3d P_est; // 结果保存到这个变量
    P_est.setZero();
    /* your code begin */
    Eigen::MatrixXd D(2 * (end_frame_id - start_frame_id), 4);
    D.setZero();
    Eigen::Matrix<double, 1, 4> pk1;
    pk1.setZero();
    Eigen::Matrix<double, 1, 4> pk2;
    pk2.setZero();
    Eigen::Matrix<double, 1, 4> pk3;
    pk3.setZero();
    for (int i = start_frame_id; i < end_frame_id; i++) {
        //R^T
        Eigen::Matrix3d Rcw = camera_pose[i].Rwc.transpose();
        //-R^T*t
        Eigen::Vector3d tcw = -Rcw * camera_pose[i].twc;
        pk1 << Rcw.block<1, 3>(0, 0), tcw(0);
        pk2 << Rcw.block<1, 3>(1, 0), tcw(1);
        pk3 << Rcw.block<1, 3>(2, 0), tcw(2);

        D.block<1, 4>(2 * (i - start_frame_id), 0) = camera_pose[i].uv(0) * pk3 - pk1;
        D.block<1, 4>(2 * (i - start_frame_id) + 1, 0) = camera_pose[i].uv(1) * pk3 - pk2;
    }
    Eigen::Matrix<double, 4, 4> DTD = D.transpose() * D;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(DTD, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixU();
    Eigen::Vector4d A = svd.singularValues();

    // std::cout << U << std::endl;
    std::cout << "U:\n"
              << U << std::endl;
    std::cout << "V:\n"
              << V << std::endl;
    std::cout << "A:\n"
              << A << std::endl;

    std::cout << "sigema3/sigema4:\n"
              << A(2) / A(3) << std::endl;

    P_est = V.block<3, 1>(0, 3) / U(3, 3); //齐次坐标归一化
    /* your code end */

    std::cout << "ground truth: \n"
              << Pw.transpose() << std::endl;
    std::cout << "your result: \n"
              << P_est.transpose() << std::endl;
    // TODO:: 请如课程讲解中提到的判断三角化结果好坏的方式，绘制奇异值比值变化曲线
    return 0;
}
