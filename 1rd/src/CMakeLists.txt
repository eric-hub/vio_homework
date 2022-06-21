cmake_minimum_required(VERSION 3.0)
project(compare)

# 为使用 sophus，需要使用find_package命令找到它
find_package(Sophus REQUIRED)

# Eigen
include_directories("/usr/include/eigen3")
add_executable(compare compare.cpp)
target_link_libraries(compare Sophus::Sophus)
