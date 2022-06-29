### 第二章

[TOC]



### 1 对非ros:生成运动数据

进入vio_data_simulation项目,执行build.sh,生成数据,画图

```bash
./build.sh
cd..
./bin/data_gen
cd python_tool
python draw_trajcory.py
```

中值积分修改如下imu.testImu如下:

```c++
std::vector<MotionData> imudata;
    LoadPose(src, imudata);

    std::ofstream save_points;
    save_points.open(dist);

    double dt = param_.imu_timestep;
    Eigen::Vector3d Pwb = init_twb_;   // position :    from  imu measurements
    Eigen::Quaterniond Qwb(init_Rwb_); // quaterniond:  from imu measurements
    Eigen::Quaterniond Qwb_0(init_Rwb_);
    Eigen::Vector3d Vw = init_velocity_; // velocity  :   from imu measurements
    Eigen::Vector3d gw(0, 0, -9.81);     // ENU frame
    Eigen::Vector3d temp_a;
    Eigen::Vector3d theta;
    for (int i = 1; i < imudata.size(); ++i) {
        MotionData imupose = imudata[i];

        //delta_q = [1 , 1/2 * thetax , 1/2 * theta_y, 1/2 * theta_z]
        // Eigen::Quaterniond dq;
        // Eigen::Vector3d dtheta_half = imupose.imu_gyro * dt / 2.0;
        // dq.w() = 1;
        // dq.x() = dtheta_half.x();
        // dq.y() = dtheta_half.y();
        // dq.z() = dtheta_half.z();
        // dq.normalize();

        /// imu 动力学模型 欧拉积分
        // Eigen::Vector3d acc_w = Qwb * (imupose.imu_acc) + gw;  // aw = Rwb * ( acc_body - acc_bias ) + gw
        // Qwb = Qwb * dq;
        // Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
        // Vw = Vw + acc_w * dt;

        // 中值积分
        MotionData imupose_0 = imudata[i - 1];
        Eigen::Quaterniond dq;
        Eigen::Vector3d dtheta_half = ((imupose_0.imu_gyro + imupose.imu_gyro) / 2.0) * dt / 2.0;
        dq.w() = 1;
        dq.x() = dtheta_half.x();
        dq.y() = dtheta_half.y();
        dq.z() = dtheta_half.z();
        dq.normalize();

        Qwb = Qwb_0 * dq;

        Eigen::Vector3d acc_w = (Qwb * (imupose.imu_acc) + gw + Qwb_0 * (imupose_0.imu_acc) + gw) / 2.0;
        Pwb = Pwb + Vw * dt + 0.5 * dt * dt * acc_w;
        Vw = Vw + acc_w * dt;

        Qwb_0 = Qwb;
```

编译与以上步骤一样。以下是结果对比

|                        欧拉法                        |                     中值法                     |
| :--------------------------------------------------: | :--------------------------------------------: |
| ![ol.png](/home/eric/vio_homework/2rd/images/ol.png) | ![](/home/eric/vio_homework/2rd/images/zz.png) |

明显看出中值法优于欧拉法

### 2.对ROS: 生成静止imu数据，用于allan标定

#### 2.1 ros数据生成

- ros编译 vio_data_simulation-ros_version

- 将生成目录修改到本项目中

  ```c++
  const std::string bag_path = home_path + "/vio_work/src/vio_data_simulation-ros_version/bag/imu.bag";
  ```

  

- 创建launch文件

  ```xml
  <launch>
      <node pkg="vio_data_simulation" type="vio_data_simulation_node" name ="vio_data_simulation_node" output="screen">
      </node>
  </launch>
  ```

- 启动项目生成bag文件

  ```bash
  roslaunch vio_data_simulation gener_alldata.launch
  ```

- 在vio_data_simulation-ros_version/bag下生成了

#### 2.2 imu_utils测试

##### 	imu_utils安装

- 下载  [imu_utils](https://github.com/gaowenliang/imu_utils) ,还需要下载 [code_utils](https://github.com/gaowenliang/code_utils) 
-  在code_utils下面找到sumpixel_test.cpp，修改#include "backward.hpp"为 #include “code_utils/backward.hpp”
- 先将imu_utils移出src,code_utils编译成功后，再将imu_utils移入src，再次进行编译

   ##### 	imu_utils使用

- 编写launch文件

  ```xml
  <launch>
      <node pkg="imu_utils" type="imu_an" name="imu_an" output="screen">
          <param name="imu_topic" type="string" value= "/imu"/>
          <param name="imu_name" type="string" value= "imu"/>
          <param name="data_save_path" type="string" value= "$(find imu_utils)/data"/>
          <param name="max_time_min" type="int" value= "120"/>
          <param name="max_cluster" type="int" value= "100"/>
      </node>
  </launch>
  ```

- 启动launch

  ```bash
  roslaunch imu_utils imu.launch
  ```

  

- 启动bag数据

  ```bash
  rosbag play -r 500 imu.bag
  ```

  

- 结果如下

  ```yaml
  type: IMU
  name: imu
  Gyr:
     unit: " rad/s"
     avg-axis:
        gyr_n: 2.1246093542450672e-01
        gyr_w: 9.8254593710355950e-04
     x-axis:
        gyr_n: 2.1186040177448104e-01
        gyr_w: 9.7205107115200121e-04
     y-axis:
        gyr_n: 2.1345478592049774e-01
        gyr_w: 9.5178211626609536e-04
     z-axis:
        gyr_n: 2.1206761857854142e-01
        gyr_w: 1.0238046238925822e-03
  Acc:
     unit: " m/s^2"
     avg-axis:
        acc_n: 2.6954952191200460e-01
        acc_w: 3.4637840369062469e-03
     x-axis:
        acc_n: 2.7201494143970562e-01
        acc_w: 3.9177076874172285e-03
     y-axis:
        acc_n: 2.6597422538148940e-01
        acc_w: 3.4939714929174324e-03
     z-axis:
        acc_n: 2.7065939891481872e-01
        acc_w: 2.9796729303840793e-03
  ```

  同一单位后，与指定的bias进行比对。转换方式：未转换单位/sqrt(200)

  

  |           |       未转换单位       |       转换后单位       | 指定误差 |
  | --------- | :--------------------: | :--------------------: | :------: |
  | gro_bias  | 9.8254593710355950e-04 |  6.94764894953218e-05  | 0.00005  |
  | gro_white | 2.1246093542450672e-01 |  0.015023256817590588  |  0.015   |
  | acc_bias  | 3.4637840369062469e-03 | 0.00024492651810621214 |  0.0005  |
  | acc_white | 2.6954952191200460e-01 |  0.019060029480957034  |  0.019   |

  从结果看，标定结果还是很准的

#### 2.3 Kalibr_allan标定

##### 	matlab安装

##### 	kalibr_all安装

- 下载地址 [kalibr_allan](https://github.com/rpng/kalibr_allan)

- 修改matlab地址

  ```tex
  修改~/vio_work/src/kalibr_allan/bagconvert/cmake目录下的FindMatlab.cmake，
  找到  find_program(MATLAB_EXE_PATH matlab　　　这一行，将他修改成
   find_program(MATLAB_EXE_PATH matlab
              PATHS /usr/local/MATLAB/R2018a/bin)
  这样可以找到matlab
  ```

  这里的“/usr/local/MATLAB/R2018a/bin”需要给成你电脑中matlab 的相应位置。之后重新编译，如果之前已经失败过，则删除build文件之后重新编辑即可

- 下来编译即可

##### 	标定

- bag数据转换为mat

  将之前生成的bag数据放到kalibr_allan的dat中，执行命令转换

  ```bash
  roscore
  另起终端，进入到data目录
  rosrun bagconvert bagconvert imu.bag imu
  ```

  **注意**:上面的topic没有/

- 执行SCRIPT_allan_matparallel，这里我复制了个SCRIPT_allan_matparallel_my进行修改，主要修改了文件名称

  ```matlab
  %% Initalization
  close all
  clear all
  
  % Read in our toolboxes
  % addpath('functions/allan_v3')
  
  % Our bag information
  %mat_path = '../data/imu_mtig700.mat';
  %mat_path = '../data/imu_tango.mat';
  mat_path = '../data/imu.mat';
  
  % IMU information (todo: move this to the yaml file)
  %update_rate = 400;
  %update_rate = 100;
  update_rate = 200;
  ```

  然后在matlab中执行

  ```matlab
  >> SCRIPT_allan_matparallel_my
  opening the mat file.
  loading timeseries.
  imu frequency of 200.00.
  sample period of 0.00500.
  calculating allan deviation.
  Elapsed time is 1611.157493 seconds.
  saving to: results_20220629T092430.mat
  ```

  再修改SCRIPT_process_results，主要也是修改了名称

  ```matlab
  %% Initalization
  close all
  clear all
  
  % Read in our toolboxes
  addpath('functions')
  addpath('functions/allan_v3')
  
  % Our bag information
  %titlestr = 'XSENS MTi-G-710';
  %mat_path = '../data/bags/results_20170908T182715.mat';
  
  %titlestr = 'Tango Yellowstone #1';
  %mat_path = '../data/bags/results_20171031T115123.mat';
  
  titlestr = 'IMU-Sensor';
  mat_path = '../data/results_20220629T092430.mat';
  ```

  然后执行

  ```matlab
  >> SCRIPT_process_results
  => opening the mat file.
  => plotting accelerometer.
  tau = 1.00 | tauid1 = 1089
  h_fit1 slope = -0.5000 | y-intercept = -3.9534
  h_fit2 slope = 0.5000 | y-intercept = -8.1766
  tau = 2.99 | tauid2 = 1201
  => plotting gyroscope.
  tau = 1.00 | tauid1 = 1089
  h_fit1 slope = -0.5000 | y-intercept = -4.1862
  h_fit2 slope = 0.5000 | y-intercept = -10.2046
  tau = 2.99 | tauid2 = 1201
  => final results
  accelerometer_noise_density = 0.01918878
  accelerometer_random_walk   = 0.00048618
  gyroscope_noise_density     = 0.01520415
  gyroscope_random_walk       = 0.00006398
  ```

- 结果如下		

|                             ACC                              |                             GRO                              |
| :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![acc](/home/eric/vio_homework/2rd/images/results_20220629T092430_accel.png) | ![gro](/home/eric/vio_homework/2rd/images/results_20220629T092430_gyro.png) |
|                   bias:0.0005 noise:0.019                    |                   bias:0.00005 noise:0.015                   |

kalibr_allan的标定结果比imu_utils更准确,推荐kalibr_allan(就是装matlab有些麻烦)

### 3 论文阅读

​	处理连续时间视觉-imu定位、建图、标定

#### 	3.1连续时间表示

​		通过李代数SE3和se3表示

##### 	3.1.1 摄像头位置转换



![20](/home/eric/vio_homework/2rd/images/20.png)

​		这里 $ b,a $ 代表从a坐标系到b坐标系

​		![21](/home/eric/vio_homework/2rd/images/21.png)	

​		代表线速度

##### 	3.1.2 C 2连续曲线SE3

​		C2不清楚什么意思，可能是平方

##### 	3.1.3 B样条累计误差函数表示

![22](/home/eric/vio_homework/2rd/images/22.png)

![23](/home/eric/vio_homework/2rd/images/23.png)

![24](/home/eric/vio_homework/2rd/images/24.png)

​		$T_{w,s}(t)$是在时间t内的位姿

##### 3.1.4 立体累计B样条

​		![26](/home/eric/vio_homework/2rd/images/26.png)

​	![27](/home/eric/vio_homework/2rd/images/27.png)

​		这里没看很懂

#### 	3.2 生成视觉惯性模型

##### 	3.2.1 参数化

​	![28](/home/eric/vio_homework/2rd/images/28.png)

​		定义像素表示模型

​	![29](/home/eric/vio_homework/2rd/images/29.png)

​		定义磁力计和加速度

##### 	3.2.1 最小化公式

​	![30](/home/eric/vio_homework/2rd/images/30.png)

#### 	3.3 工程化相机模型

![31](/home/eric/vio_homework/2rd/images/31.png)

![32](/home/eric/vio_homework/2rd/images/32.png)
