
#ifndef LIDAR_FRONT_END_FRONT_END_H
#define LIDAR_FRONT_END_FRONT_END_H

#include <deque>
#include <queue>
#include <Eigen/Dense>

#include "initial/initial_alignment.h"
#include <chrono>
#include <iostream>

#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/ndt.h>
#include <pcl/registration/gicp.h>
#include "lidar_functions/fast_gicp/include/fast_gicp/gicp/fast_gicp.hpp"
#include "lidar_functions/fast_gicp/include/fast_gicp/gicp/fast_gicp_st.hpp"
#include "lidar_functions/fast_gicp/include/fast_gicp/gicp/fast_vgicp.hpp"
#include <eigen3/Eigen/Dense>

#include <ros/ros.h>
#include <std_msgs/Header.h>
#include <std_msgs/Float32.h>
#include <std_msgs/Bool.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud.h>
#include <sensor_msgs/Image.h>
#include <sensor_msgs/image_encodings.h>
#include <nav_msgs/Path.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/PointStamped.h>
#include <visualization_msgs/Marker.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/PointCloud2.h>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <omp.h>

// #include "utility/geometry_utils.h"

class CloudData
{
public:
  // EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  //using POINT = pcl::PointXYZ;
  using POINT = pcl::PointXYZI;
  using CLOUD = pcl::PointCloud<POINT>;
  using CLOUD_PTR = CLOUD::Ptr;

public:
  CloudData()
      : cloud_ptr(new CLOUD())
  {
  }

public:
  double time = 0.0;
  CLOUD_PTR cloud_ptr;
  bool valid = false;
  std_msgs::Header header;
};

class LidarFrame
{
public:
  // EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  LidarFrame(){};
  LidarFrame(const CloudData &_cloud) : frameID(0), keylidar(false)
  {
    lidarData.point_cloud = _cloud;
    lidarData.lidar_R = Matrix3d::Identity();
    lidarData.lidar_T = Vector3d::Zero();
    lidarData.mode = 0;
    vioData.acci = Vector3d::Zero();
    vioData.accj = Vector3d::Zero();
    vioData.gyri = Vector3d::Zero();
    vioData.gyrj = Vector3d::Zero();
    vioData.Rwbi = Matrix3d::Identity();
    vioData.Pwbi = Vector3d::Zero();
    vioData.Rwbj = Matrix3d::Identity();
    vioData.Pwbj = Vector3d::Zero();
    vioData.Vbi = Vector3d::Zero();
    vioData.Vbj = Vector3d::Zero();
    vioData.ti = 0.0;
    vioData.tj = 0.0;
    vioData.dt = 0.0;
  };

  int frameID;
  bool keylidar;

  struct LidarData
  {
    ImageFrame last_image;
    ImageFrame next_image;
    CloudData point_cloud;
    Matrix3d lidar_R;
    Vector3d lidar_T;
    int mode;
  };
  struct VioData
  {
    Vector3d acci;
    Vector3d accj;
    Vector3d gyri;
    Vector3d gyrj;
    Matrix3d Rwbi;
    Vector3d Pwbi;
    Matrix3d Rwbj;
    Vector3d Pwbj;
    Vector3d Vbi;
    Vector3d Vbj;
    double ti;
    double tj;
    double dt;
  };

  LidarData lidarData;
  VioData vioData;
};

class LidarCalibration
{
public:
  LidarCalibration();
  bool CalibrationLidarExRotation(const LidarFrame &lidarframei, const LidarFrame &lidarframej, double fitness_core,
                                  Matrix3d delta_lidar, Matrix3d &calib_rlb_result,
                                  Vector3d delta_lidar_t, Vector3d &calib_tlb_result);
  void Reset();
  bool lidar_init_flag;

private:
  inline bool Lidar_align(vector<Matrix3d> &Rlidar, vector<Vector3d> &Tlidar, vector<Vector3d> &Tbody, Vector3d &Tlb, Matrix3d &Rlb, int frame_count);
  inline Matrix3d PredictRelativeR(const LidarFrame &lidarframei, const LidarFrame &lidarframej);
  inline Matrix3d PredictR(const LidarFrame &lidarframe);
  inline Vector3d PredictRelativeT(const LidarFrame &lidarframei, const LidarFrame &lidarframej);
  inline Vector3d PredictT(const LidarFrame &lidarframe);
  int frame_count;
  Vector3d temp_ypr;
  vector<Matrix3d> Rcam;
  vector<Matrix3d> Rlidar;
  vector<Matrix3d> Rc_g;
  vector<Vector3d> Tlidar;
  vector<Vector3d> Tbody;
  vector<double> FS;
  Matrix3d rlb;
};

struct LidarInitConstraint
{

  LidarInitConstraint(const Matrix3d rl_ij_,
                      const Vector3d tl_ij_, const Vector3d tb_ij_, const Matrix3d rb_ij_)

      : rl_ij(rl_ij_), tl_ij(tl_ij_), tb_ij(tb_ij_), rb_ij(rb_ij_)
  {
  }

  template <typename T>
  bool operator()(const T *RLB, const T *TLB, T *residuals) const
  {

    Eigen::Matrix<T, 3, 3> Rlij;
    Rlij << T(rl_ij(0, 0)), T(rl_ij(0, 1)), T(rl_ij(0, 2)),
        T(rl_ij(1, 0)), T(rl_ij(1, 1)), T(rl_ij(1, 2)),
        T(rl_ij(2, 0)), T(rl_ij(2, 1)), T(rl_ij(2, 2));
    Eigen::Matrix<T, 3, 3> Rbij;
    Rbij << T(rb_ij(0, 0)), T(rb_ij(0, 1)), T(rb_ij(0, 2)),
        T(rb_ij(1, 0)), T(rb_ij(1, 1)), T(rb_ij(1, 2)),
        T(rb_ij(2, 0)), T(rb_ij(2, 1)), T(rb_ij(2, 2));
    Eigen::Matrix<T, 3, 1> Tlij{T(tl_ij[0]), T(tl_ij[1]), T(tl_ij[2])};
    T Tbij[3]{T(tb_ij[0]), T(tb_ij[1]), T(tb_ij[2])};
    Eigen::Matrix<T, 3, 1> Tlb{TLB[0], TLB[1], TLB[2]};

    T RlbTij_[3];
    ceres::QuaternionRotatePoint(RLB, Tbij, RlbTij_);

    T temr[3 * 3];
    ceres::QuaternionToRotation(RLB, temr);
    Eigen::Matrix<T, 3, 3> RLB_;
    RLB_ << temr[0], temr[1], temr[2],
        temr[3], temr[4], temr[5],
        temr[6], temr[7], temr[8];

    Eigen::Matrix<T, 3, 3> RES_;
    RES_ = RLB_ * Rbij * RLB_.transpose() * Rlij.transpose();
    T temres[3 * 3];

    temres[0] = RES_(0, 0);
    temres[1] = RES_(0, 1);
    temres[2] = RES_(0, 2);
    temres[3] = RES_(1, 0);
    temres[4] = RES_(1, 1);
    temres[5] = RES_(1, 2);
    temres[6] = RES_(2, 0);
    temres[7] = RES_(2, 1);
    temres[8] = RES_(2, 2);
    T res_[4];
    ceres::RotationMatrixToQuaternion(temres, res_);

    Eigen::Matrix<T, 3, 3> I;
    I << T(1.0), T(0.0), T(0.0),
        T(0.0), T(1.0), T(0.0),
        T(0.0), T(0.0), T(1.0);

    Eigen::Matrix<T, 3, 1> RlbTij{RlbTij_[0], RlbTij_[1], RlbTij_[2]};
    Eigen::Matrix<T, 3, 1> RES;
    RES = (I - Rlij) * Tlb - Tlij + RlbTij;

    residuals[0] = RES(0, 0);
    residuals[1] = RES(1, 0);
    residuals[2] = RES(2, 0);
    residuals[3] = T(2) * res_[1];
    residuals[4] = T(2) * res_[2];
    residuals[5] = T(2) * res_[3];
    return true;
  }

  static ceres::CostFunction *makeConstraint(const Matrix3d rl_ij_, const Vector3d tl_ij_, const Vector3d tb_ij_, const Matrix3d rb_ij_)
  {
    return (new ceres::AutoDiffCostFunction<
            LidarInitConstraint, 6, 4, 3>(
        new LidarInitConstraint(rl_ij_, tl_ij_, tb_ij_, rb_ij_)));
  }

  //observation
  const Matrix3d rl_ij;
  const Vector3d tl_ij;
  const Vector3d tb_ij;
  const Matrix3d rb_ij;
};

class VelocityData
{
public:
  struct LinearVelocity
  {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };

  struct AngularVelocity
  {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };

  double time = 0.0;
  LinearVelocity linear_velocity;
  AngularVelocity angular_velocity;

public:
  bool PredictData(const LidarFrame &lidarframe);
  void TransformCoordinate(Eigen::Matrix3d rotation, Eigen::Vector3d translation);
};

class DistortionAdjust
{
public:
  void SetMotionInfo(float scan_period, VelocityData velocity_data);
  bool AdjustCloud(CloudData::CLOUD_PTR &input_cloud_ptr, CloudData::CLOUD_PTR &output_cloud_ptr);

private:
  inline Eigen::Matrix3f UpdateMatrix(float real_time);

private:
  float scan_period_;
  Eigen::Vector3f velocity_;
  Eigen::Vector3f angular_rate_;
};

void calculate_ICP_COV(CloudData::CLOUD &data_pi, CloudData::CLOUD &model_qi, Eigen::Matrix4f transform, Eigen::MatrixXd &ICP_COV);
Matrix4d PredictRelative_rt(LidarFrame &lidarframei, LidarFrame &lidarframej, const Matrix3d RBL, const Vector3d TBL);
Matrix3d Predict_r(const LidarFrame &lidarframe);
Vector3d Predict_t(const LidarFrame &lidarframe);

void RotatePoint(const Eigen::Quaternionf &q, CloudData::POINT &p);
void TransformToEnd(CloudData::CLOUD_PTR &cloud, Eigen::Quaternionf transform_es_q, Eigen::Vector3f transform_es_t, float time_factor, double min, double max);
#endif