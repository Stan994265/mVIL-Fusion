#pragma once

#include <ros/ros.h>
#include <vector>
#include <eigen3/Eigen/Dense>
#include "utility/utility.h"
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <fstream>

const double FOCAL_LENGTH = 460.0; //460
const int WINDOW_SIZE = 6;         //11 7
const int NUM_OF_CAM = 1;
const int NUM_OF_F = 1000;
//#define UNIT_SPHERE_ERROR

extern double INIT_DEPTH;
extern double MIN_PARALLAX;
extern int ESTIMATE_EXTRINSIC;

extern double ACC_N, ACC_W;
extern double GYR_N, GYR_W;

extern std::vector<Eigen::Matrix3d> RIC;
extern std::vector<Eigen::Vector3d> TIC;

extern Eigen::Matrix3d RLI;
extern Eigen::Vector3d TLI;

extern Eigen::Vector3d G;
extern Eigen::Vector3d G_DIRECTION;

extern double PBC_UX, PLB_UX;
extern double PBC_LX, PLB_LX;
extern double PBC_UY, PLB_UY;
extern double PBC_LY, PLB_LY;
extern double PBC_UZ, PLB_UZ;
extern double PBC_LZ, PLB_LZ;

extern int ANGLE_VI;

extern double BIAS_ACC_THRESHOLD;
extern double BIAS_GYR_THRESHOLD;
extern double SOLVER_TIME;
extern int NUM_ITERATIONS;
extern std::string EX_CALIB_RESULT_PATH;
extern std::string VINS_RESULT_PATH;
extern std::string IMU_TOPIC;
extern std::string LIDAR_TOPIC;
extern std::string EX_RESULT_PATH;

extern double TD;
extern double TR;
extern int ESTIMATE_TD;
extern int ROLLING_SHUTTER;
extern double ROW, COL;

extern double LeafSize;
extern int NumThreads;
extern double TransformationEpsilon;
extern double MaxCorrespondenceDistance;

extern double MinDistance;
extern double MaxDistance;

extern double LidarTimeStep;
extern int SHOW_LIDAR_CONSTRAINT;
extern int ADD_LIDAR_ICP;
extern int ADD_LPS;

void readParameters(ros::NodeHandle &n);

enum SIZE_PARAMETERIZATION
{
    SIZE_POSE = 7,
    SIZE_SPEEDBIAS = 9,
    SIZE_FEATURE = 1
};

enum StateOrder
{
    O_P = 0,
    O_R = 3,
    O_V = 6,
    O_BA = 9,
    O_BG = 12
};

enum NoiseOrder
{
    O_AN = 0,
    O_GN = 3,
    O_AW = 6,
    O_GW = 9
};
