#include "parameters.h"

double INIT_DEPTH;
double MIN_PARALLAX;
double ACC_N, ACC_W;
double GYR_N, GYR_W;

double LeafSize;
int NumThreads;
double TransformationEpsilon;
double MaxCorrespondenceDistance;

double MinDistance;
double MaxDistance;
double LidarTimeStep;

double PBC_UX, PLB_UX;
double PBC_LX, PLB_LX;
double PBC_UY, PLB_UY;
double PBC_LY, PLB_LY;
double PBC_UZ, PLB_UZ;
double PBC_LZ, PLB_LZ;

int ANGLE_VI;

std::vector<Eigen::Matrix3d> RIC;
std::vector<Eigen::Vector3d> TIC;

Eigen::Matrix3d RLI;
Eigen::Vector3d TLI;

Eigen::Vector3d G{0.0, 0.0, 9.8};
Eigen::Vector3d G_DIRECTION{0.0, 0.0, 0.0};

double BIAS_ACC_THRESHOLD;
double BIAS_GYR_THRESHOLD;
double SOLVER_TIME;
int NUM_ITERATIONS;
int ESTIMATE_EXTRINSIC;
int ESTIMATE_TD;
int ROLLING_SHUTTER;
std::string EX_CALIB_RESULT_PATH;
std::string VINS_RESULT_PATH;
std::string IMU_TOPIC;
std::string LIDAR_TOPIC;
std::string EX_RESULT_PATH;

double ROW, COL;
double TD, TR;

int SHOW_LIDAR_CONSTRAINT;
int ADD_LIDAR_ICP;
int ADD_LPS;

template <typename T>
T readParam(ros::NodeHandle &n, std::string name)
{
    T ans;
    if (n.getParam(name, ans))
    {
        ROS_INFO_STREAM("Loaded " << name << ": " << ans);
    }
    else
    {
        ROS_ERROR_STREAM("Failed to load " << name);
        n.shutdown();
    }
    return ans;
}

void readParameters(ros::NodeHandle &n)
{
    std::string config_file;
    config_file = readParam<std::string>(n, "config_file");

    std::string OUTPUT_PATH = readParam<std::string>(n, "vils_folder");

    cv::FileStorage fsSettings(config_file, cv::FileStorage::READ);
    if (!fsSettings.isOpened())
    {
        std::cerr << "ERROR: Wrong path to settings" << std::endl;
    }

    fsSettings["imu_topic"] >> IMU_TOPIC;
    fsSettings["lidar_topic"] >> LIDAR_TOPIC;

    SOLVER_TIME = fsSettings["max_solver_time"];
    NUM_ITERATIONS = fsSettings["max_num_iterations"];
    MIN_PARALLAX = fsSettings["keyframe_parallax"];
    MIN_PARALLAX = MIN_PARALLAX / FOCAL_LENGTH;

    // std::string OUTPUT_PATH;
    // fsSettings["output_path"] >> OUTPUT_PATH;
    VINS_RESULT_PATH = OUTPUT_PATH + "Frontend.txt";
    std::cout << "result path " << VINS_RESULT_PATH << std::endl;
    std::ofstream fout(VINS_RESULT_PATH, std::ios::out);
    fout.close();

    ACC_N = fsSettings["acc_n"];
    ACC_W = fsSettings["acc_w"];
    GYR_N = fsSettings["gyr_n"];
    GYR_W = fsSettings["gyr_w"];
    G.z() = fsSettings["g_norm"];

    ROW = fsSettings["image_height"];
    COL = fsSettings["image_width"];
    ROS_INFO("ROW: %f COL: %f ", ROW, COL);

    PBC_UX = fsSettings["PBC_UX"];
    PBC_LX = fsSettings["PBC_LX"];
    PBC_UY = fsSettings["PBC_UY"];
    PBC_LY = fsSettings["PBC_LY"];
    PBC_UZ = fsSettings["PBC_UZ"];
    PBC_LZ = fsSettings["PBC_LZ"];
    PLB_UX = fsSettings["PLB_UX"];
    PLB_LX = fsSettings["PLB_LX"];
    PLB_UY = fsSettings["PLB_UY"];
    PLB_LY = fsSettings["PLB_LY"];
    PLB_UZ = fsSettings["PLB_UZ"];
    PLB_LZ = fsSettings["PLB_LZ"];

    SHOW_LIDAR_CONSTRAINT = fsSettings["show_lidar_constraint"];
    ADD_LIDAR_ICP = fsSettings["add_lidar2lidar"];
    ADD_LPS = fsSettings["add_lps"];

    ANGLE_VI = fsSettings["ANGLE_VI"];

    ESTIMATE_EXTRINSIC = fsSettings["estimate_extrinsic"];

    EX_RESULT_PATH = OUTPUT_PATH + "/ex_results.txt";
    std::cout << "EX result path " << EX_RESULT_PATH << std::endl;
    std::ofstream fout_ex(EX_RESULT_PATH, std::ios::out);
    fout_ex.close();

    {
        cv::Mat cv_RLI, cv_TLI;
        fsSettings["gt_rli"] >> cv_RLI;
        fsSettings["gt_tli"] >> cv_TLI;
        Eigen::Matrix3d eigen_RLI;
        Eigen::Vector3d eigen_TLI;
        cv::cv2eigen(cv_RLI, eigen_RLI);
        cv::cv2eigen(cv_TLI, eigen_TLI);
        Eigen::Quaterniond QLI(eigen_RLI);
        eigen_RLI = QLI.normalized();
        RLI = eigen_RLI;
        TLI = eigen_TLI;
    }

    if (ESTIMATE_EXTRINSIC == 2)
    {
        ROS_WARN("have no prior about extrinsic param, calibrate extrinsic param");
        RIC.push_back(Eigen::Matrix3d::Identity());
        TIC.push_back(Eigen::Vector3d::Zero());
        EX_CALIB_RESULT_PATH = OUTPUT_PATH + "/extrinsic_parameter.csv";

        cv::Mat gd;
        fsSettings["g_direction"] >> gd;
        Eigen::Vector3d eigen_gd;
        cv::cv2eigen(gd, eigen_gd);
        G_DIRECTION = eigen_gd;
    }
    else
    {
        if (ESTIMATE_EXTRINSIC == 1)
        {
            ROS_WARN(" Optimize extrinsic param around initial guess!");
            EX_CALIB_RESULT_PATH = OUTPUT_PATH + "/extrinsic_parameter.csv";
        }
        if (ESTIMATE_EXTRINSIC == 0)
            ROS_WARN(" fix extrinsic param ");

        cv::Mat cv_R, cv_T;
        fsSettings["extrinsicRotation"] >> cv_R;
        fsSettings["extrinsicTranslation"] >> cv_T;
        Eigen::Matrix3d eigen_R;
        Eigen::Vector3d eigen_T;
        cv::cv2eigen(cv_R, eigen_R);
        cv::cv2eigen(cv_T, eigen_T);
        Eigen::Quaterniond Q(eigen_R);
        eigen_R = Q.normalized();
        RIC.push_back(eigen_R);
        TIC.push_back(eigen_T);
        ROS_INFO_STREAM("Extrinsic_R : " << std::endl
                                         << RIC[0]);
        ROS_INFO_STREAM("Extrinsic_T : " << std::endl
                                         << TIC[0].transpose());
    }

    INIT_DEPTH = 5.0;
    BIAS_ACC_THRESHOLD = 0.1;
    BIAS_GYR_THRESHOLD = 0.1;

    TD = fsSettings["td"];
    ESTIMATE_TD = fsSettings["estimate_td"];
    if (ESTIMATE_TD)
        ROS_INFO_STREAM("Unsynchronized sensors, online estimate time offset, initial td: " << TD);
    else
        ROS_INFO_STREAM("Synchronized sensors, fix time offset: " << TD);

    ROLLING_SHUTTER = fsSettings["rolling_shutter"];
    if (ROLLING_SHUTTER)
    {
        TR = fsSettings["rolling_shutter_tr"];
        ROS_INFO_STREAM("rolling shutter camera, read out time per line: " << TR);
    }
    else
    {
        TR = 0;
    }

    LeafSize = fsSettings["LeafSize"];
    NumThreads = fsSettings["NumThreads"];
    TransformationEpsilon = fsSettings["TransformationEpsilon"];
    MaxCorrespondenceDistance = fsSettings["MaxCorrespondenceDistance"];

    MinDistance = fsSettings["MinDistance"];
    MaxDistance = fsSettings["MaxDistance"];
    LidarTimeStep = fsSettings["LidarTimeStep"];

    fsSettings.release();
}
