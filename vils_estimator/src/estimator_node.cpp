#include <stdio.h>
#include <queue>
#include <deque>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
//
#include "estimator.h"
#include "parameters.h"
#include "utility/visualization.h"

#include <sensor_msgs/PointCloud2.h>

Estimator estimator;

std::condition_variable con;
double current_time = -1;
double current_img_time = -1;
queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::PointCloudConstPtr> relo_buf;

deque<CloudData> lidar_buf;
queue<nav_msgs::Odometry::ConstPtr> lps_buf;
// CloudData current_cloud;

int sum_of_wait = 0;

std::mutex m_buf;
std::mutex m_state;
std::mutex l_buf;
std::mutex m_estimator;

double latest_time;
Eigen::Vector3d tmp_P;
Eigen::Quaterniond tmp_Q;
Eigen::Vector3d tmp_V;
Eigen::Vector3d tmp_Ba;
Eigen::Vector3d tmp_Bg;
Eigen::Vector3d acc_0;
Eigen::Vector3d gyr_0;
bool init_feature = 0;
bool init_imu = 1;
double last_imu_t = 0;

double last_lidar_time = 0.0;

void predict(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    if (init_imu)
    {
        latest_time = t;
        init_imu = 0;
        return;
    }
    double dt = t - latest_time;
    latest_time = t;

    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    Eigen::Vector3d linear_acceleration{dx, dy, dz};

    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Eigen::Vector3d angular_velocity{rx, ry, rz};

    Eigen::Vector3d un_acc_0 = tmp_Q * (acc_0 - tmp_Ba) - G;

    Eigen::Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - tmp_Bg;
    tmp_Q = tmp_Q * Utility::deltaQ(un_gyr * dt);

    Eigen::Vector3d un_acc_1 = tmp_Q * (linear_acceleration - tmp_Ba) - G;

    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);

    tmp_P = tmp_P + dt * tmp_V + 0.5 * dt * dt * un_acc;
    tmp_V = tmp_V + dt * un_acc;

    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void update()
{
    TicToc t_predict;
    latest_time = current_time;
    tmp_P = estimator.Ps[WINDOW_SIZE];
    tmp_Q = estimator.Rs[WINDOW_SIZE];
    tmp_V = estimator.Vs[WINDOW_SIZE];
    tmp_Ba = estimator.Bas[WINDOW_SIZE];
    tmp_Bg = estimator.Bgs[WINDOW_SIZE];
    acc_0 = estimator.acc_0;
    gyr_0 = estimator.gyr_0;

    queue<sensor_msgs::ImuConstPtr> tmp_imu_buf = imu_buf;
    for (sensor_msgs::ImuConstPtr tmp_imu_msg; !tmp_imu_buf.empty(); tmp_imu_buf.pop())
        predict(tmp_imu_buf.front());
}

CloudData get_current_cloud(deque<CloudData> &LidarBuf, double FeatureTime)
{
    CloudData CurrentCloud;
    CurrentCloud.valid = false;

    // cout << "lidar_buf size--" << LidarBuf.size() << endl;
    while (LidarBuf.size() > 2)
    {
        LidarBuf.pop_front();
    }
    if (LidarBuf.size() == 2)
    {
        double dt_ = LidarBuf.back().time - FeatureTime;
        //  cout<<"dt="<<dt_<<endl;
        //  time: dt-: --lidar--feature--   or dt+: ---feature---lidar--
        if (dt_ < 0.0 && dt_ > -0.1)
        {
            cout << "---------------NEW--LIDAR------------------" << endl;
            cout << "dt-: --lidar--feature--" << dt_ << endl;
            CurrentCloud = LidarBuf.back();
            CurrentCloud.valid = true;
            LidarBuf.pop_front();
            return CurrentCloud;
        }
        else if (dt_ > 0.0 && dt_ < 0.1)
        {
            cout << "---------------NEW--LIDAR------------------" << endl;
            cout << "dt+: ---feature---lidar--" << dt_ << endl;
            CurrentCloud = LidarBuf.front();
            CurrentCloud.valid = true;
            LidarBuf.pop_front();
            return CurrentCloud;
        }
    }

    return CurrentCloud;
}

std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr> >
getMeasurements()
{
    std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr> > measurements;

    while (true)
    {
        // std::cout<<"Before feature_buf size="<<feature_buf.size()<<std::endl;
        // std::cout<<"Before imu_buf size="<<imu_buf.size()<<std::endl;
        // std::cout<<"Before lidar_buf size="<<lidar_buf.size()<<std::endl;

        if (imu_buf.empty() || feature_buf.empty())
            return measurements;

        if (!(imu_buf.back()->header.stamp.toSec() > feature_buf.front()->header.stamp.toSec() + estimator.td))
        {
            ROS_WARN("wait for imu, only should happen at the beginning");
            sum_of_wait++;
            return measurements;
        }

        if (!(imu_buf.front()->header.stamp.toSec() < feature_buf.front()->header.stamp.toSec() + estimator.td))
        {
            ROS_WARN("throw img, only should happen at the beginning");
            feature_buf.pop();
            continue;
        }

        //std::cout<<"BEFORE: lidar_buf size="<<lidar_buf.size()<<std::endl;
        // current_cloud.valid = false;
        // if (lidar_buf.size() > 2 && feature_buf.size() > 0)
        // {
        //     lidar_buf.pop_front();
        //     double dt_ = lidar_buf.back().time - feature_buf.front()->header.stamp.toSec() - estimator.td;
        //     //  cout<<"dt="<<dt_<<endl;
        //     //  time: dt-: --lidar--feature--   or dt+: ---feature---lidar--
        //     if (dt_ < 0.0 && dt_ > -LidarTimeStep)
        //     {
        //         cout << "---------------NEW--LIDAR------------------" << endl;
        //         cout << "dt-: --lidar--feature--" << dt_ << endl;
        //         current_cloud = lidar_buf.back();
        //         current_cloud.valid = true;
        //         lidar_buf.pop_front();
        //     }
        //     if (dt_ > 0.0 && dt_ < LidarTimeStep)
        //     {
        //         cout << "---------------NEW--LIDAR------------------" << endl;
        //         cout << "dt+: ---feature---lidar--" << dt_ << endl;
        //         current_cloud = lidar_buf.front();
        //         current_cloud.valid = true;
        //         lidar_buf.pop_front();
        //     }
        // }

        sensor_msgs::PointCloudConstPtr img_msg = feature_buf.front();
        feature_buf.pop();

        std::vector<sensor_msgs::ImuConstPtr> IMUs;
        while (imu_buf.front()->header.stamp.toSec() < img_msg->header.stamp.toSec() + estimator.td)
        {
            IMUs.emplace_back(imu_buf.front());
            imu_buf.pop();
        }
        IMUs.emplace_back(imu_buf.front());
        if (IMUs.empty())
            ROS_WARN("no imu between two image");
        measurements.emplace_back(IMUs, img_msg);

        // std::cout<<"AFTER: lidar_buf size="<<lidar_buf.size()<<std::endl;
        // std::cout<<"After feature_buf size="<<feature_buf.size()<<std::endl;
        // std::cout<<"After imu_buf size="<<imu_buf.size()<<std::endl;

        // l_buf.lock();
        // double tem_t = img_msg->header.stamp.toSec() + estimator.td;
        // current_cloud = get_current_cloud(lidar_buf, tem_t);
        // l_buf.unlock();
    }
    return measurements;
}

void lidar_callback(const sensor_msgs::PointCloud2ConstPtr &lidar_msg)
{
    CloudData cloud_data;
    cloud_data.time = lidar_msg->header.stamp.toSec();
    cloud_data.header = lidar_msg->header;
    pcl::fromROSMsg(*lidar_msg, *(cloud_data.cloud_ptr));
    l_buf.lock();
    lidar_buf.push_back(cloud_data);
    l_buf.unlock();
}

void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    if (imu_msg->header.stamp.toSec() <= last_imu_t)
    {
        ROS_WARN("imu message in disorder!");
        return;
    }
    last_imu_t = imu_msg->header.stamp.toSec();

    m_buf.lock();
    imu_buf.push(imu_msg);
    m_buf.unlock();
    con.notify_one();

    last_imu_t = imu_msg->header.stamp.toSec();

    {
        std::lock_guard<std::mutex> lg(m_state);
        predict(imu_msg);
        std_msgs::Header header = imu_msg->header;
        header.frame_id = "world";
        if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
            pubLatestOdometry(estimator.RLB.transpose(), estimator.TBL, tmp_P, tmp_Q, tmp_V, header);
    }
}

void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    if (!init_feature)
    {
        //skip the first detected feature, which doesn't contain optical flow speed
        init_feature = 1;
        return;
    }
    m_buf.lock();
    feature_buf.push(feature_msg);
    m_buf.unlock();
    con.notify_one();
}

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        m_buf.lock();
        while (!feature_buf.empty())
            feature_buf.pop();
        while (!imu_buf.empty())
            imu_buf.pop();
        m_buf.unlock();
        m_estimator.lock();
        estimator.clearState();
        estimator.setParameter();
        m_estimator.unlock();
        current_time = -1;
        last_imu_t = 0;
    }
    return;
}

void relocalization_callback(const sensor_msgs::PointCloudConstPtr &points_msg)
{
    //printf("relocalization callback! \n");
    m_buf.lock();
    relo_buf.push(points_msg);
    m_buf.unlock();
}

void LPS_callback(const nav_msgs::Odometry::ConstPtr &LPS_msg)
{
    m_buf.lock();
    lps_buf.push(LPS_msg);
    m_buf.unlock();
}

void pubEXresults(const Estimator &estimator, const std_msgs::Header &header)
{
    Eigen::Vector3d VI_r = Utility::R2ypr(estimator.ric[0]);
    Eigen::Vector3d VI_t = estimator.tic[0];
    Eigen::Vector3d LI_r = Utility::R2ypr(estimator.RLB);
    Eigen::Vector3d LI_t = estimator.TLB;
    ofstream foutC(EX_RESULT_PATH, ios::app);
    foutC.setf(ios::fixed, ios::floatfield);
    // foutC.precision(1);
    foutC << header.stamp.toSec() << " ";
    // foutC.precision(5);
    foutC << VI_r[0] << " "
          << VI_r[1] << " "
          << VI_r[2] << " "
          << VI_t[0] << " "
          << VI_t[1] << " "
          << VI_t[2] << " "
          << LI_r[0] << " "
          << LI_r[1] << " "
          << LI_r[2] << " "
          << LI_t[0] << " "
          << LI_t[1] << " "
          << LI_t[2] << endl;
    foutC.close();
}

void process_lidar()
{
    while (true)
    {
        CloudData current_cloud;
        l_buf.lock();
        // double tem_t = current_img_time + estimator.td;
        current_cloud = get_current_cloud(lidar_buf, current_img_time);
        l_buf.unlock();

        if (current_cloud.valid && last_lidar_time != current_cloud.time)
        {

            // l_buf.lock();
            m_buf.lock();
            m_state.lock();
            m_estimator.lock();
            estimator.processLidar(current_cloud, estimator.td);

            if (estimator.lidarCalibration.lidar_init_flag == false && estimator.current_lidar.keylidar != false)
            {
                pubLidarPose(estimator, current_cloud.header);
                pubLidarCloud(current_cloud, current_cloud.header);
                last_lidar_time = current_cloud.time;
                if (SHOW_LIDAR_CONSTRAINT)
                    pubLidarICPConstraintMarker(estimator, current_cloud.header);
            }
            m_estimator.unlock();
            m_state.unlock();
            m_buf.unlock();
            // l_buf.unlock();

            current_cloud.valid = false;
        }
        // m_estimator.unlock();
        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }
}

// thread: visual-inertial odometry
void process()
{
    while (true)
    {
        std::vector<std::pair<std::vector<sensor_msgs::ImuConstPtr>, sensor_msgs::PointCloudConstPtr> > measurements;
        std::unique_lock<std::mutex> lk(m_buf);
        con.wait(lk, [&]
                 { return (measurements = getMeasurements()).size() != 0; });
        lk.unlock();
        m_estimator.lock();

        for (auto &measurement : measurements)
        {
            auto img_msg = measurement.second;
            double dx = 0, dy = 0, dz = 0, rx = 0, ry = 0, rz = 0;
            for (auto &imu_msg : measurement.first)
            {
                double t = imu_msg->header.stamp.toSec();
                double img_t = img_msg->header.stamp.toSec() + estimator.td;
                if (t <= img_t)
                {
                    if (current_time < 0)
                        current_time = t;
                    double dt = t - current_time;
                    ROS_ASSERT(dt >= 0);
                    current_time = t;
                    dx = imu_msg->linear_acceleration.x;
                    dy = imu_msg->linear_acceleration.y;
                    dz = imu_msg->linear_acceleration.z;
                    rx = imu_msg->angular_velocity.x;
                    ry = imu_msg->angular_velocity.y;
                    rz = imu_msg->angular_velocity.z;
                    estimator.processIMU(dt, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
                    //printf("imu: dt:%f a: %f %f %f w: %f %f %f\n",dt, dx, dy, dz, rx, ry, rz);
                }
                else
                {
                    double dt_1 = img_t - current_time;
                    double dt_2 = t - img_t;
                    current_time = img_t;
                    ROS_ASSERT(dt_1 >= 0);
                    ROS_ASSERT(dt_2 >= 0);
                    ROS_ASSERT(dt_1 + dt_2 > 0);
                    double w1 = dt_2 / (dt_1 + dt_2);
                    double w2 = dt_1 / (dt_1 + dt_2);
                    dx = w1 * dx + w2 * imu_msg->linear_acceleration.x;
                    dy = w1 * dy + w2 * imu_msg->linear_acceleration.y;
                    dz = w1 * dz + w2 * imu_msg->linear_acceleration.z;
                    rx = w1 * rx + w2 * imu_msg->angular_velocity.x;
                    ry = w1 * ry + w2 * imu_msg->angular_velocity.y;
                    rz = w1 * rz + w2 * imu_msg->angular_velocity.z;
                    estimator.processIMU(dt_1, Vector3d(dx, dy, dz), Vector3d(rx, ry, rz));
                    //printf("dimu: dt:%f a: %f %f %f w: %f %f %f\n",dt_1, dx, dy, dz, rx, ry, rz);
                }
            }
            // set relocalization frame
            sensor_msgs::PointCloudConstPtr relo_msg = NULL;
            while (!relo_buf.empty())
            {
                relo_msg = relo_buf.front();
                relo_buf.pop();
            }

            nav_msgs::Odometry::ConstPtr lps_msg = NULL;
            while (!lps_buf.empty())
            {
                lps_msg = lps_buf.front();
                lps_buf.pop();
            }
            if (lps_msg != NULL)
            {
                estimator.LPS_call = true;
                estimator.LPS_t.x() = lps_msg->pose.pose.position.x;
                estimator.LPS_t.y() = lps_msg->pose.pose.position.y;
                estimator.LPS_t.z() = lps_msg->pose.pose.position.z;
                estimator.LPS_q.w() = lps_msg->pose.pose.orientation.w;
                estimator.LPS_q.x() = lps_msg->pose.pose.orientation.x;
                estimator.LPS_q.y() = lps_msg->pose.pose.orientation.y;
                estimator.LPS_q.z() = lps_msg->pose.pose.orientation.z;
                estimator.LPS_time = lps_msg->header.stamp.toSec();
                ROS_WARN_STREAM("GET LPS!!! time:" << fixed << estimator.LPS_time << endl);
            }

            if (relo_msg != NULL)
            {
                vector<Vector3d> match_points;
                double frame_stamp = relo_msg->header.stamp.toSec();
                for (unsigned int i = 0; i < relo_msg->points.size(); i++)
                {
                    Vector3d u_v_id;
                    u_v_id.x() = relo_msg->points[i].x;
                    u_v_id.y() = relo_msg->points[i].y;
                    u_v_id.z() = relo_msg->points[i].z;
                    match_points.push_back(u_v_id);
                }
                Vector3d relo_t(relo_msg->channels[0].values[0], relo_msg->channels[0].values[1], relo_msg->channels[0].values[2]);
                Quaterniond relo_q(relo_msg->channels[0].values[3], relo_msg->channels[0].values[4], relo_msg->channels[0].values[5], relo_msg->channels[0].values[6]);
                Matrix3d relo_r = relo_q.toRotationMatrix();
                int frame_index;
                frame_index = relo_msg->channels[0].values[7];
                estimator.setReloFrame(frame_stamp, frame_index, match_points, relo_t, relo_r);
            }

            ROS_DEBUG("processing vision data with stamp %f \n", img_msg->header.stamp.toSec());

            TicToc t_s;
            map<int, vector<pair<int, Eigen::Matrix<double, 8, 1> > > > image;
            for (unsigned int i = 0; i < img_msg->points.size(); i++)
            {
                int v = img_msg->channels[0].values[i] + 0.5;
                int feature_id = v / NUM_OF_CAM;
                int camera_id = v % NUM_OF_CAM;
                double x = img_msg->points[i].x;
                double y = img_msg->points[i].y;
                double z = img_msg->points[i].z;
                double p_u = img_msg->channels[1].values[i];
                double p_v = img_msg->channels[2].values[i];
                double velocity_x = img_msg->channels[3].values[i];
                double velocity_y = img_msg->channels[4].values[i];
                double depth = img_msg->channels[5].values[i];
                ROS_ASSERT(z == 1);
                Eigen::Matrix<double, 8, 1> xyz_uv_velocity_depth;
                xyz_uv_velocity_depth << x, y, z, p_u, p_v, velocity_x, velocity_y, depth;
                image[feature_id].emplace_back(camera_id, xyz_uv_velocity_depth);
            }

            estimator.processImage(image, img_msg->header);
            current_img_time = img_msg->header.stamp.toSec() + estimator.td;

            double whole_t = t_s.toc();
            printStatistics(estimator, whole_t);
            std_msgs::Header header = img_msg->header;
            header.frame_id = "world";
            if (estimator.solver_flag != estimator.INITIAL)
                pubOdometry(estimator, header);
            pubKeyPoses(estimator, header);
            pubCameraPose(estimator, header);
            pubPointCloud(estimator, header);
            pubTF(estimator, header);
            pubKeyframe(estimator);
            pubEXresults(estimator, header);
            if (relo_msg != NULL)
                pubRelocalization(estimator);
            //ROS_ERROR("end: %f, at %f", img_msg->header.stamp.toSec(), ros::Time::now().toSec());
        }
        m_estimator.unlock();
        m_buf.lock();
        m_state.lock();
        if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
            update();
        m_state.unlock();
        m_buf.unlock();

        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "vins_estimator");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);
    readParameters(n);
    estimator.setParameter();

#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif
    ROS_WARN("waiting for image and imu...");

    registerPub(n);

    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay());
    ros::Subscriber sub_image = n.subscribe("/feature_tracker_/feature", 10, feature_callback, ros::TransportHints().tcpNoDelay());
    // ros::Subscriber sub_restart = n.subscribe("/feature_tracker/restart", 5, restart_callback);
    // ros::Subscriber sub_relo_points = n.subscribe("/pose_graph/match_points", 5, relocalization_callback);

    ros::Subscriber sub_lidar_points = n.subscribe(LIDAR_TOPIC, 5, lidar_callback, ros::TransportHints().tcpNoDelay());
    //aft_mapped_to_init  laser_localizer  aft_mapped_to_init_high_frec
    ros::Subscriber sub_LPSdometry = n.subscribe<nav_msgs::Odometry>("/aft_mapped_to_init", 5, LPS_callback, ros::TransportHints().tcpNoDelay());
    if (!ADD_LPS)
    {
        sub_LPSdometry.shutdown();
    }

    std::thread measurement_process{process};
    std::thread lidar_process{process_lidar};

    ros::AsyncSpinner spinner(4);
    spinner.start();
    ros::waitForShutdown();
    return 0;
}