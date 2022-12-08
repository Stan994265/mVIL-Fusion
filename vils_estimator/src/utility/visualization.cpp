#include "visualization.h"

using namespace ros;
using namespace Eigen;
ros::Publisher pub_odometry, pub_latest_odometry;
ros::Publisher pub_path, pub_relo_path, pub_lidar_path;
ros::Publisher pub_point_cloud, pub_margin_cloud, pub_lidar_cloud;
ros::Publisher pub_key_poses;
ros::Publisher pub_relo_relative_pose;
ros::Publisher pub_camera_pose, pub_lidar_pose;
ros::Publisher pub_camera_pose_visual;
nav_msgs::Path path, relo_path, lidar_path;

ros::Publisher pub_keyframe_pose;
ros::Publisher pub_keyframe_point;
ros::Publisher pub_extrinsic;

ros::Publisher pub_fake_trigger, pub_fake_local_odom; //for path planning test

CameraPoseVisualization cameraposevisual(0, 1, 0, 1);
CameraPoseVisualization keyframebasevisual(0.0, 0.0, 1.0, 1.0);
static double sum_of_path = 0;
static Vector3d last_path(0.0, 0.0, 0.0);

ros::Publisher pub_ConstraintMarker;
visualization_msgs::Marker ConstraintMarker;
ros::Publisher pub_sum_mode;
visualization_msgs::Marker sum_mode_Marker;
int sum_mode1{0}, sum_mode2{0}, sum_mode3{0}, sum_mode4{0}, sum_mode5{0};

void registerPub(ros::NodeHandle &n)
{
    pub_latest_odometry = n.advertise<nav_msgs::Odometry>("imu_propagate", 10);
    pub_path = n.advertise<nav_msgs::Path>("path", 10);
    pub_relo_path = n.advertise<nav_msgs::Path>("relocalization_path", 10);
    pub_lidar_path = n.advertise<nav_msgs::Path>("lidar_path", 10);
    pub_odometry = n.advertise<nav_msgs::Odometry>("odometry", 10);
    pub_point_cloud = n.advertise<sensor_msgs::PointCloud>("point_cloud", 10);
    pub_margin_cloud = n.advertise<sensor_msgs::PointCloud>("history_cloud", 10);
    pub_lidar_cloud = n.advertise<sensor_msgs::PointCloud2>("lidar_cloud", 10);
    pub_lidar_pose = n.advertise<nav_msgs::Odometry>("lidar_pose", 10);
    pub_key_poses = n.advertise<visualization_msgs::Marker>("key_poses", 10);
    pub_camera_pose = n.advertise<nav_msgs::Odometry>("camera_pose", 1000);
    pub_camera_pose_visual = n.advertise<visualization_msgs::MarkerArray>("camera_pose_visual", 10);
    pub_keyframe_pose = n.advertise<nav_msgs::Odometry>("keyframe_pose", 10);
    pub_keyframe_point = n.advertise<sensor_msgs::PointCloud>("keyframe_point", 10);
    pub_extrinsic = n.advertise<nav_msgs::Odometry>("extrinsic", 10);
    pub_relo_relative_pose = n.advertise<nav_msgs::Odometry>("relo_relative_pose", 10);

    pub_fake_trigger = n.advertise<geometry_msgs::PoseStamped>("traj_start_trigger", 10);
    pub_fake_local_odom = n.advertise<nav_msgs::Odometry>("fake_odom", 10);

    pub_ConstraintMarker = n.advertise<visualization_msgs::Marker>("constraint_mode", 10);
    pub_sum_mode = n.advertise<visualization_msgs::Marker>("sum_mode", 10);

    cameraposevisual.setScale(1);
    cameraposevisual.setLineWidth(0.05);
    keyframebasevisual.setScale(0.1);
    keyframebasevisual.setLineWidth(0.01);
}

void pubLatestOdometry(const Eigen::Matrix3d &ex_bl_r, const Eigen::Vector3d &ex_bl_t, const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, const std_msgs::Header &header)
{
    Eigen::Quaterniond quadrotor_Q = Q;

    nav_msgs::Odometry odometry;
    odometry.header = header;
    odometry.header.frame_id = "world";
    // odometry.child_frame_id = "body";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = quadrotor_Q.x();
    odometry.pose.pose.orientation.y = quadrotor_Q.y();
    odometry.pose.pose.orientation.z = quadrotor_Q.z();
    odometry.pose.pose.orientation.w = quadrotor_Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_latest_odometry.publish(odometry);

    Eigen::Quaterniond lidar_q = quadrotor_Q * Quaterniond(ex_bl_r);
    Eigen::Vector3d lidar_p = P + quadrotor_Q.toRotationMatrix() * ex_bl_t;

    // FOR ZED-VELODYNE MVILSs
    // Eigen::Matrix3d RBL;
    // Eigen::Vector3d TBL( 0.0,  -0.09, -0.14);//bl
    // RBL <<  0.0, 1.0,   0.0,
    //         -1.0, 0.0,  0.0,
    //         0.0,  0.0,  1.0;
    // Eigen::Quaterniond lidar_q = quadrotor_Q * Quaterniond(RBL);
    // Eigen::Vector3d lidar_p = P + quadrotor_Q.toRotationMatrix() * TBL;

    static tf::TransformBroadcaster br;
    tf::Transform transform;
    tf::Quaternion transformq;
    transform.setOrigin(tf::Vector3(lidar_p.x(),
                                    lidar_p.y(),
                                    lidar_p.z()));
    transformq.setW(lidar_q.w());
    transformq.setX(lidar_q.x());
    transformq.setY(lidar_q.y());
    transformq.setZ(lidar_q.z());
    transform.setRotation(transformq);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "world", "lidar_tem"));
}

void printStatistics(const Estimator &estimator, double t)
{
    if (estimator.solver_flag != Estimator::SolverFlag::NON_LINEAR)
        return;
    printf("position: %f, %f, %f\r", estimator.Ps[WINDOW_SIZE].x(), estimator.Ps[WINDOW_SIZE].y(), estimator.Ps[WINDOW_SIZE].z());
    ROS_DEBUG_STREAM("position: " << estimator.Ps[WINDOW_SIZE].transpose());
    ROS_DEBUG_STREAM("orientation: " << estimator.Vs[WINDOW_SIZE].transpose());
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        //ROS_DEBUG("calibration result for camera %d", i);
        ROS_DEBUG_STREAM("extirnsic tic: " << estimator.tic[i].transpose());
        ROS_DEBUG_STREAM("extrinsic ric: " << Utility::R2ypr(estimator.ric[i]).transpose());
        if (ESTIMATE_EXTRINSIC)
        {
            cv::FileStorage fs(EX_CALIB_RESULT_PATH, cv::FileStorage::WRITE);
            Eigen::Matrix3d eigen_R;
            Eigen::Vector3d eigen_T;
            eigen_R = estimator.ric[i];
            eigen_T = estimator.tic[i];
            cv::Mat cv_R, cv_T;
            cv::eigen2cv(eigen_R, cv_R);
            cv::eigen2cv(eigen_T, cv_T);
            fs << "extrinsicRotation" << cv_R << "extrinsicTranslation" << cv_T;
            fs.release();
        }
    }

    static double sum_of_time = 0;
    static int sum_of_calculation = 0;
    sum_of_time += t;
    sum_of_calculation++;
    ROS_DEBUG("vo solver costs: %f ms", t);
    ROS_DEBUG("average of time %f ms", sum_of_time / sum_of_calculation);

    sum_of_path += (estimator.Ps[WINDOW_SIZE] - last_path).norm();
    last_path = estimator.Ps[WINDOW_SIZE];
    ROS_DEBUG("sum of path %f", sum_of_path);
    if (ESTIMATE_TD)
        ROS_INFO("td %f", estimator.td);
}

void pubOdometry(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {
        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]);
        odometry.pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
        odometry.pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
        odometry.pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
        odometry.pose.pose.orientation.x = tmp_Q.x();
        odometry.pose.pose.orientation.y = tmp_Q.y();
        odometry.pose.pose.orientation.z = tmp_Q.z();
        odometry.pose.pose.orientation.w = tmp_Q.w();
        odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].x();
        odometry.twist.twist.linear.y = estimator.Vs[WINDOW_SIZE].y();
        odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();
        pub_odometry.publish(odometry);

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";
        pose_stamped.pose = odometry.pose.pose;
        path.header = header;
        path.header.frame_id = "world";
        path.poses.push_back(pose_stamped);
        pub_path.publish(path);

        Vector3d correct_t;
        Vector3d correct_v;
        Quaterniond correct_q;
        correct_t = estimator.drift_correct_r * estimator.Ps[WINDOW_SIZE] + estimator.drift_correct_t;
        correct_q = estimator.drift_correct_r * estimator.Rs[WINDOW_SIZE];
        odometry.pose.pose.position.x = correct_t.x();
        odometry.pose.pose.position.y = correct_t.y();
        odometry.pose.pose.position.z = correct_t.z();
        odometry.pose.pose.orientation.x = correct_q.x();
        odometry.pose.pose.orientation.y = correct_q.y();
        odometry.pose.pose.orientation.z = correct_q.z();
        odometry.pose.pose.orientation.w = correct_q.w();

        pose_stamped.pose = odometry.pose.pose;
        relo_path.header = header;
        relo_path.header.frame_id = "world";
        relo_path.poses.push_back(pose_stamped);
        pub_relo_path.publish(relo_path);

        // write result to file
        ofstream foutC(VINS_RESULT_PATH, ios::app);
        foutC.setf(ios::fixed, ios::floatfield);
        foutC.precision(9);
        foutC << header.stamp.toSec() << " ";
        foutC.precision(5);
        foutC << estimator.Ps[WINDOW_SIZE].x() << " "
              << estimator.Ps[WINDOW_SIZE].y() << " "
              << estimator.Ps[WINDOW_SIZE].z() << " "
              << tmp_Q.x() << " "
              << tmp_Q.y() << " "
              << tmp_Q.z() << " "
              << tmp_Q.w() << endl;
        foutC.close();
    }
}

void pubKeyPoses(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.key_poses.size() == 0)
        return;
    visualization_msgs::Marker key_poses;
    key_poses.header = header;
    key_poses.header.frame_id = "world";
    key_poses.ns = "key_poses";
    key_poses.type = visualization_msgs::Marker::SPHERE_LIST;
    key_poses.action = visualization_msgs::Marker::ADD;
    key_poses.pose.orientation.w = 1.0;
    key_poses.lifetime = ros::Duration();

    //static int key_poses_id = 0;
    key_poses.id = 0; //key_poses_id++;
    key_poses.scale.x = 0.05;
    key_poses.scale.y = 0.05;
    key_poses.scale.z = 0.05;
    key_poses.color.r = 1.0;
    key_poses.color.a = 1.0;

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        geometry_msgs::Point pose_marker;
        Vector3d correct_pose;
        correct_pose = estimator.key_poses[i];
        pose_marker.x = correct_pose.x();
        pose_marker.y = correct_pose.y();
        pose_marker.z = correct_pose.z();
        key_poses.points.push_back(pose_marker);
    }
    pub_key_poses.publish(key_poses);
}

void pubLidarICPConstraintMarker(const Estimator &estimator, const std_msgs::Header &header)
{

    // visualization_msgs::Marker ConstraintMarker;
    ConstraintMarker.header = header;
    ConstraintMarker.header.frame_id = "world";
    ConstraintMarker.ns = "key_poses";
    ConstraintMarker.id = 1;
    ConstraintMarker.type = visualization_msgs::Marker::CUBE_LIST;
    ConstraintMarker.action = visualization_msgs::Marker::ADD;
    ConstraintMarker.pose.orientation.w = 1.0;
    ConstraintMarker.lifetime = ros::Duration();
    ConstraintMarker.scale.x = 0.1;
    ConstraintMarker.scale.y = 0.1;
    ConstraintMarker.scale.z = 0.1;
    std_msgs::ColorRGBA msg_;
    switch (estimator.current_lidar.lidarData.mode)
    {
    case 1:
        msg_.r = 1.0;
        msg_.g = 1.0;
        msg_.b = 1.0;
        msg_.a = 1.0;
        sum_mode1++;
        break;
    case 2: // good vio predict (green)
        msg_.g = 1.0;
        msg_.a = 1.0;
        sum_mode2++;
        break;
    case 3: //add icp constraint (red)
        msg_.r = 1.0;
        msg_.a = 1.0;
        sum_mode3++;
        break;
    case 4: //zero velocity (blue)
        msg_.b = 1.0;
        msg_.a = 1.0;
        sum_mode4++;
        break;
    case 5: //pure rotation
        msg_.r = 1.0;
        msg_.b = 1.0;
        msg_.a = 1.0;
        sum_mode5++;
        break;
    }
    cout << "sum-mode:" << sum_mode1 << " " << sum_mode2 << " " << sum_mode3 << " " << sum_mode4 << " " << sum_mode5 << endl;
    geometry_msgs::Point msg;
    Vector3d P = estimator.current_lidar.lidarData.lidar_T;
    msg.x = P[0];
    msg.y = P[1];
    // msg.z = P[2];
    msg.z = -0.5;
    ConstraintMarker.points.push_back(msg);
    ConstraintMarker.colors.push_back(msg_);
    pub_ConstraintMarker.publish(ConstraintMarker);

    sum_mode_Marker.header = header;
    sum_mode_Marker.header.frame_id = "world";
    sum_mode_Marker.ns = "sum_mode";
    sum_mode_Marker.id = 2;
    sum_mode_Marker.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
    sum_mode_Marker.action = visualization_msgs::Marker::ADD;
    sum_mode_Marker.pose.orientation.w = 1.0;
    sum_mode_Marker.lifetime = ros::Duration();
    sum_mode_Marker.scale.x = 0.5; //2
    sum_mode_Marker.scale.y = 0.5;
    sum_mode_Marker.scale.z = 0.5;
    sum_mode_Marker.color.a = 1.0;
    sum_mode_Marker.color.r = 1;
    sum_mode_Marker.color.g = 1;
    sum_mode_Marker.color.b = 1;
    sum_mode_Marker.pose.position.x = P[0];
    sum_mode_Marker.pose.position.y = P[1]; //P[1] + 15
    sum_mode_Marker.pose.position.z = P[2] + 5;
    // sum_mode_Marker.pose.position.x = 0;
    // sum_mode_Marker.pose.position.y = 0;
    // sum_mode_Marker.pose.position.z = 15;
    std::string output_text = "VIO status detected by LiDAR : \nignored (white) : " + std::to_string(sum_mode1) + "\n" + "good (green) : " + std::to_string(sum_mode2) + "\n" + "vio drift (red) : " + std::to_string(sum_mode3) + "\n" + "zero velocity (blue) : " + std::to_string(sum_mode4) + "\n" + "pure rotation (purple) : " + std::to_string(sum_mode5);
    sum_mode_Marker.text = output_text;
    pub_sum_mode.publish(sum_mode_Marker);
}

void pubCameraPose(const Estimator &estimator, const std_msgs::Header &header)
{
    int idx2 = WINDOW_SIZE - 1;

    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {
        int i = idx2;
        Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Quaterniond R = Quaterniond(estimator.Rs[i] * estimator.ric[0]);

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();

        pub_camera_pose.publish(odometry);

        cameraposevisual.reset();
        cameraposevisual.add_pose(P, R);
        cameraposevisual.publish_by(pub_camera_pose_visual, odometry.header);
    }
}

void pubLidarCloud(const CloudData cloud, const std_msgs::Header &header)
{
    sensor_msgs::PointCloud2 lidar_cloud;
    pcl::toROSMsg(*cloud.cloud_ptr, lidar_cloud);
    lidar_cloud.header = header;
    lidar_cloud.header.frame_id = "laser_link";
    pub_lidar_cloud.publish(lidar_cloud);
}

void pubLidarCloud(const Estimator &estimator, const std_msgs::Header &header)
{
    sensor_msgs::PointCloud2 lidar_cloud;
    pcl::toROSMsg(*estimator.current_lidar.lidarData.point_cloud.cloud_ptr, lidar_cloud);
    lidar_cloud.header = header;
    lidar_cloud.header.frame_id = "laser_link";
    pub_lidar_cloud.publish(lidar_cloud);
}

void pubLidarPose(const Estimator &estimator, const std_msgs::Header &header)
{

    Vector3d P = estimator.current_lidar.lidarData.lidar_T;
    Quaterniond R = Quaterniond(estimator.current_lidar.lidarData.lidar_R);

    nav_msgs::Odometry odometry;
    odometry.header = header;
    odometry.header.frame_id = "world";
    odometry.child_frame_id = "laser_link";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = R.x();
    odometry.pose.pose.orientation.y = R.y();
    odometry.pose.pose.orientation.z = R.z();
    odometry.pose.pose.orientation.w = R.w();

    pub_lidar_pose.publish(odometry);

    geometry_msgs::PoseStamped pose_stamped;
    pose_stamped.header = header;
    pose_stamped.header.frame_id = "world";
    pose_stamped.pose = odometry.pose.pose;

    pub_fake_trigger.publish(pose_stamped);
    odometry.pose.pose.position.x = 0.0;
    odometry.pose.pose.position.y = 0.0;
    odometry.pose.pose.position.z = P.z();
    pub_fake_local_odom.publish(odometry);

    lidar_path.header = header;
    lidar_path.header.frame_id = "world";
    lidar_path.poses.push_back(pose_stamped);
    pub_lidar_path.publish(lidar_path);
}

void pubPointCloud(const Estimator &estimator, const std_msgs::Header &header)
{
    sensor_msgs::PointCloud point_cloud, loop_point_cloud;
    point_cloud.header = header;
    loop_point_cloud.header = header;

    for (auto &it_per_id : estimator.f_manager.feature)
    {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        if (it_per_id.start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id.solve_flag != 1)
            continue;
        int imu_i = it_per_id.start_frame;
        Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
        Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

        geometry_msgs::Point32 p;
        p.x = w_pts_i(0);
        p.y = w_pts_i(1);
        p.z = w_pts_i(2);
        point_cloud.points.push_back(p);
    }
    pub_point_cloud.publish(point_cloud);

    // pub margined potin
    sensor_msgs::PointCloud margin_cloud;
    margin_cloud.header = header;

    for (auto &it_per_id : estimator.f_manager.feature)
    {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        //if (it_per_id->start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id->solve_flag != 1)
        //        continue;

        if (it_per_id.start_frame == 0 && it_per_id.feature_per_frame.size() <= 2 && it_per_id.solve_flag == 1)
        {
            int imu_i = it_per_id.start_frame;
            Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
            Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

            geometry_msgs::Point32 p;
            p.x = w_pts_i(0);
            p.y = w_pts_i(1);
            p.z = w_pts_i(2);
            margin_cloud.points.push_back(p);
        }
    }
    pub_margin_cloud.publish(margin_cloud);
}

void pubTF(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.solver_flag != Estimator::SolverFlag::NON_LINEAR)
        return;
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    tf::Quaternion q;
    // body frame
    Vector3d correct_t;
    Quaterniond correct_q;
    correct_t = estimator.Ps[WINDOW_SIZE];
    correct_q = estimator.Rs[WINDOW_SIZE];

    transform.setOrigin(tf::Vector3(correct_t(0),
                                    correct_t(1),
                                    correct_t(2)));
    q.setW(correct_q.w());
    q.setX(correct_q.x());
    q.setY(correct_q.y());
    q.setZ(correct_q.z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "world", "body"));

    // camera frame
    transform.setOrigin(tf::Vector3(estimator.tic[0].x(),
                                    estimator.tic[0].y(),
                                    estimator.tic[0].z()));
    q.setW(Quaterniond(estimator.ric[0]).w());
    q.setX(Quaterniond(estimator.ric[0]).x());
    q.setY(Quaterniond(estimator.ric[0]).y());
    q.setZ(Quaterniond(estimator.ric[0]).z());
    transform.setRotation(q.normalized());
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "body", "camera"));

    // lidar frame
    if (estimator.lidarCalibration.lidar_init_flag == false)
    {
        transform.setOrigin(tf::Vector3(estimator.TBL[0],
                                        estimator.TBL[1],
                                        estimator.TBL[2]));
        q.setW(Quaterniond(estimator.RLB).w());
        q.setX(Quaterniond(estimator.RLB).x());
        q.setY(Quaterniond(estimator.RLB).y());
        q.setZ(Quaterniond(estimator.RLB).z());
        transform.setRotation(q.inverse().normalized());
        br.sendTransform(tf::StampedTransform(transform, header.stamp, "body", "laser_link"));
    }
    if (estimator.lidarCalibration.lidar_init_flag == false)
    {
        transform.setOrigin(tf::Vector3(estimator.LIDAR_INIT_t.x(),
                                        estimator.LIDAR_INIT_t.y(),
                                        estimator.LIDAR_INIT_t.z()));
        q.setW(estimator.LIDAR_INIT_R.w());
        q.setX(estimator.LIDAR_INIT_R.x());
        q.setY(estimator.LIDAR_INIT_R.y());
        q.setZ(estimator.LIDAR_INIT_R.z());
        transform.setRotation(q.normalized());
        br.sendTransform(tf::StampedTransform(transform, header.stamp, "world", "world_lidar"));
    }
    nav_msgs::Odometry odometry;
    odometry.header = header;
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = estimator.tic[0].x();
    odometry.pose.pose.position.y = estimator.tic[0].y();
    odometry.pose.pose.position.z = estimator.tic[0].z();
    Quaterniond tmp_q{estimator.ric[0]};
    odometry.pose.pose.orientation.x = tmp_q.x();
    odometry.pose.pose.orientation.y = tmp_q.y();
    odometry.pose.pose.orientation.z = tmp_q.z();
    odometry.pose.pose.orientation.w = tmp_q.w();
    pub_extrinsic.publish(odometry);
}

void pubKeyframe(const Estimator &estimator)
{
    // pub camera pose, 2D-3D points of keyframe
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR && estimator.marginalization_flag == 0)
    {
        int i = WINDOW_SIZE - 2;
        //Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Vector3d P = estimator.Ps[i];
        Quaterniond R = Quaterniond(estimator.Rs[i]);

        nav_msgs::Odometry odometry;
        odometry.header = estimator.Headers[WINDOW_SIZE - 2];
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();
        //printf("time: %f t: %f %f %f r: %f %f %f %f\n", odometry.header.stamp.toSec(), P.x(), P.y(), P.z(), R.w(), R.x(), R.y(), R.z());

        pub_keyframe_pose.publish(odometry);

        sensor_msgs::PointCloud point_cloud;
        point_cloud.header = estimator.Headers[WINDOW_SIZE - 2];
        for (auto &it_per_id : estimator.f_manager.feature)
        {
            int frame_size = it_per_id.feature_per_frame.size();
            if (it_per_id.start_frame < WINDOW_SIZE - 2 && it_per_id.start_frame + frame_size - 1 >= WINDOW_SIZE - 2 && it_per_id.solve_flag == 1)
            {

                int imu_i = it_per_id.start_frame;
                Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
                Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];
                geometry_msgs::Point32 p;
                p.x = w_pts_i(0);
                p.y = w_pts_i(1);
                p.z = w_pts_i(2);
                point_cloud.points.push_back(p);

                int imu_j = WINDOW_SIZE - 2 - it_per_id.start_frame;
                sensor_msgs::ChannelFloat32 p_2d;
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.y());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.y());
                p_2d.values.push_back(it_per_id.feature_id);
                point_cloud.channels.push_back(p_2d);
            }
        }
        pub_keyframe_point.publish(point_cloud);
    }
}

void pubRelocalization(const Estimator &estimator)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(estimator.relo_frame_stamp);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = estimator.relo_relative_t.x();
    odometry.pose.pose.position.y = estimator.relo_relative_t.y();
    odometry.pose.pose.position.z = estimator.relo_relative_t.z();
    odometry.pose.pose.orientation.x = estimator.relo_relative_q.x();
    odometry.pose.pose.orientation.y = estimator.relo_relative_q.y();
    odometry.pose.pose.orientation.z = estimator.relo_relative_q.z();
    odometry.pose.pose.orientation.w = estimator.relo_relative_q.w();
    odometry.twist.twist.linear.x = estimator.relo_relative_yaw;
    odometry.twist.twist.linear.y = estimator.relo_frame_index;

    pub_relo_relative_pose.publish(odometry);
}