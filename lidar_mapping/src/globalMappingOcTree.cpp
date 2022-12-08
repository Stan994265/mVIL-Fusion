#include "global_mapping/util.h"

#include <gtsam/base/Vector.h>
#include <gtsam/geometry/Pose3.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/linear/NoiseModel.h>
#include <gtsam/nonlinear/ISAM2.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/slam/PriorFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/navigation/GPSFactor.h>

#include "fast_gicp/fast_gicp.hpp"
#include "fast_gicp/fast_gicp_st.hpp"
#include "fast_gicp/fast_vgicp.hpp"
#include <pcl/filters/approximate_voxel_grid.h>

#define LOG_FLAG
typedef pcl::PointXYZI PointT;

typedef message_filters::sync_policies::ApproximateTime<nav_msgs::Odometry, sensor_msgs::PointCloud2> slamSyncPolicy;

std::string save_directory;
SCManager scManager;
int perform_SC{0};
int save_map{0};
int onboard{0};
double floorHeight;

class LaserMapping
{
private:
    bool is_initial;
    int update_count;
    ros::Time stamp;

    Pose6D current_odom;
    Pose6D last_odom;
    Pose6D incremental_odom;
    Pose6D current_localization;
    Pose6D current_temp_localization;
    Pose6D last_localization;
    Pose6D incremental_localization;

    pcl::PointCloud<PointT>::Ptr laserCloudIn;
    pcl::PointCloud<PointT>::Ptr query;
    pcl::PointCloud<PointT>::Ptr reference;
    pcl::PointCloud<PointT>::Ptr map_data;
    pcl::PointCloud<PointT>::Ptr aftMappingCloud;

    pcl::octree::OctreePointCloudSearch<PointT>::Ptr map_octree;

    int key;
    gtsam::Pose3 odometry;                                 
    std::unique_ptr<gtsam::ISAM2> isam;                    
    gtsam::Values values;                               
    std::map<int, pcl::PointCloud<PointT>::Ptr> keyed_scans; 
    int last_closure_key;

    int floor_id;
    std::map<int, int> key_floor;
    std::vector<double> key_header;

    std::vector<pair<int, int> > sc_factor;
    std::vector<pair<int, int> > tr_factor;

    message_filters::Subscriber<nav_msgs::Odometry> *odom_sub_;       
    message_filters::Subscriber<sensor_msgs::PointCloud2> *laser_sub_; 
    message_filters::Synchronizer<slamSyncPolicy> *sync_;

    ros::Publisher pubLaserLocalizer;
    ros::Publisher pubLaserLocalizerCloud;
    ros::Publisher pubMapCloud;
    ros::Publisher pubGraphNode;
    ros::Publisher pubBetweenFactorNode;
    ros::Publisher pubSCLoopNode;
    ros::Publisher pubTRLoopNode;
    ros::Publisher pubFactorInfo;
    tf2_ros::TransformBroadcaster tfbr;

public:
    LaserMapping(ros::NodeHandle &n) : is_initial(false), update_count(0), key(0)
    {
        laserCloudIn.reset(new pcl::PointCloud<PointT>());
        query.reset(new pcl::PointCloud<PointT>());
        reference.reset(new pcl::PointCloud<PointT>());
        map_data.reset(new pcl::PointCloud<PointT>());
        aftMappingCloud.reset(new pcl::PointCloud<PointT>());

        map_octree.reset(new pcl::octree::OctreePointCloudSearch<PointT>(map_octree_resolution));
        map_octree->setInputCloud(map_data);

        ros::NodeHandle nh(n);

        odom_sub_ = new message_filters::Subscriber<nav_msgs::Odometry>(nh, "/local_odom", 5, ros::TransportHints().tcpNoDelay());       
        laser_sub_ = new message_filters::Subscriber<sensor_msgs::PointCloud2>(nh, "/local_map", 5, ros::TransportHints().tcpNoDelay());
        sync_ = new message_filters::Synchronizer<slamSyncPolicy>(slamSyncPolicy(10), *odom_sub_, *laser_sub_);
        sync_->registerCallback(boost::bind(&LaserMapping::combineCallback, this, _1, _2));

        pubLaserLocalizer = nh.advertise<nav_msgs::Odometry>("/laser_localizer", 5);
        pubLaserLocalizerCloud = nh.advertise<sensor_msgs::PointCloud2>("/laser_localizer_cloud", 5);
        pubMapCloud = nh.advertise<sensor_msgs::PointCloud2>("/map_data_cloud", 5);
        pubGraphNode = nh.advertise<visualization_msgs::Marker>("/graph_nodes", 5);
        pubBetweenFactorNode = nh.advertise<visualization_msgs::Marker>("/between_factor_nodes", 5);
        pubSCLoopNode = nh.advertise<visualization_msgs::Marker>("/sc_loop_nodes", 1);
        pubTRLoopNode = nh.advertise<visualization_msgs::Marker>("/tr_loop_nodes", 1);
        pubFactorInfo = nh.advertise<visualization_msgs::Marker>("/factor_info", 1);

        gtsam::ISAM2Params parameters;
        parameters.relinearizeSkip = 1;
        parameters.relinearizeThreshold = 0.01;
        isam.reset(new gtsam::ISAM2(parameters));

        last_closure_key = std::numeric_limits<int>::min();
    }
    ~LaserMapping()
    {
        std::string path = save_directory + "Backend.txt";
        SerialPath2File(path);

        if (save_map)
        {
            std::string map_name(save_directory + "Map.pcd");
            pcl::PCDWriter writer;
            map_data->width = 1;
            map_data->height = map_data->points.size();
            writer.write(map_name, *map_data);
        }
    }

    void combineCallback(const nav_msgs::Odometry::ConstPtr odom, const sensor_msgs::PointCloud2::ConstPtr laser)
    {

        if (!is_initial)
        {
            stamp = laser->header.stamp;
            current_odom = Pose6D(Vector3(odom->pose.pose.position.x, odom->pose.pose.position.y, odom->pose.pose.position.z),
                                  Quaternion(odom->pose.pose.orientation.w, odom->pose.pose.orientation.x, odom->pose.pose.orientation.y, odom->pose.pose.orientation.z));
            current_localization = current_odom;

            current_temp_localization = current_localization;
            pcl::fromROSMsg(*laser, *laserCloudIn);

            transformPointsCloud2Frame(laserCloudIn, query, current_localization);
            InsertPoints(query);

            gtsam::Vector3 translation(current_localization.trans().x(), current_localization.trans().y(), current_localization.trans().z());
            gtsam::Rot3 rotation(gtsam::Rot3::RzRyRx(current_localization.rot().toEuler().x(), current_localization.rot().toEuler().y(), current_localization.rot().toEuler().z()));
            gtsam::Pose3 pose(rotation, translation);

            gtsam::Vector6 noise;
            noise << 1e-9, 1e-9, 1e-9, 1e-4, 1e-4, 1e-4;
            gtsam::noiseModel::Diagonal::shared_ptr covariance = gtsam::noiseModel::Diagonal::Variances(noise);
            gtsam::NonlinearFactorGraph new_factor;
            gtsam::Values new_value;
            new_factor.add(gtsam::PriorFactor<gtsam::Pose3>(key, pose, covariance));
            new_value.insert(key, pose);

            isam->update(new_factor, new_value);
            // isam->update();

            values = isam->calculateEstimate();
            key++;

            odometry = gtsam::Pose3::identity();

#ifdef LOG_FLAG
            cout << "LaserMapping ----->  Initial done!" << endl; 
#endif

            is_initial = true;
            return;
        }

        stamp = laser->header.stamp;
        update_count++;

        last_odom = current_odom;
        current_odom = Pose6D(Vector3(odom->pose.pose.position.x, odom->pose.pose.position.y, odom->pose.pose.position.z), Quaternion(odom->pose.pose.orientation.w, odom->pose.pose.orientation.x, odom->pose.pose.orientation.y, odom->pose.pose.orientation.z));
        incremental_odom = last_odom.inv() * current_odom;
        current_localization = current_localization * incremental_odom;

        float floor{odom->pose.pose.position.z};
        if (floorHeight > 20.0)
        {
            floor_id = 0;
        }

        else if (floor <= 2.5 && floor > -5.0)
        {
            floor_id = 0;
        }

        else if (floor <= 5.0 && floor > 2.5)
        {
            floor_id = 1;
        }

        else if (floor <= 10.0 && floor > 5.0)
        {
            floor_id = 2;
        }
        else
        {
            floor_id = 3;
        }

        if (update_count == 1)
        {

#ifdef LOG_FLAG
            double start, finish;
            start = clock(); 
#endif

            gtsam::Point3 ops(odom->pose.pose.position.x, odom->pose.pose.position.y, odom->pose.pose.position.z);
            gtsam::Vector3 noise_p;
            if (floorHeight < 20.0)
            {
                noise_p << 1e9, 1e9, 250;
            }
            else
            {
                noise_p << 1e9, 1e9, 1e6;
            }
            gtsam::noiseModel::Diagonal::shared_ptr covariance_p = gtsam::noiseModel::Diagonal::Variances(noise_p);
            gtsam::NonlinearFactorGraph new_factor;
            new_factor.add(gtsam::GPSFactor(key, ops, covariance_p));

            update_count = 0;
            pcl::fromROSMsg(*laser, *laserCloudIn);

            if (onboard)
            {
                pcl::ApproximateVoxelGrid<PointT> voxelgrid;
                voxelgrid.setLeafSize(0.2, 0.2, 0.2);
                voxelgrid.setInputCloud(laserCloudIn);
                voxelgrid.filter(*laserCloudIn);
            }

            transformPointsCloud2Frame(laserCloudIn, query, current_localization); 

            ApproxNearestNeighbors(query, reference);
            pcl::PointCloud<PointT>::Ptr query_base(new pcl::PointCloud<PointT>());
            pcl::PointCloud<PointT>::Ptr reference_base(new pcl::PointCloud<PointT>());
            copyPointCloud(*laserCloudIn, *query_base);
            transformPointsCloud2Frame(reference, reference_base, current_localization.inv());

            //scan-to-map matching
            float fitnessScore;
            Pose6D poseUpdate = updateLocalization(query_base, reference_base, fitnessScore);
            // ROS_WARN_STREAM("!!!!!!!!!!!!!!!!!!!!!!! " << fitnessScore);
            last_localization = current_temp_localization;
            current_localization = current_localization * poseUpdate;
 
            current_temp_localization = current_localization;
            incremental_localization = last_localization.inv() * current_localization;

            // gtsam::NonlinearFactorGraph new_factor;
            gtsam::Values new_value;

            gtsam::Vector6 noise;
            noise << fitnessScore, fitnessScore, fitnessScore, fitnessScore, fitnessScore, fitnessScore;
            gtsam::noiseModel::Diagonal::shared_ptr covariance = gtsam::noiseModel::Diagonal::Variances(noise);
            gtsam::Pose3 new_odometry = ToGtsam(incremental_localization);
            new_factor.add(gtsam::BetweenFactor<gtsam::Pose3>(key - 1, key, new_odometry, covariance));

            gtsam::Pose3 last_pose = values.at<gtsam::Pose3>(key - 1);
            new_value.insert(key, last_pose.compose(new_odometry));

            isam->update(new_factor, new_value);
            // isam->update();

            values = isam->calculateEstimate();
            int current_pose_key = key++;

            if (check_for_loop_closure) // && isKeyFrame(new_odometry)
            {

#ifdef LOG_FLAG
                cout << "LaserMapping ----->  key = " << key << endl; 
#endif
                pcl::PointCloud<PointT>::Ptr laser_cloud_in_copy(new pcl::PointCloud<PointT>());
                copyPointCloud(*laserCloudIn, *laser_cloud_in_copy);
                keyed_scans.insert(std::pair<int, pcl::PointCloud<PointT>::Ptr>(current_pose_key, laser_cloud_in_copy));
                key_floor.insert(std::pair<int, int>(current_pose_key, floor_id));
                key_header.push_back(odom->header.stamp.toSec());

                if (perform_SC)
                {
                    scManager.makeAndSaveScancontextAndKeys(*laser_cloud_in_copy);
                }

                double temDiff = fabs(values.at<gtsam::Pose3>(current_pose_key).translation().z() - current_localization.z());
                if (temDiff > 1.0) //&& floorHeight > 10.0
                {
                    pcl::PointCloud<PointT>::Ptr regenerated_map(new pcl::PointCloud<PointT>());
                    GetMaximumLikelihoodPoints(regenerated_map);

                    map_data.reset(new pcl::PointCloud<PointT>());
                    map_octree.reset(new pcl::octree::OctreePointCloudSearch<PointT>(map_octree_resolution));
                    map_octree->setInputCloud(map_data);
                    current_localization = FromGtsam(values.at<gtsam::Pose3>(current_pose_key));
                    current_temp_localization = FromGtsam(values.at<gtsam::Pose3>(current_pose_key));

                    pcl::PointCloud<PointT>::Ptr unused(new pcl::PointCloud<PointT>);
                    InsertPoints(regenerated_map);
                }

                int closure_key;
                if (findLoopClosure(current_pose_key, closure_key))
                {
#ifdef LOG_FLAG
                    cout << "LaserMapping ----->  Loop closure detected: between " << current_pose_key << " and " << closure_key << endl;
#endif
                    pcl::PointCloud<PointT>::Ptr regenerated_map(new pcl::PointCloud<PointT>());
                    GetMaximumLikelihoodPoints(regenerated_map);

                    map_data.reset(new pcl::PointCloud<PointT>());
                    map_octree.reset(new pcl::octree::OctreePointCloudSearch<PointT>(map_octree_resolution));
                    map_octree->setInputCloud(map_data);
                    current_localization = FromGtsam(values.at<gtsam::Pose3>(current_pose_key));
                    current_temp_localization = FromGtsam(values.at<gtsam::Pose3>(current_pose_key));

                    pcl::PointCloud<PointT>::Ptr unused(new pcl::PointCloud<PointT>);
                    InsertPoints(regenerated_map);

                    tr_factor.emplace_back(make_pair(closure_key, current_pose_key));
                }
                else
                {
                    transformPointsCloud2Frame(laserCloudIn, aftMappingCloud, current_localization);
                    InsertPoints(aftMappingCloud);
                }
            }
            else
            {
                transformPointsCloud2Frame(laserCloudIn, aftMappingCloud, current_localization);
                InsertPoints(aftMappingCloud);
            }
            pubMapData();
            pubLowerFreq();

#ifdef LOG_FLAG
            finish = clock();                                                                                          
            cout << "LaserMapping ----->  Time consuming " << (finish - start) / CLOCKS_PER_SEC << " seconds" << endl; 
#endif
        }

        pubHigherFreq();

        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }

    void ScanContextThread()
    {
        ros::Rate rate(1);
        while (ros::ok())
        {
            rate.sleep();
            performSC();
        }
    }

    void performSC()
    {
        if (keyed_scans.size() < 5)
            return;
        int currentID;
        auto detectResult = scManager.detectLoopClosureID(currentID);
        int historyID = detectResult.first;
        int current_floor = key_floor[currentID];
        int history_floor = key_floor[historyID];
        if (historyID != -1 && historyID != 0 && current_floor == history_floor)
        {
            ROS_WARN_STREAM("TEST SC LOOP FROM " << historyID << " TO " << currentID);
            performSC_ICP(historyID, currentID);
        }
    }

    void performSC_ICP(int historyID, int currentID)
    {
        Pose6D pose1 = FromGtsam(values.at<gtsam::Pose3>(currentID));
        Pose6D pose2 = FromGtsam(values.at<gtsam::Pose3>(historyID));
        const pcl::PointCloud<PointT>::Ptr scan1 = keyed_scans[currentID];
        const pcl::PointCloud<PointT>::Ptr scan2 = keyed_scans[historyID];

        Pose6D distance_pose = pose2.inv() * pose1;
        Pose6D deltaTemp;
        float fitnessTemp = 0.0;
        // ROS_WARN_STREAM("key_scan size = " << keyed_scans.size());
        if (performICP(scan1, scan2, deltaTemp, fitnessTemp, distance_pose))
        {
            ROS_WARN_STREAM("performICP fitnessTemp = " << fitnessTemp);

            if (fitnessTemp < max_tolerable_fitness)
            {
                ROS_WARN_STREAM("SC LOOP FIND!!! FROM " << historyID << " TO " << currentID);
                gtsam::NonlinearFactorGraph new_factor;
                gtsam::Vector6 noise;
                noise << fitnessTemp, fitnessTemp, fitnessTemp, fitnessTemp, fitnessTemp, fitnessTemp;
                gtsam::noiseModel::Diagonal::shared_ptr covariance = gtsam::noiseModel::Diagonal::Variances(noise);

                new_factor.add(gtsam::BetweenFactor<gtsam::Pose3>(currentID, historyID, ToGtsam(deltaTemp), covariance));
                isam->update(new_factor, gtsam::Values());
                isam->update();
                values = isam->calculateEstimate();
                sc_factor.emplace_back(make_pair(historyID, currentID));
            }
        }
    }

    void GetMaximumLikelihoodPoints(pcl::PointCloud<PointT>::Ptr points)
    {
        for (const auto &keyed_pose : values)
        {
            const unsigned int key_index = keyed_pose.key;
            if (!keyed_scans.count(key_index))
                continue;

            Pose6D pose = FromGtsam(values.at<gtsam::Pose3>(key_index));
            Eigen::Matrix4d b2w;

            Eigen::Matrix<double, 3, 3> R;
            Eigen::Matrix<double, 3, 1> T;

            vector<double> temp(9, 0.0);
            pose.rot().toRotMatrix(temp);
            R(0, 0) = temp[0];
            R(0, 1) = temp[1];
            R(0, 2) = temp[2];
            R(1, 0) = temp[3];
            R(1, 1) = temp[4];
            R(1, 2) = temp[5];
            R(2, 0) = temp[6];
            R(2, 1) = temp[7];
            R(2, 2) = temp[8];

            T(0, 0) = pose.trans().x();
            T(1, 0) = pose.trans().y();
            T(2, 0) = pose.trans().z();

            b2w.block(0, 0, 3, 3) = R;
            b2w.block(0, 3, 3, 1) = T;

            pcl::PointCloud<PointT> scan_world;
            pcl::transformPointCloud(*keyed_scans[key_index], scan_world, b2w);
            *points += scan_world;
        }
    }
    bool findLoopClosure(int current_key, int &closure_key)
    {

        if (std::fabs(current_key - last_closure_key) < poses_before_reclosing)
            return false;

        const Pose6D pose1 = FromGtsam(values.at<gtsam::Pose3>(current_key));
        const pcl::PointCloud<PointT>::Ptr scan1 = keyed_scans[current_key];
        // pcl::PointCloud<PointT>::Ptr scan1_(new pcl::PointCloud<PointT>);
        // transformPointsCloud2Frame(scan1, scan1_, pose1);

        bool closed_loop = false;
        float fitnessMinResult = max_tolerable_fitness; 
        Pose6D deltaResult;                           
        int keyResut = 0;               

        for (const auto &keyed_pose : values)
        {
            const int other_key = keyed_pose.key;
            if (other_key == current_key)
                continue;

            if (std::fabs(current_key - other_key) < skip_recent_poses)
                continue;

            if (!keyed_scans.count(other_key))
                continue;

            Pose6D pose2 = FromGtsam(values.at<gtsam::Pose3>(other_key));
            Pose6D distance_pose = pose2.inv() * pose1;

            float diff_x{distance_pose.trans().x()};
            float diff_y{distance_pose.trans().y()};
            float diff = sqrt(diff_x * diff_x + diff_y * diff_y);

            int current_floor = key_floor.find(current_key)->second;
            int other_floor = key_floor.find(other_key)->second;

            if (diff < proximity_threshold && current_floor == other_floor) //distance_pose.trans().norm()
            {

                ROS_WARN_STREAM("LOOP IN FLOOR " << current_floor);
                const pcl::PointCloud<PointT>::Ptr scan2 = keyed_scans[other_key];
                // pcl::PointCloud<PointT>::Ptr scan2_(new pcl::PointCloud<PointT>);
                // transformPointsCloud2Frame(scan2, scan2_, pose2);
                Pose6D deltaTemp;
                float fitnessTemp = 0.0;

                if (performICP(scan1, scan2, deltaTemp, fitnessTemp, distance_pose))
                {
                    ROS_WARN_STREAM("performICP fitnessTemp = " << fitnessTemp);

                    if (fitnessTemp < fitnessMinResult)
                    {
                        fitnessMinResult = fitnessTemp;
                        keyResut = other_key;
                        deltaResult = deltaTemp;
                        closed_loop = true;
                    }
                }
            }
        }
        if (closed_loop)
        {
#ifdef LOG_FLAG
            cout << "LaserMapping ----->  Closure detected with fitnessScore = " << fitnessMinResult << endl;
#endif
            ROS_WARN_STREAM("Add between factor from " << current_key << " to " << keyResut);
            gtsam::NonlinearFactorGraph new_factor;
            gtsam::Vector6 noise;
            noise << fitnessMinResult, fitnessMinResult, fitnessMinResult, fitnessMinResult, fitnessMinResult, fitnessMinResult;
            gtsam::noiseModel::Diagonal::shared_ptr covariance = gtsam::noiseModel::Diagonal::Variances(noise);

            new_factor.add(gtsam::BetweenFactor<gtsam::Pose3>(current_key, keyResut, ToGtsam(deltaResult), covariance));
            isam->update(new_factor, gtsam::Values());
            isam->update();
            // isam->update();
            // values = isam->calculateEstimate();

            last_closure_key = current_key;
            closure_key = keyResut;
        }
        isam->update();
        values = isam->calculateEstimate();
        return closed_loop;
    }

    bool performICP(pcl::PointCloud<PointT>::Ptr query_base, pcl::PointCloud<PointT>::Ptr reference_base, Pose6D &delta, float &fitnessScore, Pose6D guess)
    {
        pcl::PointCloud<PointT>::Ptr aligned(new pcl::PointCloud<PointT>);
        pcl::PointCloud<PointT>::Ptr src(new pcl::PointCloud<PointT>);
        pcl::copyPointCloud(*query_base, *src);
        pcl::PointCloud<PointT>::Ptr tgt(new pcl::PointCloud<PointT>);
        pcl::copyPointCloud(*reference_base, *tgt);
        pcl::ApproximateVoxelGrid<PointT> voxelgrid;
        voxelgrid.setLeafSize(0.3, 0.3, 0.3);
        voxelgrid.setInputCloud(src);
        voxelgrid.filter(*src);
        voxelgrid.setLeafSize(0.3, 0.3, 0.3);
        voxelgrid.setInputCloud(tgt);
        voxelgrid.filter(*tgt);
        // boost::shared_ptr<fast_gicp::FastGICP<PointT, PointT> > gicp(new fast_gicp::FastGICP<PointT, PointT>());
        // gicp->setNumThreads(4);
        // gicp->setTransformationEpsilon(0.0005);
        // gicp->setMaxCorrespondenceDistance(0.8);
        // gicp->setMaximumIterations(100);
        boost::shared_ptr<fast_gicp::FastVGICP<PointT, PointT> > gicp(new fast_gicp::FastVGICP<PointT, PointT>());
        gicp->setResolution(0.5);
        gicp->setInputSource(src);
        gicp->setInputTarget(tgt);
        gicp->setNumThreads(4);
        Eigen::Matrix<float, 3, 3> R_;
        Eigen::Matrix<float, 3, 1> T_;
        vector<double> temp(9, 0.0);
        guess.rot().toRotMatrix(temp);
        R_(0, 0) = temp[0];
        R_(0, 1) = temp[1];
        R_(0, 2) = temp[2];
        R_(1, 0) = temp[3];
        R_(1, 1) = temp[4];
        R_(1, 2) = temp[5];
        R_(2, 0) = temp[6];
        R_(2, 1) = temp[7];
        R_(2, 2) = temp[8];
        T_(0, 0) = guess.trans().x();
        T_(1, 0) = guess.trans().y();
        T_(2, 0) = guess.trans().z();
        Eigen::Matrix4f init_guss(Eigen::Matrix4f::Identity());
        init_guss.block<3, 3>(0, 0) = R_;
        init_guss.block<3, 1>(0, 3) = T_;
        gicp->align(*aligned, init_guss);

        if (!gicp->hasConverged())
            return false;

        const Eigen::Matrix4f T = gicp->getFinalTransformation();
        fitnessScore = gicp->getFitnessScore();

        float tem_T = fabs(T_(0, 0) - T(0, 3)) + fabs(T_(1, 0) - T(1, 3)) + fabs(T_(2, 0) - T(2, 3));
        ROS_WARN_STREAM("DIFF T=" << tem_T);

        Pose6D poseUpdate(Vector3(T(0, 3), T(1, 3), T(2, 3)), Quaternion(T(0, 0), T(0, 1), T(0, 2), T(1, 0), T(1, 1), T(1, 2), T(2, 0), T(2, 1), T(2, 2)));
        delta = poseUpdate.inv();

        return true;
    }

    bool isKeyFrame(gtsam::Pose3 &new_odometry)
    {
        odometry = odometry.compose(new_odometry);
        if (odometry.translation().norm() > translation_threshold)
        {
            odometry = gtsam::Pose3::identity();
            return true;
        }
        return false;
    }

    gtsam::Pose3 ToGtsam(const Pose6D &pose)
    {
        gtsam::Vector3 translation(pose.trans().x(), pose.trans().y(), pose.trans().z());
        gtsam::Rot3 rotation(gtsam::Rot3::RzRyRx(pose.rot().toEuler().x(), pose.rot().toEuler().y(), pose.rot().toEuler().z()));
        return gtsam::Pose3(rotation, translation);
    }

    Pose6D FromGtsam(const gtsam::Pose3 &pose)
    {
        Vector3 v(pose.translation().x(), pose.translation().y(), pose.translation().z());
        Quaternion q(pose.rotation().roll(), pose.rotation().pitch(), pose.rotation().yaw());
        return Pose6D(v, q);
    }

    Pose6D updateLocalization(pcl::PointCloud<PointT>::Ptr query_base, pcl::PointCloud<PointT>::Ptr reference_base, float &fitnessScore)
    {

        pcl::PointCloud<PointT>::Ptr aligned(new pcl::PointCloud<PointT>);
        pcl::PointCloud<PointT>::Ptr src(new pcl::PointCloud<PointT>);
        pcl::copyPointCloud(*query_base, *src);
        pcl::PointCloud<PointT>::Ptr tgt(new pcl::PointCloud<PointT>);
        pcl::copyPointCloud(*reference_base, *tgt);

        pcl::ApproximateVoxelGrid<PointT> voxelgrid;
        voxelgrid.setLeafSize(0.3, 0.3, 0.3);
        voxelgrid.setInputCloud(src);
        voxelgrid.filter(*src);
        voxelgrid.setLeafSize(0.3, 0.3, 0.3);
        voxelgrid.setInputCloud(tgt);
        voxelgrid.filter(*tgt);

        // boost::shared_ptr<fast_gicp::FastGICP<PointT, PointT> > gicp(new fast_gicp::FastGICP<PointT, PointT>());
        // gicp->setNumThreads(4);
        // gicp->setTransformationEpsilon(0.0001);
        // gicp->setMaxCorrespondenceDistance(0.5);
        // gicp->setMaximumIterations(100);
        boost::shared_ptr<fast_gicp::FastVGICP<PointT, PointT> > gicp(new fast_gicp::FastVGICP<PointT, PointT>());
        gicp->setResolution(0.5);
        gicp->setInputSource(src);
        gicp->setInputTarget(tgt);
        gicp->setNumThreads(4);
        gicp->align(*aligned);

        const Eigen::Matrix4f T = gicp->getFinalTransformation();
        fitnessScore = gicp->getFitnessScore();

        Pose6D poseUpdate(Vector3(T(0, 3), T(1, 3), T(2, 3)),
                          Quaternion(T(0, 0), T(0, 1), T(0, 2),
                                     T(1, 0), T(1, 1), T(1, 2),
                                     T(2, 0), T(2, 1), T(2, 2)));

        return poseUpdate;
    }

    void transformPointsCloud2Frame(pcl::PointCloud<PointT>::Ptr src, pcl::PointCloud<PointT>::Ptr tgt, Pose6D t)
    {

        Eigen::Matrix<double, 3, 3> R;
        Eigen::Matrix<double, 3, 1> T;

        vector<double> temp(9, 0.0);
        t.rot().toRotMatrix(temp);
        R(0, 0) = temp[0];
        R(0, 1) = temp[1];
        R(0, 2) = temp[2];
        R(1, 0) = temp[3];
        R(1, 1) = temp[4];
        R(1, 2) = temp[5];
        R(2, 0) = temp[6];
        R(2, 1) = temp[7];
        R(2, 2) = temp[8];

        T(0, 0) = t.trans().x();
        T(1, 0) = t.trans().y();
        T(2, 0) = t.trans().z();

        Eigen::Matrix4d tf(Eigen::Matrix4d::Identity());
        tf.block(0, 0, 3, 3) = R;
        tf.block(0, 3, 3, 1) = T;

        pcl::transformPointCloud(*src, *tgt, tf);
    }
    void InsertPoints(pcl::PointCloud<PointT>::Ptr cloud)
    {
        for (int i = 0; i < cloud->points.size(); ++i)
        {
            const PointT p = cloud->points[i];
            if (!map_octree->isVoxelOccupiedAtPoint(p))
            {
                map_octree->addPointToCloud(p, map_data);
            }
        }
    }

    void ApproxNearestNeighbors(pcl::PointCloud<PointT>::Ptr que, pcl::PointCloud<PointT>::Ptr ref)
    {
        ref->clear();
        pcl::octree::OctreePointCloudSearch<PointT>::Ptr ref_octree;
        ref_octree.reset(new pcl::octree::OctreePointCloudSearch<PointT>(ref_octree_resolution));
        ref_octree->setInputCloud(ref);

        pcl::PointCloud<PointT>::Ptr copy(new pcl::PointCloud<PointT>());
        copyPointCloud(*que, *copy);
        pcl::ApproximateVoxelGrid<PointT> voxelgrid;
        voxelgrid.setLeafSize(0.3, 0.3, 0.3);
        voxelgrid.setInputCloud(copy);
        voxelgrid.filter(*copy);

        // pcl::octree::OctreePointCloudSearch<PointT>::Ptr ref_octree_tem;
        // pcl::PointCloud<PointT>::Ptr ref_tem(new pcl::PointCloud<PointT>());
        // ref_octree_tem.reset(new pcl::octree::OctreePointCloudSearch<PointT>(ref_octree_resolution));
        // ref_octree_tem->setInputCloud(ref_tem);
        // for (int i = 0; i < copy->points.size(); ++i)
        // {
        //     vector<int> result_index_tem;
        //     if (map_octree->voxelSearch(copy->points[i], result_index_tem))
        //     {
        //         for (int j = 0; j < result_index_tem.size(); j++)
        //         {
        //             if (!ref_octree_tem->isVoxelOccupiedAtPoint(map_data->points[result_index_tem[j]]))
        //             {
        //                 ref_octree_tem->addPointToCloud(map_data->points[result_index_tem[j]], ref_tem);
        //             }
        //         }
        //     }
        // }

        // double start, finish;
        // start = clock();
        for (int i = 0; i < copy->points.size(); ++i)
        {
            vector<int> result_index;
            vector<float> unused;
            map_octree->radiusSearch(copy->points[i], map_octree_radius, result_index, unused);
            for (int j = 0; j < result_index.size(); j++)
            {

                if (!ref_octree->isVoxelOccupiedAtPoint(map_data->points[result_index[j]]))
                {
                    ref_octree->addPointToCloud(map_data->points[result_index[j]], ref);
                }
            }
        }
        // finish = clock();
        // ROS_WARN_STREAM("radiusSearch cost " << (finish - start) / CLOCKS_PER_SEC << " seconds");
    }

    void pubMapData()
    {

        // pcl::PointCloud<PointT>::Ptr map_filtered(new pcl::PointCloud<PointT>());
        // pcl::copyPointCloud(*map_data, *map_filtered);
        // pcl::VoxelGrid<PointT> grid;
        // grid.setLeafSize(0.2, 0.2, 0.2);
        // grid.setInputCloud(map_filtered);
        // grid.filter(*map_filtered);

        sensor_msgs::PointCloud2 cloud;
        pcl::toROSMsg(*map_data, cloud);
        cloud.header.stamp = stamp;
        cloud.header.frame_id = fix_frame_id;
        pubMapCloud.publish(cloud);
    }
    void pubLowerFreq()
    {
        if (pubGraphNode.getNumSubscribers() > 0)
        {
            visualization_msgs::Marker m;
            m.header.frame_id = fix_frame_id;
            m.ns = fix_frame_id;
            m.id = 1;
            m.action = visualization_msgs::Marker::ADD;
            m.type = visualization_msgs::Marker::SPHERE_LIST;
            m.lifetime = ros::Duration();
            m.scale.x = 1.0;
            m.scale.y = 1.0;
            m.scale.z = 1.0;

            visualization_msgs::Marker n;
            n.header.frame_id = fix_frame_id;
            n.ns = fix_frame_id;
            n.id = 2;
            n.action = visualization_msgs::Marker::ADD;
            n.type = visualization_msgs::Marker::LINE_STRIP;
            n.lifetime = ros::Duration();
            n.scale.x = 0.3;
            n.color.r = 1.0;
            n.color.g = 1.0;
            n.color.b = 1.0;
            n.color.a = 1.0;

            for (const auto &keyed_pose : values)
            {
                geometry_msgs::Point msg;
                std_msgs::ColorRGBA msg_;
                msg.x = values.at<gtsam::Pose3>(keyed_pose.key).translation().x();
                msg.y = values.at<gtsam::Pose3>(keyed_pose.key).translation().y();
                msg.z = values.at<gtsam::Pose3>(keyed_pose.key).translation().z();
                int floorID = key_floor.find(keyed_pose.key)->second;
                switch (floorID)
                {
                case 1:
                    msg_.g = 1.0;
                    msg_.a = 1.0;
                    break;
                case 0:
                    msg_.r = 1.0;
                    msg_.a = 1.0;
                    break;
                case 2:
                    msg_.b = 1.0;
                    msg_.a = 1.0;
                    break;
                case 3:
                    msg_.r = 1.0;
                    msg_.g = 1.0;
                    msg_.b = 1.0;
                    msg_.a = 1.0;
                    break;
                }
                m.points.push_back(msg);
                m.colors.push_back(msg_);

                n.points.push_back(msg);
            }
            pubGraphNode.publish(m);
            pubBetweenFactorNode.publish(n);
        }

        if (sc_factor.size() != 0)
        {
            visualization_msgs::Marker sc_list;
            sc_list.header.frame_id = fix_frame_id;
            sc_list.ns = fix_frame_id;
            sc_list.id = 2;
            sc_list.action = visualization_msgs::Marker::ADD;
            sc_list.type = visualization_msgs::Marker::LINE_LIST;
            sc_list.lifetime = ros::Duration();
            sc_list.scale.x = 0.3;
            sc_list.color.r = 1.0;
            sc_list.color.g = 1.0;
            sc_list.color.a = 1.0;
            ROS_WARN_STREAM("SC LOOP SIZE: " << sc_factor.size());
            for (auto iter = sc_factor.begin(); iter != sc_factor.end(); iter++)
            {
                int key1 = iter->first;
                int key2 = iter->second;
                geometry_msgs::Point p1, p2;
                p1.x = values.at<gtsam::Pose3>(key1).translation().x();
                p1.y = values.at<gtsam::Pose3>(key1).translation().y();
                p1.z = values.at<gtsam::Pose3>(key1).translation().z();
                p2.x = values.at<gtsam::Pose3>(key2).translation().x();
                p2.y = values.at<gtsam::Pose3>(key2).translation().y();
                p2.z = values.at<gtsam::Pose3>(key2).translation().z();
                sc_list.points.push_back(p1);
                sc_list.points.push_back(p2);
            }
            pubSCLoopNode.publish(sc_list);
        }

        if (tr_factor.size() != 0)
        {
            visualization_msgs::Marker tr_list;
            tr_list.header.frame_id = fix_frame_id;
            tr_list.ns = fix_frame_id;
            tr_list.id = 3;
            tr_list.action = visualization_msgs::Marker::ADD;
            tr_list.type = visualization_msgs::Marker::LINE_LIST;
            tr_list.lifetime = ros::Duration();
            tr_list.scale.x = 0.3;
            tr_list.color.g = 1.0;
            tr_list.color.a = 1.0;
            ROS_WARN_STREAM("TR LOOP SIZE: " << tr_factor.size());
            for (auto iter = tr_factor.begin(); iter != tr_factor.end(); iter++)
            {
                int key1 = iter->first;
                int key2 = iter->second;
                geometry_msgs::Point p1, p2;
                p1.x = values.at<gtsam::Pose3>(key1).translation().x();
                p1.y = values.at<gtsam::Pose3>(key1).translation().y();
                p1.z = values.at<gtsam::Pose3>(key1).translation().z();
                p2.x = values.at<gtsam::Pose3>(key2).translation().x();
                p2.y = values.at<gtsam::Pose3>(key2).translation().y();
                p2.z = values.at<gtsam::Pose3>(key2).translation().z();
                tr_list.points.push_back(p1);
                tr_list.points.push_back(p2);
            }
            pubTRLoopNode.publish(tr_list);
        }
        visualization_msgs::Marker factor_info;
        factor_info.header.frame_id = "world";
        factor_info.ns = "sum_mode";
        factor_info.id = 4;
        factor_info.type = visualization_msgs::Marker::TEXT_VIEW_FACING;
        factor_info.action = visualization_msgs::Marker::ADD;
        factor_info.pose.orientation.w = 1.0;
        factor_info.lifetime = ros::Duration();
        factor_info.scale.x = 2;
        factor_info.scale.y = 2;
        factor_info.scale.z = 2;
        factor_info.color.a = 1.0;
        factor_info.color.r = 1;
        factor_info.color.g = 1;
        factor_info.color.b = 1;
        int size = key_header.size();
        factor_info.pose.position.x = values.at<gtsam::Pose3>(size - 1).translation().x();
        factor_info.pose.position.y = values.at<gtsam::Pose3>(size - 1).translation().y() - 15;
        factor_info.pose.position.z = values.at<gtsam::Pose3>(size - 1).translation().z() + 15;
        std::string output_text = "Factor graph information statistics : \ninitial factor : 1\nprior factor : " + std::to_string(size) + "\n" + "odom factor (white line) : " + std::to_string(size - 1) + "\n" + "TRLoop factor (green line) : " + std::to_string(tr_factor.size()) + "\n" + "SCLoop factor (yellow line) : " + std::to_string(sc_factor.size());
        factor_info.text = output_text;
        pubFactorInfo.publish(factor_info);
    }
    void pubHigherFreq()
    {
        sensor_msgs::PointCloud2 cloud;
        pcl::toROSMsg(*laserCloudIn, cloud);
        cloud.header.stamp = stamp;
        cloud.header.frame_id = base_frame_id;
        pubLaserLocalizerCloud.publish(cloud);

        nav_msgs::Odometry laserLocalizer;
        laserLocalizer.header.frame_id = fix_frame_id;
        laserLocalizer.child_frame_id = base_frame_id;
        laserLocalizer.header.stamp = stamp;
        laserLocalizer.pose.pose.orientation.x = current_localization.rot().x();
        laserLocalizer.pose.pose.orientation.y = current_localization.rot().y();
        laserLocalizer.pose.pose.orientation.z = current_localization.rot().z();
        laserLocalizer.pose.pose.orientation.w = current_localization.rot().u();
        laserLocalizer.pose.pose.position.x = current_localization.x();
        laserLocalizer.pose.pose.position.y = current_localization.y();
        laserLocalizer.pose.pose.position.z = current_localization.z();
        pubLaserLocalizer.publish(laserLocalizer);

        geometry_msgs::TransformStamped tf;
        geometry_msgs::Vector3 vec;
        vec.x = current_localization.x();
        vec.y = current_localization.y();
        vec.z = current_localization.z();
        geometry_msgs::Quaternion quat;
        quat.w = current_localization.rot().u();
        quat.x = current_localization.rot().x();
        quat.y = current_localization.rot().y();
        quat.z = current_localization.rot().z();
        geometry_msgs::Transform trans;
        trans.translation = vec;
        trans.rotation = quat;
        tf.transform = trans;
        tf.header.stamp = stamp;
        tf.header.frame_id = fix_frame_id;
        tf.child_frame_id = base_frame_id;
        tfbr.sendTransform(tf);
    }

    void SerialPath2File(string posePath)
    {
        ofstream foutC(posePath, ios::trunc);
        foutC.setf(ios::fixed, ios::floatfield);
        for (int i = 0; i < key_header.size(); i++)
        {
            foutC.precision(9);
            foutC << key_header[i] << " ";
            foutC.precision(5);
            Pose6D pose = FromGtsam(values.at<gtsam::Pose3>(i));
            foutC << pose.x() << " "
                  << pose.y() << " "
                  << pose.z() << " "
                  << pose.rot().x() << " "
                  << pose.rot().y() << " "
                  << pose.rot().z() << " "
                  << pose.rot().u() << endl;
        }
        foutC.close();
    }
};

int main(int argc, char **argv)
{

    ros::init(argc, argv, "GlobalMapping");
    ros::NodeHandle n;
    n.param<std::string>("save_directory", save_directory, "/");
    double sc_dist_thres, sc_max_radius;
    n.param<double>("sc_dist_thres", sc_dist_thres, 0.2);
    n.param<double>("sc_max_radius", sc_max_radius, 40.0);
    n.param<double>("floorHeight", floorHeight, 2.5);
    n.param<int>("performSC", perform_SC, 0);
    n.param<int>("save_map", save_map, 0);
    n.param<int>("onboard", onboard, 0);
    scManager.setSCdistThres(sc_dist_thres);
    scManager.setMaximumRadius(sc_max_radius);
    if (floorHeight > 10)
    {
        max_tolerable_fitness = 2.0;
        proximity_threshold = 8.0;
    }
    LaserMapping lm(n);
    std::thread SC_Thread(&LaserMapping::ScanContextThread, &lm);

    ros::AsyncSpinner spinner(2);
    spinner.start();
    ros::waitForShutdown();
    lm.~LaserMapping();
    return 0;
}
