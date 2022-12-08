#ifndef _UTIL_H_
#define _UTIL_H_

#include <ros/ros.h>

#include <sensor_msgs/PointCloud2.h>
#include <sensor_msgs/Imu.h>
#include <nav_msgs/Odometry.h>
#include <geometry_msgs/Quaternion.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/Point32.h>
#include <geometry_msgs/Pose.h>
#include <geometry_msgs/Transform.h>
#include <visualization_msgs/Marker.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl_ros/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/range_image/range_image.h>
#include <pcl/filters/filter.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/radius_outlier_removal.h>

#include <pcl/search/impl/search.hpp>
#ifndef PCL_NO_PRECOMPILE
#include <pcl/impl/instantiate.hpp>
#include <pcl/point_types.h>
PCL_INSTANTIATE(Search, PCL_POINT_TYPES)
#endif // PCL_NO_PRECOMPILE

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/octree/octree_search.h>
#include <pcl/common/common.h>

#include <pcl/features/normal_3d.h>
#include <pcl/registration/icp.h>
#include <pcl/registration/transformation_estimation_point_to_plane.h>

#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>

#include "global_mapping/Vector3.h"
#include "global_mapping/Quaternion.h"
#include "global_mapping/Pose6D.h"
#include "scancontext/Scancontext.h"

#include <vector>
#include <cmath>
#include <algorithm>
#include <queue>
#include <deque>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cfloat>
#include <iterator>
#include <sstream>
#include <string>
#include <limits>
#include <iomanip>
#include <array>
#include <thread>
#include <mutex>
#include "time.h"

using namespace std;
using namespace alg;

extern const string fix_frame_id = "/world";
extern const string odom_frame_id = "/global_odom";
extern const string base_frame_id = "/global_odom_high_freq";

extern const float map_octree_resolution = 0.05; 
extern const float map_octree_radius = 0.5;      
extern const float ref_octree_resolution = 0.3;  

extern const bool check_for_loop_closure = true;
extern const float translation_threshold = 1.0; 
extern const int poses_before_reclosing = 10;   
float max_tolerable_fitness = 1.0;
float proximity_threshold = 5;
extern const int skip_recent_poses = 10;

#endif
