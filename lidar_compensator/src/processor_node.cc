#include <pcl/io/pcd_io.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <pcl_conversions/pcl_conversions.h>
#include <tf/transform_datatypes.h>

#include <geometry_msgs/Quaternion.h> 

#include "PointProcessor.h"
#include "TicToc.h"
using namespace std;

int main(int argc, char **argv)
{

  ros::init(argc, argv, "point_processor");

  ros::NodeHandle nh("~");

  int sensor_type;
  double rad_diff;
  bool infer_start_ori;
  nh.param("sensor_type", sensor_type, 16);
  nh.param("rad_diff", rad_diff, 0.2);
  nh.param("infer_start_ori", infer_start_ori, false);
  double scan_period;
  nh.param("scan_period", scan_period, 0.1);

  PointProcessor processor; // Default sensor_type is 16

  if (sensor_type == 32)
  {
    processor = PointProcessor(-30.67f, 10.67f, 32);
  }
  else if (sensor_type == 64)
  {
    processor = PointProcessor(-24.9f, 2, 64);
  }
  else if (sensor_type == 320)
  {
    processor = PointProcessor(-25, 15, 32, true);
  }

  PointProcessorConfig config;
  config.rad_diff = rad_diff;
  config.infer_start_ori_ = infer_start_ori;
  config.scan_period = scan_period;
  processor.SetupConfig(config);
  processor.SetupRos(nh);
  ros::Rate r(100);
  // while (ros::ok())
  // {
  //   ros::spinOnce();
  //   r.sleep();
  // }
  ros::AsyncSpinner spinner(1);
  spinner.start();
  ros::waitForShutdown();

  return 0;
}