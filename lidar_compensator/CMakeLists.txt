cmake_minimum_required(VERSION 2.8.3)
project(lidar_compensator)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(catkin REQUIRED COMPONENTS
  dynamic_reconfigure
  geometry_msgs
  message_generation
  message_runtime
  nav_msgs
  pcl_conversions
  pcl_ros
  rosbag
  roscpp
  rospy
  sensor_msgs
  std_msgs
  tf
  tf_conversions
)


find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIR})
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})


catkin_package()

include_directories(
    ${catkin_INCLUDE_DIRS}
    )

add_executable(lidar_compensator
    src/PointProcessor.cc
    src/processor_node.cc
    )

target_link_libraries(lidar_compensator  ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBRARIES})
