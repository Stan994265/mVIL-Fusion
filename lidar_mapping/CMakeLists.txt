cmake_minimum_required(VERSION 3.8)
project(vils_mapping)

set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -pthread")

find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  nav_msgs
  sensor_msgs
  roscpp
  rospy
  rosbag
  std_msgs
  image_transport
  cv_bridge
  tf
)

#find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)

find_package(OpenMP REQUIRED)
find_package(GTSAM REQUIRED QUIET)
find_package(Boost REQUIRED COMPONENTS timer)
include_directories(
  include
	${catkin_INCLUDE_DIRS} 
	${PCL_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${GTSAM_INCLUDE_DIR}
)

catkin_package(
  CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs
  DEPENDS EIGEN3 PCL 
  INCLUDE_DIRS include
)


add_executable(feature src/scanRegistration.cpp)
target_link_libraries(feature ${catkin_LIBRARIES} ${PCL_LIBRARIES})

add_executable(local src/localMapping.cpp)
target_link_libraries(local ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${CERES_LIBRARIES})

include_directories(${PROJECT_SOURCE_DIR}/fast_gicp)         

add_executable(global_ikdtree 
  src/globalMappingIkdTree.cpp
  include/ikd_Tree/ikd_Tree.cpp
  include/scancontext/Scancontext.cpp
)
target_compile_options(global_ikdtree 
  PRIVATE ${OpenMP_CXX_FLAGS}
)

target_link_libraries(global_ikdtree 
  ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenMP_CXX_FLAGS}
  # gtsam ${PROJECT_SOURCE_DIR}/include/fast_gicp/libfast_gicp.a
  gtsam ${PROJECT_SOURCE_DIR}/../vils_estimator/src/lidar_functions/fast_gicp/build/libfast_gicp.a
)

add_executable(global_octree 
  src/globalMappingOcTree.cpp
  include/scancontext/Scancontext.cpp
)
target_compile_options(global_octree 
  PRIVATE ${OpenMP_CXX_FLAGS}
)
target_link_libraries(global_octree 
  ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenMP_CXX_FLAGS}
  #gtsam ${PROJECT_SOURCE_DIR}/include/fast_gicp/libfast_gicp.a
  gtsam ${PROJECT_SOURCE_DIR}/../vils_estimator/src/lidar_functions/fast_gicp/build/libfast_gicp.a
)




