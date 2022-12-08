cmake_minimum_required(VERSION 2.8.3)
project(feature_tracker_)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    sensor_msgs
    cv_bridge
    camera_model
    tf
    )
find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)

catkin_package()

# include_directories(
#     ${catkin_INCLUDE_DIRS}
#     ${PCL_INCLUDE_DIRS}
#     )
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3)
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${PCL_INCLUDE_DIRS}
)

add_executable(feature_tracker_
    src/feature_tracker_node.cpp
    src/parameters.cpp
    src/feature_tracker.cpp
    )

target_link_libraries(feature_tracker_ ${catkin_LIBRARIES} ${OpenCV_LIBS} ${PCL_LIBRARIES})
