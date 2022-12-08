cmake_minimum_required(VERSION 2.8.3)
project(vils_estimator)

# add_definitions(-msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2)
# set(CMAKE_C_FLAGS "-msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2")
# set(CMAKE_CXX_FLAGS "-msse -msse2 -msse3 -msse4 -msse4.1 -msse4.2")

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -pthread ")

find_package(catkin REQUIRED COMPONENTS
    roscpp
    std_msgs
    geometry_msgs
    nav_msgs
    tf
    cv_bridge
    )

find_package(OpenMP)
    if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()

find_package(OpenCV REQUIRED)

# message(WARNING "OpenCV_VERSION: ${OpenCV_VERSION}")

find_package(Ceres REQUIRED)

find_package(PCL REQUIRED)


include_directories(${catkin_INCLUDE_DIRS} ${CERES_INCLUDE_DIRS} ${PCL_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3)
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
)

catkin_package()

set(LIDAR_FUNCTIONS_DIR ${PROJECT_SOURCE_DIR}/src/lidar_functions)
set(FAST_GICP_DIR ${LIDAR_FUNCTIONS_DIR}/fast_gicp)
include_directories(${FAST_GICP_DIR}/include)


add_executable(vils_estimator
    src/estimator_node.cpp
    src/parameters.cpp
    src/estimator.cpp
    src/feature_manager.cpp
    src/factor/pose_local_parameterization.cpp
    src/factor/projection_factor.cpp
    src/factor/projection_td_factor.cpp
    src/factor/marginalization_factor.cpp
    src/utility/utility.cpp
    src/utility/visualization.cpp
    src/utility/CameraPoseVisualization.cpp
    src/initial/solve_5pts.cpp
    src/initial/initial_aligment.cpp
    src/initial/initial_sfm.cpp
    src/initial/initial_ex_rotation.cpp
    src/lidar_frontend.cpp
    src/lidar_backend.cpp
    )


target_link_libraries(vils_estimator ${catkin_LIBRARIES} ${OpenCV_LIBS} ${CERES_LIBRARIES} ${PCL_LIBRARIES} 
${FAST_GICP_DIR}/build/libfast_gicp.a) 

