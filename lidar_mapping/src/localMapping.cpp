// This is an advanced implementation of the algorithm described in the following paper:
//   J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time.
//     Robotics: Science and Systems Conference (RSS). Berkeley, CA, July 2014.

// Modifier: Tong Qin               qintonguav@gmail.com
// 	         Shaozu Cao 		    saozu.cao@connect.ust.hk

// Copyright 2013, Ji Zhang, Carnegie Mellon University
// Further contributions copyright (c) 2016, Southwest Research Institute
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <math.h>
#include <vector>
#include <aloam_velodyne/common.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include <eigen3/Eigen/Dense>
#include <ceres/ceres.h>
#include <mutex>
#include <queue>
#include <thread>
#include <iostream>
#include <string>

#include "lidarFactor.hpp"
#include "aloam_velodyne/common.h"
#include "aloam_velodyne/tic_toc.h"

#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <omp.h>

#define FOR_GLOBAL

typedef pcl::PointXYZI PointT;
int frameCount = 0;
double timeLaserCloudCornerLast = 0;
double timeLaserCloudSurfLast = 0;
double timeLaserCloudFullRes = 0;
double timeLaserOdometry = 0;

const int laserCloudWidth = 11;	 //21 11 7
const int laserCloudHeight = 11; //21 11 7
const int laserCloudDepth = 7;	 //11 7 5
const double CubeLength = 10.0;	 // 50.0  20.0  10.0
const double CubeHeight = 5.0;	 // 50.0  2.0

int laserCloudCenWidth = (laserCloudWidth - 1) / 2;
int laserCloudCenHeight = (laserCloudHeight - 1) / 2; //10
int laserCloudCenDepth = (laserCloudDepth - 1) / 2;	  //5
const double HalfCubeLength = 1 / 2 * CubeLength;
const double HalfCubeHeight = 1 / 2 * CubeHeight;

const int laserCloudNum = laserCloudWidth * laserCloudHeight * laserCloudDepth; //4851

int laserCloudValidInd[125];
// int laserCloudSurroundInd[125];

// input: from odom
pcl::PointCloud<PointT>::Ptr laserCloudCornerLast(new pcl::PointCloud<PointT>());
pcl::PointCloud<PointT>::Ptr laserCloudSurfLast(new pcl::PointCloud<PointT>());

// ouput: all visualble cube points
// pcl::PointCloud<PointT>::Ptr laserCloudSurround(new pcl::PointCloud<PointT>());

// surround points in map to build tree
pcl::PointCloud<PointT>::Ptr laserCloudCornerFromMap(new pcl::PointCloud<PointT>());
pcl::PointCloud<PointT>::Ptr laserCloudSurfFromMap(new pcl::PointCloud<PointT>());

//input & output: points in one frame. local --> global
pcl::PointCloud<PointT>::Ptr laserCloudFullRes(new pcl::PointCloud<PointT>());

// points in every cube
pcl::PointCloud<PointT>::Ptr laserCloudCornerArray[laserCloudNum];
pcl::PointCloud<PointT>::Ptr laserCloudSurfArray[laserCloudNum];

//kd-tree
pcl::KdTreeFLANN<PointT>::Ptr kdtreeCornerFromMap(new pcl::KdTreeFLANN<PointT>());
pcl::KdTreeFLANN<PointT>::Ptr kdtreeSurfFromMap(new pcl::KdTreeFLANN<PointT>());

double parameters[7] = {0, 0, 0, 1, 0, 0, 0};
Eigen::Map<Eigen::Quaterniond> q_w_curr(parameters);
Eigen::Map<Eigen::Vector3d> t_w_curr(parameters + 4);

// wmap_T_odom * odom_T_curr = wmap_T_curr;
// transformation between odom's world and map's world frame
Eigen::Quaterniond q_wmap_wodom(1, 0, 0, 0);
Eigen::Vector3d t_wmap_wodom(0, 0, 0);

Eigen::Quaterniond q_wodom_curr(1, 0, 0, 0);
Eigen::Vector3d t_wodom_curr(0, 0, 0);

bool first_odom{true};
Eigen::Quaterniond q_w_last(1, 0, 0, 0);
Eigen::Vector3d t_w_last(0, 0, 0);
int last_frame_count;

std::queue<sensor_msgs::PointCloud2ConstPtr> cornerLastBuf;
std::queue<sensor_msgs::PointCloud2ConstPtr> surfLastBuf;
std::queue<sensor_msgs::PointCloud2ConstPtr> fullResBuf;
std::queue<nav_msgs::Odometry::ConstPtr> odometryBuf;
std::mutex mBuf;

pcl::ApproximateVoxelGrid<PointT> downSizeFilterCorner; //ApproximateVoxelGrid
pcl::ApproximateVoxelGrid<PointT> downSizeFilterSurf;
pcl::StatisticalOutlierRemoval<PointT> OutRemFilter;

std::vector<int> pointSearchInd;
std::vector<float> pointSearchSqDis;

PointT pointOri, pointSel;

ros::Publisher pubLaserCloudSurround, pubLaserCloudMap, pubLaserCloudFullRes, pubOdomAftMapped, pubOdomAftMappedHighFrec, pubLaserAfterMappedPath;
ros::Publisher pubLaserCloudFullResLocal;

ros::Publisher LocalMap, LocalOdom;

nav_msgs::Path laserAfterMappedPath;

std::string result_path;

// set initial guess
void transformAssociateToMap()
{
	q_w_curr = q_wmap_wodom * q_wodom_curr;
	t_w_curr = q_wmap_wodom * t_wodom_curr + t_wmap_wodom;
}

void transformUpdate()
{
	q_wmap_wodom = q_w_curr * q_wodom_curr.inverse();
	t_wmap_wodom = t_w_curr - q_wmap_wodom * t_wodom_curr;
}

void pointAssociateToMap(PointT const *const pi, PointT *const po)
{
	Eigen::Vector3d point_curr(pi->x, pi->y, pi->z);
	Eigen::Vector3d point_w = q_w_curr * point_curr + t_w_curr;
	po->x = point_w.x();
	po->y = point_w.y();
	po->z = point_w.z();
	po->intensity = pi->intensity;
	//po->intensity = 1.0;
}

void pointAssociateTobeMapped_(PointT const *const pi, PointT *const po, Eigen::Quaterniond dq) // if use fake local map (pitch = roll = 0)
{
	Eigen::Vector3d point_w(pi->x, pi->y, pi->z);
	Eigen::Vector3d point_curr = dq * q_w_curr.inverse() * (point_w - t_w_curr);
	po->x = point_curr.x();
	po->y = point_curr.y();
	po->z = point_curr.z();
	po->intensity = pi->intensity;
}

void pointAssociateTobeMapped(PointT const *const pi, PointT *const po)
{
	Eigen::Vector3d point_w(pi->x, pi->y, pi->z);
	Eigen::Vector3d point_curr = q_w_curr.inverse() * (point_w - t_w_curr);
	po->x = point_curr.x();
	po->y = point_curr.y();
	po->z = point_curr.z();
	po->intensity = pi->intensity;
}

void laserCloudCornerLastHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudCornerLast2)
{
	mBuf.lock();
	cornerLastBuf.push(laserCloudCornerLast2);
	mBuf.unlock();
}

void laserCloudSurfLastHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudSurfLast2)
{
	mBuf.lock();
	surfLastBuf.push(laserCloudSurfLast2);
	mBuf.unlock();
}

void laserCloudFullResHandler(const sensor_msgs::PointCloud2ConstPtr &laserCloudFullRes2)
{
	mBuf.lock();
	fullResBuf.push(laserCloudFullRes2);
	mBuf.unlock();
}

//receive odomtry
void laserOdometryHandler(const nav_msgs::Odometry::ConstPtr &laserOdometry)
{
	mBuf.lock();
	odometryBuf.push(laserOdometry);
	mBuf.unlock();

	// high frequence publish
	Eigen::Quaterniond q_wodom_curr;
	Eigen::Vector3d t_wodom_curr;
	q_wodom_curr.x() = laserOdometry->pose.pose.orientation.x;
	q_wodom_curr.y() = laserOdometry->pose.pose.orientation.y;
	q_wodom_curr.z() = laserOdometry->pose.pose.orientation.z;
	q_wodom_curr.w() = laserOdometry->pose.pose.orientation.w;
	t_wodom_curr.x() = laserOdometry->pose.pose.position.x;
	t_wodom_curr.y() = laserOdometry->pose.pose.position.y;
	t_wodom_curr.z() = laserOdometry->pose.pose.position.z;

	Eigen::Quaterniond q_w_curr = q_wmap_wodom * q_wodom_curr;
	Eigen::Vector3d t_w_curr = q_wmap_wodom * t_wodom_curr + t_wmap_wodom;
	// Eigen::Quaterniond q_w_curr = q_wodom_curr;
	// Eigen::Vector3d t_w_curr = t_wodom_curr;

	nav_msgs::Odometry odomAftMapped;
	odomAftMapped.header.frame_id = "/world";
	odomAftMapped.child_frame_id = "/aft_mapped";
	odomAftMapped.header.stamp = laserOdometry->header.stamp;
	odomAftMapped.pose.pose.orientation.x = q_w_curr.x();
	odomAftMapped.pose.pose.orientation.y = q_w_curr.y();
	odomAftMapped.pose.pose.orientation.z = q_w_curr.z();
	odomAftMapped.pose.pose.orientation.w = q_w_curr.w();
	odomAftMapped.pose.pose.position.x = t_w_curr.x();
	odomAftMapped.pose.pose.position.y = t_w_curr.y();
	odomAftMapped.pose.pose.position.z = t_w_curr.z();
	pubOdomAftMappedHighFrec.publish(odomAftMapped);
}

void process()
{
	while (1)
	{
		while (!cornerLastBuf.empty() && !surfLastBuf.empty() &&
			   !fullResBuf.empty() && !odometryBuf.empty())
		{
			mBuf.lock();
			while (!odometryBuf.empty() && odometryBuf.front()->header.stamp.toSec() < cornerLastBuf.front()->header.stamp.toSec())
				odometryBuf.pop();
			if (odometryBuf.empty())
			{
				mBuf.unlock();
				break;
			}

			while (!surfLastBuf.empty() && surfLastBuf.front()->header.stamp.toSec() < cornerLastBuf.front()->header.stamp.toSec())
				surfLastBuf.pop();
			if (surfLastBuf.empty())
			{
				mBuf.unlock();
				break;
			}

			while (!fullResBuf.empty() && fullResBuf.front()->header.stamp.toSec() < cornerLastBuf.front()->header.stamp.toSec())
				fullResBuf.pop();
			if (fullResBuf.empty())
			{
				mBuf.unlock();
				break;
			}

			timeLaserCloudCornerLast = cornerLastBuf.front()->header.stamp.toSec();
			timeLaserCloudSurfLast = surfLastBuf.front()->header.stamp.toSec();
			timeLaserCloudFullRes = fullResBuf.front()->header.stamp.toSec();
			timeLaserOdometry = odometryBuf.front()->header.stamp.toSec();

			// if (timeLaserCloudCornerLast != timeLaserOdometry ||
			// 	timeLaserCloudSurfLast != timeLaserOdometry ||
			// 	timeLaserCloudFullRes != timeLaserOdometry)
			// {
			// 	printf("time corner %f surf %f full %f odom %f \n", timeLaserCloudCornerLast, timeLaserCloudSurfLast, timeLaserCloudFullRes, timeLaserOdometry);
			// 	printf("unsync messeage!");
			// 	// std::cout << "!!!!!!!!!!!!test" << std::endl;
			// 	mBuf.unlock();
			// 	break;
			// }

			laserCloudCornerLast->clear();
			pcl::fromROSMsg(*cornerLastBuf.front(), *laserCloudCornerLast);
			cornerLastBuf.pop();

			laserCloudSurfLast->clear();
			pcl::fromROSMsg(*surfLastBuf.front(), *laserCloudSurfLast);
			surfLastBuf.pop();

			laserCloudFullRes->clear();
			pcl::fromROSMsg(*fullResBuf.front(), *laserCloudFullRes);
			fullResBuf.pop();

			q_wodom_curr.x() = odometryBuf.front()->pose.pose.orientation.x;
			q_wodom_curr.y() = odometryBuf.front()->pose.pose.orientation.y;
			q_wodom_curr.z() = odometryBuf.front()->pose.pose.orientation.z;
			q_wodom_curr.w() = odometryBuf.front()->pose.pose.orientation.w;
			t_wodom_curr.x() = odometryBuf.front()->pose.pose.position.x;
			t_wodom_curr.y() = odometryBuf.front()->pose.pose.position.y;
			t_wodom_curr.z() = odometryBuf.front()->pose.pose.position.z;
			odometryBuf.pop();

			while (!cornerLastBuf.empty())
			{
				cornerLastBuf.pop();
				//printf("drop lidar frame in mapping for real time performance \n");
			}

			mBuf.unlock();

			TicToc t_whole;

			transformAssociateToMap();

			TicToc t_shift;
			int centerCubeI = int((t_w_curr.x() + HalfCubeLength) / CubeLength) + laserCloudCenWidth;
			int centerCubeJ = int((t_w_curr.y() + HalfCubeLength) / CubeLength) + laserCloudCenHeight;
			int centerCubeK = int((t_w_curr.z() + HalfCubeHeight) / CubeHeight) + laserCloudCenDepth;

			if (t_w_curr.x() + HalfCubeLength < 0)
				centerCubeI--;
			if (t_w_curr.y() + HalfCubeLength < 0)
				centerCubeJ--;
			if (t_w_curr.z() + HalfCubeHeight < 0)
				centerCubeK--;

			while (centerCubeI < 3)
			{
				for (int j = 0; j < laserCloudHeight; j++)
				{
					for (int k = 0; k < laserCloudDepth; k++)
					{
						int i = laserCloudWidth - 1;
						pcl::PointCloud<PointT>::Ptr laserCloudCubeCornerPointer =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						pcl::PointCloud<PointT>::Ptr laserCloudCubeSurfPointer =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						for (; i >= 1; i--)
						{
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudCornerArray[i - 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudSurfArray[i - 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						}
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeCornerPointer;
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeSurfPointer;
						laserCloudCubeCornerPointer->clear();
						laserCloudCubeSurfPointer->clear();
					}
				}

				centerCubeI++;
				laserCloudCenWidth++;
			}

			while (centerCubeI >= laserCloudWidth - 3)
			{
				for (int j = 0; j < laserCloudHeight; j++)
				{
					for (int k = 0; k < laserCloudDepth; k++)
					{
						int i = 0;
						pcl::PointCloud<PointT>::Ptr laserCloudCubeCornerPointer =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						pcl::PointCloud<PointT>::Ptr laserCloudCubeSurfPointer =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						for (; i < laserCloudWidth - 1; i++)
						{
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudCornerArray[i + 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudSurfArray[i + 1 + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						}
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeCornerPointer;
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeSurfPointer;
						laserCloudCubeCornerPointer->clear();
						laserCloudCubeSurfPointer->clear();
					}
				}

				centerCubeI--;
				laserCloudCenWidth--;
			}

			while (centerCubeJ < 3)
			{
				for (int i = 0; i < laserCloudWidth; i++)
				{
					for (int k = 0; k < laserCloudDepth; k++)
					{
						int j = laserCloudHeight - 1;
						pcl::PointCloud<PointT>::Ptr laserCloudCubeCornerPointer =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						pcl::PointCloud<PointT>::Ptr laserCloudCubeSurfPointer =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						for (; j >= 1; j--)
						{
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudCornerArray[i + laserCloudWidth * (j - 1) + laserCloudWidth * laserCloudHeight * k];
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudSurfArray[i + laserCloudWidth * (j - 1) + laserCloudWidth * laserCloudHeight * k];
						}
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeCornerPointer;
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeSurfPointer;
						laserCloudCubeCornerPointer->clear();
						laserCloudCubeSurfPointer->clear();
					}
				}

				centerCubeJ++;
				laserCloudCenHeight++;
			}

			while (centerCubeJ >= laserCloudHeight - 3)
			{
				for (int i = 0; i < laserCloudWidth; i++)
				{
					for (int k = 0; k < laserCloudDepth; k++)
					{
						int j = 0;
						pcl::PointCloud<PointT>::Ptr laserCloudCubeCornerPointer =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						pcl::PointCloud<PointT>::Ptr laserCloudCubeSurfPointer =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						for (; j < laserCloudHeight - 1; j++)
						{
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudCornerArray[i + laserCloudWidth * (j + 1) + laserCloudWidth * laserCloudHeight * k];
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudSurfArray[i + laserCloudWidth * (j + 1) + laserCloudWidth * laserCloudHeight * k];
						}
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeCornerPointer;
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeSurfPointer;
						laserCloudCubeCornerPointer->clear();
						laserCloudCubeSurfPointer->clear();
					}
				}

				centerCubeJ--;
				laserCloudCenHeight--;
			}

			while (centerCubeK < 3)
			{
				for (int i = 0; i < laserCloudWidth; i++)
				{
					for (int j = 0; j < laserCloudHeight; j++)
					{
						int k = laserCloudDepth - 1;
						pcl::PointCloud<PointT>::Ptr laserCloudCubeCornerPointer =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						pcl::PointCloud<PointT>::Ptr laserCloudCubeSurfPointer =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						for (; k >= 1; k--)
						{
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k - 1)];
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k - 1)];
						}
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeCornerPointer;
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeSurfPointer;
						laserCloudCubeCornerPointer->clear();
						laserCloudCubeSurfPointer->clear();
					}
				}

				centerCubeK++;
				laserCloudCenDepth++;
			}

			while (centerCubeK >= laserCloudDepth - 3)
			{
				for (int i = 0; i < laserCloudWidth; i++)
				{
					for (int j = 0; j < laserCloudHeight; j++)
					{
						int k = 0;
						pcl::PointCloud<PointT>::Ptr laserCloudCubeCornerPointer =
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						pcl::PointCloud<PointT>::Ptr laserCloudCubeSurfPointer =
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k];
						for (; k < laserCloudDepth - 1; k++)
						{
							laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k + 1)];
							laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
								laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * (k + 1)];
						}
						laserCloudCornerArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeCornerPointer;
						laserCloudSurfArray[i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k] =
							laserCloudCubeSurfPointer;
						laserCloudCubeCornerPointer->clear();
						laserCloudCubeSurfPointer->clear();
					}
				}

				centerCubeK--;
				laserCloudCenDepth--;
			}

			int laserCloudValidNum = 0;
			// int laserCloudSurroundNum = 0;

			for (int i = centerCubeI - 2; i <= centerCubeI + 2; i++)
			{
				for (int j = centerCubeJ - 2; j <= centerCubeJ + 2; j++)
				{
					for (int k = centerCubeK - 1; k <= centerCubeK + 1; k++)
					{
						if (i >= 0 && i < laserCloudWidth &&
							j >= 0 && j < laserCloudHeight &&
							k >= 0 && k < laserCloudDepth)
						{
							laserCloudValidInd[laserCloudValidNum] = i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k;
							laserCloudValidNum++;
							// laserCloudSurroundInd[laserCloudSurroundNum] = i + laserCloudWidth * j + laserCloudWidth * laserCloudHeight * k;
							// laserCloudSurroundNum++;
						}
					}
				}
			}

			laserCloudCornerFromMap->clear();
			laserCloudSurfFromMap->clear();
			for (int i = 0; i < laserCloudValidNum; i++)
			{
				*laserCloudCornerFromMap += *laserCloudCornerArray[laserCloudValidInd[i]];
				*laserCloudSurfFromMap += *laserCloudSurfArray[laserCloudValidInd[i]];
			}
			int laserCloudCornerFromMapNum = laserCloudCornerFromMap->points.size();
			int laserCloudSurfFromMapNum = laserCloudSurfFromMap->points.size();

			pcl::PointCloud<PointT>::Ptr laserCloudCornerStack(new pcl::PointCloud<PointT>());
			downSizeFilterCorner.setInputCloud(laserCloudCornerLast);
			downSizeFilterCorner.filter(*laserCloudCornerStack);
			// OutRemFilter.setInputCloud(laserCloudCornerStack);
			// OutRemFilter.filter(*laserCloudCornerStack);
			int laserCloudCornerStackNum = laserCloudCornerStack->points.size();

			pcl::PointCloud<PointT>::Ptr laserCloudSurfStack(new pcl::PointCloud<PointT>());
			downSizeFilterSurf.setInputCloud(laserCloudSurfLast);
			downSizeFilterSurf.filter(*laserCloudSurfStack);
			// OutRemFilter.setInputCloud(laserCloudSurfStack);
			// OutRemFilter.filter(*laserCloudSurfStack);
			int laserCloudSurfStackNum = laserCloudSurfStack->points.size();

			//printf("map prepare time %f ms\n", t_shift.toc());
			//printf("map corner num %d  surf num %d \n", laserCloudCornerFromMapNum, laserCloudSurfFromMapNum);
			if (laserCloudCornerFromMapNum > 10 && laserCloudSurfFromMapNum > 50)
			{
				TicToc t_opt;
				TicToc t_tree;
				kdtreeCornerFromMap->setInputCloud(laserCloudCornerFromMap);
				kdtreeSurfFromMap->setInputCloud(laserCloudSurfFromMap);
				//printf("build tree time %f ms \n", t_tree.toc());

				for (int iterCount = 0; iterCount < 2; iterCount++)
				{
					//ceres::LossFunction *loss_function = NULL;
					ceres::LossFunction *loss_function = new ceres::HuberLoss(0.1);
					ceres::LocalParameterization *q_parameterization =
						new ceres::EigenQuaternionParameterization();
					ceres::Problem::Options problem_options;

					ceres::Problem problem(problem_options);
					problem.AddParameterBlock(parameters, 4, q_parameterization);
					problem.AddParameterBlock(parameters + 4, 3);
					// problem.SetParameterUpperBound(parameters + 4, 2, t_wodom_curr.z() + 0.3); //1.0
					// problem.SetParameterLowerBound(parameters + 4, 2, t_wodom_curr.z() - 0.3);

					TicToc t_data;
					int corner_num = 0;

					for (int i = 0; i < laserCloudCornerStackNum; i++)
					{
						pointOri = laserCloudCornerStack->points[i];
						//double sqrtDis = pointOri.x * pointOri.x + pointOri.y * pointOri.y + pointOri.z * pointOri.z;
						pointAssociateToMap(&pointOri, &pointSel);
						kdtreeCornerFromMap->nearestKSearch(pointSel, 5, pointSearchInd, pointSearchSqDis);

						// std::vector<std::pair<float, size_t> > diff_v;
						// for(int k = 0; k < pointSearchInd.size(); k++ ){
						// 	float diff = fabs(laserCloudCornerFromMap->points[pointSearchInd[k]].intensity - pointSel.intensity);
						// 	std::pair<float, size_t> diff_p (diff, pointSearchInd[k]);
						// 	diff_v.push_back(diff_p);
						// }
						// std::sort(diff_v.begin(), diff_v.end(), std::less<>());
						// pointSearchInd.clear();
						// for(int n = 0; n < 5; n++){
						// 	pointSearchInd.push_back(diff_v[n].second);
						// }

						if (pointSearchSqDis[4] < 1.0)
						{
							std::vector<Eigen::Vector3d> nearCorners;
							Eigen::Vector3d center(0, 0, 0);
							for (int j = 0; j < 5; j++)
							{
								Eigen::Vector3d tmp(laserCloudCornerFromMap->points[pointSearchInd[j]].x,
													laserCloudCornerFromMap->points[pointSearchInd[j]].y,
													laserCloudCornerFromMap->points[pointSearchInd[j]].z);
								center = center + tmp;
								nearCorners.push_back(tmp);
							}
							center = center / 5.0;

							Eigen::Matrix3d covMat = Eigen::Matrix3d::Zero();
							for (int j = 0; j < 5; j++)
							{
								Eigen::Matrix<double, 3, 1> tmpZeroMean = nearCorners[j] - center;
								covMat = covMat + tmpZeroMean * tmpZeroMean.transpose();
							}

							Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> saes(covMat);

							// if is indeed line feature
							// note Eigen library sort eigenvalues in increasing order
							Eigen::Vector3d unit_direction = saes.eigenvectors().col(2);
							Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
							if (saes.eigenvalues()[2] > 3 * saes.eigenvalues()[1])
							{
								Eigen::Vector3d point_on_line = center;
								Eigen::Vector3d point_a, point_b;
								point_a = 0.1 * unit_direction + point_on_line;
								point_b = -0.1 * unit_direction + point_on_line;

								ceres::CostFunction *cost_function = LidarEdgeFactor::Create(curr_point, point_a, point_b, 1.0);
								problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
								corner_num++;
							}
						}
						/*
						else if(pointSearchSqDis[4] < 0.01 * sqrtDis)
						{
							Eigen::Vector3d center(0, 0, 0);
							for (int j = 0; j < 5; j++)
							{
								Eigen::Vector3d tmp(laserCloudCornerFromMap->points[pointSearchInd[j]].x,
													laserCloudCornerFromMap->points[pointSearchInd[j]].y,
													laserCloudCornerFromMap->points[pointSearchInd[j]].z);
								center = center + tmp;
							}
							center = center / 5.0;	
							Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
							ceres::CostFunction *cost_function = LidarDistanceFactor::Create(curr_point, center);
							problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
						}
						*/
					}

					int surf_num = 0;

					for (int i = 0; i < laserCloudSurfStackNum; i++)
					{
						pointOri = laserCloudSurfStack->points[i];
						//double sqrtDis = pointOri.x * pointOri.x + pointOri.y * pointOri.y + pointOri.z * pointOri.z;
						pointAssociateToMap(&pointOri, &pointSel);
						kdtreeSurfFromMap->nearestKSearch(pointSel, 10, pointSearchInd, pointSearchSqDis);

						std::vector<std::pair<float, size_t> > diff_v;
						for (int k = 0; k < pointSearchInd.size(); k++)
						{
							float diff = fabs(laserCloudSurfFromMap->points[pointSearchInd[k]].intensity - pointSel.intensity);
							std::pair<float, size_t> diff_p(diff, pointSearchInd[k]);
							diff_v.push_back(diff_p);
						}
						std::sort(diff_v.begin(), diff_v.end(), std::less<>());
						pointSearchInd.clear();
						for (int n = 0; n < 5; n++)
						{
							pointSearchInd.push_back(diff_v[n].second);
						}

						Eigen::Matrix<double, 5, 3> matA0;
						Eigen::Matrix<double, 5, 1> matB0 = -1 * Eigen::Matrix<double, 5, 1>::Ones();
						if (pointSearchSqDis[4] < 1.0)
						{

							for (int j = 0; j < 5; j++)
							{
								matA0(j, 0) = laserCloudSurfFromMap->points[pointSearchInd[j]].x;
								matA0(j, 1) = laserCloudSurfFromMap->points[pointSearchInd[j]].y;
								matA0(j, 2) = laserCloudSurfFromMap->points[pointSearchInd[j]].z;
								////printf(" pts %f %f %f ", matA0(j, 0), matA0(j, 1), matA0(j, 2));
							}
							// find the norm of plane
							Eigen::Vector3d norm = matA0.colPivHouseholderQr().solve(matB0);
							double negative_OA_dot_norm = 1 / norm.norm();
							norm.normalize();

							// Here n(pa, pb, pc) is unit norm of plane
							bool planeValid = true;
							for (int j = 0; j < 5; j++)
							{
								// if OX * n > 0.2, then plane is not fit well
								if (fabs(norm(0) * laserCloudSurfFromMap->points[pointSearchInd[j]].x +
										 norm(1) * laserCloudSurfFromMap->points[pointSearchInd[j]].y +
										 norm(2) * laserCloudSurfFromMap->points[pointSearchInd[j]].z + negative_OA_dot_norm) > 0.2)
								{
									planeValid = false;
									break;
								}
							}
							Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
							if (planeValid)
							{
								ceres::CostFunction *cost_function = LidarPlaneNormFactor::Create(curr_point, norm, negative_OA_dot_norm);
								problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
								surf_num++;
							}
						}
						/*
						else if(pointSearchSqDis[4] < 0.01 * sqrtDis)
						{
							Eigen::Vector3d center(0, 0, 0);
							for (int j = 0; j < 5; j++)
							{
								Eigen::Vector3d tmp(laserCloudSurfFromMap->points[pointSearchInd[j]].x,
													laserCloudSurfFromMap->points[pointSearchInd[j]].y,
													laserCloudSurfFromMap->points[pointSearchInd[j]].z);
								center = center + tmp;
							}
							center = center / 5.0;	
							Eigen::Vector3d curr_point(pointOri.x, pointOri.y, pointOri.z);
							ceres::CostFunction *cost_function = LidarDistanceFactor::Create(curr_point, center);
							problem.AddResidualBlock(cost_function, loss_function, parameters, parameters + 4);
						}
						*/
					}

					////printf("corner num %d used corner num %d \n", laserCloudCornerStackNum, corner_num);
					////printf("surf num %d used surf num %d \n", laserCloudSurfStackNum, surf_num);

					//printf("mapping data assosiation time %f ms \n", t_data.toc());

					TicToc t_solver;
					ceres::Solver::Options options;
					options.linear_solver_type = ceres::DENSE_SCHUR;
					options.trust_region_strategy_type = ceres::DOGLEG;
					options.max_solver_time_in_seconds = 0.05;
					// options.linear_solver_type = ceres::DENSE_QR;
					options.max_num_iterations = 4;
					options.minimizer_progress_to_stdout = false;
					options.check_gradients = false;
					options.gradient_check_relative_precision = 1e-4;
					ceres::Solver::Summary summary;
					ceres::Solve(options, &problem, &summary);
					//printf("mapping solver time %f ms \n", t_solver.toc());

					////printf("time %f \n", timeLaserOdometry);
					////printf("corner factor num %d surf factor num %d\n", corner_num, surf_num);
					////printf("result q %f %f %f %f result t %f %f %f\n", parameters[3], parameters[0], parameters[1], parameters[2],
					//	   parameters[4], parameters[5], parameters[6]);
				}
				//printf("mapping optimization time %f \n", t_opt.toc());
			}
			else
			{
				ROS_WARN("time Map corner and surf num are not enough");
			}
			transformUpdate();

			TicToc t_add;

			for (int i = 0; i < laserCloudCornerStackNum; i++)
			{
				pointAssociateToMap(&laserCloudCornerStack->points[i], &pointSel);

				int cubeI = int((pointSel.x + HalfCubeLength) / CubeLength) + laserCloudCenWidth;
				int cubeJ = int((pointSel.y + HalfCubeLength) / CubeLength) + laserCloudCenHeight;
				int cubeK = int((pointSel.z + HalfCubeHeight) / CubeHeight) + laserCloudCenDepth;

				if (pointSel.x + HalfCubeLength < 0)
					cubeI--;
				if (pointSel.y + HalfCubeLength < 0)
					cubeJ--;
				if (pointSel.z + HalfCubeHeight < 0)
					cubeK--;

				if (cubeI >= 0 && cubeI < laserCloudWidth &&
					cubeJ >= 0 && cubeJ < laserCloudHeight &&
					cubeK >= 0 && cubeK < laserCloudDepth)
				{
					int cubeInd = cubeI + laserCloudWidth * cubeJ + laserCloudWidth * laserCloudHeight * cubeK;
					laserCloudCornerArray[cubeInd]->push_back(pointSel);
				}
			}

			for (int i = 0; i < laserCloudSurfStackNum; i++)
			{
				pointAssociateToMap(&laserCloudSurfStack->points[i], &pointSel);

				int cubeI = int((pointSel.x + HalfCubeLength) / CubeLength) + laserCloudCenWidth;
				int cubeJ = int((pointSel.y + HalfCubeLength) / CubeLength) + laserCloudCenHeight;
				int cubeK = int((pointSel.z + HalfCubeHeight) / CubeHeight) + laserCloudCenDepth;

				if (pointSel.x + HalfCubeLength < 0)
					cubeI--;
				if (pointSel.y + HalfCubeLength < 0)
					cubeJ--;
				if (pointSel.z + HalfCubeHeight < 0)
					cubeK--;

				if (cubeI >= 0 && cubeI < laserCloudWidth &&
					cubeJ >= 0 && cubeJ < laserCloudHeight &&
					cubeK >= 0 && cubeK < laserCloudDepth)
				{
					int cubeInd = cubeI + laserCloudWidth * cubeJ + laserCloudWidth * laserCloudHeight * cubeK;
					laserCloudSurfArray[cubeInd]->push_back(pointSel);
				}
			}
			//printf("add points time %f ms\n", t_add.toc());

			TicToc t_filter;

			for (int i = 0; i < laserCloudValidNum; i++)
			{
				int ind = laserCloudValidInd[i];

				pcl::PointCloud<PointT>::Ptr tmpCorner(new pcl::PointCloud<PointT>());
				downSizeFilterCorner.setInputCloud(laserCloudCornerArray[ind]);
				downSizeFilterCorner.filter(*tmpCorner);

				laserCloudCornerArray[ind] = tmpCorner;

				pcl::PointCloud<PointT>::Ptr tmpSurf(new pcl::PointCloud<PointT>());
				downSizeFilterSurf.setInputCloud(laserCloudSurfArray[ind]);
				downSizeFilterSurf.filter(*tmpSurf);

				laserCloudSurfArray[ind] = tmpSurf;
			}
			//printf("filter time %f ms \n", t_filter.toc());

			TicToc t_pub;
			//publish surround map for every 5 frame
			// if (frameCount % 5 == 0) //5
			// {
			// 	laserCloudSurround->clear();
			// 	for (int i = 0; i < laserCloudSurroundNum; i++)
			// 	{
			// 		int ind = laserCloudSurroundInd[i];
			// 		*laserCloudSurround += *laserCloudCornerArray[ind];
			// 		*laserCloudSurround += *laserCloudSurfArray[ind];
			// 	}

			// 	// OutRemFilter.setInputCloud(laserCloudSurround);
			// 	// OutRemFilter.filter(*laserCloudSurround);

			// 	sensor_msgs::PointCloud2 laserCloudSurround3;
			// 	pcl::toROSMsg(*laserCloudSurround, laserCloudSurround3);
			// 	laserCloudSurround3.header.stamp = ros::Time().fromSec(timeLaserOdometry);
			// 	laserCloudSurround3.header.frame_id = "/world";
			// 	pubLaserCloudSurround.publish(laserCloudSurround3);
			// }

			bool publocal{false};
			if (first_odom)
			{
				q_w_last = q_w_curr;
				t_w_last = t_w_curr;
				// first_odom = false;
				last_frame_count = frameCount;
				if (frameCount >= 10)
				{
					first_odom = false;
					publocal = true;
				}
			}
			else
			{
				Eigen::Matrix3d R_W_c = q_w_curr.toRotationMatrix();
				Eigen::Vector3d diff = R_W_c.transpose() * (t_w_last - t_w_curr);
				int diff_c = frameCount - last_frame_count;
				// if((diff.norm() > 2.0 && diff_c > 5) || diff_c > 10){
				if (diff.norm() > 2.0 || diff_c > 30)
				{
					q_w_last = q_w_curr;
					t_w_last = t_w_curr;
					publocal = true;
					last_frame_count = frameCount;
				}
			}

			if (publocal) //20  frameCount % 5 == 0
			{
				pcl::PointCloud<PointT> laserCloudMap;

				for (int i = 0; i < laserCloudNum; i++)
				{
					laserCloudMap += *laserCloudCornerArray[i];
					laserCloudMap += *laserCloudSurfArray[i];
#ifdef FOR_GLOBAL
					laserCloudCornerArray[i].reset(new pcl::PointCloud<PointT>());
					laserCloudSurfArray[i].reset(new pcl::PointCloud<PointT>());
#endif
				}

				// if use fake local map (pitch = roll = 0)
				// Eigen::Quaterniond temq = q_w_curr;
				// Eigen::Vector3d fake_euler = q_w_curr.toRotationMatrix().eulerAngles(2, 1, 0);
				// Eigen::Matrix3d fake_r;
				// fake_r = Eigen::AngleAxisd(fake_euler[0], Eigen::Vector3d::UnitZ()) * Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(0.0, Eigen::Vector3d::UnitX());
				// Eigen::Quaterniond fake_q{fake_r};
				// Eigen::Quaterniond fake_dq = fake_q.inverse() * q_w_curr;

				sensor_msgs::PointCloud2 laserCloudMsg;
				pcl::toROSMsg(laserCloudMap, laserCloudMsg);
				laserCloudMsg.header.stamp = ros::Time().fromSec(timeLaserOdometry);
				laserCloudMsg.header.frame_id = "/world";
				pubLaserCloudMap.publish(laserCloudMsg);
#ifdef FOR_GLOBAL
				int MapNum = laserCloudMap.points.size();
				for (int i = 0; i < MapNum; i++)
				{
					pointAssociateTobeMapped(&laserCloudMap.points[i], &laserCloudMap.points[i]);
					// pointAssociateTobeMapped_(&laserCloudMap.points[i], &laserCloudMap.points[i], fake_dq); // if use fake local map (pitch = roll = 0)
				}
				sensor_msgs::PointCloud2 laserCloudMapLocal; //local_map
				pcl::toROSMsg(laserCloudMap, laserCloudMapLocal);
				laserCloudMapLocal.header.stamp = ros::Time().fromSec(timeLaserOdometry);
				laserCloudMapLocal.header.frame_id = "/world";
				LocalMap.publish(laserCloudMapLocal);

				nav_msgs::Odometry odomLocalMap;
				odomLocalMap.header.frame_id = "/world";
				odomLocalMap.child_frame_id = "/aft_mapped";
				odomLocalMap.header.stamp = ros::Time().fromSec(timeLaserOdometry);
				odomLocalMap.pose.pose.orientation.x = q_w_curr.x();
				odomLocalMap.pose.pose.orientation.y = q_w_curr.y();
				odomLocalMap.pose.pose.orientation.z = q_w_curr.z();
				odomLocalMap.pose.pose.orientation.w = q_w_curr.w();
				// if use fake local map (pitch = roll = 0)
				// odomLocalMap.pose.pose.orientation.x = fake_q.x();
				// odomLocalMap.pose.pose.orientation.y = fake_q.y();
				// odomLocalMap.pose.pose.orientation.z = fake_q.z();
				// odomLocalMap.pose.pose.orientation.w = fake_q.w();
				odomLocalMap.pose.pose.position.x = t_w_curr.x();
				odomLocalMap.pose.pose.position.y = t_w_curr.y();
				odomLocalMap.pose.pose.position.z = t_w_curr.z();
				LocalOdom.publish(odomLocalMap);
				// // ROS_WARN_STREAM("TIME LM-H-NODE: "  << std::fixed << odomLocalMap.header.stamp.toSec() << std::endl);
				q_wmap_wodom = Eigen::Quaterniond::Identity(); // better outdoor if not set ZERO, not big difference for indoor.
				t_wmap_wodom = Eigen::Vector3d::Zero();
#endif
			}

			sensor_msgs::PointCloud2 laserCloudFullRes3Local; //velodyne_cloud_registered_local
			pcl::toROSMsg(*laserCloudFullRes, laserCloudFullRes3Local);
			laserCloudFullRes3Local.header.stamp = ros::Time().fromSec(timeLaserOdometry);
			laserCloudFullRes3Local.header.frame_id = "/world";
			pubLaserCloudFullResLocal.publish(laserCloudFullRes3Local);

			int laserCloudFullResNum = laserCloudFullRes->points.size();

			for (int i = 0; i < laserCloudFullResNum; i++)
			{
				pointAssociateToMap(&laserCloudFullRes->points[i], &laserCloudFullRes->points[i]);
			}

#ifndef FOR_GLOBAL
            pcl::ApproximateVoxelGrid<PointT> voxelgrid;
            voxelgrid.setLeafSize(0.1, 0.1, 0.1);
            voxelgrid.setInputCloud(laserCloudFullRes);
            voxelgrid.filter(*laserCloudFullRes);
#endif

			sensor_msgs::PointCloud2 laserCloudFullRes3; //velodyne_cloud_registered
			pcl::toROSMsg(*laserCloudFullRes, laserCloudFullRes3);
			laserCloudFullRes3.header.stamp = ros::Time().fromSec(timeLaserOdometry);
			laserCloudFullRes3.header.frame_id = "/world";
			pubLaserCloudFullRes.publish(laserCloudFullRes3);

			//printf("mapping pub time %f ms \n", t_pub.toc());

			//printf("whole mapping time %f ms ++++++++++\n", t_whole.toc());

			nav_msgs::Odometry odomAftMapped;
			odomAftMapped.header.frame_id = "/world";
			odomAftMapped.child_frame_id = "/aft_mapped";
			odomAftMapped.header.stamp = ros::Time().fromSec(timeLaserOdometry);
			odomAftMapped.pose.pose.orientation.x = q_w_curr.x();
			odomAftMapped.pose.pose.orientation.y = q_w_curr.y();
			odomAftMapped.pose.pose.orientation.z = q_w_curr.z();
			odomAftMapped.pose.pose.orientation.w = q_w_curr.w();
			odomAftMapped.pose.pose.position.x = t_w_curr.x();
			odomAftMapped.pose.pose.position.y = t_w_curr.y();
			odomAftMapped.pose.pose.position.z = t_w_curr.z();
			pubOdomAftMapped.publish(odomAftMapped);

			{
				// std::ofstream foutC(result_path, std::ios::trunc);
				std::ofstream foutC(result_path, std::ios::app);
				foutC.setf(std::ios::fixed, std::ios::floatfield);
				foutC.precision(9);
				foutC << odomAftMapped.header.stamp.toSec() << " ";
				foutC.precision(5);
				foutC << t_w_curr.x() << " "
					  << t_w_curr.y() << " "
					  << t_w_curr.z() << " "
					  << q_w_curr.x() << " "
					  << q_w_curr.y() << " "
					  << q_w_curr.z() << " "
					  << q_w_curr.w() << std::endl;
				foutC.close();
			}

			// double roll, pitch, yaw;
			// geometry_msgs::Quaternion quat = odomAftMapped.pose.pose.orientation;
			// tf::Matrix3x3(tf::Quaternion(quat.x, quat.y, quat.z, quat.w)).getRPY(roll, pitch, yaw);
			// std::cout << "pub - roll, pitch, yaw (deg): " << rad2deg(roll) << ", " << rad2deg(pitch) << ", " << rad2deg(yaw) << std::endl;

			geometry_msgs::PoseStamped laserAfterMappedPose;
			laserAfterMappedPose.header = odomAftMapped.header;
			laserAfterMappedPose.pose = odomAftMapped.pose.pose;

			laserAfterMappedPath.header.stamp = odomAftMapped.header.stamp;
			laserAfterMappedPath.header.frame_id = "/world";
			laserAfterMappedPath.poses.push_back(laserAfterMappedPose);
			pubLaserAfterMappedPath.publish(laserAfterMappedPath);

			static tf::TransformBroadcaster br;
			tf::Transform transform;
			tf::Quaternion q;
			transform.setOrigin(tf::Vector3(t_w_curr(0),
											t_w_curr(1),
											t_w_curr(2)));
			q.setW(q_w_curr.w());
			q.setX(q_w_curr.x());
			q.setY(q_w_curr.y());
			q.setZ(q_w_curr.z());
			transform.setRotation(q);
			br.sendTransform(tf::StampedTransform(transform, odomAftMapped.header.stamp, "/world", "/aft_mapped"));

			// static tf::TransformBroadcaster br_;
			// tf::Transform transform_;
			// tf::Quaternion q_;
			// transform_.setOrigin(tf::Vector3(t_wmap_wodom(0),
			// 								 t_wmap_wodom(1),
			// 								 t_wmap_wodom(2)));
			// q_.setW(q_wmap_wodom.w());
			// q_.setX(q_wmap_wodom.x());
			// q_.setY(q_wmap_wodom.y());
			// q_.setZ(q_wmap_wodom.z());
			// transform_.setRotation(q_);
			// br_.sendTransform(tf::StampedTransform(transform_, odomAftMapped.header.stamp, "/map", "/world"));

			frameCount++;
		}
		std::chrono::milliseconds dura(2);
		std::this_thread::sleep_for(dura);
	}
}

int main(int argc, char **argv)
{
	ros::init(argc, argv, "laserMapping");
	ros::NodeHandle nh;

	float lineRes = 0;
	float planeRes = 0;
	std::string save_directory;
	nh.param<float>("mapping_line_resolution", lineRes, 0.2);
	nh.param<float>("mapping_plane_resolution", planeRes, 0.4);
	nh.param<std::string>("save_directory", save_directory, "/");
	result_path = save_directory + "Midend.txt";
	std::ofstream fout(result_path, std::ios::out);
	fout.close();

	downSizeFilterCorner.setLeafSize(lineRes, lineRes, lineRes);
	downSizeFilterSurf.setLeafSize(planeRes, planeRes, planeRes);
	OutRemFilter.setMeanK(50);
	OutRemFilter.setStddevMulThresh(0.3);
	ros::Subscriber subLaserCloudCornerLast = nh.subscribe<sensor_msgs::PointCloud2>("/laser_cloud_less_sharp", 5, laserCloudCornerLastHandler, ros::TransportHints().tcpNoDelay()); //laser_cloud_less_sharp

	ros::Subscriber subLaserCloudSurfLast = nh.subscribe<sensor_msgs::PointCloud2>("/laser_cloud_less_flat", 5, laserCloudSurfLastHandler, ros::TransportHints().tcpNoDelay()); //laser_cloud_less_flat

	ros::Subscriber subLaserOdometry = nh.subscribe<nav_msgs::Odometry>("/vils_estimator/lidar_pose", 5, laserOdometryHandler, ros::TransportHints().tcpNoDelay()); // /laser_odom_to_init  /vins_estimator/lidar_pose

	ros::Subscriber subLaserCloudFullRes = nh.subscribe<sensor_msgs::PointCloud2>("/velodyne_cloud_2", 5, laserCloudFullResHandler, ros::TransportHints().tcpNoDelay()); //velodyne_cloud_2

	// pubLaserCloudSurround = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_surround", 100);

	pubLaserCloudMap = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_map", 1);

	pubLaserCloudFullRes = nh.advertise<sensor_msgs::PointCloud2>("/velodyne_cloud_registered", 1);
	pubLaserCloudFullResLocal = nh.advertise<sensor_msgs::PointCloud2>("/velodyne_cloud_registered_local", 1);

	pubOdomAftMapped = nh.advertise<nav_msgs::Odometry>("/aft_mapped_to_init", 1);

	pubOdomAftMappedHighFrec = nh.advertise<nav_msgs::Odometry>("/aft_mapped_to_init_high_frec", 1);

	pubLaserAfterMappedPath = nh.advertise<nav_msgs::Path>("/aft_mapped_path", 1);

	LocalMap = nh.advertise<sensor_msgs::PointCloud2>("/local_map", 1);
	LocalOdom = nh.advertise<nav_msgs::Odometry>("/local_odom", 1);

	for (int i = 0; i < laserCloudNum; i++)
	{
		laserCloudCornerArray[i].reset(new pcl::PointCloud<PointT>());
		laserCloudSurfArray[i].reset(new pcl::PointCloud<PointT>());
	}

	std::thread mapping_process{process};

	// ros::spin();
	ros::AsyncSpinner spinner(4);
	spinner.start();
	ros::waitForShutdown();

	return 0;
}