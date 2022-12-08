#ifndef POINTPROCESSOR_H_
#define POINTPROCESSOR_H_

#include <pcl/filters/voxel_grid.h>

#include <pcl/features/normal_3d.h>

#include "common_ros.h"
// #include "TicToc.h"
#include "CircularBuffer.h"
#include "math_utils.h"
#include <std_msgs/Float32.h>

using namespace std;
typedef pcl::PointXYZI PointT;
typedef typename pcl::PointCloud<PointT> PointCloud;
typedef typename pcl::PointCloud<PointT>::Ptr PointCloudPtr;
typedef typename pcl::PointCloud<PointT>::ConstPtr PointCloudConstPtr;
typedef std::pair<size_t, size_t> IndexRange;

// adapted from LOAM
/** Point label options. */
enum PointLabel
{
  CORNER_SHARP = 2,      ///< sharp corner point
  CORNER_LESS_SHARP = 1, ///< less sharp corner point
  SURFACE_LESS_FLAT = 0, ///< less flat surface point
  SURFACE_FLAT = -1      ///< flat surface point
};

struct PointProcessorConfig
{
  bool deskew = false;
  double scan_period = 0.1;
  int num_scan_subregions = 8;
  int num_curvature_regions = 5;
  float surf_curv_th = 1.0;
  int max_corner_sharp = 3;
  int max_corner_less_sharp = 10 * max_corner_sharp;
  int max_surf_flat = 4;
  float less_flat_filter_size = 0.2;

  string capture_frame_id = "/laser_link"; // /laser_link  /velodyne /base_link

  double rad_diff = 0.2;

  bool infer_start_ori_ = false;
};

class PointProcessor
{

public:
  PointProcessor();
  PointProcessor(float lower_bound, float upper_bound, int num_rings, bool uneven = false);

  // WARNING: is it useful to separate Process and SetInputCloud?
  // void Process(const PointCloudConstPtr &cloud_in, PointCloud &cloud_out);
  void Process();

  void PointCloudHandler(const sensor_msgs::PointCloud2ConstPtr &raw_points_msg);

  void SetupConfig(PointProcessorConfig config)
  {
    config_ = config;
  }

  void SetupRos(ros::NodeHandle &nh);

  void SetInputCloud(const PointCloudConstPtr &cloud_in, ros::Time time_in = ros::Time::now());

  void PointToRing();
  void PointToRing(const PointCloudConstPtr &cloud_in,
                   vector<PointCloudPtr> &ring_out,
                   vector<PointCloudPtr> &intensity_out);

  inline int ElevationToRing(float rad)
  {
    double in = (RadToDeg(rad) - lower_bound_) * factor_ + 0.5;
    return int((RadToDeg(rad) - lower_bound_) * factor_ + 0.5);
  }

  void ExtractFeaturePoints();

  void PublishResults();

  // TODO: not necessary data?
  vector<PointCloudPtr> laser_scans;
  vector<PointCloudPtr> intensity_scans;
  vector<IndexRange> scan_ranges;

protected:
  ros::Time sweep_start_;
  ros::Time scan_time_;

  float lower_bound_;
  float upper_bound_;
  int num_rings_;
  float factor_;
  PointProcessorConfig config_;
  // TicToc tic_toc_;

  PointCloudConstPtr cloud_ptr_;
  PointCloud cloud_in_rings_;

  PointCloud corner_points_sharp_;
  PointCloud corner_points_less_sharp_;
  PointCloud surface_points_flat_;
  PointCloud surface_points_less_flat_;

  // the following will be assigened or resized
  vector<int> scan_ring_mask_;
  vector<pair<float, size_t> > curvature_idx_pairs_; // in subregion
  vector<PointLabel> subregion_labels_;              ///< point label buffer

  void Reset(const ros::Time &scan_time, const bool &is_new_sweep = true);
  void PrepareRing(const PointCloud &scan);
  void PrepareSubregion(const PointCloud &scan, const size_t idx_start, const size_t idx_end);
  void MaskPickedInRing(const PointCloud &scan, const size_t in_scan_idx);

  ros::Subscriber sub_raw_points_; ///< input cloud message subscriber

  ros::Publisher pub_full_cloud_;               ///< full resolution cloud message publisher
  ros::Publisher pub_corner_points_sharp_;      ///< sharp corner cloud message publisher
  ros::Publisher pub_corner_points_less_sharp_; ///< less sharp corner cloud message publisher
  ros::Publisher pub_surf_points_flat_;         ///< flat surface cloud message publisher
  ros::Publisher pub_surf_points_less_flat_;    ///< less flat surface cloud message publisher

  bool is_ros_setup_ = false;

  bool uneven_ = false;

private:
  float start_ori_, end_ori_;
  CircularBuffer<float> start_ori_buf1_{10};
  CircularBuffer<float> start_ori_buf2_{10};
};

#endif //POINTPROCESSOR_H_
