#include "../../lidar_compensator/src/PointProcessor.h"
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <ros/ros.h>
#include <sensor_msgs/PointCloud2.h>
#include <pcl_conversions/pcl_conversions.h>
#include <tf/transform_datatypes.h>
#include <geometry_msgs/Quaternion.h>

using namespace std;

double AbsRadDistance(double a, double b)
{
  return fabs(NormalizeRad(a - b));
}

PointProcessor::PointProcessor() : PointProcessor(-15.0, 15.0, 16, false) {}

PointProcessor::PointProcessor(float lower_bound, float upper_bound, int num_rings, bool uneven)
    : lower_bound_(lower_bound),
      upper_bound_(upper_bound),
      num_rings_(num_rings),
      factor_((num_rings - 1) / (upper_bound - lower_bound)),
      uneven_(uneven)
{
  PointCloud scan;
  laser_scans.clear();
  for (int i = 0; i < num_rings_; ++i)
  {
    PointCloudPtr scan(new PointCloud());
    laser_scans.push_back(scan);
  }

  intensity_scans.clear();
  for (int i = 0; i < num_rings_; ++i)
  {
    PointCloudPtr scan(new PointCloud());
    intensity_scans.push_back(scan);
  }
}

void PointProcessor::Process()
{
  PointToRing();
  ExtractFeaturePoints();
  PublishResults();
}

void PointProcessor::ExtractFeaturePoints()
{

  // tic_toc_.Tic();

  ///< i is #ring, j is #subregion, k is # in region
  for (size_t i = 0; i < num_rings_; ++i)
  {

    PointCloud::Ptr surf_points_less_flat_ptr(new PointCloud());

    size_t start_idx = scan_ranges[i].first;
    size_t end_idx = scan_ranges[i].second;

    // skip too short scans
    if (end_idx <= start_idx + 2 * config_.num_curvature_regions)
    {
      continue;
    }

    const PointCloud &scan_ring = *laser_scans[i];
    size_t scan_size = scan_ring.size();

    PrepareRing(scan_ring);

    // extract features from equally sized scan regions
    for (int j = 0; j < config_.num_scan_subregions; ++j)
    {
      // ((s+d)*N+j*(e-s-2d))/N, ((s+d)*N+(j+1)*(e-s-2d))/N-1
      size_t sp = ((0 + config_.num_curvature_regions) * (config_.num_scan_subregions - j) + (scan_size - config_.num_curvature_regions) * j) / config_.num_scan_subregions;
      size_t ep = ((0 + config_.num_curvature_regions) * (config_.num_scan_subregions - 1 - j) + (scan_size - config_.num_curvature_regions) * (j + 1)) / config_.num_scan_subregions - 1;

      // skip empty regions
      if (ep <= sp)
      {
        continue;
      }

      size_t region_size = ep - sp + 1;
      PrepareSubregion(scan_ring, sp, ep);

      // extract corner features
      int num_largest_picked = 0;
      for (size_t k = region_size; k > 0 && num_largest_picked < config_.max_corner_less_sharp;)
      {
        // k must be greater than 0
        const pair<float, size_t> &curvature_idx = curvature_idx_pairs_[--k];
        float curvature = curvature_idx.first;
        size_t idx = curvature_idx.second;
        size_t in_scan_idx = idx - 0; // scan start index is 0 for all ring scans
        size_t in_region_idx = idx - sp;

        if (scan_ring_mask_[in_scan_idx] == 0 && curvature > config_.surf_curv_th / 2)
        {
          ++num_largest_picked;
          if (num_largest_picked <= config_.max_corner_sharp)
          {
            subregion_labels_[in_region_idx] = CORNER_SHARP;
            corner_points_sharp_.push_back(scan_ring[in_scan_idx]);
          }
          else
          {
            subregion_labels_[in_region_idx] = CORNER_LESS_SHARP;
          }
          corner_points_less_sharp_.push_back(scan_ring[in_scan_idx]);

          MaskPickedInRing(scan_ring, in_scan_idx);
        }
      }

      // extract flat surface features
      int num_smallest_picked = 0;
      for (int k = 0; k < region_size && num_smallest_picked < config_.max_surf_flat; ++k)
      {
        const pair<float, size_t> &curvature_idx = curvature_idx_pairs_[k];
        float curvature = curvature_idx.first;
        size_t idx = curvature_idx.second;
        size_t in_scan_idx = idx - 0; // scan start index is 0 for all ring scans
        size_t in_region_idx = idx - sp;

        if (scan_ring_mask_[in_scan_idx] == 0 && curvature < config_.surf_curv_th / 10)
        {
          ++num_smallest_picked;
          subregion_labels_[in_region_idx] = SURFACE_FLAT;
          surface_points_flat_.push_back(scan_ring[in_scan_idx]);

          MaskPickedInRing(scan_ring, in_scan_idx);
        }
      }

      // extract less flat surface features
      for (int k = 0; k < region_size; ++k)
      {
        if (subregion_labels_[k] <= SURFACE_LESS_FLAT)
        {
          surf_points_less_flat_ptr->push_back(scan_ring[sp + k]);
        }
      }

    } /// j

    // down size less flat surface point cloud of current scan
    PointCloud surf_points_less_flat_downsampled;
    pcl::ApproximateVoxelGrid<PointT> down_size_filter;
    down_size_filter.setInputCloud(surf_points_less_flat_ptr);
    down_size_filter.setLeafSize(config_.less_flat_filter_size,
                                 config_.less_flat_filter_size,
                                 config_.less_flat_filter_size);

    if (surf_points_less_flat_ptr->empty())
    {
      continue;
    }
    else
    {

      down_size_filter.filter(surf_points_less_flat_downsampled);
      surface_points_less_flat_ += surf_points_less_flat_downsampled;
    }

  } /// i

  // size_t less_flat_cloud_size = surface_points_less_flat_.size();

  // for (int i = 0; i < less_flat_cloud_size; ++i)
  // {
  //   PointT &p = surface_points_less_flat_[i];

  //   float azi_rad = 2 * M_PI - atan2(p.y, p.x);

  //   if (azi_rad >= 2 * M_PI)
  //   {
  //     azi_rad -= 2 * M_PI;
  //   }

  //   float azi_rad_rel = azi_rad - start_ori_;
  //   if (azi_rad_rel < 0)
  //   {
  //     azi_rad_rel += 2 * M_PI;
  //   }

  //   float rel_time = config_.scan_period * azi_rad_rel / (2 * M_PI);

  // //   if (rel_time > config_.scan_period || rel_time < 0.0)
  // //   {
  // //     LOG(ERROR) << "config_.scan_period: " << config_.scan_period;
  // //     LOG(ERROR) << "rel_time error: " << rel_time;
  // //   }
  //   p.intensity = int(p.intensity) + rel_time;
  // }

  // ROS_DEBUG_STREAM("extract features time: " << tic_toc_.Toc() << " ms");
  // DLOG(INFO) << "extract features time: " << tic_toc_.Toc() << " ms";

} // ExtractFeaturePoints

void PointProcessor::PointCloudHandler(const sensor_msgs::PointCloud2ConstPtr &raw_points_msg)
{

  PointCloud laser_cloud_in;
  pcl::fromROSMsg(*raw_points_msg, laser_cloud_in);
  PointCloudConstPtr laser_cloud_in_ptr(new PointCloud(laser_cloud_in));

  SetInputCloud(laser_cloud_in_ptr, raw_points_msg->header.stamp);

  Process();
}

void PointProcessor::SetupRos(ros::NodeHandle &nh)
{
  is_ros_setup_ = true;

  sub_raw_points_ = nh.subscribe<sensor_msgs::PointCloud2>("/vils_estimator/lidar_cloud", 5, &PointProcessor::PointCloudHandler, this, ros::TransportHints().tcpNoDelay());
  pub_full_cloud_ = nh.advertise<sensor_msgs::PointCloud2>("/velodyne_cloud_2", 2);
  pub_corner_points_sharp_ = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_sharp", 2);
  pub_corner_points_less_sharp_ = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_less_sharp", 2);
  pub_surf_points_flat_ = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_flat", 2);
  pub_surf_points_less_flat_ = nh.advertise<sensor_msgs::PointCloud2>("/laser_cloud_less_flat", 2);
}

void PointProcessor::Reset(const ros::Time &scan_time, const bool &is_new_sweep)
{
  scan_time_ = scan_time;

  // clear internal cloud buffers at the beginning of a sweep
  if (is_new_sweep)
  {
    sweep_start_ = scan_time_;

    // clear cloud buffers
    cloud_in_rings_.clear();
    corner_points_sharp_.clear();
    corner_points_less_sharp_.clear();
    surface_points_flat_.clear();
    surface_points_less_flat_.clear();

    // clear scan indices vector
    scan_ranges.clear();

    for (int i = 0; i < laser_scans.size(); ++i)
    {
      laser_scans[i]->clear();
    }

    for (int i = 0; i < intensity_scans.size(); ++i)
    {
      intensity_scans[i]->clear();
    }

    //    laser_scans.clear();
    //
    //    for (int i = 0; i < num_rings_; ++i) {
    //      PointCloudPtr scan(new PointCloud());
    //      laser_scans.push_back(scan);
    //    }
  }
}

void PointProcessor::SetInputCloud(const PointCloudConstPtr &cloud_in, ros::Time time_in)
{
  Reset(time_in);
  cloud_ptr_ = cloud_in;
}

void PointProcessor::PointToRing()
{

  PointToRing(cloud_ptr_, laser_scans, intensity_scans);

  // construct sorted full resolution cloud
  size_t cloud_size = 0;
  for (int i = 0; i < num_rings_; i++)
  {
    cloud_in_rings_ += (*intensity_scans[i]);

    IndexRange range(cloud_size, 0);
    cloud_size += (*laser_scans[i]).size();
    range.second = (cloud_size > 0 ? cloud_size - 1 : 0);
    scan_ranges.push_back(range);
  }

  // ROS_WARN_STREAM("point size: " << cloud_in_rings_.size());
}

void PointProcessor::PointToRing(const PointCloudConstPtr &cloud_in,
                                 vector<PointCloudPtr> &ring_out,
                                 vector<PointCloudPtr> &intensity_out)
{
  auto &points = cloud_in->points;
  size_t cloud_size = points.size();

  // tic_toc_.Tic();

  ///>
  float startOri = -atan2(points[0].y, points[0].x);
  float endOri = -atan2(points[cloud_size - 1].y,
                        points[cloud_size - 1].x) +
                 2 * M_PI;

  if (endOri - startOri > 3 * M_PI)
  {
    endOri -= 2 * M_PI;
  }
  else if (endOri - startOri < M_PI)
  {
    endOri += 2 * M_PI;
  }
  bool halfPassed = false;
  int count = cloud_size;
  bool start_flag = false;
  ///<

  for (size_t i = 0; i < cloud_size; ++i)
  {
    PointT p = points[i];
    PointT p_with_intensity = points[i];

    // skip NaN and INF valued points
    if (!pcl_isfinite(p.x) ||
        !pcl_isfinite(p.y) ||
        !pcl_isfinite(p.z))
    {
      continue;
    }

    float dis = sqrt(p.x * p.x + p.y * p.y);
    float ele_rad = atan2(p.z, dis);
    float azi_rad = 2 * M_PI - atan2(p.y, p.x);

    ///>mine
#ifndef DEBUG_ORIGIN
    if (azi_rad >= 2 * M_PI)
    {
      azi_rad -= 2 * M_PI;
    }

    int scan_id = ElevationToRing(ele_rad);

    if (scan_id >= num_rings_ || scan_id < 0)
    {
      continue;
    }

    if (!start_flag)
    {
      start_ori_ = azi_rad;
      start_flag = true;
    }

    p.intensity = points[i].intensity;
    p_with_intensity.intensity = points[i].intensity;
#endif
    ///<mine

    ///>origin
#ifdef DEBUG_ORIGIN
    int scan_id = ElevationToRing(ele_rad);

    if (scan_id >= num_rings_ || scan_id < 0)
    {
      // DLOG(INFO) << RadToDeg(ele_rad) << ", " << scan_id;
      // DLOG(INFO) << (scan_id < 0 ? " point too low" : "point too high");
      continue;
    }

    float ori = -atan2(p.y, p.x);
    if (!halfPassed)
    {
      if (ori < startOri - M_PI / 2)
      {
        ori += 2 * M_PI;
      }
      else if (ori > startOri + M_PI * 3 / 2)
      {
        ori -= 2 * M_PI;
      }

      if (ori - startOri > M_PI)
      {
        halfPassed = true;
      }
    }
    else
    {
      ori += 2 * M_PI;

      if (ori < endOri - M_PI * 3 / 2)
      {
        ori += 2 * M_PI;
      }
      else if (ori > endOri + M_PI / 2)
      {
        ori -= 2 * M_PI;
      }
    }

    float relTime = (ori - startOri) / (endOri - startOri);
#ifndef DEBUG
    p.intensity = scan_id + config_.scan_period * relTime;
#else
    p.intensity = config_.scan_period * relTime;
#endif
#endif
    ///<origin

    ring_out[scan_id]->push_back(p);
    intensity_out[scan_id]->push_back(p_with_intensity);
  }

  // ROS_DEBUG_STREAM("ring time: " << tic_toc_.Toc() << " ms");

  // tic_toc_.Tic();
  // TODO: De-skew with rel_time, this for loop is not necessary
  // start_ori_ = NAN;

  //  for (int ring = 0; ring < num_rings_; ++ring) {
  //    if (ring_out[ring]->size() <= 0) {
  //      continue;
  //    }
  //
  //    float azi_rad = ring_out[ring]->front().intensity;
  //    if (start_ori_ != start_ori_) {
  //      start_ori_ = azi_rad;
  //    } else {
  //      start_ori_ = (RadLt(start_ori_, azi_rad) ? start_ori_ : azi_rad);
  //    }
  //    // cout << azi_rad << endl;
  //  }
  // start_ori_ = ring_out[0]->front().intensity;

  // infer right start_ori_
  //   if (config_.infer_start_ori_)
  //   {
  //     std_msgs::Float32 start_ori_msg, start_ori_inferred_msg;
  //     start_ori_msg.data = start_ori_;

  //     start_ori_buf2_.push(start_ori_);
  //     if (start_ori_buf1_.size() >= 10)
  //     {
  //       float start_ori_diff1 = start_ori_buf1_.last() - start_ori_buf1_.first();
  //       float start_ori_step1 = NormalizeRad(start_ori_diff1) / 9;

  //       float start_ori_diff2 = start_ori_buf2_.last() - start_ori_buf2_.first();
  //       float start_ori_step2 = NormalizeRad(start_ori_diff2) / 9;

  //       if (fabs(NormalizeRad(start_ori_ - start_ori_buf1_.last())) > config_.rad_diff)
  //       {
  //         start_ori_ = start_ori_buf1_.last() + start_ori_step1;
  //         start_ori_ = NormalizeRad(start_ori_);
  //         if (start_ori_ < 0)
  //         {
  //           start_ori_ += 2 * M_PI;
  //         }
  //       }
  //       if (AbsRadDistance(start_ori_step1, start_ori_step2) < 0.05 && AbsRadDistance(start_ori_buf2_[9] - start_ori_buf2_[8], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[8] - start_ori_buf2_[7], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[7] - start_ori_buf2_[6], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[6] - start_ori_buf2_[5], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[5] - start_ori_buf2_[4], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[4] - start_ori_buf2_[3], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[3] - start_ori_buf2_[2], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[2] - start_ori_buf2_[1], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[1] - start_ori_buf2_[0], start_ori_step1) < 0.05)
  //       {
  //         start_ori_ = ring_out[0]->front().intensity;
  //       }
  //     }
  //     start_ori_buf1_.push(start_ori_);

  //     start_ori_inferred_msg.data = start_ori_;
  //   } // if

  //   for (int ring = 0; ring < num_rings_; ++ring)
  //   {
  //     // points in a ring
  //     PointCloud &points_in_ring = (*ring_out[ring]);
  //     PointCloud &points_in_ring_with_intensity = (*intensity_out[ring]);
  //     size_t cloud_in_ring_size = points_in_ring.size();

  //     for (int i = 0; i < cloud_in_ring_size; ++i)
  //     {
  //       PointT &p = points_in_ring[i];
  //       PointT &p_with_intensity = points_in_ring_with_intensity[i];

  //       float azi_rad_rel = p.intensity - start_ori_;
  //       if (azi_rad_rel < 0)
  //       {
  //         azi_rad_rel += 2 * M_PI;
  //       }

  //       float rel_time = config_.scan_period * azi_rad_rel / (2 * M_PI);

  //       ///>mine
  // #ifndef DEBUG_ORIGIN

  // #ifndef DEBUG
  //       p.intensity = ring + rel_time;
  //       p_with_intensity.intensity = int(p_with_intensity.intensity) + rel_time;
  // #else
  //       p.intensity = rel_time;
  // #endif

  // #endif
  //       ///<mine
  //     } // for i
  //   }   // for ring
  // ROS_DEBUG_STREAM("reorder time: " << tic_toc_.Toc() << " ms");
}

void PointProcessor::PrepareRing(const PointCloud &scan)
{
  // const PointCloud &scan = laser_scans[idx_ring];
  size_t scan_size = scan.size();
  scan_ring_mask_.resize(scan_size);
  scan_ring_mask_.assign(scan_size, 0);
  for (size_t i = 0 + config_.num_curvature_regions; i < scan_size - config_.num_curvature_regions; ++i)
  {
    const PointT &p_prev = scan[i - 1];
    const PointT &p_curr = scan[i];
    const PointT &p_next = scan[i + 1];

    float diff_next2 = CalcSquaredDiff(p_curr, p_next);

    // about 30 cm
    if (diff_next2 > 0.1)
    {
      float depth = CalcPointDistance(p_curr);
      float depth_next = CalcPointDistance(p_next);

      if (depth > depth_next)
      {
        // to closer point
        float weighted_diff = sqrt(CalcSquaredDiff(p_next, p_curr, depth_next / depth)) / depth_next;
        // relative distance
        if (weighted_diff < 0.1)
        {
          fill_n(&scan_ring_mask_[i - 0 - config_.num_curvature_regions], config_.num_curvature_regions + 1, 1);
          continue;
        }
      }
      else
      {
        float weighted_diff = sqrt(CalcSquaredDiff(p_curr, p_next, depth / depth_next)) / depth;
        if (weighted_diff < 0.1)
        {
          fill_n(&scan_ring_mask_[i - 0 + 1], config_.num_curvature_regions + 1, 1);
          continue;
        }
      }
    }

    float diff_prev2 = CalcSquaredDiff(p_curr, p_prev);
    float dis2 = CalcSquaredPointDistance(p_curr);

    // for this point -- 1m -- 1.5cm
    if (diff_next2 > 0.0002 * dis2 && diff_prev2 > 0.0002 * dis2)
    {
      scan_ring_mask_[i - 0] = 1;
    }
  }
}

void PointProcessor::PrepareSubregion(const PointCloud &scan, const size_t idx_start, const size_t idx_end)
{

  //  cout << ">>>>>>> " << idx_ring << ", " << idx_start << ", " << idx_end << " <<<<<<<" << endl;
  //  const PointCloud &scan = laser_scans[idx_ring];
  size_t region_size = idx_end - idx_start + 1;
  curvature_idx_pairs_.resize(region_size);
  subregion_labels_.resize(region_size);
  subregion_labels_.assign(region_size, SURFACE_LESS_FLAT);

  int num_point_neighbors = 2 * config_.num_curvature_regions;

  for (size_t i = idx_start, in_region_idx = 0; i <= idx_end; ++i, ++in_region_idx)
  {
    float diff_x = -num_point_neighbors * scan[i].x;
    float diff_y = -num_point_neighbors * scan[i].y;
    float diff_z = -num_point_neighbors * scan[i].z;
    int diff_i{0};
    float diff_i_ = scan[i].intensity;

    for (int j = 1; j <= config_.num_curvature_regions; j++)
    {
      diff_x += scan[i + j].x + scan[i - j].x;
      diff_y += scan[i + j].y + scan[i - j].y;
      diff_z += scan[i + j].z + scan[i - j].z;

      float diff_a = scan[i + j].intensity / diff_i_;
      if (int(diff_a) == 1)
        diff_i = diff_i + 1;
      float diff_b = scan[i - j].intensity / diff_i_;
      if (int(diff_b) == 1)
        diff_i = diff_i + 1;
    }

    float curvature = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
    // float curvature_wi = diff_i /diff_i_ ;

    pair<float, size_t> curvature_idx_(curvature, i);
    curvature_idx_pairs_[in_region_idx] = curvature_idx_;

    if (curvature > config_.surf_curv_th / 2 && diff_i > 4)
    {
      scan_ring_mask_[i] = 1;
    }
    else if (curvature < config_.surf_curv_th / 10 && diff_i < 5)
    {
      scan_ring_mask_[i] = 1;
    }

    //    _regionCurvature[regionIdx] = diffX * diffX + diffY * diffY + diffZ * diffZ;
    //    _regionSortIndices[regionIdx] = i;
  }

  sort(curvature_idx_pairs_.begin(), curvature_idx_pairs_.end());

  //  for (const auto &pair : curvature_idx_pairs_) {
  //    cout << pair.first << " " << pair.second << endl;
  //  }
}

void PointProcessor::MaskPickedInRing(const PointCloud &scan, const size_t in_scan_idx)
{

  // const PointCloud &scan = laser_scans[idx_ring];
  scan_ring_mask_[in_scan_idx] = 1;

  for (int i = 1; i <= config_.num_curvature_regions; ++i)
  {
    /// 20cm
    if (CalcSquaredDiff(scan[in_scan_idx + i], scan[in_scan_idx + i - 1]) > 0.05)
    {
      break;
    }

    scan_ring_mask_[in_scan_idx + i] = 1;
  }

  for (int i = 1; i <= config_.num_curvature_regions; ++i)
  {
    if (CalcSquaredDiff(scan[in_scan_idx - i], scan[in_scan_idx - i + 1]) > 0.05)
    {
      break;
    }

    scan_ring_mask_[in_scan_idx - i] = 1;
  }
}

void PointProcessor::PublishResults()
{

  if (!is_ros_setup_)
  {
    return;
  }
  PublishCloudMsg(pub_full_cloud_, cloud_in_rings_, sweep_start_, "/laser_link");
  PublishCloudMsg(pub_corner_points_sharp_, corner_points_sharp_, sweep_start_, "/laser_link");
  PublishCloudMsg(pub_corner_points_less_sharp_, corner_points_less_sharp_, sweep_start_, "/laser_link");
  PublishCloudMsg(pub_surf_points_flat_, surface_points_flat_, sweep_start_, "/laser_link");
  PublishCloudMsg(pub_surf_points_less_flat_, surface_points_less_flat_, sweep_start_, "/laser_link");
}

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

  ros::AsyncSpinner spinner(1);
  spinner.start();
  ros::waitForShutdown();

  return 0;
}
