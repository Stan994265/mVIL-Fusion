#include "PointProcessor.h"
#include <pcl/filters/approximate_voxel_grid.h>

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
  PublishResults();
}

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
  std::string lidar_topic;
  nh.param<std::string>("lidar_topic", lidar_topic, "/velodyne_points");
  // subscribe to raw cloud topic  /lslidar_point_cloud  /velodyne_points  /os1_cloud_node1/points /hesai/pandar
  sub_raw_points_ = nh.subscribe<sensor_msgs::PointCloud2>(lidar_topic, 20, &PointProcessor::PointCloudHandler, this, ros::TransportHints().tcpNoDelay());
  // advertise scan registration topics
  pub_full_cloud_ = nh.advertise<sensor_msgs::PointCloud2>("/full_cloud", 20);
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

    // IndexRange range(cloud_size, 0);
    // cloud_size += (*laser_scans[i]).size();
    // range.second = (cloud_size > 0 ? cloud_size - 1 : 0);
    // scan_ranges.push_back(range);
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

    p.intensity = azi_rad;
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
  if (config_.infer_start_ori_)
  {
    std_msgs::Float32 start_ori_msg, start_ori_inferred_msg;
    start_ori_msg.data = start_ori_;

    start_ori_buf2_.push(start_ori_);
    if (start_ori_buf1_.size() >= 10)
    {
      float start_ori_diff1 = start_ori_buf1_.last() - start_ori_buf1_.first();
      float start_ori_step1 = NormalizeRad(start_ori_diff1) / 9;

      float start_ori_diff2 = start_ori_buf2_.last() - start_ori_buf2_.first();
      float start_ori_step2 = NormalizeRad(start_ori_diff2) / 9;

      if (fabs(NormalizeRad(start_ori_ - start_ori_buf1_.last())) > config_.rad_diff)
      {
        start_ori_ = start_ori_buf1_.last() + start_ori_step1;
        start_ori_ = NormalizeRad(start_ori_);
        if (start_ori_ < 0)
        {
          start_ori_ += 2 * M_PI;
        }
      }
      if (AbsRadDistance(start_ori_step1, start_ori_step2) < 0.05 && AbsRadDistance(start_ori_buf2_[9] - start_ori_buf2_[8], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[8] - start_ori_buf2_[7], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[7] - start_ori_buf2_[6], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[6] - start_ori_buf2_[5], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[5] - start_ori_buf2_[4], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[4] - start_ori_buf2_[3], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[3] - start_ori_buf2_[2], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[2] - start_ori_buf2_[1], start_ori_step1) < 0.05 && AbsRadDistance(start_ori_buf2_[1] - start_ori_buf2_[0], start_ori_step1) < 0.05)
      {
        start_ori_ = ring_out[0]->front().intensity;
      }
    }
    start_ori_buf1_.push(start_ori_);

    start_ori_inferred_msg.data = start_ori_;
  } // if

  for (int ring = 0; ring < num_rings_; ++ring)
  {
    // points in a ring
    PointCloud &points_in_ring = (*ring_out[ring]);
    PointCloud &points_in_ring_with_intensity = (*intensity_out[ring]);
    size_t cloud_in_ring_size = points_in_ring.size();

    for (int i = 0; i < cloud_in_ring_size; ++i)
    {
      PointT &p = points_in_ring[i];
      PointT &p_with_intensity = points_in_ring_with_intensity[i];

      float azi_rad_rel = p.intensity - start_ori_;
      if (azi_rad_rel < 0)
      {
        azi_rad_rel += 2 * M_PI;
      }

      float rel_time = config_.scan_period * azi_rad_rel / (2 * M_PI);

      ///>mine
#ifndef DEBUG_ORIGIN

#ifndef DEBUG
      p.intensity = ring + rel_time;
      p_with_intensity.intensity = int(p_with_intensity.intensity) + rel_time;
#else
      p.intensity = rel_time;
#endif

#endif
      ///<mine
    } // for i
  }   // for ring
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

    for (int j = 1; j <= config_.num_curvature_regions; j++)
    {
      diff_x += scan[i + j].x + scan[i - j].x;
      diff_y += scan[i + j].y + scan[i - j].y;
      diff_z += scan[i + j].z + scan[i - j].z;
    }

    float curvature = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
    pair<float, size_t> curvature_idx_(curvature, i);
    curvature_idx_pairs_[in_region_idx] = curvature_idx_;
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
}
