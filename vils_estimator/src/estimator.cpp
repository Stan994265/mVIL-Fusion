#include "estimator.h"
// #define USE_ES

Estimator::Estimator() : f_manager{Rs}
{
    ROS_INFO("init begins");
    clearState();
}

void Estimator::setParameter()
{
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
    }
    f_manager.setRic(ric);
    ProjectionFactor::sqrt_info = FOCAL_LENGTH / 2.0 * Matrix2d::Identity(); //1.5
    ProjectionTdFactor::sqrt_info = FOCAL_LENGTH / 2.0 * Matrix2d::Identity();
    td = TD;
}

void Estimator::clearState()
{
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr)
            delete pre_integrations[i];
        pre_integrations[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }

    for (auto &it : all_image_frame)
    {
        if (it.second.pre_integration != nullptr)
        {
            delete it.second.pre_integration;
            it.second.pre_integration = nullptr;
        }
    }

    solver_flag = INITIAL;
    first_imu = false,
    sum_of_back = 0;
    sum_of_front = 0;
    frame_count = 0;
    solver_flag = INITIAL;
    initial_timestamp = 0;
    all_image_frame.clear();
    td = TD;

    if (tmp_pre_integration != nullptr)
        delete tmp_pre_integration;
    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();

    f_manager.clearState();

    failure_occur = 0;
    relocalization_info = 0;

    drift_correct_r = Matrix3d::Identity();
    drift_correct_t = Vector3d::Zero();

    lidarCalibration.Reset();
}

void Estimator::processIMU(double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
    if (!first_imu)
    {
        first_imu = true;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }

    if (!pre_integrations[frame_count])
    {
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
    }
    if (frame_count != 0)
    {
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
        //if(solver_flag != NON_LINEAR)
        tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);

        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);

        int j = frame_count;
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity;
}

void Estimator::processLidar(CloudData &cloud, const double time_)
{
    current_lidar.lidarData.point_cloud = cloud;
    current_lidar.frameID = lidar_count_;

    if (solver_flag != INITIAL && ESTIMATE_EXTRINSIC != 2) //ESTIMATE_EXTRINSIC != 2 && lidar_count % 2 == 0
    {
        double tm = cloud.time; //lidar time
        // double ti = time;                                     //image time
        current_lidar.keylidar = false;
        int idl, idr;
        map<double, ImageFrame>::iterator iteri;
        map<double, ImageFrame>::iterator iterj;
        if (FindNearest2ID(Headers, tm, idl, idr))
        {
            iterj = all_image_frame.find(Headers[idr].stamp.toSec());
            iteri = all_image_frame.find(Headers[idl].stamp.toSec());
            if (iterj != all_image_frame.end() && iteri != all_image_frame.end())
            {
                cout << "FIND ID!!!!!!!!!!!:" << idl << " " << idr << endl;
                current_lidar.keylidar = true;
            }

            if (current_lidar.keylidar) //ti >= tm
            {
                current_lidar.lidarData.next_image = iterj->second;
                current_lidar.vioData.tj = iterj->first + time_;
                current_lidar.lidarData.last_image = iteri->second;
                current_lidar.vioData.ti = iteri->first + time_;
                current_lidar.vioData.dt = time_;

                current_lidar.vioData.acci = current_lidar.lidarData.last_image.pre_integration->acc_1 - Bas[idl];
                current_lidar.vioData.gyri = current_lidar.lidarData.last_image.pre_integration->gyr_1 - Bgs[idl];
                current_lidar.vioData.accj = current_lidar.lidarData.next_image.pre_integration->acc_1 - Bas[idr];
                current_lidar.vioData.gyrj = current_lidar.lidarData.next_image.pre_integration->gyr_1 - Bgs[idr];

                current_lidar.vioData.Rwbi = Rs[idl];
                current_lidar.vioData.Pwbi = Ps[idl];
                current_lidar.vioData.Vbi = Vs[idl];

                current_lidar.vioData.Rwbj = Rs[idr];
                current_lidar.vioData.Pwbj = Ps[idr];
                current_lidar.vioData.Vbj = Vs[idr];

                // step1: distortion_adjust
                CloudData::CLOUD_PTR lidar_pre(new CloudData::CLOUD);
                lidar_pre = current_lidar.lidarData.point_cloud.cloud_ptr;
                // step2: downsampling
                // pcl::ApproximateVoxelGrid<CloudData::POINT> voxelgrid;
                // voxelgrid.setLeafSize(LeafSize, LeafSize, LeafSize);
                // voxelgrid.setInputCloud(lidar_pre);
                // voxelgrid.filter(*lidar_pre);

                if (lidarCalibration.lidar_init_flag == false)
                {
                    //VERS1
                    // std::shared_ptr<VelocityData> current_velocity_data;
                    // current_velocity_data = std::make_shared<VelocityData>();
                    // current_velocity_data->PredictData(current_lidar);
                    // current_velocity_data->TransformCoordinate(RLB, TLB); //RLB.transpose(), TBL  RLB, TLB
                    // std::shared_ptr<DistortionAdjust> distortion_adjust_ptr_;
                    // distortion_adjust_ptr_ = std::make_shared<DistortionAdjust>();
                    // distortion_adjust_ptr_->SetMotionInfo(LidarTimeStep, *current_velocity_data);
                    // CloudData::CLOUD_PTR out_cloud(new CloudData::CLOUD);
                    // distortion_adjust_ptr_->AdjustCloud(lidar_pre, out_cloud);
                    // lidar_pre = out_cloud;

                    //VERS2
                    float time_factor = 1 / LidarTimeStep;
                    double ta, tb, tl, tls, tle;
                    ta = current_lidar.vioData.ti;
                    tb = current_lidar.vioData.tj;
                    tl = tm;
                    // lidar time on mid
                    tls = tl - 0.5 * LidarTimeStep;
                    tle = tl + 0.5 * LidarTimeStep;
                    // lidar time on start
                    // tls = tl;
                    // tle = tl + LidarTimeStep;

                    Quaterniond temQa(current_lidar.vioData.Rwbi);
                    Quaterniond temQb(current_lidar.vioData.Rwbj);

                    // use imu data
                    // double dts = tls - ta;
                    // double dte = tb - tle;
                    // Vector3d wa = current_lidar.vioData.gyri;
                    // Vector3d wb = current_lidar.vioData.gyrj;
                    // Quaterniond Qs{1.0, 0.5 * wa[0] * dts, 0.5 * wa[1] * dts, 0.5 * wa[2] * dts};
                    // Quaterniond Qe{1.0, -0.5 * wb[0] * dte, -0.5 * wb[1] * dte, -0.5 * wb[2] * dte};
                    // Quaterniond temQls = temQa * Qs;
                    // Quaterniond temQle = temQb * Qe;

                    // use vio data(better when zero velocity)
                    double ss = (tls - ta) / (tb - ta);
                    double se = (tle - ta) / (tb - ta);
                    Quaterniond temQls = temQa.slerp(ss, temQb);
                    Quaterniond temQle = temQa.slerp(se, temQb);

                    Matrix4d trans_b(Matrix4d::Identity()); //Tes
                    trans_b.block<3, 3>(0, 0) = (temQle.inverse() * temQls).toRotationMatrix();
                    trans_b.block<3, 1>(0, 3) = (-current_lidar.vioData.Rwbj.transpose() * current_lidar.vioData.Pwbj + current_lidar.vioData.Rwbj.transpose() * current_lidar.vioData.Pwbi) * LidarTimeStep / (tb - ta);

                    Matrix4d trans_lb(Matrix4d::Identity());
                    trans_lb.block<3, 3>(0, 0) = RLB;
                    trans_lb.block<3, 1>(0, 3) = TLB;

                    Matrix4d trans_l = trans_lb * trans_b * trans_lb.inverse();
                    Matrix3f trans_l_r = trans_l.block<3, 3>(0, 0).cast<float>();
                    Eigen::Quaternionf transform_es_q(trans_l_r);
                    Eigen::Vector3f transform_es_t(trans_l.block<3, 1>(0, 3).cast<float>());
                    TransformToEnd(lidar_pre, transform_es_q, transform_es_t, time_factor, MinDistance, MaxDistance);
                    lidar_pre->is_dense = false;
                    std::vector<int> indices;
                    pcl::removeNaNFromPointCloud(*lidar_pre, *lidar_pre, indices);
                    *cloud.cloud_ptr = *lidar_pre;
                }

                // step2: downsampling
                CloudData::CLOUD_PTR copy(new CloudData::CLOUD);
                copyPointCloud(*lidar_pre, *copy);
                pcl::ApproximateVoxelGrid<CloudData::POINT> voxelgrid;
                voxelgrid.setLeafSize(LeafSize, LeafSize, LeafSize);
                voxelgrid.setInputCloud(copy);
                voxelgrid.filter(*copy);
                current_lidar.lidarData.point_cloud.cloud_ptr = copy;

                // step3: fast-gicp
                all_lidar_frame.insert(make_pair(tm, current_lidar));
                map<double, LidarFrame>::iterator lframej;
                map<double, LidarFrame>::iterator lframei;

                if (all_lidar_frame.size() > 1)
                { //&& solver_flag == INITIAL
                    if (all_lidar_frame.size() > 2)
                    {
                        all_lidar_frame.erase(all_lidar_frame.begin());
                    }
                    lframej = all_lidar_frame.end();
                    lframej--;
                    lframei = prev(lframej);
                    CloudData::CLOUD_PTR target_cloud(new CloudData::CLOUD);
                    CloudData::CLOUD_PTR source_cloud(new CloudData::CLOUD);
                    *target_cloud = *lframei->second.lidarData.point_cloud.cloud_ptr;
                    *source_cloud = *lframej->second.lidarData.point_cloud.cloud_ptr;

                    // cout << "------------------------------------ fgicp_mt ---------------------------------" << std::endl;
                    boost::shared_ptr<fast_gicp::FastVGICP<CloudData::POINT, CloudData::POINT> > gicp(new fast_gicp::FastVGICP<CloudData::POINT, CloudData::POINT>());
                    gicp->setResolution(0.5);
                    // boost::shared_ptr<fast_gicp::FastGICP<CloudData::POINT, CloudData::POINT> > gicp(new fast_gicp::FastGICP<CloudData::POINT, CloudData::POINT>());
                    gicp->setNumThreads(NumThreads);
                    // gicp->setTransformationEpsilon(TransformationEpsilon);
                    // gicp->setMaxCorrespondenceDistance(MaxCorrespondenceDistance);

                    CloudData::CLOUD_PTR aligned(new CloudData::CLOUD);
                    gicp->setInputSource(source_cloud);
                    gicp->setInputTarget(target_cloud);
                    Vector3d init_guss_;
                    Matrix4d init_guss(Matrix4d::Identity());
                    Matrix3d tem_r;
                    Vector3d tem_t;
                    if (lidarCalibration.lidar_init_flag == false)
                    {
                        Matrix4d TRANSFORM_Lij = PredictRelative_rt(lframei->second, lframej->second, RLB.transpose(), TBL);
                        tem_r = TRANSFORM_Lij.block<3, 3>(0, 0);
                        tem_t = TRANSFORM_Lij.block<3, 1>(0, 3);
                        init_guss.block<3, 3>(0, 0) = RLB * tem_r * RLB.transpose();
                        init_guss.block<3, 1>(0, 3) = RLB * tem_r * TBL + TLB +
                                                      RLB * tem_t;
                        Matrix<float, 4, 4> guss = init_guss.cast<float>();
                        gicp->align(*aligned, guss);
                        init_guss_ = init_guss.block<3, 1>(0, 3);
                        current_lidar.lidarData.lidar_R = lframej->second.lidarData.lidar_R;
                        current_lidar.lidarData.lidar_T = lframej->second.lidarData.lidar_T;
                    }
                    else
                    {
                        gicp->align(*aligned);
                    }

                    // step4: evaluate
                    double fitness_core = gicp->getFitnessScore();
                    std::cout << "Fitness core=" << fitness_core << std::endl;

                    // step5: icp_cov(not necessary)
                    // Eigen::Matrix4f transform_m{gicp->getFinalTransformation()};
                    // Eigen::MatrixXd ICP_COV;
                    // calculate_ICP_COV(*target_cloud, *source_cloud, transform_m, ICP_COV);
                    // double trace_l{ICP_COV.trace()};
                    // cout << "ICP_COV TRACE=" << fixed << trace_l << endl;

                    Eigen::Matrix4d transformation_matrix{gicp->getFinalTransformation().cast<double>()};

                    // step6: use viodata-R for lidar(if lidar data not good), otherwise ignore this if
                    // if (lidarCalibration.lidar_init_flag == false)
                    // {
                    //     transformation_matrix.block<3, 3>(0, 0) = init_guss.block<3, 3>(0, 0);
                    // }

                    // step7: add LiDAR constrain to VIO
                    Matrix3d Rij(transformation_matrix.block<3, 3>(0, 0));
                    Vector3d Tij(transformation_matrix.block<3, 1>(0, 3));
                    if (lidarCalibration.lidar_init_flag == false) //&& trace_l > 1e-7
                    {
                        double tem_T;
                        Vector3d tem_Tij{init_guss_ - Tij};
                        tem_T = fabs(tem_Tij[0]) + fabs(tem_Tij[1]) + fabs(tem_Tij[2]);
                        cout << "DIFF T=" << tem_T << endl;
                        cout << "Ltrans =" << Tij.transpose() << endl;
                        // cout << "guess-R =" << Utility::R2ypr(init_guss.block<3, 3>(0, 0)).transpose() << endl;
                        // cout << "l-R =" << Utility::R2ypr(Rij).transpose() << endl;
                        // cout << "l-R-t =" << Utility::R2ypr(Rij.transpose()).transpose() << endl;
                        LidarICPConstraint LidarICPConstraint;
                        LidarICPConstraint.constraint_mode = 0;

                        double temYaw = Utility::R2ypr(init_guss.block<3, 3>(0, 0)).x();
                        ROS_WARN_STREAM("temYaw: " << temYaw << endl);

                        if (fitness_core < 1.0 && tem_T > 0.1) //tem_T > 0.15
                        {
                            LidarICPConstraint.constraint_mode = 3;
                        }

                        else if (fitness_core < 1.0 && tem_T <= 0.1) //tem_T > 0.15
                        {
                            LidarICPConstraint.constraint_mode = 2;
                        }

                        else if (fitness_core > 1.0) //tem_T > 0.15
                        {
                            LidarICPConstraint.constraint_mode = 1;
                        }

                        if (fabs(Tij[0]) + fabs(Tij[1]) + fabs(Tij[2]) < 0.01)
                        {
                            // double temYaw = Utility::R2ypr(init_guss.block<3, 3>(0, 0)).x();
                            // ROS_WARN_STREAM("temYaw: " << temYaw << endl);
                            if (fabs(temYaw) < 0.5)
                            {
                                LidarICPConstraint.constraint_mode = 4;
                                ROS_WARN("zero velocity!!!!!!");
                                // cout << "zero velocity!!!!!!" << endl;
                            }
                            else
                            {
                                LidarICPConstraint.constraint_mode = 5;
                                ROS_WARN("pure rotation!!!!!!");
                                // cout << "pure rotation!!!!!!" << endl;
                            }
                        }
                        current_lidar.lidarData.mode = LidarICPConstraint.constraint_mode;
                        if (!ADD_LIDAR_ICP)
                        {
                            LidarICPConstraint.constraint_mode = 0;
                        }

                        LidarICPConstraint.lidar_ta = lframei->second.lidarData.last_image.t;
                        LidarICPConstraint.lidar_tb = lframei->second.lidarData.next_image.t;
                        LidarICPConstraint.lidar_tc = lframej->second.lidarData.last_image.t;
                        LidarICPConstraint.lidar_td = lframej->second.lidarData.next_image.t;
                        LidarICPConstraint.lidar_ti = lframei->first;
                        LidarICPConstraint.lidar_tj = lframej->first;
                        // lidar_trans = transformation_matrix;
                        // lidar_cov = ICP_COV;
                        // lidar_sqrt_info = Eigen::LLT<Eigen::Matrix<double, 6, 6> >(lidar_cov.inverse()).matrixL().transpose();
                        Eigen::Matrix<double, 6, 6> tem_cov;
                        if (LidarICPConstraint.constraint_mode == 4)
                        {
                            LidarICPConstraint.lidar_trans.setIdentity();
                            LidarICPConstraint.lidar_sqrt_info = 1e12 * (tem_cov.setIdentity());
                            if (first_zv)
                            {

                                tem_zv_r = current_lidar.lidarData.lidar_R = lframei->second.lidarData.lidar_R;
                                tem_zv_t = current_lidar.lidarData.lidar_T = lframei->second.lidarData.lidar_T;
                                ROS_WARN_STREAM("First time zero!!!!!! " << tem_zv_t.transpose() << endl);
                                first_zv = false;
                                // LidarICPConstraints.push_back(LidarICPConstraint);
                                while (LidarICPConstraints.size() > 1)
                                {
                                    LidarICPConstraints.pop_front();
                                }
                            }
                            else
                            {
                                ROS_WARN_STREAM("Not First time zero!!!!!! " << tem_zv_t.transpose() << endl);
                                current_lidar.lidarData.lidar_R = tem_zv_r;
                                current_lidar.lidarData.lidar_T = tem_zv_t;
                            }
                            // LidarICPConstraints.push_back(LidarICPConstraint);
                        }
                        else if (LidarICPConstraint.constraint_mode == 3) //|| LidarICPConstraint.constraint_mode == 2
                        {
                            LidarICPConstraint.lidar_trans = EX_LB.inverse() * transformation_matrix * EX_LB;
                            // LidarICPConstraint.lidar_trans = transformation_matrix;
                            LidarICPConstraint.lidar_sqrt_info = 1 / fitness_core * 100 * (tem_cov.setIdentity());
                            // LidarICPConstraint.lidar_sqrt_info =  100 * (tem_cov.setIdentity());
                            LidarICPConstraint.lidar_sqrt_info(3, 3) = 500;
                            LidarICPConstraint.lidar_sqrt_info(4, 4) = 500;
                            LidarICPConstraint.lidar_sqrt_info(5, 5) = 500;

                            if (first_zv == false && LidarICPConstraints.size() == 1)
                            {
                                ROS_WARN_STREAM("Start Moving again!!!!!! " << tem_zv_t.transpose() << endl);
                                LidarICPConstraints.pop_front();
                                first_zv = true;
                            }
                            // current_lidar.lidarData.lidar_R = lframej->second.lidarData.lidar_R;
                            // current_lidar.lidarData.lidar_T = lframej->second.lidarData.lidar_T;

                            // LidarICPConstraints.push_back(LidarICPConstraint);
                        }

                        LidarICPConstraints.push_back(LidarICPConstraint);
                    }

                    // step8: LI initialization
                    if (solver_flag != INITIAL && lidarCalibration.lidar_init_flag == true && current_lidar.frameID > 15) //&& current_lidar.frameID > 20
                    {
#ifdef USE_ES
                        Matrix3d calib_rlb;
                        Vector3d calib_tlb;
                        if (lidarCalibration.CalibrationLidarExRotation(lframei->second, lframej->second, fitness_core,
                                                                        Rij, calib_rlb, Tij, calib_tlb))
                        {
#endif
                            // use gt
                            RLB = RLI;
                            TLB = TLI;
                            TBL = -RLB.transpose() * TLB;
#ifdef USE_ES
                            //use estimate
                            RLB = calib_rlb;
                            TLB = calib_tlb;
                            TBL = -RLB.transpose() * TLB;
#endif
                            EX_LB.block<3, 3>(0, 0) = RLB;
                            EX_LB.block<3, 1>(0, 3) = TLB;
                            Matrix3d TEM_RLC = RLB * ric[0];
                            Vector3d TEM_TLC = RLB * tic[0] + TLB;

                            ROS_WARN("initial extrinsic calib success");
                            ROS_WARN_STREAM("frame ID: " << current_lidar.frameID << endl);

                            ROS_WARN_STREAM("lidar initial extrinsic rotation: " << endl
                                                                                 << Utility::R2ypr(RLB).transpose());
                            ROS_WARN_STREAM("lidar initial extrinsic translation: " << endl
                                                                                    << TLB.transpose());
                            ROS_WARN_STREAM("cam initial extrinsic rotation: " << endl
                                                                               << Utility::R2ypr(ric[0]).transpose());
                            ROS_WARN_STREAM("cam initial extrinsic translation: " << endl
                                                                                  << tic[0].transpose());

                            ROS_WARN_STREAM("cam lidar extrinsic rotation: " << endl
                                                                             << TEM_RLC);
                            ROS_WARN_STREAM("cam lidar extrinsic translation: " << endl
                                                                                << TEM_TLC.transpose());

                            lidarCalibration.lidar_init_flag = false;
                            Matrix3d LidarInit_R = Predict_r(current_lidar);
                            Vector3d LidarInit_t = Predict_t(current_lidar);

                            LIDAR_INIT_t = LidarInit_t + LidarInit_R * TBL;
                            LIDAR_INIT_R = Quaterniond(LidarInit_R * RLB.transpose());
                            current_lidar.lidarData.lidar_T = LIDAR_INIT_t;
                            current_lidar.lidarData.lidar_R = LIDAR_INIT_R;
#ifdef USE_ES
                        }
                        else
                        {
                            RLB = calib_rlb;
                            TLB = calib_tlb;
                        }
#endif
                    }
                }
                lidar_count_++;
            }
        }
    }

    lidar_count++;
}

void Estimator::processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 8, 1> > > > &image, const std_msgs::Header &header)
{

    ROS_DEBUG("new image coming ------------------------------------------");
    ROS_DEBUG("Adding feature points %lu", image.size());

    if (f_manager.addFeatureCheckParallax(frame_count, image, td))
        marginalization_flag = MARGIN_OLD;
    else
        marginalization_flag = MARGIN_SECOND_NEW;

    ROS_DEBUG("this frame is--------------------%s", marginalization_flag ? "reject" : "accept");
    ROS_DEBUG("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_DEBUG("Solving %d", frame_count);
    ROS_DEBUG("number of feature: %d", f_manager.getFeatureCount());
    Headers[frame_count] = header;

    ImageFrame imageframe(image, header.stamp.toSec());
    imageframe.pre_integration = tmp_pre_integration;
    all_image_frame.insert(make_pair(header.stamp.toSec(), imageframe));
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};

    if (ESTIMATE_EXTRINSIC == 2)
    {
        ROS_INFO("calibrating extrinsic param, rotation movement is needed");

        if (frame_count != 0)
        {
            vector<pair<Vector3d, Vector3d> > corres = f_manager.getCorresponding(frame_count - 1, frame_count);
            Matrix3d calib_ric;
            if (initial_ex_rotation.CalibrationExRotation(corres, pre_integrations[frame_count]->delta_q, calib_ric))
            {
                // ROS_WARN("initial extrinsic rotation calib success");
                // ROS_WARN_STREAM("initial extrinsic rotation: " << endl << Utility::R2ypr(calib_ric).transpose());
                ric[0] = calib_ric;
                RIC[0] = calib_ric;
                ESTIMATE_EXTRINSIC = 1;
            }
            else
                ric[0] = calib_ric;
        }
    }

    if (solver_flag == INITIAL)
    {

        if (frame_count == WINDOW_SIZE)
        {
            bool result = false;
            if (ESTIMATE_EXTRINSIC != 2 && (header.stamp.toSec() - initial_timestamp) > 0.1)
            {

                result = initialStructure();
                initial_timestamp = header.stamp.toSec();
            }
            if (result)
            {
                solver_flag = NON_LINEAR;
                solveOdometry();
                slideWindow();
                f_manager.removeFailures();
                ROS_INFO("Initialization finish!");
                last_R = Rs[WINDOW_SIZE];
                last_P = Ps[WINDOW_SIZE];
                last_R0 = Rs[0];
                last_P0 = Ps[0];
            }
            else
            {
                slideWindow();
                RIC[0] = ric[0];
            }
        }
        else
            frame_count++;
    }
    else
    {
        TicToc t_solve;
        solveOdometry();
        ROS_DEBUG("solver costs: %fms", t_solve.toc());

        if (failureDetection())
        {
            ROS_WARN("failure detection!");
            failure_occur = 1;
            clearState();
            setParameter();
            ROS_WARN("system reboot!");
            // ESTIMATE_EXTRINSIC = 2;
            return;
        }

        TicToc t_margin;
        slideWindow();
        f_manager.removeFailures();
        ROS_DEBUG("marginalization costs: %fms", t_margin.toc());
        // prepare output of VINS
        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
    }

    // imageframe.R = Rs[10];
    // imageframe.T = Ps[10];
}

bool Estimator::initialStructure()
{
    TicToc t_sfm;
    //check imu observibility
    {
        map<double, ImageFrame>::iterator frame_it;
        Vector3d sum_g;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            sum_g += tmp_g;
        }
        Vector3d aver_g;
        aver_g = sum_g * 1.0 / ((int)all_image_frame.size() - 1);
        double var = 0;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            //cout << "frame g " << tmp_g.transpose() << endl;
        }
        var = sqrt(var / ((int)all_image_frame.size() - 1));
        //ROS_WARN("IMU variation %f!", var);
        if (var < 0.25)
        {
            ROS_INFO("IMU excitation not enouth!");
            //return false;
        }
    }

    // global sfm
    Quaterniond Q[frame_count + 1];
    Vector3d T[frame_count + 1];
    map<int, Vector3d> sfm_tracked_points;
    vector<SFMFeature> sfm_f;
    for (auto &it_per_id : f_manager.feature)
    {
        int imu_j = it_per_id.start_frame - 1;
        SFMFeature tmp_feature;
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            Vector3d pts_j = it_per_frame.point;
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()}));
        }
        sfm_f.push_back(tmp_feature);
    }
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
    if (!relativePose(relative_R, relative_T, l))
    {
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
    GlobalSFM sfm;
    // cout<<"frame_count="<<frame_count<<endl;
    if (!sfm.construct(frame_count + 1, Q, T, l,
                       relative_R, relative_T,
                       sfm_f, sfm_tracked_points))
    {
        ROS_DEBUG("global SFM failed!");
        marginalization_flag = MARGIN_OLD;
        return false;
    }

    //solve pnp for all frame
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin();

    for (int i = 0; frame_it != all_image_frame.end(); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
        if ((frame_it->first) == Headers[i].stamp.toSec())
        {
            frame_it->second.is_key_frame = true;
            //frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose();
            frame_it->second.R = Q[i].toRotationMatrix();
            frame_it->second.T = T[i];
            i++;
            continue;
        }
        if ((frame_it->first) > Headers[i].stamp.toSec())
        {
            i++;
        }
        // cout<<"pnp i="<<i<<endl;
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        // cout<<"pnp R_inital=\n"<<R_inital<<endl;
        Vector3d P_inital = -R_inital * T[i];
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec);
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id);
                if (it != sfm_tracked_points.end())
                {
                    Vector3d world_pts = it->second;
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);
        if (pts_3_vector.size() < 6)
        {
            cout << "pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_DEBUG("Not enough points for solve pnp !");
            return false;
        }
        // if (!cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        if (!cv::solvePnPRansac(pts_3_vector, pts_2_vector, K, D, rvec, t, true, cv::SOLVEPNP_EPNP))
        {
            ROS_DEBUG("solve pnp fail!");
            return false;
        }
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp, tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp);
        T_pnp = R_pnp * (-T_pnp);
        //frame_it->second.R = R_pnp * RIC[0].transpose();
        frame_it->second.R = R_pnp;
        frame_it->second.T = T_pnp;
    }
    if (visualInitialAlign())
        return true;
    else
    {
        ROS_INFO("misalign visual structure with IMU");
        return false;
    }
}

bool Estimator::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
    double *td_;
    VectorXd s;
    bool result = VisualIMUAlignment(all_image_frame, Bgs, Bas, g, x, td_, RIC[0], ric[0], TIC[0], s);
    if (!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // cout<<"frame_count="<<frame_count<<endl;
    // cout<<"all_image_frame="<<all_image_frame.size()<<endl;

    // change state
    // double scale(0.0);
    // for (int i = 0; i <= frame_count; i++){
    //     scale += s[i];
    // }
    // scale = scale/(frame_count+1);
    // cout<<"average scale="<<scale<<endl;
    double sumt(0.0);
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i].stamp.toSec()].R;
        Vector3d Pi = all_image_frame[Headers[i].stamp.toSec()].T;
        Ps[i] = Pi;
        Rs[i] = Ri;
        Ps[i] = s[i] * Ps[i] - Rs[i] * TIC[0];

        all_image_frame[Headers[i].stamp.toSec()].is_key_frame = true;
        sumt += td_[i];
    }

    Estimator::td = sumt / frame_count;
    // VectorXd dep = f_manager.getDepthVector();
    // for (int i = 0; i < dep.size(); i++)
    //     dep[i] = -1;
    // f_manager.clearDepth(dep);

    //triangulat on cam pose , no tic
    // Vector3d TIC_TMP[NUM_OF_CAM];
    // for(int i = 0; i < NUM_OF_CAM; i++)
    //     TIC_TMP[i].setZero();
    ric[0] = RIC[0];
    f_manager.setRic(ric);
    f_manager.triangulate(Ps, &(TIC[0]), &(RIC[0]));
    //f_manager.triangulate(Ps, &(TIC[0]), &(RIC[0]));

    //double s = (x.tail<1>())(0);

    // for (int i = 0; i <= WINDOW_SIZE; i++)
    // {
    //     pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    // }

    // for (int i = frame_count; i >= 0; i--)
    //     Ps[i] = s[i] * Ps[i] - Rs[i] * TIC[0];

    // for (int i = frame_count; i >= 0; i--)
    //     Ps[i] = s[i] * Ps[i] - Rs[i] * TIC[0] - (s[i] * Ps[0] - Rs[0] * TIC[0]);
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if (frame_i->second.is_key_frame)
        {
            kv++;
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }
    // for (auto &it_per_id : f_manager.feature)
    // {
    //     it_per_id.used_num = it_per_id.feature_per_frame.size();
    //     if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
    //         continue;
    //     it_per_id.estimated_depth *= s[it_per_id.endFrame()];
    // }

    Matrix3d R0 = Utility::g2R(g);
    // double yaw = Utility::R2ypr(R0 * Rs[0]).x();
    // R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!before gc0=" << g.transpose() << endl;
    g = R0 * g;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!after gc0=" << g.transpose() << endl;
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i];
        Rs[i] = rot_diff * Rs[i];
        Vs[i] = rot_diff * Vs[i];
    }
    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose());

    return true;
}

bool Estimator::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d> > corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;
            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if (average_parallax * 460 > 30 && m_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                return true;
            }
        }
    }
    return false;
}

void Estimator::solveOdometry()
{
    if (frame_count < WINDOW_SIZE)
        return;
    if (solver_flag == NON_LINEAR)
    {
        TicToc t_tri;
        f_manager.triangulate(Ps, tic, ric);
        ROS_DEBUG("triangulation costs %f", t_tri.toc());
        optimization();
    }
}

void Estimator::vector2double()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        para_SpeedBias[i][0] = Vs[i].x();
        para_SpeedBias[i][1] = Vs[i].y();
        para_SpeedBias[i][2] = Vs[i].z();

        para_SpeedBias[i][3] = Bas[i].x();
        para_SpeedBias[i][4] = Bas[i].y();
        para_SpeedBias[i][5] = Bas[i].z();

        para_SpeedBias[i][6] = Bgs[i].x();
        para_SpeedBias[i][7] = Bgs[i].y();
        para_SpeedBias[i][8] = Bgs[i].z();
    }
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        para_Feature[i][0] = dep(i);
    if (ESTIMATE_TD)
        para_Td[0][0] = td;
}

void Estimator::double2vector()
{
    Vector3d origin_R0 = Utility::R2ypr(Rs[0]);
    Vector3d origin_P0 = Ps[0];

    if (failure_occur)
    {
        origin_R0 = Utility::R2ypr(last_R0);
        origin_P0 = last_P0;
        failure_occur = 0;
    }
    Vector3d origin_R00 = Utility::R2ypr(Quaterniond(para_Pose[0][6],
                                                     para_Pose[0][3],
                                                     para_Pose[0][4],
                                                     para_Pose[0][5])
                                             .toRotationMatrix());
    double y_diff = origin_R0.x() - origin_R00.x();
    //TODO
    Matrix3d rot_diff = Utility::ypr2R(Vector3d(y_diff, 0, 0));
    if (abs(abs(origin_R0.y()) - 90) < 1.0 || abs(abs(origin_R00.y()) - 90) < 1.0)
    {
        ROS_DEBUG("euler singular point!");
        rot_diff = Rs[0] * Quaterniond(para_Pose[0][6],
                                       para_Pose[0][3],
                                       para_Pose[0][4],
                                       para_Pose[0][5])
                               .toRotationMatrix()
                               .transpose();
    }

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {

        Rs[i] = rot_diff * Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();

        Ps[i] = rot_diff * Vector3d(para_Pose[i][0] - para_Pose[0][0],
                                    para_Pose[i][1] - para_Pose[0][1],
                                    para_Pose[i][2] - para_Pose[0][2]) +
                origin_P0;

        Vs[i] = rot_diff * Vector3d(para_SpeedBias[i][0],
                                    para_SpeedBias[i][1],
                                    para_SpeedBias[i][2]);

        Bas[i] = Vector3d(para_SpeedBias[i][3],
                          para_SpeedBias[i][4],
                          para_SpeedBias[i][5]);

        Bgs[i] = Vector3d(para_SpeedBias[i][6],
                          para_SpeedBias[i][7],
                          para_SpeedBias[i][8]);
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d(para_Ex_Pose[i][0],
                          para_Ex_Pose[i][1],
                          para_Ex_Pose[i][2]);
        ric[i] = Quaterniond(para_Ex_Pose[i][6],
                             para_Ex_Pose[i][3],
                             para_Ex_Pose[i][4],
                             para_Ex_Pose[i][5])
                     .normalized()
                     .toRotationMatrix();
    }

    if (ESTIMATE_TLB == true)
    {
        TLB = Vector3d(para_Ex_Pose_Lb[0],
                       para_Ex_Pose_Lb[1],
                       para_Ex_Pose_Lb[2]);
        RLB = Quaterniond(para_Ex_Pose_Lb[6],
                          para_Ex_Pose_Lb[3],
                          para_Ex_Pose_Lb[4],
                          para_Ex_Pose_Lb[5])
                  .normalized()
                  .toRotationMatrix();
        TBL = -RLB.transpose() * TLB;
        // ROS_WARN_STREAM("lidar extrinsic rotation changed to: " << endl
        //                                                         << Utility::R2ypr(RLB).transpose());
        // ROS_WARN_STREAM("lidar extrinsic translation changed to: " << endl
        //                                                            << TLB.transpose());
        ESTIMATE_TLB = false;
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        dep(i) = para_Feature[i][0];
    f_manager.setDepth(dep);
    if (ESTIMATE_TD)
        td = para_Td[0][0];

    // relative info between two loop frame
    if (relocalization_info)
    {
        Matrix3d relo_r;
        Vector3d relo_t;
        relo_r = rot_diff * Quaterniond(relo_Pose[6], relo_Pose[3], relo_Pose[4], relo_Pose[5]).normalized().toRotationMatrix();
        relo_t = rot_diff * Vector3d(relo_Pose[0] - para_Pose[0][0],
                                     relo_Pose[1] - para_Pose[0][1],
                                     relo_Pose[2] - para_Pose[0][2]) +
                 origin_P0;
        double drift_correct_yaw;
        drift_correct_yaw = Utility::R2ypr(prev_relo_r).x() - Utility::R2ypr(relo_r).x();
        drift_correct_r = Utility::ypr2R(Vector3d(drift_correct_yaw, 0, 0));
        drift_correct_t = prev_relo_t - drift_correct_r * relo_t;
        relo_relative_t = relo_r.transpose() * (Ps[relo_frame_local_index] - relo_t);
        relo_relative_q = relo_r.transpose() * Rs[relo_frame_local_index];
        relo_relative_yaw = Utility::normalizeAngle(Utility::R2ypr(Rs[relo_frame_local_index]).x() - Utility::R2ypr(relo_r).x());
        //cout << "vins relo " << endl;
        //cout << "vins relative_t " << relo_relative_t.transpose() << endl;
        //cout << "vins relative_yaw " <<relo_relative_yaw << endl;
        relocalization_info = 0;
    }
}

bool Estimator::failureDetection()
{
    if (f_manager.last_track_num < 2)
    {
        ROS_INFO(" little feature %d", f_manager.last_track_num);
        //return true;
    }
    if (Bas[WINDOW_SIZE].norm() > 2.5)
    {
        ROS_INFO(" big IMU acc bias estimation %f", Bas[WINDOW_SIZE].norm());
        return true;
    }
    if (Bgs[WINDOW_SIZE].norm() > 1.0)
    {
        ROS_INFO(" big IMU gyr bias estimation %f", Bgs[WINDOW_SIZE].norm());
        return true;
    }
    /*
    if (tic(0) > 1)
    {
        ROS_INFO(" big extri param estimation %d", tic(0) > 1);
        return true;
    }
    */
    Vector3d tmp_P = Ps[WINDOW_SIZE];
    if ((tmp_P - last_P).norm() > 10) //5
    {
        ROS_INFO(" big translation");
        return true;
    }
    if (abs(tmp_P.z() - last_P.z()) > 1)
    {
        ROS_INFO(" big z translation");
        return true;
    }
    Matrix3d tmp_R = Rs[WINDOW_SIZE];
    Matrix3d delta_R = tmp_R.transpose() * last_R;
    Quaterniond delta_Q(delta_R);
    double delta_angle;
    delta_angle = acos(delta_Q.w()) * 2.0 / 3.14 * 180.0;
    if (delta_angle > 50)
    {
        ROS_INFO(" big delta_angle ");
        //return true;
    }
    return false;
}

void Estimator::optimization()
{
    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = new ceres::HuberLoss(1.0);
    loss_function = new ceres::CauchyLoss(1.0);
    // loss_function = new ceres::CauchyLoss(2.3849);
    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
        problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);
        // problem.SetParameterUpperBound(para_SpeedBias[i], 3, 0.5); //bas
        // problem.SetParameterLowerBound(para_SpeedBias[i], 3, -0.5);
        // problem.SetParameterUpperBound(para_SpeedBias[i], 4, 0.5);
        // problem.SetParameterLowerBound(para_SpeedBias[i], 4, -0.5);
        // problem.SetParameterUpperBound(para_SpeedBias[i], 5, 0.5);
        // problem.SetParameterLowerBound(para_SpeedBias[i], 5, -0.5);
        // problem.SetParameterUpperBound(para_SpeedBias[i], 6, 0.3);//bgs
        // problem.SetParameterLowerBound(para_SpeedBias[i], 6, -0.3);
        // problem.SetParameterUpperBound(para_SpeedBias[i], 7, 0.3);
        // problem.SetParameterLowerBound(para_SpeedBias[i], 7, -0.3);
        // problem.SetParameterUpperBound(para_SpeedBias[i], 8, 0.3);
        // problem.SetParameterLowerBound(para_SpeedBias[i], 8, -0.3);
    }
    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Ex_Pose[i], SIZE_POSE, local_parameterization);

        if (!ESTIMATE_EXTRINSIC)
        {
            ROS_DEBUG("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose[i]);
        }
        else
            ROS_DEBUG("estimate extinsic param");
    }
    if (ESTIMATE_TD)
    {
        problem.AddParameterBlock(para_Td[0], 1);
        //problem.SetParameterBlockConstant(para_Td[0]);
    }

    TicToc t_whole, t_prepare;
    vector2double();

    if (last_marginalization_info)
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
    }

    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        int j = i + 1;
        if (pre_integrations[j]->sum_dt > 10.0)
            continue;
        IMUFactor *imu_factor = new IMUFactor(pre_integrations[j]);
        problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
    }
    int f_m_cnt = 0;
    int feature_index = -1;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;

        ++feature_index;

        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;

        Vector3d pts_i = it_per_id.feature_per_frame[0].point;

        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i == imu_j)
            {
                continue;
            }
            Vector3d pts_j = it_per_frame.point;
            if (ESTIMATE_TD) //
            {

                ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                  it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td,
                                                                  it_per_id.feature_per_frame[0].uv.y(), it_per_frame.uv.y());
                problem.AddResidualBlock(f_td, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]);

                if (it_per_id.lidar_depth_flag == true)
                {
                    // para_Feature[feature_index][0] = it_per_id.estimated_depth;
                    problem.SetParameterBlockConstant(para_Feature[feature_index]);
                }

                /*
                    double **para = new double *[5];
                    para[0] = para_Pose[imu_i];
                    para[1] = para_Pose[imu_j];
                    para[2] = para_Ex_Pose[0];
                    para[3] = para_Feature[feature_index];
                    para[4] = para_Td[0];
                    f_td->check(para);
                    */
            }
            else
            {
                ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]);
                if (it_per_id.lidar_depth_flag == true)
                    problem.SetParameterBlockConstant(para_Feature[feature_index]);
            }
            f_m_cnt++;
        }
    }

    ROS_DEBUG("visual measurement count: %d", f_m_cnt);
    ROS_DEBUG("prepare for ceres: %f", t_prepare.toc());

    if (relocalization_info)
    {
        //printf("set relocalization factor! \n");
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(relo_Pose, SIZE_POSE, local_parameterization);
        int retrive_feature_index = 0;
        int feature_index = -1;
        for (auto &it_per_id : f_manager.feature)
        {
            it_per_id.used_num = it_per_id.feature_per_frame.size();
            if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
                continue;
            ++feature_index;
            int start = it_per_id.start_frame;
            if (start <= relo_frame_local_index)
            {
                while ((int)match_points[retrive_feature_index].z() < it_per_id.feature_id)
                {
                    retrive_feature_index++;
                }
                if ((int)match_points[retrive_feature_index].z() == it_per_id.feature_id)
                {
                    Vector3d pts_j = Vector3d(match_points[retrive_feature_index].x(), match_points[retrive_feature_index].y(), 1.0);
                    Vector3d pts_i = it_per_id.feature_per_frame[0].point;

                    ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                    problem.AddResidualBlock(f, loss_function, para_Pose[start], relo_Pose, para_Ex_Pose[0], para_Feature[feature_index]);
                    retrive_feature_index++;
                }
            }
        }
    }

    LidarLPSConstraint LPSmarg;
    bool NeedLPSmarg{false};
    int LPSMargID_l, LPSMargID_r;
    while (LidarLPSConstraints.size() > 7)
    {
        LidarLPSConstraints.pop_front();
    }
    if (LPS_call)
    {
        LPS_t = LPS_q.toRotationMatrix() * TLB + LPS_t;
        LPS_q = LPS_q * Quaterniond(RLB);
        LidarLPSConstraint current_lps;
        current_lps.LPSq = LPS_q;
        current_lps.LPSt = LPS_t;
        current_lps.lidar_t = LPS_time;

        LidarLPSConstraints.push_back(current_lps);

        LPS_call = false;
        // }
        std::list<LidarLPSConstraint>::iterator iter_;
        for (iter_ = LidarLPSConstraints.begin(); iter_ != LidarLPSConstraints.end(); ++iter_)
        {
            // problem.SetParameterBlockConstant(para_Pose[0]);
            // cout << "LidarLPSConstraints size:" << LidarLPSConstraints.size() << endl;
            LidarLPSConstraint tem_ = *iter_;
            int id_l{0}, id_r{0};
            if (FindNearest2ID(Headers, tem_.lidar_t, id_l, id_r))
            {
                double check_time = Headers[id_r].stamp.toSec() - Headers[id_l].stamp.toSec();
                cout << "LPS time check: " << check_time << endl;
                if (id_l == 0 && check_time < 0.2)
                {
                    LPSmarg = tem_;
                    LPSMargID_l = id_l;
                    LPSMargID_r = id_r;
                    NeedLPSmarg = true;
                }
                // problem.SetParameterBlockConstant(para_Pose[id_l]);
                if (check_time < 0.2)
                {
                    ceres::CostFunction *cost_function = LPSConstraint::makeConstraint(Headers[id_l].stamp.toSec(), Headers[id_r].stamp.toSec(), tem_.lidar_t, tem_.LPSq, tem_.LPSt, 0.1);
                    problem.AddResidualBlock(cost_function, loss_function, para_Pose[id_l], para_Pose[id_r]);
                }
            }
        }
    }

    // Quaterniond rlb_q(RLB);
    // para_Ex_Pose_Lb[0] = TLB[0];
    // para_Ex_Pose_Lb[1] = TLB[1];
    // para_Ex_Pose_Lb[2] = TLB[2];
    // para_Ex_Pose_Lb[6] = rlb_q.w();
    // para_Ex_Pose_Lb[3] = rlb_q.x();
    // para_Ex_Pose_Lb[4] = rlb_q.y();
    // para_Ex_Pose_Lb[5] = rlb_q.z();
    // ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
    // problem.AddParameterBlock(para_Ex_Pose_Lb, SIZE_POSE, local_parameterization);
    // problem.SetParameterBlockConstant(para_Ex_Pose_Lb);

    LidarICPConstraint ICPmarg;
    bool NeedICPmarg{false};
    int ICPMargID_a, ICPMargID_b, ICPMargID_c, ICPMargID_d;
    std::list<LidarICPConstraint>::iterator iter;

    while (LidarICPConstraints.size() > 5)
    {
        LidarICPConstraints.pop_front();
    }

    for (iter = LidarICPConstraints.begin(); iter != LidarICPConstraints.end(); ++iter)
    {
        int id_a{0}, id_b{0}, id_c{0}, id_d{0};
        LidarICPConstraint tem = *iter;
        if (tem.constraint_mode == 4)
        {
            // for (int i = 0; i < WINDOW_SIZE + 1; i++)
            // {
            //     para_SpeedBias[i][0] = para_SpeedBias[i][1] = para_SpeedBias[i][2] = 0.0;
            //     problem.SetParameterBlockConstant(para_SpeedBias[i]);
            //     problem.SetParameterBlockConstant(para_Pose[i]);
            // }

            ROS_WARN("ZERO VELOCITY!!!!!");
            // para_SpeedBias[WINDOW_SIZE][0] = para_SpeedBias[WINDOW_SIZE][1] = para_SpeedBias[WINDOW_SIZE][2] = 0.0;
            // problem.SetParameterBlockConstant(para_SpeedBias[WINDOW_SIZE]);
            // problem.SetParameterBlockConstant(para_Pose[WINDOW_SIZE]);

            para_SpeedBias[WINDOW_SIZE - 1][0] = para_SpeedBias[WINDOW_SIZE - 1][1] = para_SpeedBias[WINDOW_SIZE - 1][2] = 0.0;
            problem.SetParameterBlockConstant(para_SpeedBias[WINDOW_SIZE - 1]);
            problem.SetParameterBlockConstant(para_Pose[WINDOW_SIZE - 1]);

            // ceres::CostFunction *cost_function = LidarICPConstraint_b::makeConstraint_b(tem.lidar_ta, tem.lidar_tb, tem.lidar_tc, tem.lidar_td,
            //                                                                             tem.lidar_ti, tem.lidar_tj, tem.lidar_trans, tem.lidar_sqrt_info);
            // problem.AddResidualBlock(cost_function, NULL, para_Pose[3], para_Pose[4], para_Pose[5], para_Pose[6]);
        }
        else if (tem.constraint_mode == 3 && FindWindowsID(Headers,
                               tem.lidar_ta, tem.lidar_tb, tem.lidar_tc, tem.lidar_td,
                               id_a, id_b, id_c, id_d))
        {
            cout << "Windows ID find!!" << id_a << " " << id_b << " " << id_c << " " << id_d << " " << endl;
            if (id_a == 0)
            {
                ICPmarg = tem;
                ICPMargID_a = id_a;
                ICPMargID_b = id_b;
                ICPMargID_c = id_c;
                ICPMargID_d = id_d;
                NeedICPmarg = true;
            }

            ceres::CostFunction *cost_function = LidarICPConstraint_b::makeConstraint_b(tem.lidar_ta, tem.lidar_tb, tem.lidar_tc, tem.lidar_td,
                                                                                        tem.lidar_ti, tem.lidar_tj, tem.lidar_trans, tem.lidar_sqrt_info);

            // problem.AddResidualBlock(cost_function, NULL, para_Pose[id_a], para_Pose[id_b], para_Pose[id_c], para_Pose[id_d], para_Ex_Pose_Lb);
            problem.AddResidualBlock(cost_function, loss_function, para_Pose[id_a], para_Pose[id_b], para_Pose[id_c], para_Pose[id_d]);
            // ESTIMATE_TLB = false;
        }
    }

    ceres::Solver::Options options;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = NUM_ITERATIONS;
    //options.use_explicit_schur_complement = true;
    //options.minimizer_progress_to_stdout = true;
    //options.use_nonmonotonic_steps = true;
    // if (marginalization_flag == MARGIN_OLD)
    //     options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0;
    // else
    options.max_solver_time_in_seconds = SOLVER_TIME;
    TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // cout << summary.BriefReport() << endl;
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    ROS_DEBUG("solver costs: %f", t_solver.toc());

    double2vector();

    cov_estimate = false;
    if (cov_estimate)
    {
        ceres::Covariance::Options cov_options;
        cov_options.algorithm_type = ceres::DENSE_SVD;
        cov_options.min_reciprocal_condition_number = 10e-14;
        cov_options.num_threads = 4;
        cov_options.null_space_rank = -1;
        ceres::Covariance covariance(cov_options);
        std::vector<std::pair<const double *, const double *> > covariance_blocks;

        // covariance_blocks.push_back(std::make_pair(para_Pose[0], para_Pose[0]));
        // covariance_blocks.push_back(std::make_pair(para_SpeedBias[3], para_SpeedBias[3]));
        // covariance_blocks.push_back(std::make_pair(para_Ex_Pose[0], para_Ex_Pose[0]));
        // covariance.Compute(covariance_blocks, &problem);

        // Eigen::Matrix<double, 6, 6, Eigen::RowMajor> cov_pose = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Zero();
        // covariance.GetCovarianceBlockInTangentSpace(para_Pose[0],para_Pose[0], cov_pose.data());

        // Eigen::Matrix<double, 6, 6, Eigen::RowMajor> cov_ex = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Zero();
        // covariance.GetCovarianceBlockInTangentSpace(para_Ex_Pose[0],para_Ex_Pose[0], cov_ex.data());

        // Eigen::Matrix<double, 9, 9, Eigen::RowMajor> cov_speedbias = Eigen::Matrix<double, 9, 9, Eigen::RowMajor>::Zero();
        // covariance.GetCovarianceBlock(para_SpeedBias[3],para_SpeedBias[3], cov_speedbias.data());

        // ROS_WARN_STREAM("POSE COV TEST!!!" << endl << cov_pose);
        // ROS_WARN_STREAM("SPEED-BIAS COV TEST!!!" << endl << cov_speedbias);
        // ROS_WARN_STREAM("EX COV TEST!!!" << endl << cov_ex);

        // Eigen::JacobiSVD<Eigen::MatrixXd> svd1(cov_pose.inverse(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        // std::cout << "!!!!!!!!!!!!TEST--POSE--HESSIAN!!!!!!!!!!\n"<< svd1.singularValues() << std::endl;
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd2(cov_speedbias.inverse(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        // std::cout << "!!!!!!!!!!!!TEST--SPEED-BIAS--HESSIAN!!!!!!!!!!\n"<< svd2.singularValues() << std::endl;
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd3(cov_ex.inverse(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        // std::cout << "!!!!!!!!!!!!TEST--EX--HESSIAN!!!!!!!!!!\n"<< svd3.singularValues() << std::endl;

        for (int i = 0; i < WINDOW_SIZE; i++)
        {
            for (int j = 0; j < WINDOW_SIZE; j++)
            {
                covariance_blocks.push_back(std::make_pair(para_Pose[i], para_Pose[j]));
            }
        }
        covariance.Compute(covariance_blocks, &problem);
        Eigen::Matrix<double, 42, 42, Eigen::RowMajor> H_pose = Eigen::Matrix<double, 42, 42, Eigen::RowMajor>::Zero();

        for (int i = 0; i < WINDOW_SIZE; i++)
        {
            for (int j = 0; j < WINDOW_SIZE; j++)
            {
                Eigen::Matrix<double, 6, 6, Eigen::RowMajor> cov_pose = Eigen::Matrix<double, 6, 6, Eigen::RowMajor>::Zero();
                covariance.GetCovarianceBlockInTangentSpace(para_Pose[i], para_Pose[j], cov_pose.data());
                // ROS_WARN_STREAM("POSE COV TEST!!!" << endl << cov_pose);
                H_pose.block(i * 6, j * 6, 6, 6) = cov_pose;
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(H_pose.inverse(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        std::cout << "!!!!!!!!!!!!TEST--POSE--HESSIAN!!!!!!!!!!\n"
                  << svd.singularValues() << std::endl;
    }

    TicToc t_whole_marginalization;
    if (marginalization_flag == MARGIN_OLD)
    {
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        vector2double();

        if (last_marginalization_info)
        {
            vector<int> drop_set;
            for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
            {
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] ||
                    last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);
            }
            // construct new marginlization_factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                           last_marginalization_parameter_blocks,
                                                                           drop_set);

            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        {
            if (NeedICPmarg)
            {
                ROS_WARN("ADD ICPMarg!!!!!");
                ceres::CostFunction *cost_function = LidarICPConstraint_b::makeConstraint_b(ICPmarg.lidar_ta, ICPmarg.lidar_tb, ICPmarg.lidar_tc, ICPmarg.lidar_td,
                                                                                            ICPmarg.lidar_ti, ICPmarg.lidar_tj, ICPmarg.lidar_trans, ICPmarg.lidar_sqrt_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(cost_function, loss_function,
                                                                               vector<double *>{para_Pose[0], para_Pose[ICPMargID_b], para_Pose[ICPMargID_c], para_Pose[ICPMargID_d]},
                                                                               vector<int>{0});
                marginalization_info->addResidualBlockInfo(residual_block_info);
                NeedICPmarg = false;
            }
        }

        {
            if (NeedLPSmarg)
            {
                ROS_WARN("ADD LPSMarg!!!!!");
                ceres::CostFunction *cost_function = LPSConstraint::makeConstraint(Headers[LPSMargID_l].stamp.toSec(), Headers[LPSMargID_r].stamp.toSec(),
                                                                                   LPSmarg.lidar_t, LPSmarg.LPSq, LPSmarg.LPSt, 0.1);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(cost_function, loss_function,
                                                                               vector<double *>{para_Pose[0], para_Pose[LPSMargID_r]},
                                                                               vector<int>{0});
                marginalization_info->addResidualBlockInfo(residual_block_info);
                NeedICPmarg = false;
            }
        }

        {
            if (pre_integrations[1]->sum_dt < 10.0)
            {
                IMUFactor *imu_factor = new IMUFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                               vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                               vector<int>{0, 1});
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }

        {
            int feature_index = -1;
            for (auto &it_per_id : f_manager.feature)
            {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                if (!(it_per_id.used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
                    continue;

                ++feature_index;

                int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
                if (imu_i != 0)
                    continue;

                Vector3d pts_i = it_per_id.feature_per_frame[0].point;

                for (auto &it_per_frame : it_per_id.feature_per_frame)
                {
                    imu_j++;
                    if (imu_i == imu_j)
                        continue;

                    Vector3d pts_j = it_per_frame.point;
                    if (ESTIMATE_TD)
                    {
                        ProjectionTdFactor *f_td = new ProjectionTdFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td,
                                                                          it_per_id.feature_per_frame[0].uv.y(), it_per_frame.uv.y());
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f_td, loss_function,
                                                                                       vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]},
                                                                                       vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    else
                    {
                        ProjectionFactor *f = new ProjectionFactor(pts_i, pts_j);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                                                                                       vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index]},
                                                                                       vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                }
            }
        }

        TicToc t_pre_margin;
        marginalization_info->preMarginalize();
        ROS_DEBUG("pre marginalization %f ms", t_pre_margin.toc());

        TicToc t_margin;
        marginalization_info->marginalize();

        ROS_DEBUG("marginalization %f ms", t_margin.toc());
        std::unordered_map<long, double *> addr_shift;
        for (int i = 1; i <= WINDOW_SIZE; i++)
        {
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
            addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
        }
        for (int i = 0; i < NUM_OF_CAM; i++)
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
        if (ESTIMATE_TD)
        {
            addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
        }
        vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);

        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;
        last_marginalization_parameter_blocks = parameter_blocks;
    }
    else
    {
        if (last_marginalization_info &&
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();
            vector2double();
            if (last_marginalization_info)
            {
                vector<int> drop_set;
                for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
                {
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }
                // construct new marginlization_factor
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks,
                                                                               drop_set);

                marginalization_info->addResidualBlockInfo(residual_block_info);
            }

            TicToc t_pre_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->preMarginalize();
            ROS_DEBUG("end pre marginalization, %f ms", t_pre_margin.toc());

            TicToc t_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->marginalize();
            ROS_DEBUG("end marginalization, %f ms", t_margin.toc());

            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                }
                else
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                }
            }
            for (int i = 0; i < NUM_OF_CAM; i++)
                addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];
            if (ESTIMATE_TD)
            {
                addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
            }

            vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = parameter_blocks;
        }
    }
    ROS_DEBUG("whole marginalization costs: %f", t_whole_marginalization.toc());

    ROS_DEBUG("whole time for ceres: %f", t_whole.toc());
}

void Estimator::slideWindow()
{
    TicToc t_margin;
    if (marginalization_flag == MARGIN_OLD)
    {
        double t_0 = Headers[0].stamp.toSec();

        back_R0 = Rs[0];
        back_P0 = Ps[0];
        if (frame_count == WINDOW_SIZE)
        {
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                Rs[i].swap(Rs[i + 1]);

                std::swap(pre_integrations[i], pre_integrations[i + 1]);

                dt_buf[i].swap(dt_buf[i + 1]);
                linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                Headers[i] = Headers[i + 1];
                Ps[i].swap(Ps[i + 1]);
                Vs[i].swap(Vs[i + 1]);
                Bas[i].swap(Bas[i + 1]);
                Bgs[i].swap(Bgs[i + 1]);
            }
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];
            Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
            Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];

            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();

            if (true || solver_flag == INITIAL)
            {
                // double t_0 = Headers[0].stamp.toSec();
                map<double, ImageFrame>::iterator it_0;
                it_0 = all_image_frame.find(t_0);

                delete it_0->second.pre_integration;
                it_0->second.pre_integration = nullptr;

                for (map<double, ImageFrame>::iterator it = all_image_frame.begin(); it != it_0; ++it)
                {
                    if (it->second.pre_integration)
                        delete it->second.pre_integration;
                    it->second.pre_integration = nullptr;
                }

                all_image_frame.erase(all_image_frame.begin(), it_0);
                all_image_frame.erase(t_0);
            }
            slideWindowOld();
        }
    }
    else
    {
        if (frame_count == WINDOW_SIZE)
        {

            double t_0 = Headers[WINDOW_SIZE - 1].stamp.toSec();
            for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
            {
                double tmp_dt = dt_buf[frame_count][i];
                Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];
                pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity);
                dt_buf[frame_count - 1].push_back(tmp_dt);
                linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
            }

            Headers[frame_count - 1] = Headers[frame_count];
            Ps[frame_count - 1] = Ps[frame_count];
            Vs[frame_count - 1] = Vs[frame_count];
            Rs[frame_count - 1] = Rs[frame_count];
            Bas[frame_count - 1] = Bas[frame_count];
            Bgs[frame_count - 1] = Bgs[frame_count];
            delete pre_integrations[WINDOW_SIZE];
            pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

            dt_buf[WINDOW_SIZE].clear();
            linear_acceleration_buf[WINDOW_SIZE].clear();
            angular_velocity_buf[WINDOW_SIZE].clear();
            if (solver_flag == INITIAL)
            {
                all_image_frame.erase(t_0);
            }
            slideWindowNew();
        }
    }
}

// real marginalization is removed in solve_ceres()
void Estimator::slideWindowNew()
{
    sum_of_front++;
    f_manager.removeFront(frame_count);
}
// real marginalization is removed in solve_ceres()
void Estimator::slideWindowOld()
{
    sum_of_back++;

    bool shift_depth = solver_flag == NON_LINEAR ? true : false;
    if (shift_depth)
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        R0 = back_R0 * ric[0];
        R1 = Rs[0] * ric[0];
        P0 = back_P0 + back_R0 * tic[0];
        P1 = Ps[0] + Rs[0] * tic[0];
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
    }
    else
        f_manager.removeBack();
}

void Estimator::setReloFrame(double _frame_stamp, int _frame_index, vector<Vector3d> &_match_points, Vector3d _relo_t, Matrix3d _relo_r)
{
    relo_frame_stamp = _frame_stamp;
    relo_frame_index = _frame_index;
    match_points.clear();
    match_points = _match_points;
    prev_relo_t = _relo_t;
    prev_relo_r = _relo_r;
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        if (relo_frame_stamp == Headers[i].stamp.toSec())
        {
            relo_frame_local_index = i;
            relocalization_info = 1;
            for (int j = 0; j < SIZE_POSE; j++)
                relo_Pose[j] = para_Pose[i][j];
        }
    }
}