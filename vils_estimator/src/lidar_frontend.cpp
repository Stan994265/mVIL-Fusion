#include "lidar_frontend.h"

LidarCalibration::LidarCalibration()
{
    Reset();
    Rcam.push_back(Matrix3d::Identity());
    Rc_g.push_back(Matrix3d::Identity());
    Rlidar.push_back(Matrix3d::Identity());
    Tlidar.push_back(Vector3d::Zero());
    Tbody.push_back(Vector3d::Zero());
    FS.push_back(0.0);
    frame_count = 0;
    temp_ypr.setZero();
    rlb = Matrix3d::Identity();
    lidar_init_flag = true;
}

void LidarCalibration::Reset()
{
    frame_count = 0;
    temp_ypr.setZero();
    Rcam.clear();
    Rc_g.clear();
    Rlidar.clear();
    Tlidar.clear();
    Tbody.clear();
    FS.clear();
    rlb = Matrix3d::Identity();
    lidar_init_flag = true;
}

bool LidarCalibration::CalibrationLidarExRotation(const LidarFrame &lidarframei, const LidarFrame &lidarframej, double fitness_core,
                                                  Matrix3d delta_lidar_r, Matrix3d &calib_rlb_result,
                                                  Vector3d delta_lidar_t, Vector3d &calib_tlb_result)
{
    frame_count++;
    Rcam.push_back(PredictRelativeR(lidarframei, lidarframej)); //Rcij
    Rlidar.push_back(delta_lidar_r);                            //qij
    Rc_g.push_back(rlb.inverse() * delta_lidar_r * rlb);
    Tlidar.push_back(delta_lidar_t);                             //Tlidarij
    Tbody.push_back(PredictRelativeT(lidarframei, lidarframej)); //Tbodyij
    FS.push_back(fitness_core);                                  //weight
    cout << "!!!!!!!!!!!test: tlidar=" << delta_lidar_t.transpose() << "  tbody=" << Tbody.back().transpose() << endl;

    //temp_ypr
    // for (int i = 3; i <= frame_count; i++)
    // {
    //     temp_ypr += Utility::R2ypr(Rlidar[i]).cwiseAbs();
    // }
    // cout << "tem_pyr=" << temp_ypr.transpose() << endl;
    // cout << "init lidar frame_count=" << frame_count << endl;
    // if (frame_count < 3 || temp_ypr[0] < 5.0 || temp_ypr[1] < 5.0 || temp_ypr[2] < 5.0)
    // {
    //     cout << "Lidar rotation not enough!! " << endl;
    //     return false;
    // }

    if (frame_count >= 200)
    {
        cout << "Too long time, restart!!! " << endl;
        Reset();
        return false;
    }

    Eigen::MatrixXd A(frame_count * 4, 4);
    A.setZero();
    int sum_ok = 0;
    for (int i = 1; i <= frame_count; i++)
    {
        Quaterniond r1(Rcam[i]);
        Quaterniond r2(Rc_g[i]);

        double angular_distance = 180 / M_PI * r1.angularDistance(r2);
        ROS_DEBUG(
            "%d %f", i, angular_distance);
        cout << "angular_distance" << i << "=" << angular_distance << endl;

        double huber = angular_distance > 3.0 ? 3.0 / angular_distance : 1.0;
        double fs;
        if (LeafSize > 0.21)
        {
            fs = FS[i] > 1.0 ? 1 / FS[i] : 1.0; //outdoor
        }
        if (LeafSize < 0.16)
        {
            fs = FS[i] > 0.05 ? 0.0 : 1.0; //indoor
        }

        huber = huber * fs;
        //double huber = angular_distance > 3.0 ? 0.5 : 1.0;
        // double huber;
        // if (angular_distance>=20.0) huber = 0.0;
        // if (angular_distance>=10.0 && angular_distance<20) huber = 0.1;
        // if (angular_distance>5.0 && angular_distance<10) huber = 0.5;
        // if (angular_distance<=5.0) huber = 1.0;

        ++sum_ok;
        Matrix4d L, R;
        Quaterniond R_ab(Rcam[i]);
        double w = R_ab.w();
        Vector3d q = R_ab.vec();
        L.block<3, 3>(0, 0) = w * Matrix3d::Identity() + Utility::skewSymmetric(q);
        L.block<3, 1>(0, 3) = q;
        L.block<1, 3>(3, 0) = -q.transpose();
        L(3, 3) = w;
        Quaterniond R_ij(Rlidar[i]);
        w = R_ij.w();
        q = R_ij.vec();
        R.block<3, 3>(0, 0) = w * Matrix3d::Identity() - Utility::skewSymmetric(q);
        R.block<3, 1>(0, 3) = q;
        R.block<1, 3>(3, 0) = -q.transpose();
        R(3, 3) = w;

        A.block<4, 4>((i - 1) * 4, 0) = huber * (L - R);
    }

    JacobiSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);
    Matrix<double, 4, 1> x = svd.matrixV().transpose().row(3);
    Quaterniond estimated_R(x);
    rlb = estimated_R.toRotationMatrix().transpose();

    cout << "svd.singularValues()=\n"
         << svd.singularValues() << endl;
    cout << "RLB=\n"
         << Utility::R2ypr(rlb).transpose() << endl;
    calib_rlb_result = rlb;
    //if (frame_count >= 10 && temp_ypr[0]>30 && temp_ypr[1]>30 && temp_ypr[2]>30) //&& temp_ypr[0]>30 && temp_ypr[1]>30 && temp_ypr[2]>30
    Vector3d rlb_cov;
    rlb_cov = svd.singularValues().tail<3>();
    if (frame_count >= 30 && rlb_cov(1) > 0.15)
    {
        calib_rlb_result = rlb;
        //calib_rlb_result = Matrix3d::Identity();
        if (Lidar_align(Rlidar, Tlidar, Tbody, calib_tlb_result, calib_rlb_result, frame_count))
            return true;
        else
            return false;
    }
    else
        return false;
}

Matrix3d LidarCalibration::PredictRelativeR(const LidarFrame &lidarframei, const LidarFrame &lidarframej)
{
    Matrix3d Rci_cj;
    Matrix3d Rcj = PredictR(lidarframej);
    Matrix3d Rci = PredictR(lidarframei);
    Rci_cj = Rci.inverse() * Rcj;
    return Rci_cj;
}

// Matrix3d LidarCalibration::PredictRelativeR(const LidarFrame &lidarframei, const LidarFrame &lidarframej)
// {
//     Matrix3d Rci_cj;
//     Matrix3d Rci = lidarframei.vioData.Rwbi;
//     Matrix3d Rcj = lidarframej.vioData.Rwbj;
//     Rci_cj = Rci.inverse() * Rcj;
//     return Rci_cj;
// }

Matrix3d LidarCalibration::PredictR(const LidarFrame &lidarframe)
{
    Quaterniond qc0ca(lidarframe.vioData.Rwbi);
    Quaterniond qc0cb(lidarframe.vioData.Rwbj);
    double ta, tl, tb;
    ta = lidarframe.vioData.ti;
    tb = lidarframe.vioData.tj;
    tl = lidarframe.lidarData.point_cloud.time;
    double t = (tl - ta) / (tb - ta);
    // cout << "t=" << t << endl;
    Quaterniond q;
    if (t > 0)
    {
        q = qc0ca.slerp(t, qc0cb);
        q.normalized();
    }
    else
    {
        q = qc0ca;
        q.normalized();
    }

    return q.toRotationMatrix();
}

Vector3d LidarCalibration::PredictRelativeT(const LidarFrame &lidarframei, const LidarFrame &lidarframej)
{
    Vector3d Tbi_bj;
    Vector3d Tbj = PredictT(lidarframej);
    Vector3d Tbi = PredictT(lidarframei);
    Matrix3d Rbi = PredictR(lidarframei);
    // cout << "!!!!!!!test: Tbi=" << Tbi.transpose() << "  Tbj=" << Tbj.transpose() << endl;
    Tbi_bj = Rbi.transpose() * Tbj - Rbi.transpose() * Tbi;
    // cout << "!!!!!!!test:Tbij=" << Tbi_bj.transpose() << endl;
    return Tbi_bj;
}

// Vector3d LidarCalibration::PredictRelativeT(const LidarFrame &lidarframei, const LidarFrame &lidarframej)
// {
//     Vector3d Tbi_bj;
//     Vector3d Tbi = lidarframei.vioData.Pwbi;
//     Vector3d Tbj = lidarframej.vioData.Pwbj;
//     Matrix3d Rbi = lidarframei.vioData.Rwbi;
//     // cout << "!!!!!!!test: Tbi=" << Tbi.transpose() << "  Tbj=" << Tbj.transpose() << endl;
//     Tbi_bj = Rbi.transpose() * Tbj - Rbi.transpose() * Tbi;
//     // cout << "!!!!!!!test:Tbij=" << Tbi_bj.transpose() << endl;
//     return Tbi_bj;
// }

Vector3d LidarCalibration::PredictT(const LidarFrame &lidarframe)
{
    Vector3d Pba(lidarframe.vioData.Pwbi);
    Vector3d Va(lidarframe.vioData.Vbi);
    Vector3d Vb(lidarframe.vioData.Vbj);
    double ta, tl, tb;
    ta = lidarframe.vioData.ti;
    tb = lidarframe.vioData.tj;
    tl = lidarframe.lidarData.point_cloud.time;
    double dt = tl - ta;
    Vector3d Pk;
    if (dt >= 0)
    {
        Vector3d a = (Vb - Va) / (tb - ta);
        Pk = Pba + Va * dt + 0.5 * a * dt * dt;
    }
    else
    {
        Pk = Pba;
    }
    return Pk;
}

bool LidarCalibration::Lidar_align(vector<Matrix3d> &Rlidar, vector<Vector3d> &Tlidar, vector<Vector3d> &Tbody,
                                   Vector3d &Tlb, Matrix3d &Rlb, int frame_count)
{

    //ldlt solver
    // MatrixXd A{3,3};
    // VectorXd b{3};
    // Vector3d tlb;
    // A.setZero();
    // b.setZero();

    // for (int i = 1; i <= frame_count; i++)
    // {
    //     Matrix3d rl_ij = Rlidar[i];
    //     Vector3d tl_ij = Tlidar[i];
    //     Vector3d tb_ij = Tbody[i];
    //     Matrix3d tmp_A(3, 3);
    //     Vector3d tmp_b(3);
    //     tmp_A = Matrix3d::Identity() - rl_ij;
    //     tmp_b = tl_ij - Rlb * tb_ij;
    //     Matrix3d r_A = tmp_A.transpose() * tmp_A;
    //     Vector3d r_b = tmp_A.transpose() * tmp_b;
    //     A += r_A;
    //     b += r_b;
    // }
    // Tlb = A.ldlt().solve(b);

    ceres::Problem problem;

    Quaterniond rlb_q(Rlb);
    double RLB_Q[4];
    RLB_Q[0] = rlb_q.w();
    RLB_Q[1] = rlb_q.x();
    RLB_Q[2] = rlb_q.y();
    RLB_Q[3] = rlb_q.z();
    // RLB_Q[0] = 1.0;
    // RLB_Q[1] = 0.0;
    // RLB_Q[2] = 0.0;
    // RLB_Q[3] = 0.0;
    ceres::LocalParameterization *local_parameterization = new ceres::QuaternionParameterization();
    problem.AddParameterBlock(RLB_Q, 4, local_parameterization);

    double TLB[3];
    TLB[0] = 0.0;
    TLB[1] = 0.0;
    TLB[2] = 0.0;
    problem.AddParameterBlock(TLB, 3);
    problem.SetParameterUpperBound(TLB, 0, PLB_UX);
    problem.SetParameterLowerBound(TLB, 0, PLB_LX);
    problem.SetParameterUpperBound(TLB, 1, PLB_UY);
    problem.SetParameterLowerBound(TLB, 1, PLB_LY);
    problem.SetParameterUpperBound(TLB, 2, PLB_UZ);
    problem.SetParameterLowerBound(TLB, 2, PLB_LZ);

    ceres::LossFunction *loss_function;
    loss_function = new ceres::HuberLoss(1.0);
    // loss_function = new ceres::CauchyLoss(1.0);
    // loss_function = new ceres::CauchyLoss(2.3849);

    for (int i = 1; i <= frame_count; i++)
    {
        Matrix3d rl_ij = Rlidar[i];
        Vector3d tl_ij = Tlidar[i];
        Vector3d tb_ij = Tbody[i];
        Matrix3d rb_ij = Rcam[i];
        ceres::CostFunction *cost_function = LidarInitConstraint::makeConstraint(rl_ij, tl_ij, tb_ij, rb_ij);
        problem.AddResidualBlock(cost_function, loss_function, RLB_Q, TLB);
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR; //DENSE_NORMAL_CHOLESKY  DENSE_SCHUR  DENSE_QR
    // options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG; //LEVENBERG_MARQUARDT  DOGLEG
    options.max_solver_time_in_seconds = 0.1;
    options.max_num_iterations = 60;
    //options.line_search_direction_type = ceres::BFGS;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout << summary.BriefReport() << endl;
    if (summary.final_cost > 1.0)
    {
        return false;
    }
    Quaterniond RLB_;
    RLB_.w() = RLB_Q[0];
    RLB_.x() = RLB_Q[1];
    RLB_.y() = RLB_Q[2];
    RLB_.z() = RLB_Q[3];
    RLB_.normalized();
    Rlb = RLB_.toRotationMatrix();
    Tlb[0] = TLB[0];
    Tlb[1] = TLB[1];
    Tlb[2] = TLB[2];

    cout << "Rlb_refine=" << Utility::R2ypr(Rlb).transpose() << endl;
    cout << "Tlb=" << Tlb.transpose() << endl;

    return true;
}

bool VelocityData::PredictData(const LidarFrame &lidarframe)
{
    Matrix3d Ra(lidarframe.vioData.Rwbi.transpose());
    Matrix3d Rb(lidarframe.vioData.Rwbj.transpose());
    Vector3d Va(Ra * lidarframe.vioData.Vbi);
    Vector3d Vb(Rb * lidarframe.vioData.Vbj);
    Vector3d Aa(lidarframe.vioData.gyri);
    Vector3d Ab(lidarframe.vioData.gyrj);
    double ta, tl, tb;
    ta = lidarframe.vioData.ti;
    tb = lidarframe.vioData.tj;
    tl = lidarframe.lidarData.point_cloud.time;
    double dt = tl - ta;
    Vector3d Vk;
    Vector3d Ak;
    if (dt >= 0)
    {
        Vector3d a = (Vb - Va) / (tb - ta);
        Vk = Va + a * dt;
        Vector3d a_ = (Ab - Aa) / (tb - ta);
        Ak = Aa + a_ * dt;
    }
    else
    {
        Vk = Ra * Va;
        Ak = Aa;
    }

    linear_velocity.x = Vk[0];
    linear_velocity.y = Vk[1];
    linear_velocity.z = Vk[2];
    angular_velocity.x = Ak[0];
    angular_velocity.y = Ak[1];
    angular_velocity.z = Ak[2];

    return true;
}

void VelocityData::TransformCoordinate(Eigen::Matrix3d rotation, Eigen::Vector3d translation)
{
    //Eigen::Matrix4d matrix = transform_matrix.cast<double>();
    Eigen::Matrix3d t_R(rotation);
    Eigen::Vector3d w(angular_velocity.x, angular_velocity.y, angular_velocity.z);
    Eigen::Vector3d v(linear_velocity.x, linear_velocity.y, linear_velocity.z);
    w = t_R * w;

    // v = t_R * v;

    Eigen::Vector3d r(translation);
    Eigen::Vector3d delta_v;
    delta_v(0) = w(1) * r(2) - w(2) * r(1);
    delta_v(1) = w(2) * r(0) - w(0) * r(2);
    delta_v(2) = w(0) * r(1) - w(1) * r(0);
    v = t_R * v + delta_v;

    angular_velocity.x = w(0);
    angular_velocity.y = w(1);
    angular_velocity.z = w(2);
    linear_velocity.x = v(0);
    linear_velocity.y = v(1);
    linear_velocity.z = v(2);
}

void DistortionAdjust::SetMotionInfo(float scan_period, VelocityData velocity_data)
{
    scan_period_ = scan_period;
    velocity_ << velocity_data.linear_velocity.x, velocity_data.linear_velocity.y, velocity_data.linear_velocity.z;
    angular_rate_ << velocity_data.angular_velocity.x, velocity_data.angular_velocity.y, velocity_data.angular_velocity.z;
}

bool DistortionAdjust::AdjustCloud(CloudData::CLOUD_PTR &input_cloud_ptr, CloudData::CLOUD_PTR &output_cloud_ptr)
{
    CloudData::CLOUD_PTR origin_cloud_ptr(new CloudData::CLOUD(*input_cloud_ptr));
    output_cloud_ptr->points.clear();

    float orientation_space = 2.0 * M_PI;
    float delete_space = 5.0 * M_PI / 180.0;
    float start_orientation = atan2(origin_cloud_ptr->points[0].y, origin_cloud_ptr->points[0].x);

    Eigen::AngleAxisf t_V(start_orientation, Eigen::Vector3f::UnitZ());
    Eigen::Matrix3f rotate_matrix = t_V.matrix();
    Eigen::Matrix4f transform_matrix = Eigen::Matrix4f::Identity();
    transform_matrix.block<3, 3>(0, 0) = rotate_matrix.inverse();
    pcl::transformPointCloud(*origin_cloud_ptr, *origin_cloud_ptr, transform_matrix);

    velocity_ = rotate_matrix * velocity_;
    angular_rate_ = rotate_matrix * angular_rate_;

    for (size_t point_index = 1; point_index < origin_cloud_ptr->points.size(); ++point_index)
    {
        float orientation = atan2(origin_cloud_ptr->points[point_index].y, origin_cloud_ptr->points[point_index].x);
        if (orientation < 0.0)
            orientation += 2.0 * M_PI;

        if (orientation < delete_space || 2.0 * M_PI - orientation < delete_space)
            continue;

        float real_time = fabs(orientation) / orientation_space * scan_period_ - scan_period_ / 2.0;

        Eigen::Vector3f origin_point(origin_cloud_ptr->points[point_index].x,
                                     origin_cloud_ptr->points[point_index].y,
                                     origin_cloud_ptr->points[point_index].z);

        Eigen::Matrix3f current_matrix = UpdateMatrix(real_time);
        Eigen::Vector3f rotated_point = current_matrix * origin_point;
        Eigen::Vector3f adjusted_point = rotated_point + velocity_ * real_time;
        CloudData::POINT point;
        point.x = adjusted_point(0);
        point.y = adjusted_point(1);
        point.z = adjusted_point(2);
        output_cloud_ptr->points.push_back(point);
    }

    pcl::transformPointCloud(*output_cloud_ptr, *output_cloud_ptr, transform_matrix.inverse());
    return true;
}

Eigen::Matrix3f DistortionAdjust::UpdateMatrix(float real_time)
{
    Eigen::Vector3f angle = angular_rate_ * real_time;
    Eigen::AngleAxisf t_Vz(angle(2), Eigen::Vector3f::UnitZ());
    Eigen::AngleAxisf t_Vy(angle(1), Eigen::Vector3f::UnitY());
    Eigen::AngleAxisf t_Vx(angle(0), Eigen::Vector3f::UnitX());
    Eigen::AngleAxisf t_V;
    t_V = t_Vz * t_Vy * t_Vx;
    return t_V.matrix();
}

void calculate_ICP_COV(CloudData::CLOUD &data_pi, CloudData::CLOUD &model_qi, Eigen::Matrix4f transform, Eigen::MatrixXd &ICP_COV)
{

    double Tx = transform(0, 3);
    double Ty = transform(1, 3);
    double Tz = transform(2, 3);

    // Eigen::Matrix3d R = transform.block<3, 3>(0, 0).cast<double>();
    // Eigen::Vector3d nn = R.col(0);
    // Eigen::Vector3d oo = R.col(1);
    // Eigen::Vector3d aa = R.col(2);

    // double yaw = atan2(nn(1), nn(0));
    // double pitch = atan2(-nn(2), nn(0) * cos(yaw) + nn(1) * sin(yaw));
    // double roll = atan2(aa(0) * sin(yaw) - aa(1) * cos(yaw), -oo(0) * sin(yaw) + oo(1) * cos(yaw));

    double roll = atan2f(transform(2, 1), transform(2, 2));
    double pitch = asinf(-transform(2, 0));
    double yaw = atan2f(transform(1, 0), transform(0, 0));

    double x, y, z, a, b, c;
    x = Tx;
    y = Ty;
    z = Tz;
    a = yaw;
    b = pitch;
    c = roll; // important // According to the rotation matrix I used and after verification, it is Yaw Pitch ROLL = [a,b,c]== [R] matrix used in the MatLab also :)

    /* Flushing out in the form of XYZ ABC */
    //std::cout << "\nPrinting out [x, y, z, a, b, c] =  " <<x<<"    "<<y<<"    "<<z<<"    "<<a<<"    "<<b<<"    "<<c<<std::endl;

    //Matrix initialization
    Eigen::MatrixXd d2J_dX2(6, 6);
    d2J_dX2 = Eigen::MatrixXd::Zero(6, 6);

    /****  Calculating d2J_dX2  ****/
    for (size_t s = 0; s < data_pi.points.size(); ++s)
    {
        double pix = data_pi.points[s].x;
        double piy = data_pi.points[s].y;
        double piz = data_pi.points[s].z;
        double qix = model_qi.points[s].x;
        double qiy = model_qi.points[s].y;
        double qiz = model_qi.points[s].z;

        /************************************************************

d2J_dX2 -- X is the [R|T] in the form of [x, y, z, a, b, c]
x, y, z is the translation part 
a, b, c is the rotation part in Euler format
[x, y, z, a, b, c] is acquired from the Transformation Matrix returned by ICP.

Now d2J_dX2 is a 6x6 matrix of the form

d2J_dx2                   
d2J_dxdy    d2J_dy2
d2J_dxdz    d2J_dydz    d2J_dz2
d2J_dxda    d2J_dyda    d2J_dzda   d2J_da2
d2J_dxdb    d2J_dydb    d2J_dzdb   d2J_dadb   d2J_db2
d2J_dxdc    d2J_dydc    d2J_dzdc   d2J_dadc   d2J_dbdc   d2J_dc2


*************************************************************/

        double d2J_dx2, d2J_dydx, d2J_dzdx, d2J_dadx, d2J_dbdx, d2J_dcdx,
            d2J_dxdy, d2J_dy2, d2J_dzdy, d2J_dady, d2J_dbdy, d2J_dcdy,
            d2J_dxdz, d2J_dydz, d2J_dz2, d2J_dadz, d2J_dbdz, d2J_dcdz,
            d2J_dxda, d2J_dyda, d2J_dzda, d2J_da2, d2J_dbda, d2J_dcda,
            d2J_dxdb, d2J_dydb, d2J_dzdb, d2J_dadb, d2J_db2, d2J_dcdb,
            d2J_dxdc, d2J_dydc, d2J_dzdc, d2J_dadc, d2J_dbdc, d2J_dc2;

        // These terms are generated from the provided Matlab scipts. We just have to copy
        // the expressions from the matlab output with two very simple changes.
        // The first one being the the sqaure of a number 'a' is shown as a^2 in matlab,
        // which is converted to pow(a,2) in the below expressions.
        // The second change is to add ';' at the end of each expression :)
        // In this way, matlab can be used to generate these terms for various objective functions of ICP
        // and they can simply be copied to the C++ files and with appropriate changes to ICP estimation,
        // its covariance can be easily estimated.

        d2J_dx2 =

            2;

        d2J_dy2 =

            2;

        d2J_dz2 =

            2;

        d2J_dydx =

            0;

        d2J_dxdy =

            0;

        d2J_dzdx =

            0;

        d2J_dxdz =

            0;

        d2J_dydz =

            0;

        d2J_dzdy =

            0;

        d2J_da2 =

            (piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + pix * cos(a) * cos(b)) * (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) - (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) + (piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) * (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) - (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b));

        d2J_db2 =

            (pix * cos(b) + piz * cos(c) * sin(b) + piy * sin(b) * sin(c)) * (2 * pix * cos(b) + 2 * piz * cos(c) * sin(b) + 2 * piy * sin(b) * sin(c)) - (2 * piz * cos(b) * cos(c) - 2 * pix * sin(b) + 2 * piy * cos(b) * sin(c)) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) - (2 * pix * cos(a) * cos(b) + 2 * piz * cos(a) * cos(c) * sin(b) + 2 * piy * cos(a) * sin(b) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + (piz * cos(a) * cos(b) * cos(c) - pix * cos(a) * sin(b) + piy * cos(a) * cos(b) * sin(c)) * (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) - (2 * pix * cos(b) * sin(a) + 2 * piz * cos(c) * sin(a) * sin(b) + 2 * piy * sin(a) * sin(b) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) + (piz * cos(b) * cos(c) * sin(a) - pix * sin(a) * sin(b) + piy * cos(b) * sin(a) * sin(c)) * (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c));

        d2J_dc2 =

            (piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) + (piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) - (2 * piz * cos(b) * cos(c) + 2 * piy * cos(b) * sin(c)) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) + (2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) - 2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b))) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + (piy * cos(b) * cos(c) - piz * cos(b) * sin(c)) * (2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c)) - (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b))) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dxda =

            2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) - 2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * pix * cos(b) * sin(a);

        d2J_dadx =

            2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) - 2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * pix * cos(b) * sin(a);

        d2J_dyda =

            2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b);

        d2J_dady =

            2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b);

        d2J_dzda =

            0;

        d2J_dadz =

            0;

        d2J_dxdb =

            2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c);

        d2J_dbdx =

            2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c);

        d2J_dydb =

            2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c);

        d2J_dbdy =

            2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c);

        d2J_dzdb =

            -2 * pix * cos(b) - 2 * piz * cos(c) * sin(b) - 2 * piy * sin(b) * sin(c);

        d2J_dbdz =

            -2 * pix * cos(b) - 2 * piz * cos(c) * sin(b) - 2 * piy * sin(b) * sin(c);

        d2J_dxdc =

            2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c));

        d2J_dcdx =

            2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c));

        d2J_dydc =

            -2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) - 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c));

        d2J_dcdy =

            -2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) - 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c));

        d2J_dzdc =

            2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c);

        d2J_dcdz =

            2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c);

        d2J_dadb =

            (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c)) * (piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + pix * cos(a) * cos(b)) - (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) * (piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) + (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) - (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b));

        d2J_dbda =

            (piz * cos(b) * cos(c) * sin(a) - pix * sin(a) * sin(b) + piy * cos(b) * sin(a) * sin(c)) * (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) - (piz * cos(a) * cos(b) * cos(c) - pix * cos(a) * sin(b) + piy * cos(a) * cos(b) * sin(c)) * (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) + (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) - (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b));

        d2J_dbdc =

            (2 * piy * cos(a) * cos(b) * cos(c) - 2 * piz * cos(a) * cos(b) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + (2 * piy * cos(b) * cos(c) * sin(a) - 2 * piz * cos(b) * sin(a) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) - (2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c)) * (pix * cos(b) + piz * cos(c) * sin(b) + piy * sin(b) * sin(c)) - (2 * piy * cos(c) * sin(b) - 2 * piz * sin(b) * sin(c)) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) + (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (piz * cos(a) * cos(b) * cos(c) - pix * cos(a) * sin(b) + piy * cos(a) * cos(b) * sin(c)) - (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (piz * cos(b) * cos(c) * sin(a) - pix * sin(a) * sin(b) + piy * cos(b) * sin(a) * sin(c));

        d2J_dcdb =

            (2 * piy * cos(a) * cos(b) * cos(c) - 2 * piz * cos(a) * cos(b) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + (2 * piy * cos(b) * cos(c) * sin(a) - 2 * piz * cos(b) * sin(a) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) - (piy * cos(b) * cos(c) - piz * cos(b) * sin(c)) * (2 * pix * cos(b) + 2 * piz * cos(c) * sin(b) + 2 * piy * sin(b) * sin(c)) - (2 * piy * cos(c) * sin(b) - 2 * piz * sin(b) * sin(c)) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) + (piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) - (piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c));

        d2J_dcda =

            (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) - (piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) - (piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) + (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dadc =

            (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) - (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) - (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) * (piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + pix * cos(a) * cos(b)) + (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        Eigen::MatrixXd d2J_dX2_temp(6, 6);

        d2J_dX2_temp << d2J_dx2, d2J_dydx, d2J_dzdx, d2J_dadx, d2J_dbdx, d2J_dcdx,
            d2J_dxdy, d2J_dy2, d2J_dzdy, d2J_dady, d2J_dbdy, d2J_dcdy,
            d2J_dxdz, d2J_dydz, d2J_dz2, d2J_dadz, d2J_dbdz, d2J_dcdz,
            d2J_dxda, d2J_dyda, d2J_dzda, d2J_da2, d2J_dbda, d2J_dcda,
            d2J_dxdb, d2J_dydb, d2J_dzdb, d2J_dadb, d2J_db2, d2J_dcdb,
            d2J_dxdc, d2J_dydc, d2J_dzdc, d2J_dadc, d2J_dbdc, d2J_dc2;

        d2J_dX2 = d2J_dX2 + d2J_dX2_temp;

    } // End of the FOR loop!!!

    //std::cout << "\n**************\n Successfully Computed d2J_dX2 \n**************\n" << std::endl;

    // Now its time to calculate d2J_dZdX , where Z are the measurements Pi and Qi, X = [x,y,z,a,b,c]

    // n is the number of correspondences
    int n = data_pi.points.size();

    /*  Here we check if the number of correspondences between the source and the target point clouds are greater than 200.
      if yes, we only take the first 200 correspondences to calculate the covariance matrix
      You can try increasing it but if its too high, the system may run out of memory and give an exception saying

     terminate called after throwing an instance of 'std::bad_alloc'
     what():  std::bad_alloc
     Aborted (core dumped)

  */

    if (n > 50)
        n = 50; ////////////****************************IMPORTANT CHANGE********but may not affect*****************/////////////////////////////////////////

    //std::cout << "\nNumber of Correspondences used for ICP's covariance estimation = " << n << std::endl;

    Eigen::MatrixXd d2J_dZdX(6, 6 * n);

    for (int k = 0; k < n; ++k) // row
    {
        //here the current correspondences are loaded into Pi and Qi
        double pix = data_pi.points[k].x;
        double piy = data_pi.points[k].y;
        double piz = data_pi.points[k].z;
        double qix = model_qi.points[k].x;
        double qiy = model_qi.points[k].y;
        double qiz = model_qi.points[k].z;

        Eigen::MatrixXd d2J_dZdX_temp(6, 6);

        double d2J_dpix_dx, d2J_dpiy_dx, d2J_dpiz_dx, d2J_dqix_dx, d2J_dqiy_dx, d2J_dqiz_dx,
            d2J_dpix_dy, d2J_dpiy_dy, d2J_dpiz_dy, d2J_dqix_dy, d2J_dqiy_dy, d2J_dqiz_dy,
            d2J_dpix_dz, d2J_dpiy_dz, d2J_dpiz_dz, d2J_dqix_dz, d2J_dqiy_dz, d2J_dqiz_dz,
            d2J_dpix_da, d2J_dpiy_da, d2J_dpiz_da, d2J_dqix_da, d2J_dqiy_da, d2J_dqiz_da,
            d2J_dpix_db, d2J_dpiy_db, d2J_dpiz_db, d2J_dqix_db, d2J_dqiy_db, d2J_dqiz_db,
            d2J_dpix_dc, d2J_dpiy_dc, d2J_dpiz_dc, d2J_dqix_dc, d2J_dqiy_dc, d2J_dqiz_dc;

        d2J_dpix_dx =

            2 * cos(a) * cos(b);

        d2J_dpix_dy =

            2 * cos(b) * sin(a);

        d2J_dpix_dz =

            -2 * sin(b);

        d2J_dpix_da =

            cos(b) * sin(a) * (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) - cos(a) * cos(b) * (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) - 2 * cos(b) * sin(a) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + 2 * cos(a) * cos(b) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dpix_db =

            sin(b) * (2 * pix * cos(b) + 2 * piz * cos(c) * sin(b) + 2 * piy * sin(b) * sin(c)) - 2 * cos(b) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) + cos(a) * cos(b) * (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) - 2 * sin(a) * sin(b) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) + cos(b) * sin(a) * (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c)) - 2 * cos(a) * sin(b) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b));

        d2J_dpix_dc =

            cos(a) * cos(b) * (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) - sin(b) * (2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c)) - cos(b) * sin(a) * (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)));

        d2J_dpiy_dx =

            2 * cos(a) * sin(b) * sin(c) - 2 * cos(c) * sin(a);

        d2J_dpiy_dy =

            2 * cos(a) * cos(c) + 2 * sin(a) * sin(b) * sin(c);

        d2J_dpiy_dz =

            2 * cos(b) * sin(c);

        d2J_dpiy_da =

            (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) * (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) + (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) * (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) - (2 * cos(a) * cos(c) + 2 * sin(a) * sin(b) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) - (2 * cos(c) * sin(a) - 2 * cos(a) * sin(b) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dpiy_db =

            (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) * (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c)) - (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) * (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) - 2 * sin(b) * sin(c) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) - cos(b) * sin(c) * (2 * pix * cos(b) + 2 * piz * cos(c) * sin(b) + 2 * piy * sin(b) * sin(c)) + 2 * cos(a) * cos(b) * sin(c) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + 2 * cos(b) * sin(a) * sin(c) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dpiy_dc =

            (2 * sin(a) * sin(c) + 2 * cos(a) * cos(c) * sin(b)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) - (2 * cos(a) * sin(c) - 2 * cos(c) * sin(a) * sin(b)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) - (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) * (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) - (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) * (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) + 2 * cos(b) * cos(c) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) + cos(b) * sin(c) * (2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c));

        d2J_dpiz_dx =

            2 * sin(a) * sin(c) + 2 * cos(a) * cos(c) * sin(b);

        d2J_dpiz_dy =

            2 * cos(c) * sin(a) * sin(b) - 2 * cos(a) * sin(c);

        d2J_dpiz_dz =

            2 * cos(b) * cos(c);

        d2J_dpiz_da =

            (2 * cos(a) * sin(c) - 2 * cos(c) * sin(a) * sin(b)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) - (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) * (2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a)) - (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) * (2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + 2 * pix * cos(a) * cos(b)) + (2 * sin(a) * sin(c) + 2 * cos(a) * cos(c) * sin(b)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dpiz_db =

            (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) * (2 * piz * cos(a) * cos(b) * cos(c) - 2 * pix * cos(a) * sin(b) + 2 * piy * cos(a) * cos(b) * sin(c)) - (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) * (2 * piz * cos(b) * cos(c) * sin(a) - 2 * pix * sin(a) * sin(b) + 2 * piy * cos(b) * sin(a) * sin(c)) - 2 * cos(c) * sin(b) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c)) - cos(b) * cos(c) * (2 * pix * cos(b) + 2 * piz * cos(c) * sin(b) + 2 * piy * sin(b) * sin(c)) + 2 * cos(a) * cos(b) * cos(c) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) + 2 * cos(b) * cos(c) * sin(a) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a));

        d2J_dpiz_dc =

            (2 * cos(c) * sin(a) - 2 * cos(a) * sin(b) * sin(c)) * (x - qix - piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) + piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + pix * cos(a) * cos(b)) - (2 * cos(a) * cos(c) + 2 * sin(a) * sin(b) * sin(c)) * (y - qiy + piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + pix * cos(b) * sin(a)) + (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) * (2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) + 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c))) + (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) * (2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c))) + cos(b) * cos(c) * (2 * piy * cos(b) * cos(c) - 2 * piz * cos(b) * sin(c)) - 2 * cos(b) * sin(c) * (z - qiz - pix * sin(b) + piz * cos(b) * cos(c) + piy * cos(b) * sin(c));

        d2J_dqix_dx =

            -2;

        d2J_dqix_dy =

            0;

        d2J_dqix_dz =

            0;

        d2J_dqix_da =

            2 * piy * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c)) - 2 * piz * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * pix * cos(b) * sin(a);

        d2J_dqix_db =

            2 * pix * cos(a) * sin(b) - 2 * piz * cos(a) * cos(b) * cos(c) - 2 * piy * cos(a) * cos(b) * sin(c);

        d2J_dqix_dc =

            -2 * piy * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * piz * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c));

        d2J_dqiy_dx =

            0;

        d2J_dqiy_dy =

            -2;

        d2J_dqiy_dz =

            0;

        d2J_dqiy_da =

            2 * piy * (cos(c) * sin(a) - cos(a) * sin(b) * sin(c)) - 2 * piz * (sin(a) * sin(c) + cos(a) * cos(c) * sin(b)) - 2 * pix * cos(a) * cos(b);

        d2J_dqiy_db =

            2 * pix * sin(a) * sin(b) - 2 * piz * cos(b) * cos(c) * sin(a) - 2 * piy * cos(b) * sin(a) * sin(c);

        d2J_dqiy_dc =

            2 * piy * (cos(a) * sin(c) - cos(c) * sin(a) * sin(b)) + 2 * piz * (cos(a) * cos(c) + sin(a) * sin(b) * sin(c));

        d2J_dqiz_dx =

            0;

        d2J_dqiz_dy =

            0;

        d2J_dqiz_dz =

            -2;

        d2J_dqiz_da =

            0;

        d2J_dqiz_db =

            2 * pix * cos(b) + 2 * piz * cos(c) * sin(b) + 2 * piy * sin(b) * sin(c);

        d2J_dqiz_dc =

            2 * piz * cos(b) * sin(c) - 2 * piy * cos(b) * cos(c);

        d2J_dZdX_temp << d2J_dpix_dx, d2J_dpiy_dx, d2J_dpiz_dx, d2J_dqix_dx, d2J_dqiy_dx, d2J_dqiz_dx,
            d2J_dpix_dy, d2J_dpiy_dy, d2J_dpiz_dy, d2J_dqix_dy, d2J_dqiy_dy, d2J_dqiz_dy,
            d2J_dpix_dz, d2J_dpiy_dz, d2J_dpiz_dz, d2J_dqix_dz, d2J_dqiy_dz, d2J_dqiz_dz,
            d2J_dpix_da, d2J_dpiy_da, d2J_dpiz_da, d2J_dqix_da, d2J_dqiy_da, d2J_dqiz_da,
            d2J_dpix_db, d2J_dpiy_db, d2J_dpiz_db, d2J_dqix_db, d2J_dqiy_db, d2J_dqiz_db,
            d2J_dpix_dc, d2J_dpiy_dc, d2J_dpiz_dc, d2J_dqix_dc, d2J_dqiy_dc, d2J_dqiz_dc;

        d2J_dZdX.block<6, 6>(0, 6 * k) = d2J_dZdX_temp;
    }

    //By reaching here both the matrices d2J_dX2 and d2J_dZdX are calculated and lets print those values out;

    //std::cout << "\n Finally here are the results \n\n" << "d2J_dX2 = \n " << d2J_dX2 <<std::endl;
    //std::cout << "\n\n\n" << "d2J_dZdX = \n " << d2J_dZdX <<std::endl;

    //By reaching here both the matrices d2J_dX2 and d2J_dZdX are calculated and lets print those values out;

    //std::cout << "\n Finally here are the two matrices \n\n" << "d2J_dX2 = \n " << d2J_dX2 <<std::endl;
    //std::cout << "\n\n\n" << "d2J_dZdX = \n " << d2J_dZdX <<std::endl;

    /**************************************
     *
     * Here we create the matrix cov(z) as mentioned in Section 3.3 in the paper, "Covariance of ICP with 3D Point to Point and Point to Plane Error Metrics"
     *
     * ************************************/

    Eigen::MatrixXd cov_z(6 * n, 6 * n);
    cov_z = 0.01 * Eigen::MatrixXd::Identity(6 * n, 6 * n);

    ICP_COV = d2J_dX2.inverse() * d2J_dZdX * cov_z * d2J_dZdX.transpose() * d2J_dX2.inverse();

    //std::cout << "\n\n********************** \n\n" << "ICP_COV = \n" << ICP_COV <<"\n*******************\n\n"<< std::endl;

    //std::cout << "\nSuccessfully Computed the ICP's Covariance !!!\n" << std::endl;
}

Matrix4d PredictRelative_rt(LidarFrame &lidarframei, LidarFrame &lidarframej, const Matrix3d RBL, const Vector3d TBL)
{
    Matrix3d Rbi_bj;
    Matrix3d Rbj = Predict_r(lidarframej);
    Matrix3d Rbi = Predict_r(lidarframei);
    Rbi_bj = Rbi.inverse() * Rbj;
    Vector3d Tbi_bj;
    Vector3d Tbj = Predict_t(lidarframej);
    Vector3d Tbi = Predict_t(lidarframei);
    Tbi_bj = Rbi.transpose() * Tbj - Rbi.transpose() * Tbi;
    Matrix4d Lij(Matrix4d::Identity());
    Lij.block<3, 3>(0, 0) = Rbi_bj;
    Lij.block<3, 1>(0, 3) = Tbi_bj;
    lidarframei.lidarData.lidar_R = Rbi * RBL;
    lidarframej.lidarData.lidar_R = Rbj * RBL;
    lidarframei.lidarData.lidar_T = Tbi + Rbi * TBL;
    lidarframej.lidarData.lidar_T = Tbj + Rbj * TBL;
    return Lij;
}

Matrix3d Predict_r(const LidarFrame &lidarframe)
{
    Quaterniond qc0ca(lidarframe.vioData.Rwbi);
    Quaterniond qc0cb(lidarframe.vioData.Rwbj);
    double ta, tl, tb;
    ta = lidarframe.vioData.ti;
    tb = lidarframe.vioData.tj;
    tl = lidarframe.lidarData.point_cloud.time;
    double t = (tl - ta) / (tb - ta);
    cout << "t=" << t << endl;
    Quaterniond q;
    if (t > 0)
    {
        q = qc0ca.slerp(t, qc0cb);
        q.normalized();
    }
    else
    {
        q = qc0ca;
        q.normalized();
    }

    return q.toRotationMatrix();
}

Vector3d Predict_t(const LidarFrame &lidarframe)
{
    Vector3d Pba(lidarframe.vioData.Pwbi);
    Vector3d Va(lidarframe.vioData.Vbi);
    Vector3d Vb(lidarframe.vioData.Vbj);
    double ta, tl, tb;
    ta = lidarframe.vioData.ti;
    tb = lidarframe.vioData.tj;
    tl = lidarframe.lidarData.point_cloud.time;
    double dt = tl - ta;
    Vector3d Pk;
    if (dt >= 0)
    {
        Vector3d a = (Vb - Va) / (tb - ta);
        Pk = Pba + Va * dt + 0.5 * a * dt * dt;
    }
    else
    {
        Pk = Pba;
    }
    return Pk;
}

void RotatePoint(const Eigen::Quaternionf &q, CloudData::POINT &p)
{
    Eigen::Vector3f vec, vec_out;
    vec.x() = p.x;
    vec.y() = p.y;
    vec.z() = p.z;
    vec_out = q * vec;
    p.x = vec_out.x();
    p.y = vec_out.y();
    p.z = vec_out.z();
}

void TransformToEnd(CloudData::CLOUD_PTR &cloud, Eigen::Quaternionf transform_es_q, Eigen::Vector3f transform_es_t, float time_factor, double min, double max)
{
    size_t cloud_size = cloud->points.size();
#pragma omp parallel for num_threads(2)
    for (size_t i = 0; i < cloud_size; i++)
    {
        CloudData::POINT &point = cloud->points[i];
        double distance = sqrt(point.x * point.x + point.y * point.y);
        float s = time_factor * (point.intensity - int(point.intensity));

        if (s < 0 || s > 1.001 || distance < min || distance > max) //s < 0.005 || s > 0.995
        {
            point.x = std::numeric_limits<float>::quiet_NaN();
            point.y = std::numeric_limits<float>::quiet_NaN();
            point.z = std::numeric_limits<float>::quiet_NaN();
            // cout << "POINT SET 2 NaN!!!!!!!!!!!! " << endl;
            continue;
        }
        else
        {
            point.x -= s * transform_es_t.x();
            point.y -= s * transform_es_t.y();
            point.z -= s * transform_es_t.z();
        }
        // cout << "trans2end scale = " << s << endl;
        Eigen::Quaternionf q_id, q_s, q_e;
        q_e = transform_es_q;
        q_id.setIdentity();
        q_s = q_id.slerp(s, q_e);
        // RotatePoint(q_s.inverse().normalized(), point);
        RotatePoint(q_s.conjugate().normalized(), point);

        //Trans2End
        RotatePoint(q_e, point);
        point.x += transform_es_t.x();
        point.y += transform_es_t.y();
        point.z += transform_es_t.z();
        point.intensity = int(point.intensity);
        // cout << "point intensity = " << point.intensity << endl;
    }
}
