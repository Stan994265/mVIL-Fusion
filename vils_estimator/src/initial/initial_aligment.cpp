#include "initial_alignment.h"

void solveGyroscopeBias(map<double, ImageFrame> &all_image_frame, Vector3d *Bgs)
{
    Matrix3d A;
    Vector3d b;
    Vector3d delta_bg;
    A.setZero();
    b.setZero();
    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    // for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    // {
    //     frame_j = next(frame_i);
    //     MatrixXd tmp_A(3, 3);
    //     tmp_A.setZero();
    //     VectorXd tmp_b(3);
    //     tmp_b.setZero();
    //     Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
    //     tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
    //     tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
    //     A += tmp_A.transpose() * tmp_A;
    //     b += tmp_A.transpose() * tmp_b;

    // }
    // delta_bg = A.ldlt().solve(b);
    // ROS_WARN_STREAM("gyroscope bias initial calibration " << delta_bg.transpose());

    // for (int i = 0; i <= WINDOW_SIZE; i++)
    //     Bgs[i] += delta_bg;

    // for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end( ); frame_i++)
    // {
    //     frame_j = next(frame_i);
    //     frame_j->second.pre_integration->repropagate(Vector3d::Zero(), Bgs[0]);
    // }
    vector<Vector3d> dbg;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++)
    {
        frame_j = next(frame_i);
        MatrixXd tmp_A(3, 3);
        tmp_A.setZero();
        VectorXd tmp_b(3);
        tmp_b.setZero();
        Eigen::Quaterniond q_ij(frame_i->second.R.transpose() * frame_j->second.R);
        tmp_A = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        tmp_b = 2 * (frame_j->second.pre_integration->delta_q.inverse() * q_ij).vec();
        A = tmp_A.transpose() * tmp_A;
        b = tmp_A.transpose() * tmp_b;
        delta_bg = A.ldlt().solve(b);
        dbg.push_back(delta_bg);
        frame_j->second.pre_integration->init_refine_delta_pvq_bgs(delta_bg);
    }

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {

        cout << "Bgs" << i << "=" << Bgs[i].transpose() << endl;
        Bgs[i] += dbg[i];
        cout << "dbg" << i << "=" << dbg[i].transpose() << endl;
    }
}

MatrixXd TangentBasis(Vector3d &g0)
{
    Vector3d b, c;
    Vector3d a = g0.normalized();
    Vector3d tmp(0, 0, 1);
    if (a == tmp)
        tmp << 1, 0, 0;
    b = (tmp - a * (a.transpose() * tmp)).normalized();
    c = a.cross(b);
    MatrixXd bc(3, 2);
    bc.block<3, 1>(0, 0) = b;
    bc.block<3, 1>(0, 1) = c;
    return bc;
}

void RefineGravity(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    Vector3d g0 = g.normalized() * G.norm();
    Vector3d lx, ly;
    //VectorXd x;
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 2 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    for (int k = 0; k < 4; k++)
    {
        MatrixXd lxly(3, 2);
        lxly = TangentBasis(g0);
        int i = 0;
        for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
        {
            frame_j = next(frame_i);

            MatrixXd tmp_A(6, 9);
            tmp_A.setZero();
            VectorXd tmp_b(6);
            tmp_b.setZero();

            double dt = frame_j->second.pre_integration->sum_dt;

            tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
            tmp_A.block<3, 2>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity() * lxly;
            tmp_A.block<3, 1>(0, 8) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
            tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0] - frame_i->second.R.transpose() * dt * dt / 2 * g0;

            tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
            tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
            tmp_A.block<3, 2>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity() * lxly;
            tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v - frame_i->second.R.transpose() * dt * Matrix3d::Identity() * g0;

            Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
            //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
            //MatrixXd cov_inv = cov.inverse();
            cov_inv.setIdentity();

            MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
            VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

            A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
            b.segment<6>(i * 3) += r_b.head<6>();

            A.bottomRightCorner<3, 3>() += r_A.bottomRightCorner<3, 3>();
            b.tail<3>() += r_b.tail<3>();

            A.block<6, 3>(i * 3, n_state - 3) += r_A.topRightCorner<6, 3>();
            A.block<3, 6>(n_state - 3, i * 3) += r_A.bottomLeftCorner<3, 6>();
        }
        A = A * 1000.0;
        b = b * 1000.0;
        x = A.ldlt().solve(b);
        VectorXd dg = x.segment<2>(n_state - 3);
        g0 = (g0 + lxly * dg).normalized() * G.norm();
        //double s = x(n_state - 1);
    }
    g = g0;
}

bool LinearAlignment(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x)
{
    int all_frame_count = all_image_frame.size();
    int n_state = all_frame_count * 3 + 3 + 1;

    MatrixXd A{n_state, n_state};
    A.setZero();
    VectorXd b{n_state};
    b.setZero();

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        MatrixXd tmp_A(6, 10);
        tmp_A.setZero();
        VectorXd tmp_b(6);
        tmp_b.setZero();

        double dt = frame_j->second.pre_integration->sum_dt;

        tmp_A.block<3, 3>(0, 0) = -dt * Matrix3d::Identity();
        tmp_A.block<3, 3>(0, 6) = frame_i->second.R.transpose() * dt * dt / 2 * Matrix3d::Identity();
        tmp_A.block<3, 1>(0, 9) = frame_i->second.R.transpose() * (frame_j->second.T - frame_i->second.T) / 100.0;
        tmp_b.block<3, 1>(0, 0) = frame_j->second.pre_integration->delta_p + frame_i->second.R.transpose() * frame_j->second.R * TIC[0] - TIC[0];
        //cout << "delta_p   " << frame_j->second.pre_integration->delta_p.transpose() << endl;
        tmp_A.block<3, 3>(3, 0) = -Matrix3d::Identity();
        tmp_A.block<3, 3>(3, 3) = frame_i->second.R.transpose() * frame_j->second.R;
        tmp_A.block<3, 3>(3, 6) = frame_i->second.R.transpose() * dt * Matrix3d::Identity();
        tmp_b.block<3, 1>(3, 0) = frame_j->second.pre_integration->delta_v;
        //cout << "delta_v   " << frame_j->second.pre_integration->delta_v.transpose() << endl;

        Matrix<double, 6, 6> cov_inv = Matrix<double, 6, 6>::Zero();
        //cov.block<6, 6>(0, 0) = IMU_cov[i + 1];
        //MatrixXd cov_inv = cov.inverse();
        cov_inv.setIdentity();

        MatrixXd r_A = tmp_A.transpose() * cov_inv * tmp_A;
        VectorXd r_b = tmp_A.transpose() * cov_inv * tmp_b;

        A.block<6, 6>(i * 3, i * 3) += r_A.topLeftCorner<6, 6>();
        b.segment<6>(i * 3) += r_b.head<6>();

        A.bottomRightCorner<4, 4>() += r_A.bottomRightCorner<4, 4>();
        b.tail<4>() += r_b.tail<4>();

        A.block<6, 4>(i * 3, n_state - 4) += r_A.topRightCorner<6, 4>();
        A.block<4, 6>(n_state - 4, i * 3) += r_A.bottomLeftCorner<4, 6>();
    }
    A = A * 1000.0;
    b = b * 1000.0;
    x = A.ldlt().solve(b);
    double s = x(n_state - 1) / 100.0;
    ROS_DEBUG("estimated scale: %f", s);
    g = x.segment<3>(n_state - 4);
    ROS_DEBUG_STREAM(" result g     " << g.norm() << " " << g.transpose());
    if (fabs(g.norm() - G.norm()) > 1.0 || s < 0)
    {
        return false;
    }
    cout << "before refine g" << g.norm() << " " << g.transpose() << endl;
    RefineGravity(all_image_frame, g, x);
    s = (x.tail<1>())(0) / 100.0;
    (x.tail<1>())(0) = s;
    ROS_DEBUG_STREAM(" refine     " << g.norm() << " " << g.transpose());
    if (s < 0.0)
        return false;
    else
        return true;
}

bool Estimate_ric_td_bg(std::map<double, ImageFrame> &all_image_frame, Eigen::Matrix3d &RIC, double *init_dt, Vector3d *Bgs, Eigen::Matrix3d ric)
{

    ceres::Problem problem;
    int frame_num(all_image_frame.size());
    Quaterniond ric_q(RIC);
    double RIC_Q[4];
    RIC_Q[0] = ric_q.w();
    RIC_Q[1] = ric_q.x();
    RIC_Q[2] = ric_q.y();
    RIC_Q[3] = ric_q.z();
    ceres::LocalParameterization *local_parameterization = new ceres::QuaternionParameterization();
    problem.AddParameterBlock(RIC_Q, 4, local_parameterization);

    if (ESTIMATE_EXTRINSIC == 0)
    {
        problem.SetParameterBlockConstant(RIC_Q);
    }

    double BGS[frame_num][3];
    double Td[frame_num][1];
    for (int i = 0; i < frame_num; i++)
    {
        BGS[i][0] = 0.0;
        BGS[i][1] = 0.0;
        BGS[i][2] = 0.0;
        Td[i][0] = 0.0;
        problem.AddParameterBlock(BGS[i], 3);
        problem.AddParameterBlock(Td[i], 1);
        problem.SetParameterUpperBound(Td[i], 0, 0.1);
        problem.SetParameterLowerBound(Td[i], 0, -0.1);
        problem.SetParameterUpperBound(BGS[i], 0, 0.1);
        problem.SetParameterLowerBound(BGS[i], 0, -0.1);
        problem.SetParameterUpperBound(BGS[i], 1, 0.1);
        problem.SetParameterLowerBound(BGS[i], 1, -0.1);
        problem.SetParameterUpperBound(BGS[i], 2, 0.1);
        problem.SetParameterLowerBound(BGS[i], 2, -0.1);
    }

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int i = 0;

    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, i++)
    {
        frame_j = next(frame_i);

        Vector3d Wj(frame_j->second.pre_integration->gyr_1);
        Vector3d Wi(frame_j->second.pre_integration->gyr_0);
        Quaterniond qcjc0(frame_j->second.R.transpose());
        Quaterniond qc0ci(frame_i->second.R);
        Quaterniond q_bibj(frame_j->second.pre_integration->delta_q);
        Matrix3d J = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_R, O_BG);
        ceres::CostFunction *cost_function = InitRotationConstraint::makeConstraint_r(Wj, Wi, qcjc0, qc0ci, q_bibj, J);
        problem.AddResidualBlock(cost_function, NULL, BGS[i], RIC_Q, Td[i], Td[i + 1]);
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.max_solver_time_in_seconds = 0.1;
    options.max_num_iterations = 50;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.minimizer_progress_to_stdout = true;
    // options.num_threads = 8;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout << summary.BriefReport() << endl;
    if (summary.final_cost > 1e-5)
    {
        return false;
    }

    Quaterniond RIC_;
    RIC_.w() = RIC_Q[0];
    RIC_.x() = RIC_Q[1];
    RIC_.y() = RIC_Q[2];
    RIC_.z() = RIC_Q[3];
    RIC_.normalized();
    RIC = RIC_.toRotationMatrix();

    // Matrix3d tem(RIC.transpose()*ric);
    // Vector3d tem_(Utility::R2ypr(tem));
    // double _tem(abs(tem_[0])+abs(tem_[1])+abs(tem_[2]));
    // if(_tem>3.0) return false;

    int j = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, j++)
    {
        double tt;
        tt = Td[j][0];
        cout << "Td=" << tt << endl;
        frame_j = next(frame_i);
        Vector3d Bgs_{BGS[j][0], BGS[j][1], BGS[j][2]};
        Bgs[j] = Bgs_;
        frame_j->second.pre_integration->init_refine_delta_pvq_bgs(Bgs_);
        cout << "bgs=" << Bgs[j].transpose() << endl;
        Vector3d Wi(frame_j->second.pre_integration->gyr_0);
        Quaterniond Qr{1.0, 0.5 * Wi[0] * tt, 0.5 * Wi[1] * tt, 0.5 * Wi[2] * tt};
        Quaterniond qclck(frame_i->second.R);
        Quaterniond qrci(RIC.transpose());
        qclck = qclck * qrci * Qr;
        frame_i->second.R = qclck.toRotationMatrix();
    }

    return true;
}

bool Estimate_vel_g_s_tic(map<double, ImageFrame> &all_image_frame, Vector3d &g, VectorXd &x, Vector3d &TIC, VectorXd &S, Vector3d *Bas)
{

    ceres::Problem problem;
    int frame_num(all_image_frame.size());
    double g_normal[3];
    g_normal[0] = G_DIRECTION[0];
    g_normal[1] = G_DIRECTION[1];
    g_normal[2] = G_DIRECTION[2];
    ceres::LocalParameterization *local_parameterization = new ceres::IdentityParameterization(3);
    problem.AddParameterBlock(g_normal, 3, local_parameterization);

    double s[frame_num][1];
    double bas[frame_num][3];

    double velocity[frame_num][3];
    for (int i = 0; i < frame_num; i++)
    {
        velocity[i][0] = 0.0;
        velocity[i][1] = 0.0;
        velocity[i][2] = 0.0;
        s[i][0] = 0.0;
        bas[i][0] = 0.0;
        bas[i][1] = 0.0;
        bas[i][2] = 0.0;
        problem.AddParameterBlock(velocity[i], 3);
        problem.AddParameterBlock(s[i], 1);
        problem.SetParameterLowerBound(s[i], 0, 0.0);
        problem.AddParameterBlock(bas[i], 3);
        problem.SetParameterUpperBound(bas[i], 0, 0.2);
        problem.SetParameterLowerBound(bas[i], 0, -0.2);
        problem.SetParameterUpperBound(bas[i], 1, 0.2);
        problem.SetParameterLowerBound(bas[i], 1, -0.2);
        problem.SetParameterUpperBound(bas[i], 2, 0.2);
        problem.SetParameterLowerBound(bas[i], 2, -0.2);
    }

    double pbc[3];

    if (ESTIMATE_EXTRINSIC == 0)
    {
        pbc[0] = TIC[0];
        pbc[1] = TIC[1];
        pbc[2] = TIC[2];
        problem.AddParameterBlock(pbc, 3);
        problem.SetParameterBlockConstant(pbc);
    }
    else if (ESTIMATE_EXTRINSIC == 1)
    {
        pbc[0] = TIC[0];
        pbc[1] = TIC[1];
        pbc[2] = TIC[2];
        problem.AddParameterBlock(pbc, 3);
    }
    else
    {
        pbc[0] = 0.0;
        pbc[1] = 0.0;
        pbc[2] = 0.0;
        problem.AddParameterBlock(pbc, 3);
        problem.SetParameterUpperBound(pbc, 0, PBC_UX);
        problem.SetParameterLowerBound(pbc, 0, PBC_LX);
        problem.SetParameterUpperBound(pbc, 1, PBC_UY);
        problem.SetParameterLowerBound(pbc, 1, PBC_LY);
        problem.SetParameterUpperBound(pbc, 2, PBC_UZ);
        problem.SetParameterLowerBound(pbc, 2, PBC_LZ);
    }

    map<double, ImageFrame>::iterator frame_i;
    map<double, ImageFrame>::iterator frame_j;
    int k = 0;
    for (frame_i = all_image_frame.begin(); next(frame_i) != all_image_frame.end(); frame_i++, k++)
    {
        frame_j = next(frame_i);
        Vector3d deltap(frame_j->second.pre_integration->delta_p);
        Vector3d deltav(frame_j->second.pre_integration->delta_v);
        double dt = frame_j->second.pre_integration->sum_dt;
        cout << "sumt=" << dt << endl;

        Matrix3d rbic0(frame_i->second.R.transpose());
        Matrix3d rcobj(frame_j->second.R);

        // Matrix3d rbic0(RIC[0]*frame_i->second.R.transpose());
        // Matrix3d rcobj(frame_j->second.R*RIC[0].transpose());

        double gnorm(G.z());
        Vector3d pc0cj(frame_j->second.T);
        Vector3d pc0ci(frame_i->second.T);
        Matrix3d jp = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_P, O_BA);
        Matrix3d jv = frame_j->second.pre_integration->jacobian.template block<3, 3>(O_V, O_BA);
        //cout<<"pcj-pci="<<(pc0cj-pc0ci).transpose()<<endl;
        ceres::CostFunction *cost_function = InitTranslationConstraint::makeConstraint_t(deltap, deltav, dt, rbic0, rcobj, gnorm, pc0cj, pc0ci, jp, jv);

        problem.AddResidualBlock(cost_function, NULL, velocity[k], velocity[k + 1], pbc, g_normal, s[k], s[k + 1], bas[k]);
    }
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.max_solver_time_in_seconds = 0.1;
    options.max_num_iterations = 100;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    cout << summary.BriefReport() << endl;
    if (summary.final_cost > 5e-3) //2e-3 1e-5
    {
        return false;
    }

    g[0] = g_normal[0];
    g[1] = g_normal[1];
    g[2] = g_normal[2];
    g = G.z() * g.normalized();
    //g[0] = G.z()* g_normal[0];g[1] = G.z()* g_normal[1];g[2] = G.z()* g_normal[2];
    TIC[0] = pbc[0];
    TIC[1] = pbc[1];
    TIC[2] = pbc[2];

    int state_num = frame_num * 3 + 3 + 3;
    x.resize(state_num);
    x.setZero();
    S.resize(frame_num);
    S.setZero();
    // for (int j = 0; j < frame_num ;j++){
    //     S[j]=s[j][0];
    //     Vector3d Bas_{bas[j][0],bas[j][1],bas[j][2]};
    //     if(S[j]<0.0||abs(Bas_[0])>3.0||abs(Bas_[1])>3.0||abs(Bas_[2])>3.0)  //||abs(Bas_[0])>3.0||abs(Bas_[1])>3.0||abs(Bas_[2])>3.0
    //         return false;

    // }

    for (int j = 0; j < frame_num; j++)
    {
        x[3 * j] = velocity[j][0];
        x[3 * j + 1] = velocity[j][1];
        x[3 * j + 2] = velocity[j][2];
        S[j] = s[j][0];
        Vector3d Bas_{bas[j][0], bas[j][1], bas[j][2]};
        Bas[j] = Bas_;
        frame_j->second.pre_integration->init_refine_delta_pvq_bas(Bas_);
        cout << "bas=" << Bas[j].transpose() << endl;
    }

    x[state_num - 3] = g_normal[0];
    x[state_num - 2] = g_normal[1];
    x[state_num - 1] = g_normal[2];

    x[state_num - 6] = pbc[0];
    x[state_num - 5] = pbc[1];
    x[state_num - 4] = pbc[2];

    cout << "g=" << g.transpose() << endl;
    cout << "pbc=" << pbc[0] << " " << pbc[1] << " " << pbc[2] << " " << endl;
    cout << "scale=" << S.transpose() << endl;

    return true;
}

bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d *Bgs, Vector3d *Bas, Vector3d &g, VectorXd &x, double *init_dt, Eigen::Matrix3d &RIC, Eigen::Matrix3d ric, Vector3d &TIC, VectorXd &S)
{
    // solveGyroscopeBias(all_image_frame, Bgs);
    cout << "*************START Estimate_ric_td_bg***************" << endl;
    cout << "BEFORE SLOVE RIC=\n"
         << Utility::R2ypr(RIC).transpose() << endl;
    if (!Estimate_ric_td_bg(all_image_frame, RIC, init_dt, Bgs, ric))
        return false;
    cout << "AFTER SLOVE RIC=\n"
         << Utility::R2ypr(RIC).transpose() << endl;

    cout << "***************START Estimate_vel_g_s_tic***************" << endl;
    if (Estimate_vel_g_s_tic(all_image_frame, g, x, TIC, S, Bas))
    {
        cout << "***************************END**********************************" << endl;
        return true;
    }
    else
        return false;

    // for (int j = 0; j < all_image_frame.size() ;j++){
    //     if (S[j]<0.0)
    //         return false;
    //     else
    //         continue;
    // }

    // return true;

    // if(LinearAlignment(all_image_frame, g, x))
    //     return true;
    // else
    //     return false;
}
