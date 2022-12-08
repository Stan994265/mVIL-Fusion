#ifndef LIDAR_BACK_END_FRONT_END_H
#define LIDAR_BACK_END_FRONT_END_H

#include "lidar_frontend.h"
struct LidarICPConstraint
{
    double lidar_ta;
    double lidar_tb;
    double lidar_tc;
    double lidar_td;
    double lidar_ti;
    double lidar_tj;
    // Eigen::MatrixXd lidar_cov;
    Eigen::Matrix<double, 6, 6> lidar_sqrt_info;
    Eigen::Matrix4d lidar_trans;
    int constraint_mode;
};

struct LidarLPSConstraint
{
    // double lidar_tl;
    // double lidar_tr;
    double lidar_t;
    Eigen::Quaterniond LPSq;
    Eigen::Vector3d LPSt;
};

bool FindWindowsID(const std_msgs::Header Headers[(WINDOW_SIZE + 1)],
                   const double ta, const double tb, const double tc, const double td,
                   int &id_a, int &id_b, int &id_c, int &id_d);
bool FindNearest2ID(const std_msgs::Header Headers[(WINDOW_SIZE + 1)],
                    const double tl,
                    int &id_a, int &id_b);

struct LPSConstraint
{

    LPSConstraint(const double tl_, const double tr_, const double tk_,
                  const Quaterniond trans_q_, const Vector3d trans_t_, const double sqrt_info_)

        : tl(tl_), tr(tr_), tk(tk_), trans_q(trans_q_), trans_t(trans_t_), sqrt_info(sqrt_info_)
    {
    }

    template <typename T>
    bool operator()(const T *POSEa, const T *POSEb, T *residuals) const
    {
        Eigen::Quaternion<T> Qa{POSEa[6], POSEa[3], POSEa[4], POSEa[5]};
        Eigen::Quaternion<T> Qb{POSEb[6], POSEb[3], POSEb[4], POSEb[5]};
        Eigen::Quaternion<T> Qi;
        double t_i = (tk - tl) / (tr - tl);
        Qi = Qa.slerp(T(t_i), Qb);
        Qi.normalized();

        Eigen::Quaternion<T> Q1{T(trans_q.w()), T(trans_q.x()), T(trans_q.y()), T(trans_q.z())};
        Eigen::Quaternion<T> Q12;
        Q12 = Qi.inverse() * Q1;

        residuals[0] = T(2) * Q12.x() / T(0.01);
        residuals[1] = T(2) * Q12.y() / T(0.01);
        residuals[2] = T(2) * Q12.z() / T(0.01);

        // residuals[3] = T(2) * Q12.x() / T(0.01);
        // residuals[4] = T(2) * Q12.y() / T(0.01);
        // residuals[5] = T(2) * Q12.z() / T(0.01);

        // Eigen::Matrix<T, 3, 1> Pa{POSEa[0], POSEa[1], POSEa[2]};
        // Eigen::Matrix<T, 3, 1> Pb{POSEb[0], POSEb[1], POSEb[2]};
        // Eigen::Matrix<T, 3, 1> Pi;
        // Pi = Pa + (Pb - Pa) / T(tr - tl) * T(tk - tl);

        // Eigen::Matrix<T, 3, 1> transltion{T{trans_t[0]}, T{trans_t[1]}, T{trans_t[2]}};
        // Eigen::Matrix<T, 3, 1> dff_t;
        // dff_t = Qi.inverse() * (transltion - Pi);

        // residuals[0] = dff_t(0, 0) / T(0.1);
        // residuals[1] = T(0.0); //Altitude 10000
        // residuals[2] = dff_t(2, 0) / T(0.1);
        return true;
    }

    static ceres::CostFunction *makeConstraint(const double tl_, const double tr_, const double tk_,
                                               const Quaterniond trans_q_, const Vector3d trans_t_, const double sqrt_info_)
    {

        return (new ceres::AutoDiffCostFunction<
                LPSConstraint, 3, 7, 7>(
            new LPSConstraint(tl_, tr_, tk_, trans_q_, trans_t_, sqrt_info_)));
    }

    //observation
    const double tl, tr, tk, sqrt_info;
    const Eigen::Quaterniond trans_q;
    const Eigen::Vector3d trans_t;
};

struct LidarICPConstraint_b
{

    LidarICPConstraint_b(const double ta_, const double tb_, const double tc_, const double td_,
                         const double ti_, const double tj_, const Matrix4d trans_, const Matrix<double, 6, 6> sqrt_info_)

        : ta(ta_), tb(tb_), tc(tc_), td(td_), ti(ti_), tj(tj_), trans(trans_), sqrt_info(sqrt_info_)
    {
    }

    template <typename T>
    bool operator()(const T *POSEa, const T *POSEb, const T *POSEc, const T *POSEd, T *residuals) const
    {
        Eigen::Quaternion<T> Qa{POSEa[6], POSEa[3], POSEa[4], POSEa[5]};
        Eigen::Matrix<T, 3, 1> Pa{POSEa[0], POSEa[1], POSEa[2]};
        Eigen::Quaternion<T> Qb{POSEb[6], POSEb[3], POSEb[4], POSEb[5]};
        Eigen::Matrix<T, 3, 1> Pb{POSEb[0], POSEb[1], POSEb[2]};
        Eigen::Quaternion<T> Qc{POSEc[6], POSEc[3], POSEc[4], POSEc[5]};
        Eigen::Matrix<T, 3, 1> Pc{POSEc[0], POSEc[1], POSEc[2]};
        Eigen::Quaternion<T> Qd{POSEd[6], POSEd[3], POSEd[4], POSEd[5]};
        Eigen::Matrix<T, 3, 1> Pd{POSEd[0], POSEd[1], POSEd[2]};

        Eigen::Quaternion<T> Qi;
        Eigen::Quaternion<T> Qj;
        double t_i = (ti - ta) / (tb - ta);
        double t_j = (tj - tc) / (td - tc);
        Qi = Qa.slerp(T(t_i), Qb);
        Qj = Qc.slerp(T(t_j), Qd);
        Qi.normalized();
        Qj.normalized();

        Eigen::Matrix<T, 3, 1> Pi; //pw_bi
        Eigen::Matrix<T, 3, 1> Pj;
        Pi = Pa + (Pb - Pa) / T(tb - ta) * T(ti - ta);
        Pj = Pc + (Pd - Pc) / T(td - tc) * T(tj - tc);
        //0
        // Eigen::Quaternion<T> temQ;
        // temQ = Qi * Qj.inverse();

        // Eigen::Matrix<T, 3, 1> PIJ;
        // PIJ = Pi - temQ * Pj;

        // T TIJ[3]{T(PIJ(0, 0)), T(PIJ(1, 0)), T(PIJ(2, 0))};
        // T TIJ_[3]{T(trans(0, 3)), T(trans(1, 3)), T(trans(2, 3))}; //TBIJ

        // residuals[0] = (TIJ[0] - TIJ_[0]) * T(sqrt_info(0, 0));
        // residuals[1] = (TIJ[1] - TIJ_[1]) * T(sqrt_info(0, 0));
        // residuals[2] = (TIJ[2] - TIJ_[2]) * T(sqrt_info(0, 0));
        //1
        Eigen::Quaternion<T> temQ;
        temQ = Qj.inverse() * Qi;

        Eigen::Matrix<T, 3, 1> temPIJ;
        temPIJ = Qi.inverse() * (Pj - Pi);

        Eigen::Matrix<T, 3, 1> PIJ{T(trans(0, 3)), T(trans(1, 3)), T(trans(2, 3))};
        Eigen::Matrix<T, 3, 1> RES;
        RES = temQ * (PIJ - temPIJ);

        residuals[0] = (RES(0, 0)) * T(sqrt_info(0, 0));
        residuals[1] = T(0.0);
        residuals[2] = (RES(2, 0)) * T(sqrt_info(0, 0));
        //2
        // Eigen::Matrix<T, 3, 1> PIJ{T(trans(0, 3)), T(trans(1, 3)), T(trans(2, 3))};
        // Eigen::Matrix<T, 3, 1> RES;
        // RES = Pi - Pj - PIJ;

        // residuals[0] = (RES(0, 0)) * T(sqrt_info(0, 0));
        // residuals[1] = (RES(1, 0)) * T(sqrt_info(0, 0));
        // residuals[2] = (RES(2, 0)) * T(sqrt_info(0, 0));

        return true;
    }

    static ceres::CostFunction *makeConstraint_b(const double ta_, const double tb_, const double tc_, const double td_,
                                                 const double ti_, const double tj_, const Matrix4d trans_, const Matrix<double, 6, 6> sqrt_info_)
    {

        return (new ceres::AutoDiffCostFunction<
                LidarICPConstraint_b, 3, 7, 7, 7, 7>(
            new LidarICPConstraint_b(ta_, tb_, tc_, td_, ti_, tj_, trans_, sqrt_info_)));
    }

    //observation
    const double ta, tb, tc, td, ti, tj;
    const Matrix4d trans;
    const Matrix<double, 6, 6> sqrt_info;
};

#endif