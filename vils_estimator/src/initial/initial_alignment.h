#pragma once
#include <eigen3/Eigen/Dense>
#include <iostream>
#include "../factor/imu_factor.h"
#include "../utility/utility.h"
#include <ros/ros.h>
#include <map>
#include "../feature_manager.h"
#include <ceres/rotation.h>

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl_conversions/pcl_conversions.h>

using namespace Eigen;
using namespace std;

class ImageFrame
{
public:
    ImageFrame(){};
    ImageFrame(const map<int, vector<pair<int, Eigen::Matrix<double, 8, 1> > > > &_points, double _t) : t{_t}, is_key_frame{false}
    {
        points = _points;
    };
    map<int, vector<pair<int, Eigen::Matrix<double, 8, 1> > > > points;
    double t;
    Matrix3d R;
    Vector3d T;
    IntegrationBase *pre_integration;
    bool is_key_frame;
};

//bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d* Bgs, Vector3d &g, VectorXd &x);
bool VisualIMUAlignment(map<double, ImageFrame> &all_image_frame, Vector3d *Bgs, Vector3d *Bas,
                        Vector3d &g, VectorXd &x, double *init_dt, Eigen::Matrix3d &RIC, Eigen::Matrix3d ric, Vector3d &TIC, VectorXd &S);

struct InitRotationConstraint
{

    InitRotationConstraint(const Vector3d wj_, const Vector3d wi_,
                           const Quaterniond qcjc0_, const Quaterniond qc0ci_,
                           const Quaterniond qbibj_, const Matrix3d j_dq_dbg_)

        : wj(wj_), wi(wi_), qcjc0(qcjc0_), qc0ci(qc0ci_), qbibj(qbibj_), j_dq_dbg(j_dq_dbg_)
    {
    }

    template <typename T>
    bool operator()(const T *BGS, const T *RIC_Q, const T *Tdi, const T *Tdj, T *residuals) const
    {

        Eigen::Quaternion<T> Ql{T(1.0), -0.5 * T(wj[0]) * Tdj[0], -0.5 * T(wj[1]) * Tdj[0], -0.5 * T(wj[2]) * Tdj[0]};
        Eigen::Quaternion<T> Qr{T(1.0), 0.5 * T(wi[0]) * Tdi[0], 0.5 * T(wi[1]) * Tdi[0], 0.5 * T(wi[2]) * Tdi[0]};
        Eigen::Quaternion<T> Qbc{RIC_Q[0], RIC_Q[1], RIC_Q[2], RIC_Q[3]};
        Eigen::Quaternion<T> Qcb{RIC_Q[0], -RIC_Q[1], -RIC_Q[2], -RIC_Q[3]};
        Eigen::Quaternion<T> Qcj{T(qcjc0.w()), T(qcjc0.x()), T(qcjc0.y()), T(qcjc0.z())};
        Eigen::Quaternion<T> Qci{T(qc0ci.w()), T(qc0ci.x()), T(qc0ci.y()), T(qc0ci.z())};
        Eigen::Quaternion<T> Qbij{T(qbibj.w()), T(qbibj.x()), T(qbibj.y()), T(qbibj.z())};

        Eigen::Matrix<T, 3, 3> J;
        J << T(j_dq_dbg(0, 0)), T(j_dq_dbg(0, 1)), T(j_dq_dbg(0, 2)),
            T(j_dq_dbg(1, 0)), T(j_dq_dbg(1, 1)), T(j_dq_dbg(1, 2)),
            T(j_dq_dbg(2, 0)), T(j_dq_dbg(2, 1)), T(j_dq_dbg(2, 2));
        Eigen::Matrix<T, 3, 1> Bgs{BGS[0], BGS[1], BGS[2]};
        Eigen::Matrix<T, 3, 1> Qjbg_v;
        Qjbg_v = T(0.5) * J * Bgs;

        Eigen::Quaternion<T> Qjbg{T(1.0), Qjbg_v(0, 0), Qjbg_v(1, 0), Qjbg_v(2, 0)};

        Eigen::Quaternion<T> Qresiduals;

        //Qresiduals = Qbc*Ql*Qcj*Qci*Qr*Qcb*Qbij*Qjbg;
        Qresiduals = Ql * Qbc * Qcj * Qci * Qcb * Qr * Qbij * Qjbg;
        //Qresiduals.normalized();
        residuals[0] = T(2) * Qresiduals.x();
        residuals[1] = T(2) * Qresiduals.y();
        residuals[2] = T(2) * Qresiduals.z();

        return true;
    }

    static ceres::CostFunction *makeConstraint_r(const Vector3d wj_, const Vector3d wi_,
                                                 const Quaterniond qcjc0_, const Quaterniond qc0ci_,
                                                 const Quaterniond qbibj_, const Matrix3d j_dq_dbg_)
    {

        return (new ceres::AutoDiffCostFunction<
                InitRotationConstraint, 3, 3, 4, 1, 1>(
            new InitRotationConstraint(wj_, wi_, qcjc0_, qc0ci_, qbibj_, j_dq_dbg_)));
    }

    //observation
    const Vector3d wj;
    const Vector3d wi;
    const Quaterniond qcjc0;
    const Quaterniond qc0ci;
    const Quaterniond qbibj;
    const Matrix3d j_dq_dbg;
};

struct InitTranslationConstraint
{

    InitTranslationConstraint(const Vector3d deltap_, const Vector3d deltav_, const double dt_,
                              const Matrix3d rbic0_, const Matrix3d rcobj_, const double gnorm_,
                              const Vector3d pc0cj_, const Vector3d pc0ci_,
                              const Matrix3d j_dp_dba_, const Matrix3d j_dv_dba_)

        : deltap(deltap_), deltav(deltav_), dt(dt_), rbic0(rbic0_), rcobj(rcobj_), gnorm(gnorm_),
          pc0cj(pc0cj_), pc0ci(pc0ci_), j_dp_dba(j_dp_dba_), j_dv_dba(j_dv_dba_)
    {
    }

    template <typename T>
    bool operator()(const T *VI, const T *VJ, const T *PBC, const T *GN, const T *SI, const T *SJ, const T *BAS, T *residuals) const
    {

        Eigen::Matrix<T, 3, 1> Vi{VI[0], VI[1], VI[2]};
        Eigen::Matrix<T, 3, 1> Vj{VJ[0], VJ[1], VJ[2]};
        Eigen::Matrix<T, 3, 1> Pbc{PBC[0], PBC[1], PBC[2]};
        //Eigen::Matrix <T,3,1>  Pcb; Pcb=-Pbc;
        Eigen::Matrix<T, 3, 1> Gc0{T(gnorm) * GN[0], T(gnorm) * GN[1], T(gnorm) * GN[2]};
        Eigen::Matrix<T, 3, 1> Dp{T(deltap[0]), T(deltap[1]), T(deltap[2])};
        Eigen::Matrix<T, 3, 1> Dv{T(deltav[0]), T(deltav[1]), T(deltav[2])};
        Eigen::Matrix<T, 3, 1> Pcj{T(pc0cj[0]), T(pc0cj[1]), T(pc0cj[2])};
        Eigen::Matrix<T, 3, 1> Pci{T(pc0ci[0]), T(pc0ci[1]), T(pc0ci[2])};
        Eigen::Matrix<T, 3, 3> Rbic0;
        Rbic0 << T(rbic0(0, 0)), T(rbic0(0, 1)), T(rbic0(0, 2)),
            T(rbic0(1, 0)), T(rbic0(1, 1)), T(rbic0(1, 2)),
            T(rbic0(2, 0)), T(rbic0(2, 1)), T(rbic0(2, 2));
        Eigen::Matrix<T, 3, 3> Rc0bi;
        Rc0bi = Rbic0.transpose();

        Eigen::Matrix<T, 3, 3> Rcobj;
        Rcobj << T(rcobj(0, 0)), T(rcobj(0, 1)), T(rcobj(0, 2)),
            T(rcobj(1, 0)), T(rcobj(1, 1)), T(rcobj(1, 2)),
            T(rcobj(2, 0)), T(rcobj(2, 1)), T(rcobj(2, 2));

        Eigen::Matrix<T, 3, 1> RES0_2;
        Eigen::Matrix<T, 3, 1> RES3_5;

        Eigen::Matrix<T, 3, 3> JP;
        JP << T(j_dp_dba(0, 0)), T(j_dp_dba(0, 1)), T(j_dp_dba(0, 2)),
            T(j_dp_dba(1, 0)), T(j_dp_dba(1, 1)), T(j_dp_dba(1, 2)),
            T(j_dp_dba(2, 0)), T(j_dp_dba(2, 1)), T(j_dp_dba(2, 2));
        Eigen::Matrix<T, 3, 3> JV;
        JV << T(j_dv_dba(0, 0)), T(j_dv_dba(0, 1)), T(j_dv_dba(0, 2)),
            T(j_dv_dba(1, 0)), T(j_dv_dba(1, 1)), T(j_dv_dba(1, 2)),
            T(j_dv_dba(2, 0)), T(j_dv_dba(2, 1)), T(j_dv_dba(2, 2));
        Eigen::Matrix<T, 3, 1> Bas{BAS[0], BAS[1], BAS[2]};

        RES0_2 = Dp + JP * Bas - Pbc + Rbic0 * Rcobj * Pbc - Rbic0 * (SJ[0] * Pcj - SI[0] * Pci) + Vi * T(dt) - T(0.5) * Rbic0 * Gc0 * T(dt) * T(dt);
        RES3_5 = Dv + JV * Bas - Rbic0 * (Rcobj * Vj - Rc0bi * Vi + Gc0 * T(dt));

        // RES0_2=Dp+JP*Bas-Rbic0*SJ[0]*Pcj+Rbic0*SI[0]*Pci+Rbic0*Vi*T(dt)-Rbic0*T(0.5)*Gc0*T(dt)*T(dt)+Rbic0*Rcobj*Pbc-Pbc;
        // RES3_5=Dv+JV*Bas-Rbic0*Rcobj*Vj+Rbic0*Gc0*T(dt)+Vi;

        // RES0_2=RES0_2*T(1000);
        // RES3_5=RES3_5*T(1000);

        residuals[0] = RES0_2(0, 0);
        residuals[1] = RES0_2(1, 0);
        residuals[2] = RES0_2(2, 0);
        residuals[3] = RES3_5(0, 0);
        residuals[4] = RES3_5(1, 0);
        residuals[5] = RES3_5(2, 0);

        return true;
    }

    static ceres::CostFunction *makeConstraint_t(const Vector3d deltap_, const Vector3d deltav_, const double dt_,
                                                 const Matrix3d rbic0_, const Matrix3d rcobj_, const double gnorm_,
                                                 const Vector3d pc0cj_, const Vector3d pc0ci_,
                                                 const Matrix3d j_dp_dba_, const Matrix3d j_dv_dba_)
    {

        return (new ceres::AutoDiffCostFunction<
                InitTranslationConstraint, 6, 3, 3, 3, 3, 1, 1, 3>(
            new InitTranslationConstraint(deltap_, deltav_, dt_, rbic0_, rcobj_, gnorm_, pc0cj_, pc0ci_, j_dp_dba_, j_dv_dba_)));
    }

    //observation
    const Vector3d deltap;
    const Vector3d deltav;
    const double dt;
    const Matrix3d rbic0;
    const Matrix3d rcobj;
    const double gnorm;
    const Vector3d pc0cj;
    const Vector3d pc0ci;
    const Matrix3d j_dp_dba;
    const Matrix3d j_dv_dba;
};
