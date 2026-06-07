/**
* This file is part of SUPER
*
* Copyright 2025 Yunfan REN, MaRS Lab, University of Hong Kong, <mars.hku.hk>
* Developed by Yunfan REN <renyf at connect dot hku dot hk>
* for more information see <https://github.com/hku-mars/SUPER>.
* If you use this code, please cite the respective publications as
* listed on the above website.
*
* SUPER is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SUPER is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with SUPER. If not, see <http://www.gnu.org/licenses/>.
*/

/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#pragma once

#include <utils/header/type_utils.hpp>
#include <utils/geometry/quickhull.h>
#include <Eigen/Eigen>


namespace geometry_utils {
    using super_utils::Mat3f;
    using super_utils::Vec3f;
    using super_utils::Vec4f;
    using super_utils::vec_Vec3f;
    using super_utils::Mat3Df;
    using super_utils::MatD4f;
    using super_utils::vec_E;
    using super_utils::PolyhedronH;

    ///============ 2023-06-30: add by yunfan ============///
    // 使用简化的一维 Point-Mass(质点) 运动模型，为一段长度为 total_dis 的路径分配时间。
    // 输入：
    // a_max: 最大加速度
    // v_max: 最大速度
    // v0:    该段路径的参考初速度，用于与前一段路径做速度连续衔接
    // total_dis: 整段路径的总长度
    // cur_dis:   剩余距离，详见super_planner
    // 输出：
    // t:   从该段起点走到 cur_dis 所需的累计时间
    // vel: 到达 cur_dis 时的速度
    // 计算思路：
    // 1) 先根据总路程 total_dis 判断整段速度剖面属于哪一类：
    //    仅加速、加速后立刻减速、加速+匀速+减速；
    // 2) 再根据当前点 cur_dis 落在该速度剖面的哪个阶段，反推出该点的时间和速度。
    static void simplePMTimeAllocator(const double &a_max, const double &v_max,
                                const double &v0,
                                const double &total_dis,
                                const double &cur_dis, double &t, double &vel) {
        // Helper lambda functions
        // 已知加速度和时间，计算走过的距离 s = 1/2 * a * t^2。
        auto calc_dis = [](double a, double t) { return 0.5 * a * t * t; };
        // 已知加速度和距离，反解“从静止匀加速走完这段距离所需时间”。从该公式 s = 1/2 * a * t^2 中反解求出t 
        auto calc_time = [](double a, double cur_dis) { return sqrt(2 * cur_dis / a); };
        // 一元二次方程求根公式。用于求解速度分配里出现的一元二次方程，返回其中一个实根。
        auto solve_quadratic = [](double a, double b, double c) {
            double delta = b * b - 4 * a * c;
            return (-b + sqrt(delta)) / (2 * a);
        };

        // Precompute reusable values
        // 从 0 以 a_max 加速到 v_max 所需时间和距离。
        const double t_to_v_max = v_max / a_max;
        const double dis_to_v_max = calc_dis(a_max, t_to_v_max);

        // 从 0 以 a_max 加速到 v0 所需时间和距离。
        const double t_to_v0 = v0 / a_max;
        const double dis_to_v0 = calc_dis(a_max, t_to_v0);

        // 从 v_max 以 a_max 减速到 v0 所需时间和距离。
        const double dec_time = (v_max - v0) / a_max;
        const double dec_dis = 0.5 * (v_max + v0) * dec_time;

        // Case 1: Only acceleration to v0
        // 情况 1：总路程太短，短到连完整“升到某个更高速度再减速”都不需要，
        // 直接按单段匀加速模型分配时间即可。
        if (total_dis <= dis_to_v0) {
            t = calc_time(a_max, cur_dis);
            vel = a_max * t;
            return;
        }

        // Case 2: Acceleration to v_max, then deceleration
        // 情况 2：总路程中等，能先加速，但不足以支撑完整的“加速到 v_max + 匀速巡航 + 减速”；
        // 因此速度是“加速到某个峰值速度，再立刻减速”的三角速度模型。
        // 先完全忽略初始速度v0，使用以 a_max 为加速度的匀加速模型，此时加速时间和减速时间相等，均为t_acc
        if (total_dis <= dis_to_v_max + dec_dis) {
            // 这里通过解二次方程求出：加速段结束时间 t_acc，加速段的初速度为 v0
            // 然后根据 t_acc 该段距离 dis_acc、以及峰值速度 cur_v_max。
            // 但是这个地方的求解乱七八糟的感觉是错的，b 应该是0才对
            // 正确的方程建模应为 初速度为0，加速度为a_max的加速段 + 末速度为v0，减速度为a_max的减速段 
            const double a = 2 * a_max;
            const double b = -(a_max - v0);
            const double c = -(v0 * v0 / a_max + 2 * total_dis);
            const double t_acc = solve_quadratic(a, b, c);
            // t_dec 表示从峰值速度 cur_v_max 减速到末端速度 v0 所需的时间。
            // t_to_v0 表示 从 0 加速到 v0 的时间，
            const double t_dec = t_acc - t_to_v0;
            const double dis_acc = calc_dis(a_max, t_acc);
            const double cur_v_max = a_max * t_acc;

            // 若当前点还在加速段内，就直接按匀加速公式反解时间和速度。
            if (cur_dis <= dis_acc) {
                t = calc_time(a_max, cur_dis);
                vel = a_max * t;
            } else {
                // 否则当前点已经进入减速段，需要再解一次二次方程，求减速段内所经历的时间 t2。
                const double remaining_dis = cur_dis - dis_acc;
                const double t2 = solve_quadratic(-a_max, 2 * cur_v_max, -2 * remaining_dis);
                t = t_acc + t2;
                vel = cur_v_max - t2 * a_max;
            }
            return;
        }

        // Case 3: Acceleration + constant speed + deceleration
        // 情况 3：总路程足够长，可以形成完整的三段式模式：
        // 加速到 v_max -> 匀速巡航 -> 再减速到末端速度。
        // Case 3.1: During acceleration phase
        if (cur_dis < dis_to_v_max) {  
            // 当前点落在加速段。
            t = calc_time(a_max, cur_dis);
            vel = a_max * t;
            return;
        }
        
        // Case 3.2: During constant speed phase
        if (cur_dis < total_dis - dec_dis) {  
            // 当前点落在匀速段：总时间 = 加速时间 + 匀速阶段已走时间。
            const double remaining_dis = cur_dis - dis_to_v_max;
            const double t_const = remaining_dis / v_max;
            t = t_to_v_max + t_const;
            vel = v_max;
            return;
        }

        // Case 3.3: During deceleration phase
        // 当前点落在减速段：先扣掉前面的加速段和匀速段距离，再反解减速阶段内的时间。
        const double const_phase_dis = total_dis - dec_dis - dis_to_v_max;
        const double remaining_dis = cur_dis - dis_to_v_max - const_phase_dis;
        const double t_dec = solve_quadratic(-a_max, 2 * v_max, -2 * remaining_dis);
        t = t_to_v_max + const_phase_dis / v_max + t_dec;
        vel = v_max - t_dec * a_max;
    }

    ///============ 2023-06-30: add by yunfan ============///
    double DistancePointEllipse(double e0, double e1, double y0, double y1, double& x0, double& x1);

    double
    DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2, double& x0, double& x1,
                           double& x2);


    ///============ 2023-06-13: add by Gene ============///
    template <typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 1> quaternion_to_yrp(const Eigen::Quaternion<Scalar_t>& q_);

    ///============ 2023-06-13: add by Yunfan ============///
    Vec4f translatePlane(const Vec4f& plane, const Vec3f& translation);

    ///============ 2023-05-23: add by Yunfan ============///
    void normalizeNextYaw(const double& last_yaw, double& yaw);

    ///============ 2023-3-12: add by Yunfan ============///
    void convertFlatOutputToAttAndOmg(const Vec3f& p,
                                      const Vec3f& v,
                                      const Vec3f& a,
                                      const Vec3f& j,
                                      const double& yaw,
                                      const double& yaw_dot,
                                      Vec3f& rpy,
                                      Vec3f& omg,
                                      double& aT
    );


    ///============ 2022-12-10: add by Yunfan ============///
    bool pointInsidePolytope(const Vec3f& point, const PolyhedronH& polytope,
                             double margin = 1e-6);

    ///============ 2022-12-5: add by Yunfan ============///
    double pointLineSegmentDistance(const Vec3f& p, const Vec3f& a, const Vec3f& b);

    double computePathLength(const vec_E<Vec3f>& path);

    ///============ 2022-11-18 ====================================================================================
    int inline GetIntersection(float fDst1, float fDst2, Vec3f P1, Vec3f P2, Vec3f& Hit);

    int inline InBox(Vec3f Hit, Vec3f B1, Vec3f B2, const int Axis);

    //The box in this article is Axis-Aligned and so can be defined by only two 3D points:
    // B1 - the smallest values of X, Y, Z
    //        B2 - the largest values of X, Y, Z
    // returns true if line (L1, L2) intersects with the box (B1, B2)
    // returns intersection point in Hit
    int lineIntersectBox(Vec3f L1, Vec3f L2, Vec3f B1, Vec3f B2, Vec3f& Hit);

    Vec3f lineBoxIntersectPoint(const Vec3f& pt, const Vec3f& pos,
                                const Vec3f& box_min, const Vec3f& box_max);

    ///================================================================================================


    Eigen::Matrix3d RotationFromVec3(const Eigen::Vector3d& v);

    // 通过三点获得一个平面
    void FromPointsToPlane(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3,
                           Eigen::Vector4d& hPoly);

    void getFovCheckPlane(const Eigen::Matrix3d R, const Eigen::Vector3d t, Eigen::MatrixX4d& fov_planes,
                          std::vector<Eigen::Matrix3d>& fov_pts);

    void GetFovPlanes(const Eigen::Matrix3d R, const Eigen::Vector3d t, Eigen::MatrixX4d& fov_planes,
                      std::vector<Eigen::Matrix3d>& fov_pts);


    double findInteriorDist(const Eigen::MatrixX4d& hPoly,
                            Eigen::Vector3d& interior);

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    bool findInterior(const Eigen::MatrixX4d& hPoly,
                      Eigen::Vector3d& interior);

    bool overlap(const Eigen::MatrixX4d& hPoly0,
                 const Eigen::MatrixX4d& hPoly1,
                 const double eps = 1.0e-6);

    struct filterLess {
        bool operator()(const Eigen::Vector3d& l,
                        const Eigen::Vector3d& r) const {
            return l(0) < r(0) ||
            (l(0) == r(0) &&
                (l(1) < r(1) ||
                    (l(1) == r(1) &&
                        l(2) < r(2))));
        }
    };

    void filterVs(const Eigen::Matrix3Xd& rV,
                  const double& epsilon,
                  Eigen::Matrix3Xd& fV);

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    void enumerateVs(const Eigen::MatrixX4d& hPoly,
                     const Eigen::Vector3d& inner,
                     Eigen::Matrix3Xd& vPoly,
                     const double epsilon = 1.0e-6);

    // Each row of hPoly is defined by h0, h1, h2, h3 as
    // h0*x + h1*y + h2*z + h3 <= 0
    // proposed epsilon is 1.0e-6
    bool enumerateVs(const Eigen::MatrixX4d& hPoly,
                     Eigen::Matrix3Xd& vPoly,
                     const double epsilon = 1.0e-6);


    template <typename Scalar_t>
    Scalar_t toRad(const Scalar_t& x);

    template <typename Scalar_t>
    Scalar_t toDeg(const Scalar_t& x);

    template <typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 3> rotx(Scalar_t t);

    template <typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 3> roty(Scalar_t t);

    template <typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 3> rotz(Scalar_t t);

    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 3> ypr_to_R(const Eigen::DenseBase<Derived>& ypr);

    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 3>
    vec_to_R(const Eigen::MatrixBase<Derived>& v1, const Eigen::MatrixBase<Derived>& v2);

    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 1> R_to_ypr(const Eigen::DenseBase<Derived>& R);

    template <typename Derived>
    Eigen::Quaternion<typename Derived::Scalar> ypr_to_quaternion(const Eigen::DenseBase<Derived>& ypr);

    template <typename Scalar_t>
    Eigen::Matrix<Scalar_t, 3, 1> quaternion_to_ypr(const Eigen::Quaternion<Scalar_t>& q_);

    template <typename Scalar_t>
    Scalar_t get_yaw_from_quaternion(const Eigen::Quaternion<Scalar_t>& q);

    template <typename Scalar_t>
    Eigen::Quaternion<Scalar_t> yaw_to_quaternion(Scalar_t yaw);

    template <typename Scalar_t>
    Scalar_t normalize_angle(Scalar_t a);

    template <typename Scalar_t>
    Scalar_t angle_add(Scalar_t a, Scalar_t b);

    template <typename Scalar_t>
    Scalar_t yaw_add(Scalar_t a, Scalar_t b);

    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 3> get_skew_symmetric(const Eigen::DenseBase<Derived>& v);

    template <typename Derived>
    Eigen::Matrix<typename Derived::Scalar, 3, 1> from_skew_symmetric(const Eigen::DenseBase<Derived>& M);
}
