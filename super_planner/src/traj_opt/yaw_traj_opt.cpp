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

#include "traj_opt/yaw_traj_opt.h"
#include <utils/optimization/polynomial_interpolation.h>
using namespace geometry_utils;
namespace traj_opt {
    using namespace color_text;
    /**
     * 作用：根据整段运动的总时长，为偏航轨迹的分段插值分配每一段持续时间，
     * 让起始转向和结束转向都保留足够时间满足最大允许角速度约束。

     * 输入：
     * `duration` 表示整条偏航轨迹可使用的总时长，这里通常与位置轨迹的总时长一致。
     * `times` 是输出用的时间分配数组，函数会重写其中内容，使每个元素对应一段偏航多项式的持续时间。

     * 输出：
     * 函数没有显式返回值；它通过修改 `times` 输出时间分段结果，
     * 这个结果会被后续偏航航点插值直接使用，从而影响整条偏航轨迹的形状与可行性。
     */
    void YawTrajOpt::getYawTimeAllocation(const double &duration, VecDf &times) const {
        // 以“转过 pi 弧度所需时间”作为标准分段长度，这个时间由偏航角速度上限决定。
        double interp_dt = M_PI / yaw_dot_max_;
        if (duration < interp_dt * 2) {
            // 总时长过短时，不再插入中间航点，整段偏航轨迹只保留一个时间段。
            /// if the duration less than 2 interp, then no waypoint need.
            times.resize(1);
            times[0] = duration;
        } else {
            // 总时长足够时，首尾各预留一个标准转向时间，中间剩余部分再均匀拆成若干段。
            /// use ceil to make the seg num non-zero and, make sure the last seg
            /// have time to turn to the goal yaw.
            
            // interp_num 表示中间需要插入多少个过渡时间段；这里向上取整，避免中间段数量变成 0，
            // 同时保证尾段仍然能保留一个标准转向时间。
            int interp_num = ceil((duration - 2 * interp_dt) / interp_dt);
            // interp_t 表示扣除首尾两段固定时间后，中间每一段平均分到的时长。
            double interp_t = (duration - 2 * interp_dt) / (interp_num);
            // 总段数 = 首段 1 个 + 中间过渡段 interp_num 个 + 末段 1 个。
            times.resize(2 + interp_num);
            for (int i = 0; i < interp_num; i++) {
                times(i + 1) = interp_t;
            }
            // 第 0 段和最后一段固定为标准转向时间，中间各段保持均分。
            times(times.size() - 1) = interp_dt;
            times(0) = interp_dt;
        }
        // 如果恰好分成三段且中间段过短，就退化成两段等长分配，避免出现过小的中间时间片。
        if (times.size() == 3 && times(1) < times(0) / 3) {
            times.resize(2);
            times.setConstant(duration / 2);
        }
    }

    /**
     * 作用：根据已经生成的位置轨迹，提取偏航插值所需的中间朝向航点，
     * 使后续偏航多项式轨迹整体朝向位置轨迹的前进方向。
     * 输入：
     * `init_state` 表示偏航轨迹的起点状态，其中第 0 个分量是起始偏航角，后续分量是起始角速度等边界条件。
     * `goal_state` 表示偏航轨迹的终点状态，其中第 0 个分量是目标偏航角，函数会在末尾对这个目标角做连续化修正。
     * `way_pts` 是输出用的中间偏航航点数组，函数会根据位置轨迹重新写入其中内容。
     * `times` 表示偏航轨迹各时间段的持续时间，函数会按这些时间分界点去位置轨迹上取样。
     * `pos_traj` 表示已经确定好的位置轨迹，函数通过它在相邻时刻的位置变化方向推导当前应朝向的偏航角。
     * 输出：
     * 函数没有显式返回值；它通过修改 `way_pts` 产出中间偏航航点，
     * 并通过修正 `goal_state` 中的终点偏航角，给后续多项式插值提供连续、可用的偏航约束。
     */
    void YawTrajOpt::getYawWaypointAllocation(const Vec4f &init_state, Vec4f &goal_state, VecDf &way_pts, VecDf &times,
                                              const Trajectory &pos_traj) {
        double eval_t = 0;
        vec_Vec3f debug1, debug2;
        // 从起始偏航角出发，逐个生成后续中间航点，并保证相邻偏航角在数值上连续。
        double last_yaw = init_state(0);
        way_pts.resize(times.size() - 1);
        const double pos_traj_duration = pos_traj.getTotalDuration();

        // 根据位置轨迹在每个时刻 eval_t（eval_t 是按照 yaw_dot_max 进行分配后的时间间隔） 获取位置，并朝前看一小段时间（例如 0.5s）获取前方位置，
        // 通过两者的位移方向推导当前时刻应当面对的偏航角。
        for (long int i = 0; i < times.size() - 1; i++) {
            eval_t += times(i);
            double cur_yaw;
            Vec3f pt_i = pos_traj.getPos(eval_t);
            Vec3f pt_g;
            // 在当前时间点朝前查看一小段位置变化，用这段位移方向来估计此处应当面对的朝向。
            // 若已经接近位置轨迹末尾，则改为使用末端前后两个点来估计最后的前进方向。
            if (eval_t + 0.5 >= pos_traj_duration) {
                pt_g = pos_traj.getPos(pos_traj_duration);
                pt_i = pos_traj.getPos(pos_traj_duration - 0.5 > 0 ? pos_traj_duration - 0.5 : 0);
            } else {
                pt_g = pos_traj.getPos(eval_t + 0.5);
            }


            Vec3f dir = pt_g - pt_i;
            if (dir.norm() > 0.1) {
                // 当这一小段位移足够明显时，直接把平面运动方向转换成偏航角，
                // 再把新偏航角调整到与上一时刻尽量接近，避免角度在 pi 附近发生跳变。
                cur_yaw = atan2(dir.y(), dir.x());
                // atan2 返回的角度会落在 [-pi, pi]，因此真实朝向连续跨过边界时，
                // 数值上可能出现“假跳变”。
                // 例如上一时刻是 3.13 rad，而当前时刻算出 -3.13 rad，
                // 物理上这两个朝向几乎相同，但数值上相差接近 2pi。
                // 这里通过加减 2pi 的整数倍，把当前偏航角改写成与上一时刻最接近的等价表示，
                // 例如把 -3.13 调整成 3.15，使后续插值处理的是连续角度序列。
                normalizeNextYaw(last_yaw, cur_yaw);
            } else {
                // 当局部位移过小、方向不可靠时，沿用上一时刻的偏航角，避免在接近静止时产生抖动。
//                    print(fg(color::indian_red),
//                          " -- [SUPER] Yaw planning failed, the goal yaw is too close to the current yaw.\n");
                cur_yaw = last_yaw;
            }
            way_pts(i) = cur_yaw;
            last_yaw = cur_yaw;
        }
        // 最后再把终点偏航角调整到与前一朝向连续，避免终点约束与中间航点之间出现 2pi 跳变。
        if (way_pts.size() == 0) {
            geometry_utils::normalizeNextYaw(init_state[0], goal_state[0]);
        } else {
            geometry_utils::normalizeNextYaw(way_pts(way_pts.size() - 1), goal_state[0]);
        }
//            print("Remain dis = {}.\n", (cur_drone_state_.position - gi_.mission_waypoints.back()).norm());
//            print("Yaw way_pts init = {}\n", init_state(0));
//            for (long unsigned int i = 0; i < way_pts.size(); i++) {
//                print("Yaw way_pts[{}] = {}\n", i, way_pts[i]);
//            }
//            print("Yaw way_pts goal = {}\n", goal_state(0));
//            cout << times.transpose() << endl;
    }

    YawTrajOpt::YawTrajOpt(const double &_yaw_dot_max) : yaw_dot_max_(_yaw_dot_max) {
    }

    /**
     * 作用：根据已经确定好的位置轨迹生成一条偏航轨迹，
     * 使朝向整体跟随位置轨迹的前进方向，并满足给定的起终边界条件与插值阶次要求。
     * 输入：
     * `istate_in` 表示偏航轨迹的起点边界状态，其中包含起始偏航角以及角速度等导数约束。
     * `gstate_in` 表示偏航轨迹的终点边界状态，其中包含目标偏航角以及角速度等导数约束。
     * `pos_traj` 表示已经生成好的位置轨迹，函数会用它的总时长和局部运动方向来构造偏航轨迹。
     * `out_traj` 是输出轨迹对象，函数成功时会把生成好的偏航轨迹写入这里。
     * `order` 表示偏航插值使用的多项式阶次，不同阶次对应不同数量的边界导数约束。
     * `free_start` 表示是否放松起点偏航角约束；若为真，起点朝向将由位置轨迹起始方向自动推导。
     * `free_goal` 表示是否放松终点偏航角约束；若为真，终点朝向将由位置轨迹末端方向自动推导。
     * 输出：
     * 函数返回 `true` 表示成功生成偏航轨迹，返回 `false` 表示阶次不支持等情况。
     * 同时它会修改 `out_traj`，并让输出偏航轨迹的起始世界时间与输入位置轨迹保持一致。
     */
    bool YawTrajOpt::optimize(const Vec4f &istate_in,
                              const Vec4f &gstate_in,
                              const Trajectory &pos_traj,
                              Trajectory &out_traj,
                              const int & order,
                              const bool &free_start,
                              const bool &free_goal) {
        free_goal_ = free_goal;
        Vec4f init_state = istate_in;
        Vec4f goal_state = gstate_in;
        double pos_traj_dur = pos_traj.getTotalDuration();
        // 当起点朝向不强制固定时，从位置轨迹起点附近的运动方向估计一个更自然的起始偏航角。
        if (free_start) {
            Vec3f pt_i = pos_traj.getPos(0);

            double t_g = pos_traj_dur > 0.5 > 0 ? 0.5 : pos_traj_dur;
            Vec3f pt_g = pos_traj.getPos(t_g);
            Vec3f dir = pt_g - pt_i;
            while (dir.norm() < 0.5 && t_g < pos_traj_dur) {
                t_g += 0.1;
                pt_g = pos_traj.getPos(t_g);
                dir = pt_g - pt_i;
            }
            init_state(0) = atan2(dir.y(), dir.x());
        }
        // 当终点朝向不强制固定时，从位置轨迹末端附近的运动方向估计终点偏航角。
        if (free_goal_) {
            Vec3f pt_g = pos_traj.getPos(pos_traj_dur);
            double t_i = pos_traj_dur - 0.5 > 0 ? pos_traj_dur - 0.5 : 0;
            Vec3f pt_i = pos_traj.getPos(t_i);
            Vec3f dir = pt_g - pt_i;
            while (dir.norm() < 0.5 && t_i > 0) {
                t_i -= 0.1;
                pt_i = pos_traj.getPos(t_i);
                dir = pt_g - pt_i;
            }
            goal_state(0) = atan2(dir.y(), dir.x());
        }

        // 先根据位置轨迹总时长生成偏航时间分段，再在这些分段边界上提取中间朝向航点。
        VecDf times;
        getYawTimeAllocation(pos_traj_dur, times);
        VecDf way_pts;
        getYawWaypointAllocation(init_state, goal_state, way_pts, times, pos_traj);
        Trajectory yaw_traj;
        // 根据外部指定的阶次选择不同的多项式插值模型。
        // 三阶使用偏航角和角速度边界，五阶再加入角加速度，七阶则使用完整四维边界状态。
        switch (order) {
            case 3: {
                const super_utils::Vec2f init_state3 = init_state.head(2);
                const super_utils::Vec2f goal_state3 = goal_state.head(2);
                yaw_traj = poly_interpo::minimumAccInterpolation<1>(init_state3,
                                                                    goal_state3,
                                                                    way_pts,
                                                                    times);
                break;
            }
            case 5: {
                Vec3f init_state3 = init_state.head(3);
                Vec3f goal_state3 = goal_state.head(3);
                yaw_traj = poly_interpo::minimumJerkInterpolation<1>(init_state3,
                                                                     goal_state3,
                                                                     way_pts,
                                                                     times);
                break;
            }
            case 7: {
                yaw_traj = poly_interpo::minimumSnapInterpolation<1>(init_state,
                                                                     goal_state,
                                                                     way_pts,
                                                                     times);
                break;
            }
            default: {
                cout << "Unsupported order for yaw trajectory optimization." << endl;
                return false;
            }
        }

//            for (double eval_t = 0; eval_t < yaw_traj.getTotalDuration(); eval_t += 0.01) {
//                cout << yaw_traj.getPos(eval_t)[0] << " ";
//            }
//            cout << ";" << endl;
//            for (double eval_t = 0; eval_t < yaw_traj.getTotalDuration(); eval_t += 0.01) {
//                cout << yaw_traj.getVel(eval_t)[0] << " ";
//            }
//            cout << endl;
//
//            yaw_traj.printProfile();
        // 轨迹生成后检查最大偏航角速度是否明显超过限制；当前实现仅打印告警，不直接返回失败。
        double max_yaw_rate = yaw_traj.getMaxVelRate();
        if (max_yaw_rate > yaw_dot_max_ + 2.0) {
            cout << YELLOW << " Yaw rate too large, " << max_yaw_rate << RESET << endl;
//                return false;
        }
        // 输出生成好的偏航轨迹，并继承位置轨迹的起始世界时间，保证二者时间基准一致。
        out_traj = yaw_traj;
        out_traj.start_WT = pos_traj.start_WT;
        return true;
    }
}
