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

#include <super_core/super_planner.h>
#include <memory>
#include <super_utils/scope_timer.hpp>
#include <fmt/color.h>

using namespace super_utils;

namespace super_planner {
    SuperPlanner::SuperPlanner
            (const std::string &cfg_path,
             const ros_interface::RosInterface::Ptr &ros_ptr,
             const rog_map::ROGMapROS::Ptr &map_ptr
            ) : cfg_(Config(cfg_path)), ros_ptr_(ros_ptr), map_ptr_(map_ptr) {

        ros_ptr_->setResolution(cfg_.resolution);
        ros_ptr_->setVisualizationEn(cfg_.visualization_en);
        exp_traj_opt_ = std::make_shared<traj_opt::ExpTrajOpt>(cfg_.exp_traj_cfg, ros_ptr_);
        back_traj_opt_ = std::make_shared<traj_opt::BackupTrajOpt>(cfg_.back_traj_cfg, ros_ptr_);
        yaw_traj_opt_ = std::make_shared<traj_opt::YawTrajOpt>(cfg_.yaw_dot_max);
        const auto &rog_map_cfg = map_ptr_->getMapConfig();
        astar_ptr_ = std::make_shared<path_search::Astar>(cfg_path, ros_ptr_, map_ptr_);
        cg_ptr_ = std::make_shared<CorridorGenerator>(ros_ptr_, map_ptr_, cfg_.corridor_bound_dis,
                                                      cfg_.corridor_line_max_length,
                                                      cfg_.resolution, rog_map_cfg.virtual_ground_height,
                                                      rog_map_cfg.virtual_ceil_height,
                                                      cfg_.robot_r,
                                                      cfg_.obs_skip_num,
                                                      cfg_.iris_iter_num);
        cg_ptr_->SetLineNeighborList(cfg_.seed_line_neighbour);


        time_consuming_.resize(8);

        robot_state_.rcv = false;
        planner_process_start_WT_ = ros_ptr_->getSimTime();
        fov_checker_ = std::make_shared<FOVChecker>(FOVType::OMNI,
                                                    -1.0,
                                                    -35.0,
                                                    35.0);

        const int neighbor_step = floor(cfg_.robot_r / cfg_.resolution);
        astar_ptr_->setFineInfNeighbors(neighbor_step);
    }

    RET_CODE
    SuperPlanner::PlanFromRest(const Vec3f &goal_p,
                               const double &goal_yaw,
                               const bool &new_goal) {
        std::lock_guard<std::mutex> guard(replan_lock_);
        latest_replan.reset();
        latest_replan.setGoal(goal_p, goal_yaw, robot_state_);
        if (robot_state_.rcv == false) {
            ros_ptr_->warn(" -- [SUPER] in [PlanFromRest]: No odom, force return.");
            latest_replan.setRetCode(SUPER_RET_CODE::SUPER_NO_ODOM);
            return FAILED;
        }
        gi_.goal_p = goal_p;
        gi_.goal_yaw = goal_yaw;
        gi_.new_goal = new_goal;
        gi_.goal_valid = true;
        vec_Vec3f viz_pts{goal_p, robot_state_.p};

        {
            TimeConsuming t_viz("viz goal path", false);
            ros_ptr_->vizGoalPath(viz_pts);
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }


        /// 1) First, shift the start_point to free space.
        Vec3f local_star_pt;
        if (!map_ptr_->getNearestCellNot(GridType::OCCUPIED, robot_state_.p, local_star_pt, 3.0)) {
            ros_ptr_->error(
                    " -- [SUPER] in [PlanFromRest] Local start point is deeply occupied, which should not happened.");
            latest_replan.setRetCode(SUPER_RET_CODE::SUPER_NO_START_POINT);
            return FAILED;
        }
        latest_replan.setLocalStartP(local_star_pt);

        /// 2) Generate Exp traj
        ExpTraj exp_traj_info;
        BackupTraj back_traj_info;
        last_exp_traj_info_.setEmpty();
        local_start_p_ = local_star_pt;
        RET_CODE exp_ret_code = generateExpTraj(last_exp_traj_info_, exp_traj_info);
        //GenerateRestToRestExpTraj(local_star_pt, exp_traj_info);
        if (exp_ret_code == FAILED) {
            ros_ptr_->warn(" -- [SUPER] in [PlanFromRest] GenerateExpTrajectory failed with {}.",
                           RET_CODE_STR[exp_ret_code].c_str());
            return FAILED;
        } else {
            ros_ptr_->info(" -- [SUPER] in [PlanFromRest] GenerateExpTrajectory SUCCESS.");
        }

        back_traj_info.setEmpty();
        RET_CODE back_ret_code = generateBackupTrajectory(exp_traj_info, back_traj_info);;

        if (back_ret_code == SUCCESS) {
            if (cfg_.print_log) {
                ros_ptr_->info(" -- [SUPER] in [PlanFromRest] generateBackupTrajectory SUCCESS.");
            }

            cmd_traj_info_.setTrajectory(exp_traj_info, back_traj_info);
            last_exp_traj_info_ = exp_traj_info;
            robot_on_backup_traj_ = false;
            gi_.new_goal = false;

            // For visualization
            {
                TimeConsuming t_viz("viz goal VisualizeCommitTrajectory", false);
                ros_ptr_->vizCommittedTraj(cmd_traj_info_.posTraj(), cmd_traj_info_.getBackupTrajStartTT());
                time_consuming_[VISUALIZATION] += t_viz.stop();
                latest_replan.setRetCode(SUPER_RET_CODE::SUPER_SUCCESS_WITH_BACKUP);
            }

            return SUCCESS;
        } else if (back_ret_code == FINISH || back_ret_code == NO_NEED) {
            if (cfg_.print_log) {
                ros_ptr_->info(" -- [SUPER] in [PlanFromRest] generateBackupTrajectory Finish or NO_NEED.");
            }
            robot_on_backup_traj_ = false;
            cmd_traj_info_.setTrajectory(exp_traj_info);
            last_exp_traj_info_ = exp_traj_info;
            gi_.new_goal = false;

            // For visualization
            TimeConsuming t_viz("viz goal VisualizeCommitTrajectory", false);
            {
                ros_ptr_->vizCommittedTraj(cmd_traj_info_.posTraj(), -1);
                time_consuming_[VISUALIZATION] += t_viz.stop();
            }
            latest_replan.setRetCode(SUPER_RET_CODE::SUPER_SUCCESS_NO_BACKUP);
            return SUCCESS;
        }
        ros_ptr_->warn(" -- [SUPER] in [PlanFromRest] generateBackupTrajectory return [{}], force return",
                       RET_CODE_STR[back_ret_code].c_str());
        return FAILED;
    }


    RET_CODE
    SuperPlanner::ReplanOnce(const Vec3f &goal_p,
                             const double &goal_yaw,
                             const bool &new_goal) {
        TimeConsuming replan_total_t("ReplanOnce", false);
        std::lock_guard<std::mutex> guard(replan_lock_);

        gi_.goal_p = goal_p;
        gi_.goal_yaw = goal_yaw;
        gi_.new_goal = new_goal;
        gi_.goal_valid = true;
        latest_replan.reset();
        latest_replan.setGoal(goal_p, goal_yaw, robot_state_);

        vec_Vec3f viz_pts{goal_p, robot_state_.p};

        {
            TimeConsuming t_viz("tviz", false);
            ros_ptr_->vizGoalPath(viz_pts);
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }


        /// 1) Replan EXP traj
        ExpTraj exp_traj_info;
        TimeConsuming t_exp("t_exp", false);
        RET_CODE exp_ret_code = generateExpTraj(last_exp_traj_info_, exp_traj_info);
        time_consuming_[GENERATE_EXP_TRAJ] = t_exp.stop();

        if (exp_ret_code == FAILED) {
            ros_ptr_->warn(" -- [SUPER] in [ReplanOnce]: GenerateExpTrajectory failed, force return");
            return FAILED;
        } else if (exp_ret_code == NEW_TRAJ) {
            if (cfg_.print_log) {
                ros_ptr_->info(" -- [SUPER] in [ReplanOnce]: Last epx traj end, switch to new traj.");
            }
            return NEW_TRAJ;
        } else if (exp_ret_code == EMER) {
            ros_ptr_->warn(" -- [SUPER] in [ReplanOnce]: Replan failed, switch to emer.");
            return EMER;
        } else if (exp_ret_code == SUCCESS) {
            if (cfg_.print_log) {
                ros_ptr_->info(" -- [SUPER] in [ReplanOnce]: Replan a new exp traj success.");
            }
        } else if (exp_ret_code == NO_NEED) {
            if (cfg_.print_log)
                ros_ptr_->info(" -- [SUPER] in [ReplanOnce]: No need to replan a new exp traj, use last one.");
        }

        {
            TimeConsuming t_viz("tviz", false);
            ros_ptr_->vizYawTraj(exp_traj_info.posTraj(), exp_traj_info.yawTraj());
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }


        BackupTraj back_traj_info;
        // 2）生成back轨迹
        TimeConsuming t_back("t_back", false);
        RET_CODE back_ret_code = generateBackupTrajectory(exp_traj_info, back_traj_info);
        time_consuming_[GENERATE_BACK_TRAJ] = t_back.stop();

        {
            ft += time_consuming_[EPX_TRAJ_FRONTEND] + time_consuming_[BACK_TRAJ_FRONTEND];
            ft_cnt++;
            bt += time_consuming_[BACK_TRAJ_OPT] + time_consuming_[EXP_TRAJ_OPT];
            bt_cnt++;
        }

        double replan_dt = replan_total_t.stop();
        if (replan_dt > cfg_.replan_forward_dt * 0.9) {
            ros_ptr_->warn(" -- [SUPER] in [ReplanOnce]: Replan overtime, check parameters, replan dt = {}.", replan_dt);
            return FAILED;
        }

        if (back_ret_code == SUCCESS) {
            cmd_traj_info_.setTrajectory(exp_traj_info, back_traj_info);
            last_exp_traj_info_ = exp_traj_info;
            robot_on_backup_traj_ = false;
            gi_.new_goal = false;

            {
                // For visualization
                TimeConsuming t_viz("tviz", false);
                ros_ptr_->vizCommittedTraj(cmd_traj_info_.posTraj(), cmd_traj_info_.getBackupTrajStartTT());
                time_consuming_[VISUALIZATION] += t_viz.stop();
            }

            latest_replan.setRetCode(SUPER_SUCCESS_WITH_BACKUP);
            if (cfg_.print_log)
                ros_ptr_->info(" -- [SUPER] in [ReplanOnce]: Replan a new back traj success, all replan success.");
            return SUCCESS;
        } else if (back_ret_code == NO_NEED) {
            // 这次生成backup轨迹的点没有意义,
            robot_on_backup_traj_ = false;
            last_exp_traj_info_ = exp_traj_info;
            gi_.new_goal = false;


            {
                TimeConsuming t_viz("tviz", false);
                ros_ptr_->vizCommittedTraj(cmd_traj_info_.posTraj(), -1);
                time_consuming_[VISUALIZATION] += t_viz.stop();

            }

            if (cfg_.print_log)
                ros_ptr_->info(" -- [SUPER] in [ReplanOnce]: No need back traj success, all replan success.");
            latest_replan.setRetCode(SUPER_SUCCESS_NO_BACKUP);
            return SUCCESS;
        } else if (back_ret_code == FINISH) {
            // Which means the exp traj is all in known free, no need for backup traj
            cmd_traj_info_.setTrajectory(exp_traj_info);
            last_exp_traj_info_ = exp_traj_info;
            robot_on_backup_traj_ = false;
            gi_.new_goal = false;

            {
                TimeConsuming t_viz("tviz", false);
                ros_ptr_->vizCommittedTraj(cmd_traj_info_.posTraj(), -1);
                time_consuming_[VISUALIZATION] += t_viz.stop();
            }

            if (cfg_.print_log)
                ros_ptr_->info(" -- [SUPER] in [ReplanOnce]: No need back traj success, all replan success.");
            latest_replan.setRetCode(SUPER_SUCCESS_NO_BACKUP);
            return SUCCESS;
        }
        ros_ptr_->warn(" -- [SUPER] in [ReplanOnce]: generateBackupTrajectory return {}, replan Failed return",
                       RET_CODE_STR[back_ret_code].c_str());
        return FAILED;
    }

    void SuperPlanner::getOneHeartbeatTime(double &start_WT_pos, bool &traj_finish) {
        double eval_t = (ros_ptr_->getSimTime() - cmd_traj_info_.getStartWallTime());
        traj_finish = false;
        double total_dur = cmd_traj_info_.getTotalDuration();
        if (eval_t > total_dur) {
            traj_finish = true;
            eval_t = total_dur;
        }
        start_WT_pos = cmd_traj_info_.getStartWallTime();
        if (cmd_traj_info_.backupTrajAvilibale() && eval_t > cmd_traj_info_.getBackupTrajStartTT()) {
            robot_on_backup_traj_ = true;
        } else {
            robot_on_backup_traj_ = false;
        }
    }

    Trajectory SuperPlanner::getCommittedPositionTrajectory() {
        return cmd_traj_info_.posTraj();
    }

    Trajectory SuperPlanner::getCommittedYawTrajectory() {
        return cmd_traj_info_.yawTraj();
    }


    void SuperPlanner::getOneCommandFromTraj(StatePVAJ &pvaj,
                                             double &yaw,
                                             double &yaw_dot,
                                             bool &on_backup_traj,
                                             bool &traj_finish) {
        cmd_traj_info_.lock();
        const double &cur_t = ros_ptr_->getSimTime();
        const double &cmd_start_WT = cmd_traj_info_.getStartWallTime();
//        const bool &backup_avilibale = cmd_traj_info_.backupTrajAvilibale();
//        const double &backup_start_TT = cmd_traj_info_.getBackupTrajStartTT();
        const double &total_dur = cmd_traj_info_.getTotalDuration();

        traj_finish = (cur_t - cmd_start_WT) > total_dur;
        const double &eval_t = traj_finish ? total_dur : (cur_t - cmd_start_WT);

//        bool last_round_robot_on_backup_traj = robot_on_backup_traj_;
        robot_on_backup_traj_ = cmd_traj_info_.isTTOnBackupTraj(eval_t);
        on_backup_traj = robot_on_backup_traj_;

        pvaj = cmd_traj_info_.posTraj().getState(eval_t);


        /// Get Yaw planning
        static double last_yaw = robot_state_.yaw;

        yaw = cmd_traj_info_.getYaw((eval_t))[0];
        yaw_dot = cmd_traj_info_.getYawRate((eval_t))[0];

        if (std::isnan(yaw)) {
            yaw = last_yaw;
            yaw_dot = 0;
        } else {
            last_yaw = yaw;
        }
        if (std::isnan(yaw_dot)) {
            yaw_dot = 0;
        }

//        if (last_round_robot_on_backup_traj != robot_on_backup_traj_) {
//            if (last_round_robot_on_backup_traj) {
//                ros_ptr_->info(" -- [CMD] Emergency Stop End ========================");
//            } else {
//                ros_ptr_->info(" -- [CMD] Emergency Stop Start ========================");
//            }
//        }

//        double cur_yaw = geometry_utils::get_yaw_from_quaternion(robot_state_.q);
        cmd_traj_info_.unlock();
    }


    void SuperPlanner::getModuleTimeConsuming(vector<double> &time) {
        time = time_consuming_;
        std::fill(time_consuming_.begin(), time_consuming_.end(), 0);
    }


    RET_CODE SuperPlanner::generateExpTraj(ExpTraj &last_exp_traj_info, ExpTraj &out_exp_traj_info) {
        /* 1) Log the exp traj frontend time*/
        // 记录 exp 轨迹前端阶段的耗时，前端主要包括引导路径和安全走廊的构造。
        TimeConsuming t_exp_frontend("t_exp_frontend", false);

        // use hot init or not, just prepare a guide path, a guide t, init and fina state and sfc for exp traj opt
        // 位置轨迹优化的起点/终点状态约束，列通常对应 P/V/A/J。
        StatePVAJ pos_init_state, pos_fina_state;
        // Safe Flight Corridor，后续位置优化会被约束在这组走廊多面体内。
        PolytopeVec sfc;
        // 引导路径，既可来自旧轨迹的可复用段，也可来自新的路径搜索结果。
        vec_Vec3f guide_path;
        // the guide_stamp saves a TT
        // guide_path 中各点对应的轨迹时间戳 TT，供热启动和时间分配使用。
        vector<double> guide_stamp;
        // 引导路径末端速度，用于给后续新搜索路径做时间分配时保持速度连续性。
        double guide_path_end_vel{0.0};

        // 按规划视野长度 / 地图分辨率粗略估算可能的路径点数，并额外预留约 20% 余量。
        // reserve 只分配容量、不改变当前 size，目的是减少后续 push_back 时的动态扩容开销。
        int reserve_size = cfg_.planning_horizon / cfg_.resolution * 1.2;
        guide_path.reserve(reserve_size);
        guide_stamp.reserve(reserve_size);

        // yaw 优化的起点状态默认取当前机器人 yaw；终点先置零，后续再按目标决定是否约束。
        Vec4f init_yaw{robot_state_.yaw, 0, 0, 0};
        Vec4f fina_yaw{0, 0, 0, 0};


        // alias for last_exp_traj_info
        // 这三个轨迹变量分别指向：
        // 1) current commited trajectory 中的位置轨迹
        // 2) current commited trajectory 中的偏航轨迹
        // 3) 上一轮生成的探索轨迹
        Trajectory guide_pos_traj, guide_yaw_traj, last_exp_traj;

        // record the wall time (WT) and the trajectory time (TT) at the start of the replan.
        // WT 是世界/仿真时间；TT 是某条轨迹自己的相对时间坐标。
        const double replan_process_start_WT = ros_ptr_->getSimTime();
        // replan_process_start_TT: 当前时刻落在旧轨迹上的 TT(Trajectory Time)
        // replan_state_TT: 真正用于新轨迹接管的 TT(Trajectory Time)，通常会前瞻一个 replan_forward_dt
        double replan_process_start_TT, replan_state_TT;

        /* 2) Check last exp traj */
        if (last_exp_traj_info.empty()) {
            /* 2.1) Perform rest2rest exp traj generation */
            // just skip the first part of the guide trajectory
            pos_init_state.setZero();
            // 设置轨迹起点
            pos_init_state.col(0) = local_start_p_;
            replan_process_start_TT = -1;
            replan_state_TT = -1;
        } else {
            guide_pos_traj = cmd_traj_info_.posTraj(); // last_exp_traj;
            guide_yaw_traj = cmd_traj_info_.yawTraj(); //last_exp_traj_info.exp_yaw_traj;
            last_exp_traj = last_exp_traj_info.posTraj();
            
            // 当前这一刻落在旧轨迹内部的哪个轨迹时间点上。
            replan_process_start_TT = replan_process_start_WT - last_exp_traj.start_WT;
            // 新轨迹并不从“现在”立刻接管，而是从向前看一个 replan_forward_dt 后的时刻开始接管。
            // 预留一个前瞻时间窗口，给规划、通信和控制衔接留余量。
            replan_state_TT = replan_process_start_TT + cfg_.replan_forward_dt;
            /* 2.2) Perform collision check on last exp traj*/
            // 沿旧 commited 轨迹采样后保留下来的“时间-位置”序列。
            // 后续一方面用它判断哪些旧轨迹点仍可安全复用，另一方面用它来构造新 guide_path 的前半段。
            vector<TimePosPair> last_exp_traj_time_pos;
            // 与 last_exp_traj_time_pos 一一对应的速度模长序列。
            // 后续在确定 guide_path 末端速度、以及给新搜索路径做时间分配时会用到。
            vector<double> last_exp_traj_vel;


            // check early exit condition
            // 1) if the replan state is beyond the last cmd traj, return NO_NEED
            // 如果“新轨迹的接管时刻”已经超过当前 commited 轨迹的末尾，就没法再从这条 commited 轨迹上接管了。
            if (replan_state_TT >= cmd_traj_info_.getTotalDuration()) {
                out_exp_traj_info = last_exp_traj_info;

                // 如果机器人此时已经在 backup 轨迹上，继续等待可能不安全，直接报失败，交给更保守的流程处理。
                if (robot_on_backup_traj_) {
                    // 仅在开启日志时打印告警信息。
                    if (cfg_.print_log)
                        ros_ptr_->warn(
                                " -- [SUPER] Replan, emergency stop, return FAILED and wait for plan form rest.");
                    return FAILED;
                }

                // 正常情况下这里只是说明“不需要基于当前轨迹继续重规划”，打印提示后返回 NO_NEED。
                if (cfg_.print_log) {
                    ros_ptr_->warn(
                            " -- [generateExpTraj] replan_state_TT >= cmd_traj_info_.pos_traj.getTotalDuration(), return NONEED and wait for plan form rest.");
                }
                return NO_NEED;
            }

            // 有旧探索轨迹时，再进一步检查若干“可以直接提前退出”的条件。
            if (!last_exp_traj_info.empty()) {
                // 如果接管时刻已经落到旧 exp 轨迹末尾之外，旧 exp 也没有继续复用的价值了。
                if (replan_state_TT >= last_exp_traj.getTotalDuration()) {
                    out_exp_traj_info = last_exp_traj_info;
                    // 仅在开启日志时打印告警信息。
                    if (cfg_.print_log)
                        ros_ptr_->warn(
                                " -- [generateExpTraj] replan_state_TT >= last_exp_traj.getTotalDuration(), return NONEED and wait for plan form rest.");
                    // 如果当前还处在 backup 上，则继续拖延不安全，直接失败；否则仅表示这次无需重规划。
                    if (robot_on_backup_traj_) {
                        // 仅在开启日志时打印告警信息。
                        if (cfg_.print_log)
                            ros_ptr_->warn(
                                    " -- [SUPER] Replan, emergency stop, return FAILED and wait for plan form rest.");
                        return FAILED;
                    } else {
                        return NO_NEED;
                    }
                }

                /// 1) Check a series of early termination conditions.
                // 目标没变，且旧 exp 只对应一个 corridor，并且已经确认连到目标，说明这条轨迹已经足够简单直接，无需再重规划。
                if (!gi_.new_goal && last_exp_traj_info.getSFCSize() == 1 && last_exp_traj_info.connectedToGoal()) {
                    // 仅在开启日志时打印告警信息。
                    if (cfg_.print_log) {
                        ros_ptr_->warn(
                                " -- [SUPER] Replan, last exp have only one corridor and connected to goal return NONEED.");
                    }

                    out_exp_traj_info = last_exp_traj_info;
                    // 如果此时在 backup 上，则不再接受“无事发生”，而是直接按失败处理。
                    if (robot_on_backup_traj_) {
                        // 仅在开启日志时打印告警信息。
                        if (cfg_.print_log)
                            ros_ptr_->warn(
                                    " -- [SUPER] Replan, emergency stop, return FAILED and wait for plan form rest.");
                        return FAILED;
                    } else {
                        return NO_NEED;
                    }
                }

                // 目标没变，且在接管时刻旧轨迹已经非常接近目标，就没必要再花代价重规划一条新 exp。
                if (!gi_.new_goal &&
                    (gi_.goal_p - last_exp_traj.getPos(replan_state_TT)).norm() < cfg_.resolution * 3) {
                    // Return if the traj close to goal
                    out_exp_traj_info = last_exp_traj_info;
                    out_exp_traj_info.setGoalConnectedFlag(true);

                    ros_ptr_->warn(" -- [SUPER] Replan, close to goal and return NONEED.");
                    // 同样地，如果当前落在 backup 上，则宁可判失败也不继续维持这种状态。
                    if (robot_on_backup_traj_) {
                        ros_ptr_->warn(
                                " -- [SUPER] Replan, emergency stop, return FAILED and wait for plan form rest.");
                        return FAILED;
                    } else {
                        return NO_NEED;
                    }
                }
            }

            /// Ready for replan.
            out_exp_traj_info.setGoalConnectedFlag(false);

            // * 2) Check if in backup trajectory. While in backup trajectory,
            // *    the guide trajectory should be a part of cmd trajectory.
            // TODO: Why cannot directly replan on cmd traj? 241121

            // * 3) Perform collision check on the guide trajectory.
            // TODO 0929 critical change for hot init.
            // 从“新轨迹准备接管的时刻”开始向后采样检查，而不是从当前时刻开始。
            double eval_t = replan_state_TT; //replan_process_start_TT;
            // 引导位置轨迹的总时长，后续采样循环不能超过这个上界。
            double guide_pos_traj_total_time = guide_pos_traj.getTotalDuration();

            // temp_pt: 当前采样点；last_sample_pt: 上一个被接受的采样点。
            Vec3f temp_pt, last_sample_pt;
            // 清空上一轮记录的“时间-位置”采样结果，准备重新收集可复用的旧轨迹点。
            last_exp_traj_time_pos.clear();
            // 先默认旧轨迹后半段整体已知安全；若后续采样发现碰撞，再置为 false。
            last_exp_traj_info.setWholeTrajKnownFreeFlag(true);
            // 记录采样起点处的位置
            last_sample_pt = guide_pos_traj.getPos(eval_t);
            // 采样时刻先向前推进一个固定步长，后续循环从下一个采样点开始检查。
            eval_t += cfg_.sample_traj_dt;


            // * 4) 记录replan点在evaluated_pts上的id
            // 记录“第一个位于 replan_state_TT 之后的有效采样点”在采样点数组中的下标；-1 表示尚未记录到。
            int replan_id = -1;
            // 沿 guide_pos_traj 继续向后离散采样，直到轨迹末尾。guide_pos_traj 是当前 commited 轨迹中的位置部分
            // 采样的时间不会超过 guide_pos_traj 的总时长
            for (; eval_t < guide_pos_traj_total_time; eval_t += cfg_.sample_traj_dt) {
                temp_pt = guide_pos_traj.getPos(eval_t);
                // 如果与上一个被接受的采样点距离过近，就跳过，避免采样点过密。
                if ((temp_pt - last_sample_pt).norm() < cfg_.resolution * 0.8) {
                    continue;
                }

                // 在 inflated map 上检查当前采样点的占据状态。
                rog_map::GridType temp_grid = map_ptr_->getInfGridType(temp_pt);

                // 一旦发现点落入障碍物或超出地图范围，说明旧轨迹后半段不再整体安全，停止继续复用。
                if (temp_grid == rog_map::GridType::OCCUPIED || temp_grid == rog_map::GridType::OUT_OF_MAP) {
                    last_exp_traj_info.setWholeTrajKnownFreeFlag(false);
                    break;
                }

                // 第一次遇到“严格位于接管时刻之后”的有效采样点时，记录它在采样数组中的位置。
                if (eval_t > replan_state_TT && replan_id == -1) {
                    replan_id = last_exp_traj_time_pos.size();
                }
                // 保存可复用的旧轨迹采样点：(轨迹时间, 位置)。
                last_exp_traj_time_pos.emplace_back(eval_t, temp_pt);
                // 同时保存该点对应的速度模长，后续时间分配会用到。
                last_exp_traj_vel.emplace_back(guide_pos_traj.getVel(eval_t).norm());
                last_sample_pt = temp_pt;
            }


            // * 6) Decide where to split the original exp trajecory and re-plan a new one with an A*,
            // *    If the whole trajectory is free,  the whole trajectory should be receding and if not, or a new goal
            // *    is given, we should only receiding a small distance and replan new trajectory ASAP
            // 如果旧轨迹后半段并不完全安全，或者目标已经变化，那么就不要长时间沿用旧轨迹，而是只保留一小段必要的旧前缀，
            // 然后尽快重新搜索新路径并切换到新轨迹。也就是“尽快开始新的重规划结果接管”。
           
            // split_dis 表示这次最多保留多少“旧轨迹的可复用前缀距离”。
            // 默认只保留 cfg_.receding_dis 这么长的一小段，然后尽快切到新的路径搜索结果；
            // 但如果旧轨迹后半段整体已知安全、目标也没变，则说明没有必要激进重规划，
            // 于是把 split_dis 设为无穷大，等价于尽量长地复用旧轨迹，仅做 receding-horizon 式平移。
            double split_dis = cfg_.receding_dis;
            if (last_exp_traj_info.wholeTrajKnownFree() && !gi_.new_goal && cfg_.receding_dis > 0.0) {
                split_dis = std::numeric_limits<double>::max();
            }


            // * 7）Begin replan process, first get the replan state from the committed trajectory.
            if (!guide_pos_traj.getState(replan_state_TT, pos_init_state)) {
                ros_ptr_->warn(" -- [SUPER] Invalid traj or eval t");
                return FAILED;
            }

            // * Generate guide path with time stamp, for hot trajectory initialization
            // * the guide stamp is time from the replan start t
            guide_stamp.clear();
            guide_path.clear();
            
            // guide_path_end_vel 有时等于当前机器人速度，有时等于旧轨迹可复用前缀末端的速度；
            // 它统一表示的是“当前 guide_path 最后一个点的速度参考”。

            // 如果不打算保留旧轨迹前缀，或者前面根本没有采到可复用的旧轨迹点，
            // 那么 guide_path 就只从当前重规划起点开始，后续完全依赖新的路径搜索结果来补全。
            if (split_dis <= 0 || last_exp_traj_time_pos.empty()) {
                /// No need receding, just path search.
                guide_path.push_back(pos_init_state.col(0));
                guide_stamp.push_back(0.0);
                last_exp_traj_time_pos.clear();
                last_exp_traj_time_pos.emplace_back(replan_state_TT, pos_init_state.col(0));
                // 此时guide_path只从当前重规划起点开始，不存在复用路径，故将当前机器人速度设置为guide_path_end_vel参考速度
                guide_path_end_vel = robot_state_.v.norm();
            } else {
                // 否则尝试复用一段旧轨迹前缀，但需要先从尾部开始裁剪掉“过远”或“已不安全”的点。
                temp_pt = last_exp_traj_time_pos.back().second;
                // * 8) Pop all evaluated pts after the sampled point.
                // 只要当前末尾点处于膨胀障碍内，或它离重规划起点的距离超过允许复用的 split_dis，
                // 就持续从尾部删除，直到剩下的一段旧轨迹既安全又足够靠近。
                while (map_ptr_->isOccupiedInflate(temp_pt) ||
                       (temp_pt - pos_init_state.col(0)).norm() > split_dis) {
                    last_exp_traj_time_pos.pop_back();
                    last_exp_traj_vel.pop_back();
                    // 如果连一个可复用点都删没了，说明这次旧轨迹前缀完全不可用，只能退回到“从当前起点重新开始”。
                    if (last_exp_traj_time_pos.empty()) {
                        ros_ptr_->warn(" -- [SUPER] WARN, all traj is collide in INF2");
                        break;
                    }
                    temp_pt = last_exp_traj_time_pos.back().second;
                }
                // 若裁剪后仍保留了一段旧轨迹，就把它按顺序写入 guide_path / guide_stamp，
                // 作为后续热启动引导路径的前半段。
                // 将last_exp_traj_time写入guide_path
                if (!last_exp_traj_time_pos.empty()) {
                    for (long unsigned int i = 0; i < last_exp_traj_time_pos.size(); i++) {
                        guide_path.push_back(last_exp_traj_time_pos[i].second);
                        guide_stamp.push_back(last_exp_traj_time_pos[i].first - last_exp_traj_time_pos.front().first);
                        guide_path_end_vel = last_exp_traj_vel[i];
                    }
                } else {
                    // 裁剪后已经没有任何旧轨迹可复用时，回退到最保守方案：guide_path 只保留当前重规划起点。
                    guide_path.push_back(pos_init_state.col(0));
                    guide_stamp.push_back(0.0);
                    last_exp_traj_time_pos.emplace_back(replan_state_TT, pos_init_state.col(0));
                    guide_path_end_vel = robot_state_.v.norm();
                }
            }
        } // 判断轨迹是否为空的else情况截止在此处

        // second, geometry part of the guide path
        ///=================The Second Part of Guide Path ================================================

        double guide_path_length = geometry_utils::computePathLength(guide_path);
        double temp_horizon = cfg_.planning_horizon - guide_path_length;

        vector<int> path_passed_waypoint_id;
        vec_Vec3f inside_poly_goals;
        vector<int> sfc_waypoint_ids;

        // 检查当前 guide_path 是否已经以“本次重规划的起点 pos_init_state.col(0)”作为第一个点：
        // 1) 若 guide_path 为空，说明前面还没有任何可用引导点；
        // 2) 若 guide_path 非空但首点与重规划起点不一致，说明引导路径的起点不正确。
        // 这两种情况下都需要把当前重规划起点补到 guide_path 开头，并同步补一个 0.0 的时间戳。
        if (guide_path.empty() ||
            ((guide_path.front() - pos_init_state.col(0)).norm() > 1e-2)) {
            guide_path.insert(guide_path.begin(), pos_init_state.col(0));
            guide_stamp.insert(guide_stamp.begin(), 0.0);
        }

        // 这里开始构造 guide_path 的第二部分：几何搜索段。
        // 前面构造的 guide_path 主要来自旧轨迹可复用前缀；如果它还没有用完整个 planning_horizon，
        // 就继续从当前 guide_path 末端朝目标补一段新的几何路径。
        // temp_horizon 表示在扣除已复用前缀长度之后，还剩多少搜索预算。
        if (temp_horizon > cfg_.resolution * 2) {
            // start point TT + exp_traj start_WT
            // double path_search_start_point_WT = guide_stamp.back() + guide_pos_traj.start_WT;
            
            // 如果 guide_path 当前末端已经足够接近目标点，就没必要再做一次 A* 搜索，
            // 直接把目标点追加到 guide_path 末尾，并按“距离 / 最大速度”粗略补一个时间戳即可。
            if ((guide_path.back() - gi_.goal_p).norm() < cfg_.resolution * 5) {
                guide_stamp.push_back(guide_stamp.back() +
                                      (guide_path.back() - gi_.goal_p).norm() / cfg_.exp_traj_cfg.max_vel);
                guide_path.push_back(gi_.goal_p);
                // NO NEED
            } else {
                vec_Vec3f new_path;
                // 如果目标太远，就先把目标投影到 planning horizon 内，但后续放弃了这个选项
                // project goal within the planning horizon
//                const Vec3f dir = (gi_.goal_p - robot_state_.p).normalized();
//                const double dis2goal = (gi_.goal_p - robot_state_.p).norm();
//                Vec3f cadi_p = gi_.goal_p;
//                if(dis2goal > cfg_.planning_horizon) {
//                    double proj_l = cfg_.planning_horizon;
//                    Vec3f cadi_p = robot_state_.p + dir * proj_l;
//                    int max_iter = 100;
//                    while(map_ptr_->isOccupiedInflate(cadi_p) && max_iter-- > 0) {
//                        if(map_ptr_->getNearestInfCellNot(OCCUPIED, cadi_p, cadi_p, 1.0)) {
//                            break;
//                        }
//                        proj_l -= 2.0;
//                        if(proj_l < 1){
//                            ros_ptr_->warn(" -- [SUPER] Project goal failed");
//                            gi_.goal_valid = false;
//                            return FAILED;
//                        }
//                        cadi_p = robot_state_.p + dir * proj_l;
//                    }
//                    if(max_iter <= 0) {
//                        ros_ptr_->warn(" -- [SUPER] Project goal failed");
//                        gi_.goal_valid = false;
//                        return FAILED;
//                    }
//                }

                // 这里调用的不是 A* 底层接口本身，而是 SuperPlanner 对前端搜索做的一层封装：
                // 1) 先用 escapePathSearch() 检查当前起点是否需要从 prob map 中“逃逸”到可安全规划的位置；
                // 2) 再用 pointToPointPathSearch() 从当前 guide_path 末端朝目标点搜索几何路径；
                // 3) 默认先在 inf map 上搜索，失败后再回退到 prob map 重试。
                // 返回的 new_path 是“用于补全 guide_path 后半段”的离散几何路径，而不是最终优化后的轨迹。
                if (!PathSearch(guide_path.back(), gi_.goal_p, temp_horizon, new_path)) {
                    ros_ptr_->warn(" -- [SUPER] PathSearch for new path failed");
                    return FAILED;
                }
                if (new_path.size() < 2) {
                    ros_ptr_->warn(" -- [SUPER] PathSearch for new path failed");
                    return FAILED;
                }

                // compute total dis
                // backward compute dis for all points
                double total_dis{0.0};
                // dis[i]表示索引i代表的剩余路径长度
                vector<double> dis(new_path.size());
                Vec3f last_p = new_path.back();
                for (int i = new_path.size() - 2; i >= 0; i--) {
                    auto d = (new_path[i] - last_p).norm();
                    total_dis += d;
                    dis[i+1] = total_dis;
                    last_p = new_path[i];
                }
                // 处理guide_path终点与new_path起点之间的距离
                total_dis += (new_path.front() - guide_path.back()).norm();
                dis[0] = total_dis;

//                for (int i = 0; i < dis.size(); i++) {
//                    cout << dis[i] << " ";
//                }
//                cout << endl;
                vector<double> stamps(new_path.size(), 0);
                vector<double> dt(new_path.size(), 0);
                double last_stamp = 0;

                // 这个时间分配的策略真的是10，先计算出dis[i]表示的剩余路径长度，然后从最后一个点开始
                // 倒着建立运动学模型，原本的减速段编程了加速度，加速度变成了减速段，防御性编程这一块
                for (int i = dis.size() - 1; i >= 0; i--) {
                    double vel;
                    geometry_utils::simplePMTimeAllocator(cfg_.exp_traj_cfg.max_acc, cfg_.exp_traj_cfg.max_vel,
                                                          guide_path_end_vel,
                                                          total_dis,
                                                          dis[i], stamps[i], vel);
                    dt[i] = stamps[i] - last_stamp;
                    last_stamp = stamps[i];
                }
                double time_stamp = guide_stamp.back();

//                for (int i = 0; i < stamps.size(); i++) {
//                    cout << stamps[i] << " ";
//                }
//                cout << endl;
//
//                for (int i = 0; i < dt.size(); i++) {
//                    cout << dt[i] << " ";
//                }
//                cout << endl;

                for (long unsigned int i = 1; i < new_path.size(); i++) {
                    double t = dt[i];
                    time_stamp += t;
                    guide_path.emplace_back(new_path[i]);
                    guide_stamp.emplace_back(time_stamp);
                }
            }
        }

        // 判断当前前端构造出的 guide_path 是否已经“连到目标”：
        // 这里只比较 x-y 平面上的末端距离，若 guide_path 末端与目标点的水平距离小于 2 个分辨率，
        // 就认为这条引导路径已经到达目标附近。该标志后续会影响：
        // 1) ExpTraj 是否被视为已连接到目标；
        // 2) 后续重规划是否可以提前退出；
        // 3) 终点 yaw 是否需要严格对齐目标 yaw。
        const bool connected_goal = (guide_path.back().head(2) - gi_.goal_p.head(2)).norm() < cfg_.resolution * 2;
        out_exp_traj_info.setGoalConnectedFlag(connected_goal);



        // 这里开始将“离散 guide_path”转换为“后端优化可用的安全走廊 SFC”。
        // 前面得到的 guide_path 只是几何引导路径；后端位置优化真正需要的是一串连续的多面体约束。
        sfc.clear();
        {
            TimeConsuming t_viz("tviz", false);
            ros_ptr_->vizFrontendPath(guide_path);
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }
        // 先把 shifted_sfc_start_pt_ 设为一个明显无效的哨兵值。
        // 若 SearchPolytopeOnPath 在构造走廊时发现 guide_path 前段点落在膨胀障碍中，
        // 它会把真正可用的“后移起点”写回这个变量；否则该值保持无效，后续逻辑会自动忽略它。
        shifted_sfc_start_pt_ = Vec3f(9999,9999,9999);
        // 沿 guide_path 生成 Safe Flight Corridor。
        bool bool_ret_code = cg_ptr_->SearchPolytopeOnPath(guide_path, sfc, shifted_sfc_start_pt_, cfg_.use_fov_cut);

        // 若无法沿当前 guide_path 生成连续安全走廊，则后端优化没有可用空间约束，只能直接失败返回。
        if (!bool_ret_code) {
            ros_ptr_->warn(" -- [SUPER] SearchPolytopeOnPath for new path failed");
            return FAILED;
        }
        {
            TimeConsuming t_viz("tviz", false);
            ros_ptr_->vizExpSfc(sfc);
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }



        // 至此，exp 轨迹“前端”阶段结束：guide_path、guide_stamp、SFC 都已准备完成。
        time_consuming_[EPX_TRAJ_FRONTEND] = t_exp_frontend.stop();

        // 设置位置后端优化的终点状态 pos_fina_state：
        // 注意：guide_path.back() 不一定等于真实目标 gi_.goal_p，
        // 因为前端 PathSearch 可能只到达 planning horizon（REACH_HORIZON），而不一定真正到达目标（REACH_GOAL）。
        // 因此后端优化默认先以“当前引导路径末端”为终点，而不是直接假设终点一定是 gi_.goal_p。

        pos_fina_state.setZero();
        pos_fina_state.col(0) = guide_path.back();
        // 若启用了 goal_vel_en，且目标还比较远，则不给终点设为静止，
        // 而是给一个朝目标方向的前进速度，使轨迹在尚未真正到目标时保持向前推进。
        if (cfg_.goal_vel_en && (gi_.goal_p - robot_state_.p).norm() > cfg_.planning_horizon / 2) {
            pos_fina_state.col(1) = (gi_.goal_p - robot_state_.p).normalized() * cfg_.exp_traj_cfg.max_vel / 2;
        }
        // 若当前终点已经足够接近真实目标点，则把终点位置精确钉到 goal_p，
        // 同时将终点速度清零，表示这次轨迹应当以“到达目标并停住”为终点条件。
        if ((pos_fina_state.col(0) - gi_.goal_p).norm() < cfg_.resolution * 2) {
            pos_fina_state.col(1).setZero();
            pos_fina_state.col(0) = gi_.goal_p;
        }

        // optimize and update exp traj
        // ===== 后端位置轨迹优化阶段 =====
        // 使用前端准备好的起终点状态、guide_path、guide_stamp 和 SFC，
        bool temp_ret;
        // 优化生成一条新的位置轨迹后缀 out_traj。
        Trajectory out_traj;
        TimeConsuming t_exp_opt("t_exp_opt", false);
        auto original_sfc = sfc;
        temp_ret = exp_traj_opt_->optimize(pos_init_state,
                                           pos_fina_state,
                                           guide_path,
                                           guide_stamp,
                                           sfc,
                                           out_traj);
        time_consuming_[EXP_TRAJ_OPT] = t_exp_opt.stop();
        {
            // 记录优化器内部使用的初始化值，便于日志分析和调试复现。
            VecDf init_ts;
            vec_Vec3f init_ps;
            exp_traj_opt_->getInitValue(init_ts, init_ps);
            latest_replan.setExpCondition(init_ts, init_ps, pos_init_state, pos_fina_state, sfc);
        }
        // 若位置后端优化失败，则本次 exp 轨迹生成直接失败。
        if (!temp_ret) {
            ros_ptr_->warn(" -- [SUPER] OptimizationExpTrajInPolytopes for new path failed");
            return FAILED;
        }
        // 检查整次重规划耗时是否超过允许的 前瞻接管时间窗口 。
        // 若已超时，则说明新轨迹来不及安全接管，只能判失败。
        double replan_total_t = (ros_ptr_->getSimTime() - replan_process_start_WT);
        if (replan_total_t > cfg_.replan_forward_dt) {
            ros_ptr_->warn(" -- [SUPER] Replan over time({})!!!! Return FAILED", replan_total_t);
            return FAILED;
        }

        {
            TimeConsuming t_viz("tviz", false);
            ros_ptr_->vizExpTraj(out_traj);
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }

        // ===== 拼接最终位置 exp 轨迹：旧前缀 + 新后缀 =====
        //
        // 为什么 guide_path 已经完成了路径的拼接，此处又要进行路径的拼接呢？
        // 因为 guide_path 的注释表示的很清楚，guide_path 的拼接是从 replan_process_start_TT 开始的
        // 而并不是从 replan_process_start_WT 开始的，所以我们需要补全[replan_process_start_TT, replan_state_TT]这一段时间内的轨迹
        // 这便是轨迹的 warm start 策略
        //
        // 新轨迹的起点定义为本次重规划开始时刻。
        // replan_process_start_WT 是在该函数被调用时获取的世界时间
        double new_traj_WT = replan_process_start_WT;

        replan_process_start_TT = replan_process_start_WT - guide_pos_traj.start_WT;
        Trajectory temp_exp_traj;
        // 若之前已有旧轨迹，则先从当前 committed 位置轨迹中截取
        // [replan_process_start_TT, replan_state_TT] 这一小段旧前缀，用于平滑衔接到新位置后缀。
        if (!last_exp_traj_info_.empty() &&
            !guide_pos_traj.getPartialTrajectoryByTime(replan_process_start_TT, replan_state_TT,
                                                       temp_exp_traj)) {
            ros_ptr_->error(" -- [SUPER] in [generateExpTraj]: getPartialTrajectoryByTime failed, force return");
            return FAILED;
        }
        out_exp_traj_info.setSFC(sfc);
        // 将旧前缀与新优化得到的位置后缀拼接成完整 exp 位置轨迹。
        temp_exp_traj = temp_exp_traj + out_traj;
        temp_exp_traj.start_WT = new_traj_WT; //last_exp_traj_info.replan_start_WT ;

        // ===== 准备 yaw 轨迹优化的起点状态 =====
        // 若之前存在旧 yaw 轨迹，则从接管时刻 replan_state_TT 读取 yaw 状态，执行 yaw 轨迹的 warm_start
        if (!last_exp_traj_info.empty()) {
            StatePVAJ yaw_replan_state;
            if (!guide_yaw_traj.getState(replan_state_TT, yaw_replan_state)) {
                ros_ptr_->warn(" -- [SUPER] Invalid traj or eval t");
                return FAILED;
            }
            init_yaw = yaw_replan_state.row(0);
        }


        // ===== 生成 yaw 后缀并与旧 yaw 前缀拼接 =====
        // 默认 yaw 终点自由；只有在启用了 goal_yaw、目标 yaw 有效、且轨迹确实连到目标时，
        // 才强约束终点 yaw 对齐目标 yaw。
        bool free_end{true};
        if (cfg_.goal_yaw_en && !std::isnan(gi_.goal_yaw) && connected_goal) {
            free_end = false;
            fina_yaw[0] = gi_.goal_yaw;
        }
        Trajectory new_traj, old_traj;

        // 针对新位置后缀 out_traj 生成对应的新 yaw 轨迹
        // 这里 yaw 优化必须同时输入位置轨迹 out_traj，而不能只给起终 yaw：
        // 1) yaw 轨迹的总时长需要与位置轨迹总时长保持一致，且对应路径点之间的时间也应该与 out_traj 相同
        // 2) yaw 中间航向点是根据位置轨迹在相邻时刻的位置变化方向自动推导的，
        //    例如通过 pos_traj.getPos(t + dt) - pos_traj.getPos(t) 的方向来计算该时刻应朝向的 yaw；
        // 3) 因此这里生成的 yaw 轨迹，本质上是朝向位置轨迹的前进方向。
        if (!yaw_traj_opt_->optimize(init_yaw, fina_yaw, out_traj, new_traj, 3, false, free_end)) {
            ros_ptr_->error(" -- [SUPER] in [generateExpTraj]: YawTrajOpt failed, force return");
            return FAILED;
        }
        // 若存在旧 yaw 轨迹，则同样截取 [replan_process_start_TT, replan_state_TT] 的旧 yaw 前缀。
        if (!last_exp_traj_info.empty()) {
            if (!guide_yaw_traj.getPartialTrajectoryByTime(replan_process_start_TT, replan_state_TT,
                                                           old_traj)) {
                ros_ptr_->error(" -- [SUPER] in [generateExpTraj]: getPartialTrajectoryByTime failed, force return");
                return FAILED;
            }
        }

        // 旧 yaw 前缀 + 新 yaw 后缀 = 完整 exp yaw 轨迹。
        const auto temp_yaw_traj = old_traj + new_traj;

        // check if part of the exp on last backup
        double on_backup_end_TT{-1}, on_backup_start_TT{-1};
        // 若本次接管时刻已经落入 committed backup 轨迹区间，
        // 则记录新 exp 轨迹前半段中有多少时间实际上仍属于旧 backup 轨迹。
        if (!last_exp_traj_info.empty() && replan_state_TT > cmd_traj_info_.getBackupTrajStartTT()) {
            // on_backup_start_TT = backup开始时刻 - 新exp起点对应的旧轨迹时刻
            // on_backup_end_TT   = 新exp接管时刻 - 新exp起点对应的旧轨迹时刻
            // 这两个变量计算的目的其实是想让后续模块知道：这条新的 exp 轨迹前缀里，有一小段其实继承自旧 backup 轨迹
            on_backup_start_TT = cmd_traj_info_.getBackupTrajStartTT() - replan_process_start_TT;
            on_backup_end_TT = replan_state_TT - replan_process_start_TT;
        }

        // 将位置轨迹、yaw 轨迹及 backup 区间信息写入最终输出的 ExpTraj。
        out_exp_traj_info.setTrajectory(new_traj_WT, temp_exp_traj, temp_yaw_traj, on_backup_start_TT,
                                        on_backup_end_TT);

        // 同步记录到 latest_replan，供调试和日志查看。
        latest_replan.setExpYawTraj(temp_yaw_traj);
        latest_replan.setExpTraj(temp_exp_traj);

        return SUCCESS;
    }

    RET_CODE SuperPlanner::generateBackupTrajectory(ExpTraj &ref_exp_traj, BackupTraj &back_traj_info) {
        drone_state_mutex_.lock();
        back_traj_info.setRobotPos(robot_state_.p);
        drone_state_mutex_.unlock();
        TimeConsuming t_back_frontend("t_back_frontend", false);
        double total_dur = ref_exp_traj.getTotalDuration();
        double start_t = ros_ptr_->getSimTime() - ref_exp_traj.getStartWallTime();


        if (start_t > total_dur - 0.01) {
            if (cfg_.print_log) {
                ros_ptr_->info(" -- [SUPER] in [generateBackupTrajectory]: start_t > total_dur, return NO_NEED");
            }
            return NO_NEED;
        }

        Vec3f temp_point;
        double out_t;
        bool all_traj_visible{true};
        // 同时记录每一个点的刹车时间和刹车距离
        vector<double> min_stop_dis;
        vector<TimePosPair> eval_ps;
        Vec3f temp_vel;

        // 记录当前时刻到最远时刻的所有可视部分
        Vec3f last_pos = ref_exp_traj.getPos(start_t);
        for (out_t = start_t; out_t < total_dur; out_t += cfg_.sample_traj_dt) {
            temp_point = ref_exp_traj.getPos(out_t);
            if ((last_pos - temp_point).norm() < cfg_.resolution * 0.8) {
                continue;
            }
            last_pos = temp_point;
            temp_vel = ref_exp_traj.getVel(out_t);
            // Compute initial
            double v_norm = temp_vel.norm();
            min_stop_dis.push_back(v_norm * v_norm / 2.0 / cfg_.exp_traj_cfg.max_acc);
            eval_ps.push_back(std::pair<double, Vec3f>(out_t, temp_point));
            const double min_dis =
                    cfg_.sensing_horizon > 0 ? std::min(cfg_.sensing_horizon, cfg_.safe_corridor_line_max_length)
                                             : cfg_.safe_corridor_line_max_length;
            if (!map_ptr_->isLineFree(back_traj_info.getRobotPos(),
                                      temp_point,
                                      min_dis,
                                      cfg_.seed_line_neighbour)) {
                all_traj_visible = false;
                break;
            }
        }

        if (all_traj_visible) {
            back_traj_info.setEmpty();
            {
                double dur = ref_exp_traj.getTotalDuration();
                Vec3f seed_pt = ref_exp_traj.getPos(dur);
                Line line{back_traj_info.getRobotPos(), seed_pt};
                Polytope temp_poly;
                if (cg_ptr_->GeneratePolytopeFromLine(line, temp_poly)) {
                    back_traj_info.setSFC(temp_poly);
                    {
                        TimeConsuming t_viz("tviz", false);
                        ros_ptr_->vizBackupSfc(temp_poly);
                        time_consuming_[VISUALIZATION] += t_viz.stop();
                    }
                }
            }
            return FINISH;
        }
        Vec3f invisible_p = eval_ps.back().second;
        while (out_t > start_t) {
            out_t -= cfg_.sample_traj_dt;
            Vec3f out_p = ref_exp_traj.getPos(out_t);
            if ((out_p - invisible_p).norm() > cfg_.robot_r) {
                break;
            }
        }

        double seed_point_t = std::max(start_t, out_t);

        // TODO check this logic, comment on Dec. 13
        // if
        // 1) last exp traj has a backup traj
        // 2) last backup WT is larger than this term
        // 3) last exp is collision free
        // if (ref_exp_traj.back_traj_start_TT > 0 &&
        // seed_point_t < ref_exp_traj.back_traj_start_TT) {
        // return NO_NEED;
        // }


        Vec3f seed_point = ref_exp_traj.getPos(seed_point_t);

        Vec3f shifted_robot_p = shifted_sfc_start_pt_.norm()> 999?robot_state_.p:shifted_sfc_start_pt_;
        if (!map_ptr_->getNearestCellNot(GridType::OCCUPIED, shifted_robot_p, shifted_robot_p, 3.0)) {
            ros_ptr_->error(
                    " -- [SUPER] in [PlanFromRest] Local start point is deeply occupied, which should not happened.");
            latest_replan.setRetCode(SUPER_RET_CODE::SUPER_NO_START_POINT);
            return FAILED;
        }

        Line line{shifted_robot_p, seed_point};
        Polytope temp_poly;
        if (!cg_ptr_->GeneratePolytopeFromLine(line, temp_poly)) {
            ros_ptr_->warn(" -- [SUPER] GeneratePolytopeFromLine failed, force return");
            return FAILED;
        }
        Eigen::Vector3d inner;
        Eigen::Matrix3Xd vPoly;
        if (!geometry_utils::findInterior(temp_poly.GetPlanes(), inner)) {
            ros_ptr_->warn(" -- [SUPER] Cannot generate feasible backup sfc, force return");
            vec_Vec3f seed{back_traj_info.getRobotPos(), seed_point};
            return FAILED;
        }

        if (cfg_.use_fov_cut) {
            if (!fov_checker_->cutPolyByFov(robot_state_.p, robot_state_.q, seed_point,
                                            temp_poly)) {
                ros_ptr_->warn(" -- [SUPER] cutPolyByFov failed, force return");
                return FAILED;
            }
        }
        // cut by sensing horizon
        if (cfg_.sensing_horizon > 0 &&
            !fov_checker_->cutPolyBySensingHorizon(robot_state_.p, seed_point, cfg_.sensing_horizon,
                                                   temp_poly)) {
            ros_ptr_->warn(" -- [SUPER] cutPolyBySensingHorizon failed, force return");
            vec_Vec3f seed{back_traj_info.getRobotPos(), seed_point};
            return FAILED;
        }

        back_traj_info.setSFC(temp_poly);

        {
            TimeConsuming t_viz("tviz", false);
            ros_ptr_->vizBackupSfc(temp_poly);
            time_consuming_[VISUALIZATION] += t_viz.stop();
        }

//        Vec3f out_p = temp_point;
//        double t_R = 0.0;
        double eval_t = eval_ps.back().first + cfg_.sample_traj_dt;
        last_pos = eval_ps.back().second;
        while (temp_poly.PointIsInside(eval_ps.back().second) && eval_t < total_dur) {
            Vec3f cur_pos = ref_exp_traj.getPos(eval_t);

            if ((cur_pos - last_pos).norm() < cfg_.resolution * 0.8) {
                eval_t += cfg_.sample_traj_dt;
                continue;
            }
            temp_vel = ref_exp_traj.getVel(out_t);
            double v_norm = temp_vel.norm();
            min_stop_dis.push_back(v_norm * v_norm / 2.0 / cfg_.exp_traj_cfg.max_acc);
            eval_ps.emplace_back(eval_t, cur_pos);
            last_pos = cur_pos;
            eval_t += cfg_.sample_traj_dt;
        }
        eval_ps.pop_back();
        seed_point = eval_ps.back().second;
        seed_point_t = eval_ps.back().first;

        //        bool use_new{true};
        //        if (use_new) {
        double t0 = ros_ptr_->getSimTime() -
                    ref_exp_traj.getStartWallTime() + 0.01;
        double te = seed_point_t;
        //            cout << "t0: " << t0 << endl;
        //            cout << "te: " << te << endl;
        //            cout << "exp_traj_dur: " << ref_exp_traj.optimized_exp_traj.getTotalDuration() << endl;
        double vel_e_n = ref_exp_traj.getVel(te).norm();
        double heu_ts = std::max((t0 + te) / 2, te - vel_e_n / cfg_.back_traj_cfg.max_acc);
        double heu_dur = te - heu_ts;
        Vec3f heu_p = seed_point;
        time_consuming_[BACK_TRAJ_FRONTEND] = t_back_frontend.stop();
        TimeConsuming t_back_opt("t_back_opt", false);
        double opt_ts = heu_ts;
        Trajectory temp_pos_traj;
        auto sfc0 = back_traj_info.getSFC();
        bool temp_ret = back_traj_opt_->optimize(ref_exp_traj.posTraj(),
                                                 t0,
                                                 te,
                                                 heu_ts,
                                                 heu_p,
                                                 heu_dur,
                                                 back_traj_info.getSFC(),
                                                 temp_pos_traj,
                                                 opt_ts);
        time_consuming_[BACK_TRAJ_OPT] = t_back_opt.stop();

        {
            double init_ts;
            VecDf init_times;
            vec_Vec3f init_ps;
            back_traj_opt_->getInitValue(init_ts, init_times, init_ps);
            latest_replan.setBackupCondition(init_ts, init_times, init_ps,
                                             t0, te,
                                             back_traj_info.getSFC());
            Trajectory traj;
            double out_ts;
            back_traj_opt_->optimize(ref_exp_traj.posTraj(),
                                     t0,
                                     te,
                                     init_ts,
                                     sfc0,
                                     init_times,
                                     init_ps,
                                     traj,
                                     out_ts
            );

        }

        if (!temp_ret) {
            ros_ptr_->warn(" -- [SUPER] OptimizationBakTrajInPolytopes failed, force return");
            back_traj_info.setEmpty();
            return OPT_FAILED;
        } else {
            Vec4f yaw_init_vec = ref_exp_traj.getYawState(opt_ts).row(0);
            Vec4f yaw_goal{0, 0, 0, 0};
            bool free_end{true};
            if (cfg_.goal_yaw_en) {
                if (!std::isnan(gi_.goal_yaw)) {
                    free_end = false;
                    yaw_goal[0] = gi_.goal_yaw;
                }
            }
            Trajectory temp_yaw_traj;
            if (!yaw_traj_opt_->optimize(yaw_init_vec, yaw_goal, temp_pos_traj,
                                         temp_yaw_traj, 3, false, free_end)) {
                ros_ptr_->error(" -- [SUPER] in [generateBackupTrajectory] YawTrajOpt FAILD.");
                return OPT_FAILED;
            }


            if (opt_ts < t0) {
                ros_ptr_->error(" -- [SUPER] opt_ts {} < t0 {}", opt_ts, t0);
                return OPT_FAILED;
            }
            double new_ts_WT = ref_exp_traj.getStartWallTime() + opt_ts;
            const auto &committed_ts_WT = cmd_traj_info_.getBackupTrajStartTT();
            if (committed_ts_WT < cmd_traj_info_.getTotalDuration() && new_ts_WT < committed_ts_WT) {
                ros_ptr_->error(" -- [SUPER] new_ts_WT {} < committed_ts_WT {}", new_ts_WT, committed_ts_WT);
                return OPT_FAILED;
            }


            {
                TimeConsuming t_viz("tviz", false);
                ros_ptr_->vizBackupTraj(temp_pos_traj);
                time_consuming_[VISUALIZATION] += t_viz.stop();
            }

            back_traj_info.setTrajectory(new_ts_WT, opt_ts, temp_pos_traj, temp_yaw_traj);
            latest_replan.setBackupTraj(temp_pos_traj);
            latest_replan.setBackupYawTraj(temp_yaw_traj);
            return SUCCESS;
        }
        ros_ptr_->warn(" -- [SUPER] Cannot find backup traj start point.");
        return FAILED;
    }

    int SuperPlanner::getNearestFurtherGoalPoint(const vec_E<Vec3f> &goals, const Vec3f &start_pt) {
        if (goals.size() == 1) {
            return 0;
        }
        Vec3f a = start_pt, b;
        int min_id = 0;
        double min_dis = 1e10;
        for (long unsigned int i = 0; i < goals.size() - 1; i++) {
            b = goals[i];
            double dis = geometry_utils::pointLineSegmentDistance(start_pt, a, b);
            if (dis < min_dis) {
                min_dis = dis;
                min_id = i;
            }
            a = b;
        }
        return min_id;
    }

    bool
    SuperPlanner::PathSearch(const Vec3f &start_pt, const Vec3f &goal,
                             const double &searching_horizon,
                             vec_Vec3f &path) {
        using namespace path_search;
        if (searching_horizon <= 0.0) {
            ros_ptr_->error(" -- [SUPER] Goal waypoints empty or searching horizon negative, force return.");
            return false;
        }

        // 1) check and shift pts
        // 		For start point, must be collision free
        rog_map::GridType start_type;
        start_type = map_ptr_->getGridType(start_pt);

        /// If the start_pt is obstacle in prob map, just shift it to the nearest free point.
        if (start_type == rog_map::GridType::OCCUPIED ||
            start_type == rog_map::GridType::OUT_OF_MAP) {
            ros_ptr_->warn(
                    " -- [SUPER] The start point in obstacle, this should not happen since the start point should be shift before pathsearch.");
            return false;
        }
        vec_E<Vec3f> start_point_escape_path;

        int flag_es = ON_PROB_MAP | (cfg_.frontend_in_known_free ? UNKNOWN_AS_OCCUPIED : UNKNOWN_AS_FREE);
        vec_Vec3f out_path;
        RET_CODE ret_es = astar_ptr_->escapePathSearch(start_pt, flag_es, out_path);
        if (ret_es != NO_NEED) {
            if (ret_es != REACH_HORIZON && ret_es != REACH_GOAL) {
                ros_ptr_->error(
                        " -- [SUPER] Escape path search failed with [{}], force return.",
                        RET_CODE_STR[ret_es].c_str());
                return false;
            } else {
                start_point_escape_path = out_path;
            }
        }

        Vec3f shifted_start_pt = start_pt;

        if (!start_point_escape_path.empty()) {
            shifted_start_pt = start_point_escape_path.back();
        }

        Vec3f temp_goal_point, temp_start_point;
        temp_start_point = shifted_start_pt;
        double temp_plannning_horizon = searching_horizon;
        //            int start_id = getNearestFurtherGoalPoint(goal_waypoints, start_pt);

        int flag = ON_INF_MAP | (cfg_.frontend_in_known_free ? UNKNOWN_AS_OCCUPIED : UNKNOWN_AS_FREE) | DONT_USE_INF_NEIGHBOR;

        RET_CODE ret_code = astar_ptr_->pointToPointPathSearch(temp_start_point, goal, flag, temp_plannning_horizon,
                                                               path);

        if(ret_code == INIT_ERROR){
            gi_.goal_valid = false;
            return false;
        }
        //add may23, if failed on inf map, use prob map try again

        if (ret_code == NO_PATH) {
            flag = ON_PROB_MAP | (cfg_.frontend_in_known_free ? UNKNOWN_AS_OCCUPIED : UNKNOWN_AS_FREE) |
                   USE_INF_NEIGHBOR;
            fmt::print(fg(fmt::color::indian_red) | fmt::emphasis::bold,
                       " -- [Astar] Path search failed on inf map, try again on prob map.\n");
            ret_code = astar_ptr_->pointToPointPathSearch(temp_start_point, goal, flag, temp_plannning_horizon,
                                                          path);
            if (ret_code == SUCCESS || ret_code == REACH_HORIZON || ret_code == REACH_GOAL) {
                fmt::print(fg(fmt::color::lime_green) | fmt::emphasis::bold,
                           " -- [Astar] Path search on prob map success.\n");
            } else {
                fmt::print(fg(fmt::color::indian_red) | fmt::emphasis::bold,
                           " -- [Astar] Path search failed on prob map still failed.\n");
            }
        }
        if (ret_code != REACH_HORIZON && ret_code != REACH_GOAL) {
            ros_ptr_->error(
                    " -- [SUPER] Path search failed with [{}], force return.\n", RET_CODE_STR[ret_code].c_str());
            return false;
        }
        if (!start_point_escape_path.empty()) {
            path.insert(path.begin(), start_point_escape_path.begin(),
                        start_point_escape_path.end());
        }

        if (path.empty()) {
            ros_ptr_->warn(
                    " -- [SUPER] Path search failed with empty segments, force return.");
            return false;
        }
        path.insert(path.begin(), start_pt);
        if (ret_code == REACH_GOAL) {
            path.push_back(goal);
        }
        return true;
    }


    void SuperPlanner::getRobotState(rog_map::RobotState &out) {
        robot_state_ = map_ptr_->getRobotState();
        out = robot_state_;
    }
}
