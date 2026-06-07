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

#include <super_core/corridor_generator.h>

using namespace super_utils;

namespace super_planner {

    CorridorGenerator::CorridorGenerator(const ros_interface::RosInterface::Ptr &ros_ptr,
                                         const rog_map::ROGMapROS::Ptr &map_ptr, const double bound_dis,
                                         const double seed_line_max_dis, const double min_overlap_threshold,
                                         const double virtual_groud_height, const double virtual_ceil_height,
                                         const double robot_r, const int box_search_skip_num, const int iris_iter_num)
            : ros_ptr_(ros_ptr), map_ptr_(map_ptr) {
        ciri_ = std::make_shared<CIRI>(ros_ptr_);
        ciri_->setupParams(robot_r, iris_iter_num);
        bound_dis_ = bound_dis;
        seed_line_max_length_ = seed_line_max_dis;
        min_overlap_threshold_ = min_overlap_threshold;
        robot_r_ = robot_r;
        box_search_skip_num_ = box_search_skip_num;
        iris_iter_num_ = iris_iter_num;
        virtual_ceil_height_ = virtual_ceil_height - robot_r;
        virtual_groud_height_ = virtual_groud_height + robot_r;
//        failed_traj_log.open(DEBUG_FILE_DIR("sfc.csv"), std::ios::out | std::ios::trunc);
    }


    void CorridorGenerator::SetLineNeighborList(const vec_E<Vec3i> &_line_seed_neighbor_list) {
        this->line_seed_neighbor_list = _line_seed_neighbor_list;
    }

    // 沿一条离散路径逐段生成连续的安全走廊（Safe Flight Corridor, SFC）。
    //
    // 整体功能：
    // 1) 以 path 上的点为参考，尽量把多段连续路径点合并成一条较长的 seed line；
    // 2) 对每条 seed line 调用 CIRI/凸分解流程生成一个 polytope；
    // 3) 检查相邻 polytope 是否有足够重叠，必要时在连接点额外补一个“点走廊”；
    // 4) 最终输出一串首尾连续、可供后端轨迹优化使用的多面体走廊 sfcs。
    //
    // 输入：
    // - path: 离散引导路径，通常来自 A* 搜索或旧轨迹复用后的 guide_path。
    // - sfcs: 输出参数，保存生成的走廊序列；函数开始时会先清空。
    // - shifted_start_pt: 输出参数。当 path 前几个点已经落入膨胀障碍时，
    //   该变量会被写成第一个真正可用于生成走廊的路径点。
    // - cut_first_poly: 预留开关，当前实现里未实际使用。
    //
    // 输出：
    // - 返回 true 表示成功生成了一串非空且连续的 polytope；
    // - 返回 false 表示路径为空、某段凸分解失败、或相邻走廊无法形成连续重叠。
    bool
    CorridorGenerator::SearchPolytopeOnPath(const vec_Vec3f &path, PolytopeVec &sfcs,
                                            Vec3f &shifted_start_pt,
                                            bool cut_first_poly) {
        // https://whimsical.com/flow-3TASJFwe1dASYYY2xHEmze
        // password: wtr
        //	TimeConsuming t___("SearchPolytopeOnPath");
        sfcs.clear();
        if (path.empty()) {
            return false;
        }

        // seed_lines: 按路径分段后得到的“种子线段”，每条线段对应一个候选 polytope。
        vector<Line> seed_lines;
        // first_id / second_id: 当前正在处理的路径段在 path 中的起止索引。
        int first_id, second_id;
        // overlap: 当前新 polytope 与上一个 polytope 的交集；
        // interior_pt / interior_depth: 用于衡量交集是否“足够厚”，从而保证走廊连续可穿行。
        Polytope overlap;
        Vec3f interior_pt;
        double interior_depth;
        // temp_poly: 当前 seed line 生成的主 polytope；
        // temp_poly_fix_p: 当相邻两个主 polytope 重叠不足时，在连接点额外补的“修复 polytope”。
        Polytope temp_poly, temp_poly_fix_p;
        // max_loop / cnt_loop: 防御性循环计数，避免异常情况下死循环。
        int max_loop = 1000;
        int cnt_loop = 0;
        first_id = 0;

        // 若路径起点附近已经落入膨胀障碍，则先跳过这些点，
        // 后续从第一个安全点开始正式生成 corridor。
        while(first_id < path.size() && map_ptr_->isOccupiedInflate(path[first_id])) {
            first_id++;
        }

        // 同时记录 shifted_start_pt，供上层知道真正的 corridor 起点被沿路径向后移动了。
        if(first_id!=0){
            shifted_start_pt = path[first_id];
            double dis = (path[first_id] - path[0]).norm() * 1.2;
            // 如果前面确实跳过了一段被膨胀障碍覆盖的路径，
            // 则用一个“空走廊 polytope”把原起点到第一个安全点先包起来，
            GenerateEmptyPolytope(path[0], dis, temp_poly);
            sfcs.emplace_back(temp_poly);
        }

        while (cnt_loop++ < max_loop) {
            second_id = first_id;
            // 从 path[first_id] 开始，尽可能往后扩展一条更长的 seed line。
            // 只要 first_id -> j 这整段线仍然是 free 的，就继续向后吃更多点；
            // 一旦不再 free，就在上一个可行点附近截断。
            for (int j = first_id + 1; j < path.size(); j++) {
                bool reach_segment = false;
                // 直观理解：
                // 这里在问“能否把 path[first_id] 和 path[j] 之间的整段路径，
                // 视为一条可以直接拿来生成 corridor 的 seed line”。
                // 即把一条连续线段展开成它穿过的所有体素格子，逐格检查。
                //
                // 这个判断同时包含两层约束：
                // 1) 这条候选线段不能太长，长度上限由 seed_line_max_length_ 控制；
                // 2) 这条线段附近要留有足够余量，line_seed_neighbor_list 可以理解为
                //    围绕线段做的一圈离散“加粗/膨胀”检查，避免只是一条无厚度中心线可通。
                //
                // 一旦这个判断失败，就说明再继续把 seed line 往 path[j] 延长已经不合适，
                // 当前段应该在更早的位置截断。
                if (!map_ptr_->isLineFree(path[first_id], path[j], seed_line_max_length_,
                                          line_seed_neighbor_list)) {
                    reach_segment = true;
                }
                if (reach_segment) {
                    second_id = j - 1;
                    // 再额外保守退一格，避免 seed line 贴着碰撞边界，
                    // 从而让后续凸分解更稳定一些。
                    if (second_id - 1 > first_id) {
                        second_id -= 1;
                    }
                    break;
                }
                // 当 segment 仍然能向后扩展时，更新 second_id
                second_id = j;
            }

            // 若 first_id 这一点甚至无法和更远的点组成可行线段，
            // 至少强行让 seed line 覆盖到下一个路径点，避免停滞。
            if (second_id == first_id && second_id + 1 < path.size()) {
                second_id += 1;
            }

            seed_lines.emplace_back(path[first_id], path[second_id]);
            // 正常情况下 seed line 长度应被 isLineFree(..., seed_line_max_length_, ...)
            // 这一约束控制住；这里再做一次防御性检查。
            if ((path[first_id] - path[second_id]).norm() > seed_line_max_length_ * 1.5) {
                fmt::print("first: {}\n second: {}\n seed line max: {}\n", path[first_id].transpose(),
                           path[second_id].transpose(), seed_line_max_length_);
                throw std::runtime_error("seed line too long");
                return false;
            }
            // 对当前处理的 seed line 生成一个主 polytope。
            if (!GeneratePolytopeFromLine(seed_lines.back(), temp_poly)) {
                cout << YELLOW << " -- [SUPER] GeneratePolytopeFromLine failed." << RESET << endl;
                return false;
            }

// viz for debug
//            ros_ptr_->vizCiriPolytope(temp_poly, "debug");
//            usleep(10000);

            // 这一整段逻辑的目标是“维护 corridor 序列的连接质量”：
            // 当前已经新生成了一个主 polytope(temp_poly)，但它不能只单独可行，
            // 还必须与前面已经接受到 sfcs 中的走廊序列连续衔接。
            //
            // 因此这里会分两类情况处理：
            // 1) 如果 temp_poly 与上一段走廊的交集太薄，说明两段主走廊直接连接不可靠，
            //    就在连接点额外插入一个点 polytope 作为“过渡走廊”；
            // 2) 如果 temp_poly 与上一段已经连接良好，则进一步检查中间某段是否冗余，
            //    必要时删除冗余 corridor，以避免走廊过密。
            //
            // 简单说：这一段既负责“补连接”，也负责“去冗余”。
            if (!sfcs.empty()) {
                // 计算当前 polytope 与上一个已接受 polytope 的交集深度。
                // 深度越大，说明两段 corridor 的衔接越“厚”，连续性越好。
                overlap = sfcs.back().CrossWith(temp_poly);
                interior_depth = geometry_utils::findInteriorDist(overlap.GetPlanes(), interior_pt);
                temp_poly.overlap_depth_with_last_one = interior_depth;
                temp_poly.interior_pt_with_last_one = interior_pt;
                if (interior_depth < min_overlap_threshold_) {
                    // 若相邻两段主走廊重叠过薄，则在连接点 path[first_id]
                    // 额外生成一个点 polytope 作为“过渡走廊”。
                    if (!GeneratePolytopeFromPoint(path[first_id], temp_poly_fix_p)) {
                        cout << YELLOW << " -- [SUPER] GeneratePolytopeFromPoint failed." << RESET << endl;
                        return false;
                    }
                    // 先检查“上一段主走廊”和“点走廊”的重叠是否足够。
                    overlap = sfcs.back().CrossWith(temp_poly_fix_p);
                    interior_depth = geometry_utils::findInteriorDist(overlap.GetPlanes(), interior_pt);
                    if (interior_depth <= 0.01) {
                        ros_ptr_->warn(
                                " -- [SUPER] Cannot find continuous corridor on path, overlap only {}, force return.",
                                interior_depth);
// viz for debug
//                        ros_ptr_->vizCiriPointCloud(latest_pc);
//                        usleep(100000);
//                        exit(-1);
                        return false;
                    }
                    temp_poly_fix_p.overlap_depth_with_last_one = interior_depth;
                    temp_poly_fix_p.interior_pt_with_last_one = interior_pt;
                    sfcs.push_back(temp_poly_fix_p);
                    // 再检查“点走廊”和“当前主走廊”的重叠是否足够。
                    overlap = sfcs.back().CrossWith(temp_poly);
                    interior_depth = geometry_utils::findInteriorDist(overlap.GetPlanes(), interior_pt);
                    if (interior_depth <= 0.01) {
                        ros_ptr_->warn(
                                " -- [SUPER] Cannot find continuous corridor on path, overlap only {}, force return.",
                                interior_depth);
                        // viz for debug
//                        ros_ptr_->vizCiriPointCloud(latest_pc);
//                        usleep(100000);
//                        exit(-1);
                        return false;
                    }
                } else {
                    // 如果当前 polytope 与“上一个上一个 polytope”的重叠反而更大，
                    // 说明中间那个 polytope 可能是冗余的，可以将其删除以降低走廊密度。
                    int temp_id = sfcs.size() - 2;
                    if (temp_id > 0) {
                        overlap = sfcs[temp_id].CrossWith(temp_poly);
                        interior_depth = geometry_utils::findInteriorDist(overlap.GetPlanes(), interior_pt);
                        if (interior_depth > sfcs[temp_id + 1].overlap_depth_with_last_one * 0.25) {
                            temp_poly.overlap_depth_with_last_one = interior_depth;
                            temp_poly.interior_pt_with_last_one = interior_pt;
                            sfcs.pop_back();
                        }
                    }
                }
            }

            // 接受当前主 polytope，作为 corridor 的下一段。
            sfcs.push_back(temp_poly);
            if (second_id == path.size() - 1) {
                // 已覆盖到 path 的最后一个点，走廊生成结束。
                break;
            }
            // 从当前段末端继续向后处理下一段路径。
            first_id = second_id;
        }
        // Delete last polytope if the second last one contains the last point


        if (cnt_loop >= max_loop) {
            cout << YELLOW << " -- [SUPER] Reach max iteration, failed." << RESET << endl;
            return false;
        }

        if (sfcs.empty()) {
            return false;
        }

        // 能走到这里，说明已经成功构造出一串非空 corridor。
        return true;
    }


    // 根据一个 seed point / seed line 的两个端点，构造其对应的局部搜索包围盒。
    //
    // 功能：
    // - 先取 p1、p2 的逐维最小/最大值，得到覆盖该点或线段的轴对齐包围盒；
    // - 再在 xyz 三个方向各扩张 bound_dis_，为后续 boxSearch / 凸分解预留搜索空间。
    //
    // 何时调用：
    // - 在 GeneratePolytopeFromPoint() 中，用于围绕单个连接点建立局部搜索区域；
    // - 在 GeneratePolytopeFromLine() 中，用于围绕当前 seed line 建立局部搜索区域。
    void CorridorGenerator::getSeedBBox(const Vec3f &p1, const Vec3f &p2, Vec3f &box_min, Vec3f &box_max) {
        box_min = p1.cwiseMin(p2);
        box_max = p1.cwiseMax(p2);
        box_min -= Vec3f(bound_dis_, bound_dis_, bound_dis_);
        box_max += Vec3f(bound_dis_, bound_dis_, bound_dis_);
//        box_min.z() = std::max(box_min.z(), virtual_groud_height_);
//        box_max.z() = std::min(box_max.z(), virtual_ceil_height_);
    }

    // 围绕单个点生成一个 polytope。
    //
    // 功能：
    // - 以 pt 为 seed，在其附近先取一个局部包围盒；
    // - 收集该包围盒内的障碍物点；
    // - 若附近无障碍，则直接返回一个轴对齐盒状 polytope；
    // - 若附近有障碍，则调用 CIRI/凸分解，在障碍约束下生成围绕该点的可行 polytope。
    //
    // 何时调用：
    // - 在 SearchPolytopeOnPath() 中，当相邻两个主 corridor 的重叠不够时，
    //   会在连接点处插入一个“点走廊”作为过渡 corridor。
    bool CorridorGenerator::GeneratePolytopeFromPoint(const Vec3f &pt, Polytope &polytope) {
        Eigen::Vector3d box_max, box_min;
        vec_E<Vec3f> pc;
        getSeedBBox(pt, pt, box_min, box_max);
        // TODO the box did not consider the robot_r
        map_ptr_->boundBoxByLocalMap(box_min, box_max);
        map_ptr_->boxSearch(box_min, box_max, OCCUPIED, pc);
        box_min.z() += robot_r_;
        box_max.z() -= robot_r_;
        MatD4f planes;
        Eigen::Vector3d a = pt, b = pt;
        Eigen::Matrix<double, 6, 4> bd = Eigen::Matrix<double, 6, 4>::Zero();
        bd(0, 0) = 1.0;
        bd(1, 0) = -1.0;
        bd(2, 1) = 1.0;
        bd(3, 1) = -1.0;
        bd(4, 2) = 1.0;
        bd(5, 2) = -1.0;
        bd(0, 3) = -box_max.x();
        bd(1, 3) = box_min.x();
        bd(2, 3) = -box_max.y();
        bd(3, 3) = box_min.y();
        bd(4, 3) = -box_max.z();
        bd(5, 3) = box_min.z();
        // 将vector放到Eigen里，准备开始分解
        if (pc.empty()) {
            // 障碍物点云为空，直接返回一个方块
            // Ax + By + Cz + D = 0
            planes.resize(6, 4);
            planes.row(0) << 1, 0, 0, -box_max.x();
            planes.row(1) << 0, 1, 0, -box_max.y();
            planes.row(2) << 0, 0, 1, -box_max.z();
            planes.row(3) << -1, 0, 0, box_min.x();
            planes.row(4) << 0, -1, 0, box_min.y();
            planes.row(5) << 0, 0, -1, box_min.z();
            polytope.SetPlanes(planes);
            polytope.SetSeedLine(Line{pt, pt});
            return true;
        }
        latest_pc.insert(latest_pc.end(), pc.begin(), pc.end());
        Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>> pp(pc[0].data(), 3, pc.size());
        rog_map::TimeConsuming tc("emvp", false);
        RET_CODE success = ciri_->comvexDecomposition(bd, pp, a, b);
        double dt = tc.stop();
        if (success == SUCCESS) {
            ciri_cnt++;
            ciri_t += dt;
            ciri_->getPolytope(polytope);
            polytope.SetSeedLine(Line{pt, pt});
            return true;
        } else {
            cout << YELLOW << " -- [SUPER] CSpaceFiri failed." << RESET << endl;
            cout << YELLOW << "\t box_min =" << box_min.transpose() << endl;
            cout << YELLOW << "\t box_max = " << box_max.transpose() << endl;
            cout << YELLOW << "\t seed pt = " << pt.transpose() << endl;
            polytope.Reset();
            return false;
        }

    }

    // 直接生成一个不考虑障碍物的“空 polytope”，本质上就是以 pt 为中心的轴对齐盒。
    //
    // 功能：
    // - 不做 boxSearch，也不做凸分解；
    // - 仅按给定半径 dis 构造一个立方体/长方体 polytope。
    //
    // 何时调用：
    // - 在 SearchPolytopeOnPath() 中，如果路径起始一小段已经落入膨胀障碍，
    //   会先补一个这种过渡用的空 corridor，把原始起点和第一个可用走廊段连接起来。
    bool CorridorGenerator::GenerateEmptyPolytope(const super_utils::Vec3f &pt,
                                                  const double & dis,
                                                  Polytope & polytope){
        Eigen::Vector3d box_max, box_min;
        box_min = pt;
        box_max = pt;
        box_min -= Vec3f(dis, dis, dis);
        box_max += Vec3f(dis, dis, dis);
        MatD4f planes;
        Eigen::Matrix<double, 6, 4> bd = Eigen::Matrix<double, 6, 4>::Zero();
        bd(0, 0) = 1.0;
        bd(1, 0) = -1.0;
        bd(2, 1) = 1.0;
        bd(3, 1) = -1.0;
        bd(4, 2) = 1.0;
        bd(5, 2) = -1.0;
        bd(0, 3) = -box_max.x();
        bd(1, 3) = box_min.x();
        bd(2, 3) = -box_max.y();
        bd(3, 3) = box_min.y();
        bd(4, 3) = -box_max.z();
        bd(5, 3) = box_min.z();
        // 障碍物点云为空，直接返回一个方块
        // Ax + By + Cz + D = 0
        planes.resize(6, 4);
        planes.row(0) << 1, 0, 0, -box_max.x();
        planes.row(1) << 0, 1, 0, -box_max.y();
        planes.row(2) << 0, 0, 1, -box_max.z();
        planes.row(3) << -1, 0, 0, box_min.x();
        planes.row(4) << 0, -1, 0, box_min.y();
        planes.row(5) << 0, 0, -1, box_min.z();
        polytope.SetPlanes(planes);
        polytope.SetSeedLine(Line{pt, pt});
        return true;
    }

    // 围绕一条 seed line 生成一个 polytope，这是 corridor 生成流程中的主力函数。
    //
    // 功能：
    // - 以 line.first -> line.second 这条 seed line 为核心，在其周围建立局部包围盒；
    // - 收集该局部区域内的障碍物点；
    // - 若无障碍，则直接返回一个覆盖 seed line 的轴对齐盒状 polytope；
    // - 若有障碍，则调用 CIRI/凸分解，在障碍约束下生成包裹这条 seed line 的可行 corridor。
    //
    // 何时调用：
    // - 在 SearchPolytopeOnPath() 中，每当从路径上挑出一条可行的 seed line 后，
    //   就会调用本函数生成当前这一个“主 corridor 段”。
    bool CorridorGenerator::GeneratePolytopeFromLine(Line &line, Polytope &polytope) {
        Eigen::Vector3d box_max, box_min;
        vec_E<Vec3f> pc, pts{line.first, line.second};
        getSeedBBox(line.first, line.second, box_min, box_max);
        map_ptr_->boundBoxByLocalMap(box_min, box_max);
        map_ptr_->boxSearch(box_min, box_max, OCCUPIED, pc);
        box_min.z() += robot_r_;
        box_max.z() -= robot_r_;
        MatD4f planes;
        Eigen::Vector3d a = line.first, b = line.second;
        Eigen::Matrix<double, 6, 4> bd = Eigen::Matrix<double, 6, 4>::Zero();
        bd(0, 0) = 1.0;
        bd(1, 0) = -1.0;
        bd(2, 1) = 1.0;
        bd(3, 1) = -1.0;
        bd(4, 2) = 1.0;
        bd(5, 2) = -1.0;
        bd(0, 3) = -box_max.x();
        bd(1, 3) = box_min.x();
        bd(2, 3) = -box_max.y();
        bd(3, 3) = box_min.y();
        bd(4, 3) = -box_max.z();
        bd(5, 3) = box_min.z();
        // 将vector放到Eigen里，准备开始分解
        if (pc.empty()) {
            // 障碍物点云为空，直接返回一个方块
            // Ax + By + Cz + D = 0g
            planes.resize(6, 4);
            planes.row(0) << 1, 0, 0, -box_max.x();
            planes.row(1) << 0, 1, 0, -box_max.y();
            planes.row(2) << 0, 0, 1, -box_max.z();
            planes.row(3) << -1, 0, 0, box_min.x();
            planes.row(4) << 0, -1, 0, box_min.y();
            planes.row(5) << 0, 0, -1, box_min.z();
            polytope.SetPlanes(planes);
            polytope.SetSeedLine(line);
            return true;
        }
        // save to latest pc
        latest_pc.insert(latest_pc.end(), pc.begin(), pc.end());
        Eigen::Map<const Eigen::Matrix<double, 3, -1, Eigen::ColMajor>> pp(pc[0].data(), 3, pc.size());
        rog_map::TimeConsuming tc("emvp", false);
        RET_CODE success = ciri_->comvexDecomposition(bd, pp, a, b);
        double dt = tc.stop();
        if (success == SUCCESS) {
            ciri_cnt++;
            ciri_t += dt;
            ciri_->getPolytope(polytope);
            polytope.SetSeedLine(line);
            return true;
        } else {
            polytope.Reset();
            cout << YELLOW << "\t box_min = " << box_min.transpose() << RESET << endl;
            cout << YELLOW << "\t box_max =" << box_max.transpose() << RESET << endl;
            cout << YELLOW << "\t seed line =" << line.first.transpose() << " --> " << line.second.transpose()
                 << RESET << endl;

//            failed_traj_log << 889900 << endl;
//            failed_traj_log << bd << endl;
//            failed_traj_log << 0 << endl;
//            failed_traj_log << pp << endl;
//            failed_traj_log << 0 << endl;
//            failed_traj_log << a.transpose() << endl;
//            failed_traj_log << 0 << endl;
//            failed_traj_log << b.transpose() << endl;
//            failed_traj_log << 0 << endl;
//            failed_traj_log << robot_r_ << endl;
//            failed_traj_log << 0 << endl;
//            failed_traj_log << iris_iter_num_ << endl;
            return false;
        }

    }
}
