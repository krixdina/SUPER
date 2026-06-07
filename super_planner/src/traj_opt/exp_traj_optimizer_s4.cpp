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

#include <traj_opt/exp_traj_optimizer_s4.h>
#include <utils/optimization/lbfgs.h>
#include <ros_interface/ros_interface.hpp>

#define POS_IDX 1
#define VEL_IDX 2
#define ACC_IDX 3
#define JER_IDX 4
#define ATT_IDX 5
#define OMG_IDX 6
#define THR_IDX 7

using namespace traj_opt;
using namespace color_text;
using namespace super_utils;
using namespace math_utils;
using namespace optimization_utils;

using Vec8f = Eigen::Matrix<double, 8, 1>;
using Mat83f = Eigen::Matrix<double, 8, 3>;

/*
 * 功能:
 *   在给定分段多项式轨迹后，对整条轨迹按时间采样并累计各类约束惩罚函数，同时计算这些惩罚函数
 *   对分段时间 T 和多项式系数 coeffs 的梯度，供外层 L-BFGS 目标函数继续向优化变量回传。
 *
 * 输入:
 *   T: 每段轨迹的持续时间。
 *   coeffs: MINCO 生成的 8 阶多项式系数，每段对应 8x3 系数块。
 *   hIdx: 每段轨迹对应使用哪个 H 表示走廊多面体的索引。
 *   hPolys: H 表示的安全走廊集合，用于位置约束惩罚。
 *   waypoint_attractor: 相邻走廊重叠区域提取出的吸引点，用于把段间连接点拉向可行中间区域。
 *   waypoint_attractor_dead_d: 吸引点的免罚半径，落在该距离内不产生吸引惩罚。
 *   smoothFactor: 平滑 L1 罚函数的平滑系数。
 *   integralResolution: 每段轨迹的积分离散分辨率。
 *   magnitudeBounds: 动力学量上下界，依次包含速度、加速度、jerk、角速度和推力范围。
 *   penaltyWeights: 各类约束罚项的权重，控制位置、速度、加速度、jerk、吸引点、角速度、推力约束的相对强度。
 *   flatMap: 微分平坦映射，用于在速度/加速度/jerk 与推力/姿态/角速度之间做正反向梯度传播。
 *
 * 输出:
 *   cost: 累加后的约束罚代价，会直接叠加到外层总代价中。
 *   gradT: 约束罚代价对每段时间 T 的偏导。
 *   gradC: 约束罚代价对多项式系数 coeffs 的偏导。
 *   pena_log: 记录各类约束在本次评估中的最大违反量，便于日志和调试观察。
 */
void ExpTrajOpt::constraintsFunctional(const VecDf &T,
                                       const MatD3f &coeffs,
                                       const VecDi &hIdx,
                                       const PolyhedraH &hPolys,
                                       const Mat3Df &waypoint_attractor,
                                       const VecDf &waypoint_attractor_dead_d,
                                       const double &smoothFactor,
                                       const int &integralResolution,
                                       const VecDf &magnitudeBounds,
                                       const VecDf &penaltyWeights,
                                       flatness::FlatnessMap &flatMap,
        // outputs
                                       double &cost,
                                       VecDf &gradT,
                                       MatD3f &gradC,
                                       VecDf &pena_log) {
    // 读取约束边界与权重，并预先整理出平方或区间中心等便于后续重复使用的量。
    /* 1) define some varible alias*/
    const auto &vmax = magnitudeBounds[0];
    const auto &amax = magnitudeBounds[1];
    const auto &jmax = magnitudeBounds[2];
    const auto &omgmax = magnitudeBounds[3];
    const auto &accthrmin = magnitudeBounds[4];
    const auto &accthrmax = magnitudeBounds[5];

    const auto &vmaxSqr = vmax * vmax;
    const auto &amaxSqr = amax * amax;
    const auto &jmaxSqr = jmax * jmax;
    const auto &omgmaxSqr = omgmax * omgmax;

    const auto &thrustMean = 0.5 * (accthrmax + accthrmin);
    const auto &thrustRadi = 0.5 * std::abs(accthrmax - accthrmin);
    const auto &thrustSqrRadi = thrustRadi * thrustRadi;

    const auto &weightPos = penaltyWeights[0];
    const auto &weightVel = penaltyWeights[1];
    const auto &weightAcc = penaltyWeights[2];
    const auto &weightJer = penaltyWeights[3];
    const auto &weightAtt = penaltyWeights[4];
    const auto &weightOmg = penaltyWeights[5];
    const auto &weightAccThr = penaltyWeights[6];

    const auto &piece_num = T.size();

    const double integralFrac = 1.0 / integralResolution;
    VecDf max_pena(8);
    max_pena.setZero();

    // 对每一段轨迹做数值积分。每个采样点先通过多项式系数恢复 pos/vel/acc/jer/snap，再累计所有启用的约束罚项。
    /* 2) add integral cost */

    for (int i = 0; i < piece_num; i++) {
        const Mat83f &c = coeffs.block<8, 3>(i * 8, 0);
        const auto &step = T(i) * integralFrac;
        for (int j = 0; j <= integralResolution; j++) {
            double s1 = j * step;
            double s2 = s1 * s1;
            double s3 = s2 * s1;
            double s4 = s2 * s2;
            double s5 = s4 * s1;
            double s6 = s4 * s2;
            double s7 = s4 * s3;
            Vec8f beta0, beta1, beta2, beta3, beta4;
            beta0 << 1.0, s1, s2, s3, s4, s5, s6, s7;
            beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4, 6.0 * s5, 7.0 * s6;
            beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3, 30.0 * s4, 42.0 * s5;
            beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2, 120.0 * s3, 210.0 * s4;
            beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * s1, 360.0 * s2, 840.0 * s3;
            //beta5 << 0.0, 0.0, 0.0, 0., 0.0, 120.0, 720.0 * s1, 2520.0 * s2;

            const Vec3f pos = c.transpose() * beta0;
            const Vec3f vel = c.transpose() * beta1;
            const Vec3f acc = c.transpose() * beta2;
            const Vec3f jer = c.transpose() * beta3;
            const Vec3f sna = c.transpose() * beta4;

            double tmp_cost{0.0};
            Vec3f gradPos{0, 0, 0}, gradVel{0, 0, 0}, gradAcc{0, 0, 0}, gradJer{0, 0, 0};

            // 位置约束: 检查采样点是否穿出当前段所属安全走廊的任一半空间边界。
            /* 2.1  For position cost */
            const auto &L = hIdx(i);
            const auto &K = hPolys[L].rows();
            if (weightPos > 0) {
                for (int k = 0; k < K; k++) {
                    const Vec3f outerNormal = hPolys[L].block<1, 3>(k, 0);
                    const double violaPos = outerNormal.dot(pos) + hPolys[L](k, 3);
                    if (violaPos > max_pena(POS_IDX)) max_pena(POS_IDX) = violaPos;
                    double violaPosPena, violaPosPenaD;
                    if (gcopter::smoothedL1(violaPos, smoothFactor, violaPosPena, violaPosPenaD)) {
                        gradPos += weightPos * violaPosPenaD * outerNormal;
                        tmp_cost += weightPos * violaPosPena;
                    }
                }
            }

            // 吸引点约束: 仅在段首/段尾节点处生效，鼓励连接点落入相邻走廊重叠区附近，减少切换走廊时的不可行风险。
            /* 2.2  For attract point cost  */
            if (weightAtt > 0.0) {
                // 只在内部连接点检查吸引点约束:
                // - 当前段的起点(j == 0, i != 0)对应前一对走廊的重叠区
                // - 当前段的终点(j == integralResolution, i != piece_num - 1)对应当前与下一段走廊的重叠区
                const auto is_waypoint = (j == 0) && (i != 0);
                const auto is_end = ((j == integralResolution) && (i != piece_num - 1));
                // waypoint_attractor 的第 k 列表示第 k 个相邻走廊重叠区域的吸引点，
                // 因此段尾使用 i，段首使用 i - 1 来索引对应的重叠区。
                const auto idx = is_end ? i : i - 1;

                if (is_waypoint || is_end) {
                    Vec3f p_a = pos - waypoint_attractor.col(idx);
                    // 若连接点落在吸引点 dead zone 半径内则不罚，超出后按平方距离超量进行平滑 L1 惩罚。
                    const auto &violaAtt =
                            p_a.squaredNorm() - waypoint_attractor_dead_d(idx) * waypoint_attractor_dead_d(idx);
                    double violaAttPena, violaAttPenaD;
                    if (violaAtt > max_pena(ATT_IDX)) max_pena(ATT_IDX) = violaAtt;
                    if (gcopter::smoothedL1(violaAtt, smoothFactor, violaAttPena, violaAttPenaD)) {
                        // 对位置的梯度方向就是把当前连接点往吸引点中心拉回去。
                        gradPos += weightAtt * violaAttPenaD * 2.0 * p_a;
                        tmp_cost += weightAtt * violaAttPena;
                    }
                }
            }

            // 运动学约束: 直接对速度、加速度、jerk 的平方范数超限部分施加平滑惩罚。
            /* 2.3 For vel cost  */
            const auto &violaVel = vel.squaredNorm() - vmaxSqr;
            double violaVelPena, violaVelPenaD;
            if (weightVel > 0 && gcopter::smoothedL1(violaVel, smoothFactor, violaVelPena, violaVelPenaD)) {
                gradVel += weightVel * violaVelPenaD * 2.0 * vel;
                tmp_cost += weightVel * violaVelPena;
                if (violaVel > max_pena(VEL_IDX)) max_pena(VEL_IDX) = violaVel;
            }

            /* 2.4 For acc cost  */
            const auto &violaAcc = acc.squaredNorm() - amaxSqr;
            double violaAccPena, violaAccPenaD;
            if (weightAcc > 0 && gcopter::smoothedL1(violaAcc, smoothFactor, violaAccPena, violaAccPenaD)) {
                gradAcc += weightAcc * violaAccPenaD * 2.0 * acc;
                tmp_cost += weightAcc * violaAccPena;
                if (violaAcc > max_pena(ACC_IDX)) max_pena(ACC_IDX) = violaAcc;
            }

            /* 2.5 For acc cost  */
            const auto &violaJer = jer.squaredNorm() - jmaxSqr;
            double violaJerPena, violaJerPenaD;
            if (weightJer > 0 && gcopter::smoothedL1(violaJer, smoothFactor, violaJerPena, violaJerPenaD)) {
                gradJer += weightJer * violaJerPenaD * 2.0 * jer;
                tmp_cost += weightJer * violaJerPena;
                if (violaJer > max_pena(JER_IDX)) max_pena(JER_IDX) = violaJer;
            }

            Vec3f totalGradPos{0.0, 0.0, 0.0}, totalGradVel{0.0, 0.0, 0.0},
                    totalGradAcc{0.0, 0.0, 0.0}, totalGradJer{0.0, 0.0, 0.0};

            // 若启用动力学映射，则先通过微分平坦性把速度/加速度/jerk 映射到推力与角速度，
            // 在推力/角速度空间加罚，再把梯度反传回位置导数空间。
            /* 2.6  For omg amd thr cost  */
            if (weightOmg > 0 && weightAccThr > 0) {
                double thr;
                Vec4f quat;
                Vec3f omg;
                flatMap.forward(vel, acc, jer, 0.0, 0.0, thr, quat, omg);
                const auto &violaOmg = omg.squaredNorm() - omgmaxSqr;
                const auto &violaThrust = (thr - thrustMean) * (thr - thrustMean) - thrustSqrRadi;

                /* 2.6.1  For omg cost  */
                double violaOmgPena, violaOmgPenaD;
                Vec3f gradOmg{0, 0, 0};
                if (weightOmg > 0 && gcopter::smoothedL1(violaOmg, smoothFactor, violaOmgPena, violaOmgPenaD)) {
                    gradOmg += weightOmg * violaOmgPenaD * 2.0 * omg;
                    tmp_cost += weightOmg * violaOmgPena;
                    if (violaOmg > max_pena(OMG_IDX)) max_pena(OMG_IDX) = violaOmg;
                }

                /* 2.6.2  For thr cost  */
                // 推力约束以允许区间的中心 thrustMean 为参考，
                // 当当前推力 thr 超出区间半径 thrustSqrRadi 后，对超量部分施加平滑 L1 惩罚。
                double violaThrustPena, violaThrustPenaD;
                double gradThr{0.0};
                if (weightAccThr > 0 &&
                    gcopter::smoothedL1(violaThrust, smoothFactor, violaThrustPena, violaThrustPenaD)) {
                    // 这里得到的是罚项对标量推力 thr 的梯度，后面还需要通过 flatMap.backward
                    // 把它继续反传回平坦输出对应的 pos/vel/acc/jer 导数空间。
                    gradThr += weightAccThr * violaThrustPenaD * 2.0 * (thr - thrustMean);
                    tmp_cost += weightAccThr * violaThrustPena;
                    if (violaThrust > max_pena(THR_IDX)) max_pena(THR_IDX) = violaThrust;
                }
                // 将位置/速度/加速度/jerk 本身的梯度，以及推力/角速度空间中的附加梯度，
                // 统一通过微分平坦映射反传，得到后续对多项式系数和时间积分时真正使用的总梯度。
                double totalGradPsi{0.0}, totalGradPsiD{0.0};
                flatMap.backward(gradPos, gradVel, gradAcc, gradJer, gradThr, Vec4f(0, 0, 0, 0), gradOmg,
                                 totalGradPos, totalGradVel, totalGradAcc, totalGradJer,
                                 totalGradPsi, totalGradPsiD);
            } else {
                totalGradPos = gradPos;
                totalGradVel = gradVel;
                totalGradAcc = gradAcc;
                totalGradJer = gradJer;
            }

            // 用梯形积分把当前采样点的代价和梯度累计到整段轨迹，对系数梯度直接按基函数回传，
            // 对时间梯度则结合链式法则考虑采样时刻随段时长缩放带来的影响。
            const auto node = (j == 0 || j == integralResolution) ? 0.5 : 1.0;
            const double alpha = j * integralFrac;
            gradC.block<8, 3>(i * 8, 0) += (beta0 * totalGradPos.transpose() +
                                            beta1 * totalGradVel.transpose() +
                                            beta2 * totalGradAcc.transpose() +
                                            beta3 * totalGradJer.transpose()) *
                                           node * step;
            gradT(i) += (totalGradPos.dot(vel) +
                         totalGradVel.dot(acc) +
                         totalGradAcc.dot(jer) +
                         totalGradJer.dot(sna)) *
                        alpha * node * step +
                        node * integralFrac * tmp_cost;
            cost += node * step * tmp_cost;
        }
    }

    // 仅输出各类约束本轮评估的最大违反值，供外层日志记录和调参诊断使用。
    /* 3) log all violations */
    pena_log.tail(7) = max_pena.tail(7);
}


/*
 * @ brief: This is the callback function of the L-BFGS solver
 *
 */
/*
 * 功能:
 *   作为 L-BFGS 的目标函数回调，把优化变量 x 还原为分段时间和中间控制点，构造 MINCO 轨迹后，
 *   统一计算能量项、各类约束罚项以及总时间代价，并把总代价对优化变量的梯度写回 g。
 *
 * 输入:
 *   ptr: 指向 OptimizationVariables 的上下文指针，里面保存了当前优化问题的维度、走廊、权重、MINCO 实例等状态。
 *   x: L-BFGS 当前迭代的优化变量向量，前半部分是时间变量 tau，后半部分是空间变量 xi。
 *      tau 不是直接时间，而是经过映射后的无约束时间参数；xi 的含义取决于 pos_constraint_type，
 *      可能是直接展开的中间点坐标，也可能是走廊多面体中的参数化变量。
 *
 * 输出:
 *   返回值: 当前这组优化变量对应的总代价，包含轨迹能量、约束罚项和总时间线性代价。
 *   g: 输出总代价对 x 的梯度，供 L-BFGS 使用。
 *   副作用: 会更新 obj.iter_num、obj.minco 内部参数，以及 obj.penalty_log 中的各类代价/违反量日志。
 */
double ExpTrajOpt::costFunctional(void *ptr,
                                  const VecDf &x,
                                  VecDf &g) {
    // 从 L-BFGS 透传的上下文中取出本次代价评估所需的所有问题配置和缓存对象。
    /* 1) Decode the pointer */
    OptimizationVariables &obj = *static_cast<OptimizationVariables *>(ptr);
    const auto &dimTau = obj.temporalDim;
    const auto &dimXi = obj.spatialDim;
    const auto &weightT = obj.rho;
    const auto &vPolyIdx = obj.vPolyIdx;
    const auto &vPolytopes = obj.vPolytopes;
    const auto &hPolyIdx = obj.hPolyIdx;
    const auto &hPolytopes = obj.hPolytopes;
    const auto &waypoint_attractor = obj.waypoint_attractor;
    const auto &waypoint_attractor_dead_d = obj.waypoint_attractor_dead_d;
    const auto &smooth_eps = obj.smooth_eps;
    const auto &integral_res = obj.integral_res;
    const auto &magnitudeBounds = obj.magnitudeBounds;
    const auto &penaltyWeights = obj.penaltyWeights;
    const auto &block_energy_cost = obj.block_energy_cost;

    auto &quadrotor_flatness = obj.quadrotor_flatness;

    obj.iter_num++;
    const auto &pos_constraint_type = obj.pos_constraint_type;

    const Eigen::Map<const VecDf> tau(x.data(), dimTau);
    const Eigen::Map<const VecDf> xi(x.data() + dimTau, dimXi);
    Eigen::Map<VecDf> gradTau(g.data(), dimTau);
    Eigen::Map<VecDf> gradXi(g.data() + dimTau, dimXi);

    // 将无约束优化变量重新映射回真实轨迹参数:
    // tau -> 每段实际时间 times，xi -> 中间控制点 points。
    // 位置变量既支持直接优化控制点坐标，也支持在安全走廊内部做参数化优化。
    /* 2) Reconstruct the optimization varibles */

    Mat3Df points;
    VecDf times;
    gcopter::forwardMapTauToT(tau, times);
    switch (pos_constraint_type) {
        case 1: {
            VecDf xi_e = xi;
            points = Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic>>(xi_e.data(), 3, xi_e.size() / 3);
            break;
        }
        default: {
            gcopter::forwardP(xi, vPolyIdx, vPolytopes, points);
            break;
        }
    }

    // 基于 points 和 times 设置 MINCO 轨迹，并先计算平滑能量项及其对多项式系数/时间的偏导。
    // 若 block_energy_cost 为真，则跳过能量项，只保留约束罚项与时间代价。
    /* 3) Compute the energy const and gradient */
    double cost{0};
    obj.minco.setParameters(points, times);
    MatD3f partialGradByCoeffs(8 * times.size(), 3);
    VecDf partialGradByTimes(times.size());
    partialGradByCoeffs.setZero();
    partialGradByTimes.setZero();
    if (!block_energy_cost) {
        obj.minco.getEnergy(cost);
        obj.minco.getEnergyPartialGradByCoeffs(partialGradByCoeffs);
        obj.minco.getEnergyPartialGradByTimes(partialGradByTimes);
    }
    obj.penalty_log(0) = cost;

    // 在当前 MINCO 系数轨迹上累计位置、动力学、吸引点、角速度、推力等约束罚项，
    // 并把这些罚项对时间和多项式系数的梯度继续加到 partialGradByTimes / partialGradByCoeffs 中。
    /* 4) Compute the constrain cost and gradient  */
    constraintsFunctional(times, obj.minco.getCoeffs(),
                          hPolyIdx, hPolytopes,
                          waypoint_attractor, waypoint_attractor_dead_d,
                          smooth_eps, integral_res,
                          magnitudeBounds, penaltyWeights,
                          quadrotor_flatness,
                          cost, partialGradByTimes, partialGradByCoeffs, obj.penalty_log);

    // MINCO 内部优化实际求导的对象是多项式系数和分段时间，
    // 这里先把梯度从“系数/时间空间”反传到“中间点/时间空间”，再叠加总时间正则项 rho * sum(times) 的梯度。
    /* 5) Propagate the gradient from CT to PT */
    Mat3Df gradByPoints;
    VecDf gradByTimes;
    obj.minco.propogateGrad(partialGradByCoeffs, partialGradByTimes,
                            gradByPoints, gradByTimes);
    cost += weightT * times.sum();
    gradByTimes.array() += weightT;

    // 最后把梯度映射回 L-BFGS 真正优化的变量:
    // gradByTimes 反传到 tau，gradByPoints 根据位置参数化方式反传到 xi。
    // 若 xi 使用走廊参数化，还会额外经过 norm restriction 层，保证变量保持在有效参数域附近。
    /* 6) Propagate the gradient from PT to optimization varibles*/
    gcopter::propagateGradientTToTau(tau, gradByTimes, gradTau);
    switch (pos_constraint_type) {
        case 1: {
            // WAYPOINT 模式:
            // 此时 xi 的物理含义就是“真实的中间控制点坐标”本身，
            // 整个向量按 [p1x, p1y, p1z, p2x, p2y, p2z, ...] 的方式直接展开。
            // 前向阶段只是把 xi 重新解释成 points 矩阵，并没有经过额外的几何参数化映射，
            // 因而总代价对 points 的梯度与总代价对 xi 的梯度是一一对应的。
            // 位置变量 xi 直接就是中间控制点坐标的展开形式，
            // 因此点空间梯度 gradByPoints 只需按相同内存布局展平成向量即可写回 gradXi。
            MatDf gp = gradByPoints;
            gradXi = Eigen::Map<VecDf>(gp.data(), gp.size());
            break;
        }
        default: {
            // CORRIDOR 模式:
            // 此时 xi 不是实际三维坐标，而是每个中间点在对应安全走廊内部的参数化变量。
            // 前向构点时会先取出 xi 的一段参数，经过 normalized() 和平方权重后，
            // 再与该走廊多面体的顶点表示组合，生成真正的空间点 points。
            // 因此这里优化的不是“点的位置本身”，而是“决定该点落在走廊何处的内部参数”。
            //
            // 由于前向关系是 points = f(xi) 的非线性映射，gradByPoints 不能直接赋给 gradXi，
            // 必须通过 backwardGradP 按链式法则把 dJ/dpoints 反传成 dJ/dxi。
            gcopter::backwardGradP(xi, vPolyIdx, vPolytopes, gradByPoints, gradXi);
            // 这里再额外加入参数空间的范数限制:
            // 因为 forwardP 内部会先对 xi 做归一化，同一个 points 往往可由不同尺度的 xi 表示。
            // 若不约束 xi 的模长，优化器可能把参数放得很大而几乎不改变实际点位置，造成数值不稳定。
            // normRetrictionLayer 会在 ||xi|| > 1 时增加罚项，并把对应梯度累加到 gradXi 中。
            gcopter::normRetrictionLayer(xi, vPolyIdx, vPolytopes, cost, gradXi);
            break;
        }
    }
    return cost;
}

static void truncateToSixDecimals(double &num) {
    num = std::trunc(num * 1e6) / 1e6; // 直接截断，无四舍五入
}

/*
 * @ brief: This function pre-process the corridor
 *
 */
/*
 * 功能:
 *   对输入的 H 表示安全走廊序列做预处理，生成后续轨迹优化直接使用的几类几何缓存。
 *   当前实现主要做四件事:
 *   1) 对每个原始走廊 hPolytopes[i] 调用 enumerateVs(...)，把半空间交形式(H-representation)
 *      转成顶点集合形式(V-representation)；
 *   2) 再把该顶点集合改写成“第 1 列存一个基准顶点，其余列存相对该基准顶点的位移”的形式，
 *      存入 opt_vars.vPolytopes，供 CORRIDOR 模式的 forwardP/backwardGradP 做位置参数化；
 *   3) 对每一对相邻走廊，将两者的 H 约束按行直接拼接，得到交集区域的 H 表示
 *      opt_vars.hOverlapPolytopes[i]；
 *   4) 基于这个相邻交集区域，估计一组段间吸引约束参数:
 *      - 吸引点免罚半径 waypoint_attractor_dead_d(i) 取 findInteriorDist(curIH, interior) / 2；
 *      - 按当前这个函数的实现，吸引点 waypoint_attractor.col(i) 取的是当前走廊 curIV 的顶点均值，
 *        即当前走廊 V 表示的几何中心近似，而不是相邻重叠区重新枚举后的顶点均值。
 *
 *   这些结果会被 setupProblemAndCheck() 用来初始化中间路径点，并在 CORRIDOR 模式下作为位置参数化、
 *   走廊切换吸引约束和索引映射的几何基础。
 *
 * 输入:
 *   无显式形参。函数依赖成员变量 opt_vars.hPolytopes 作为输入走廊序列，
 *   其中每个元素是一个分段安全走廊的 H 表示半空间约束。
 *
 * 输出:
 *   返回值:
 *     true 表示所有走廊都成功完成顶点枚举和相邻重叠区预处理；
 *     false 表示某个走廊无法完成顶点枚举，当前优化问题的几何准备失败。
 *   副作用:
 *     会清空并重建 opt_vars.vPolytopes；
 *     会写入每对相邻走廊的交集 H 表示 opt_vars.hOverlapPolytopes；
 *     会写入段间吸引点 opt_vars.waypoint_attractor；
 *     会写入吸引点免罚半径 opt_vars.waypoint_attractor_dead_d。
 */
bool ExpTrajOpt::processCorridor() {
    const long sizeCorridor = static_cast<long>(opt_vars.hPolytopes.size() - 1);

    // 为几何缓存分配空间。
    // 当前实现会在每次循环中向 vPolytopes 连续压入两次由 curIV 构造的相对顶点表示，
    // 最后再补一次最后一个走廊的相对顶点表示，因此预留 2 * sizeCorridor + 1 个槽位。
    opt_vars.vPolytopes.clear();
    opt_vars.vPolytopes.reserve(2 * sizeCorridor + 1);

    long nv;
    PolyhedronH curIH;
    PolyhedronV curIV, curIOB;
    opt_vars.waypoint_attractor.resize(3, sizeCorridor);
    opt_vars.waypoint_attractor_dead_d.resize(sizeCorridor);
    opt_vars.hOverlapPolytopes.resize(sizeCorridor);

    // 逐段处理前 sizeCorridor 个走廊，并同步为每一对相邻走廊构造交集区域信息。
    for (long i = 0; i < sizeCorridor; i++) {
        // 第 1 步: 将当前走廊从 H 表示枚举成 V 表示。
        // enumerateVs 内部会先调用 findInterior 找到一个严格内部点 inner，
        // 再基于 dual/quickhull 构造所有顶点，并做去重，最终得到当前走廊的顶点集合 curIV。
        if (!geometry_utils::enumerateVs(opt_vars.hPolytopes[i], curIV)) {
            cout << YELLOW << " -- [SUPER] in [ GcopterExpS4::processCorridor]: Failed to enumerate corridor Vs." << RESET
                 << endl;
            return false;
        }
        nv = curIV.cols();
        curIOB.resize(3, nv);
        // 第 2 步: 将绝对顶点坐标 curIV 转成 “base point + relative offsets” 的表示。
        // 具体为:
        //   curIOB.col(0)      = 第一个顶点 v0
        //   curIOB.col(j > 0)  = vj - v0
        // 这种写法正是 forwardP/backwardGradP 使用的顶点参数化输入格式。
        curIOB.col(0) = curIV.col(0);
        curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
        opt_vars.vPolytopes.push_back(curIOB);

        // 第 3 步: 构造当前走廊与下一走廊交集区域的 H 表示。
        // 每一行约束都满足 h0 * x + h1 * y + h2 * z + h3 <= 0。
        // 将两组约束按行直接拼接后，新的点集必须同时满足“当前走廊约束”和“下一走廊约束”，
        // 因而这个拼接结果正好就是两者交集的 H-representation。
        curIH.resize(opt_vars.hPolytopes[i].rows() + opt_vars.hPolytopes[i + 1].rows(), 4);
        curIH.topRows(opt_vars.hPolytopes[i].rows()) = opt_vars.hPolytopes[i];
        curIH.bottomRows(opt_vars.hPolytopes[i + 1].rows()) = opt_vars.hPolytopes[i + 1];
        opt_vars.hOverlapPolytopes[i] = curIH;
        Vec3f interior;

        // 第 4 步: 计算交集区域的“内部距离尺度”。
        // findInteriorDist 会对每个平面法向量做归一化，然后求解一个线性规划，
        // 寻找能够放入该 H 多面体内部的最大内接球半径，同时返回对应的内部点 interior。
        // 这里的 dis 取该内部距离的一半，作为后续 attractor 罚项的 dead zone 半径。
        const double &dis = geometry_utils::findInteriorDist(curIH, interior) / 2;

        // 第 5 步: 记录段间吸引点和免罚半径。
        opt_vars.waypoint_attractor.col(i) = curIV.colwise().mean();
        opt_vars.waypoint_attractor_dead_d(i) = dis;

        // 按当前实现，再次把同一个 curIV 转成相对顶点表示后压入 vPolytopes。
        // 因此本函数得到的 vPolytopes 是按当前代码路径构造出来的缓存结果，而不是这里重新枚举出的交集 V 表示。
        nv = curIV.cols();
        curIOB.resize(3, nv);
        curIOB.col(0) = curIV.col(0);
        curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
        opt_vars.vPolytopes.push_back(curIOB);
    }

    // 最后一段走廊没有“与下一段的交集区域”，因此这里只补它自身的 V 表示及其相对顶点形式。
    if (!geometry_utils::enumerateVs(opt_vars.hPolytopes.back(), curIV)) {
        cout << YELLOW << " -- [SUPER] in [ GcopterExpS4::processCorridor]: Failed to enumerate corridor Vs." <<
             RESET << endl;
        return false;
    }

    nv = curIV.cols();
    curIOB.resize(3, nv);
    curIOB.col(0) = curIV.col(0);
    curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
    opt_vars.vPolytopes.push_back(curIOB);
    return true;
}

/**
 * 功能:
 * 在非默认初始化分支中，为整条安全走廊建立优化所需的几何缓存，并结合引导轨迹生成热启动的中间点与分段时间。
 * 与 processCorridor 的区别是：processCorridor 主要完成走廊与相邻交叠区的几何预处理，
 * 随后再由 defaultInitialization 根据初始路径长度补出时间和中间点；
 * 这里除了做几何预处理，还会把相邻走廊交叠区最大内接圆的圆心作为段间吸引点，opt_vars.guide_t（这些离散位置对应的时间戳）为这些段间吸引点分配时间
 * 然后直接写入 opt_vars.points（优化的中间控制点初值）与 opt_vars.times（各段时长初值），用于热启动。
 *
 * 输入:
 * 无显式形参。函数依赖成员状态中的 opt_vars.hPolytopes（按轨迹段组织的安全走廊半空间约束）、
 * opt_vars.guide_path（引导轨迹采样点）、opt_vars.guide_t（引导轨迹采样时间）以及预先分配好的 opt_vars.points 和 opt_vars.times。
 *
 * 输出:
 * 返回 bool，表示几何缓存和热启动初值是否构建成功。成功时会更新 opt_vars.vPolytopes（供位置约束使用的顶点缓存）、
 * opt_vars.hOverlapPolytopes（相邻走廊交叠区的半空间表示）、opt_vars.waypoint_attractor（段间吸引点）、
 * opt_vars.waypoint_attractor_dead_d（吸引点免罚半径）、opt_vars.points 和 opt_vars.times；失败时返回 false，调用方会终止后续问题构建。
 */
bool ExpTrajOpt::processCorridorWithGuideTraj() {
    // 为每段走廊及每对相邻走廊的交叠区预留顶点缓存空间。
    // 这里沿用“每轮压入当前走廊一次、交叠区一次，最后再补最后一段走廊一次”的存储布局。
    // * 1) allocate memory for vertex
    const int sizeCorridor = static_cast<int>(opt_vars.hPolytopes.size() - 1);

    opt_vars.vPolytopes.clear();
    opt_vars.vPolytopes.reserve(2 * sizeCorridor + 1);

    long nv;
    PolyhedronH curIH;
    PolyhedronV curIV, curIOB;
    opt_vars.waypoint_attractor.resize(3, sizeCorridor);
    opt_vars.hOverlapPolytopes.resize(sizeCorridor);
    opt_vars.waypoint_attractor_dead_d.resize(sizeCorridor);
    // 逐段枚举走廊顶点，并同步构造与下一段走廊的交叠区域。
    // 当前走廊的顶点会被改写为“首个顶点 + 其余顶点相对首点的偏移量”，供后续优化变量映射使用。
    // * 2) Process the corridor
    for (int i = 0; i < sizeCorridor; i++) {
        // 第 1 步: 将当前走廊从 H 表示枚举成 V 表示。
        // enumerateVs 内部会先调用 findInterior 找到一个严格内部点 inner，
        // 再基于 dual/quickhull 构造所有顶点，并做去重，最终得到当前走廊的顶点集合 curIV。
        // * 2.1) Get current vertex
        if (!geometry_utils::enumerateVs(opt_vars.hPolytopes[i], curIV)) {
            cout << YELLOW << " -- [SUPER] in [ GcopterExpS4::processCorridor]: Failed to enumerate corridor Vs."
                 << RESET << endl;

            return false;
        }
        // * 2.3) Conver the vertex to the frame of the first point
        nv = curIV.cols();
        curIOB.resize(3, nv);

        // 第 2 步: 将绝对顶点坐标 curIV 转成 “base point + relative offsets” 的表示。
        // 具体为:
        //   curIOB.col(0)      = 第一个顶点 v0
        //   curIOB.col(j > 0)  = vj - v0
        // 这种写法正是 forwardP/backwardGradP 使用的顶点参数化输入格式。
        // *    Save the position of the first point
        curIOB.col(0) = curIV.col(0);
        // *    Use the relative position of the rest vertex.
        curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
        // *    save the i-th corridor's vertex
        opt_vars.vPolytopes.push_back(curIOB);
        
        // 第 3 步: 构造当前走廊与下一走廊交集区域的 H 表示。
        // 每一行约束都满足 h0 * x + h1 * y + h2 * z + h3 <= 0。
        // 将两组约束按行直接拼接后，新的点集必须同时满足“当前走廊约束”和“下一走廊约束”，
        // 因而这个拼接结果正好就是两者交集的 H-representation。
        // * 2.4) Find the overlap corridor
        curIH.resize(opt_vars.hPolytopes[i].rows() + opt_vars.hPolytopes[i + 1].rows(), 4);
        curIH.topRows(opt_vars.hPolytopes[i].rows()) = opt_vars.hPolytopes[i];
        curIH.bottomRows(opt_vars.hPolytopes[i + 1].rows()) = opt_vars.hPolytopes[i + 1];
        opt_vars.hOverlapPolytopes[i] = curIH;

        // 重叠走廊区域最大内接球的圆心
        Vec3f interior;
        
        // 第 4 步: 计算交集区域的“内部距离尺度”。
        // findInteriorDist 会对每个平面法向量做归一化，然后求解一个线性规划，
        // 寻找能够放入该 H 多面体内部的最大内接球半径，同时返回对应的内部点 interior。
        // 这里的 dis 取该内部距离的一半，作为后续 attractor 罚项的 dead zone 半径。
        const double dis = geometry_utils::findInteriorDist(curIH, interior) / 2;
        // 先确认相邻两段走廊的交叠区确实可用。
        // 如果内部距离为负或无穷大，说明这片同时满足两段约束的区域不存在、退化，
        // 或者线性规划没有给出稳定结果，此时后续就不能再把它作为段间连接区域使用。
        if (dis < 0.0 || std::isinf(dis)) {

            cout << YELLOW << " -- [SUPER] in [ GcopterExpS4::processCorridor]: Failed findInteriorDist Vs." <<
                 RESET << endl;
            return false;
        }
        // 基于上一步已经求出的内部点，把交叠区从半空间形式枚举成顶点集合。
        // 这里之所以显式枚举交叠区顶点，是因为当前函数希望后面压入 vPolytopes 的第二份缓存
        // 真正对应“当前走廊与下一走廊都允许经过的连接区域”，而不是重复使用当前单段走廊的顶点。
        geometry_utils::enumerateVs(curIH, interior, curIV);
        const double test_sum = curIV.sum();
        // 再做一次数值稳定性检查，避免枚举出的交叠区顶点中出现 NaN 或 Inf，
        // 否则后续无论是把内部点记为段间吸引点，还是把交叠区顶点送入位置约束参数化，都会污染优化初值。
        if (std::isnan(test_sum) || std::isinf(test_sum)) {
            return false;
        }

        // 第 5 步: 记录段间吸引点和免罚半径。
        opt_vars.waypoint_attractor.col(i) = interior;
        opt_vars.waypoint_attractor_dead_d(i) = dis;
        nv = curIV.cols();
        curIOB.resize(3, nv);
        curIOB.col(0) = curIV.col(0);
        curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
        opt_vars.vPolytopes.push_back(curIOB);
    }

    // 依据引导轨迹为每个段间吸引点匹配最近的参考采样时刻，
    // 再据此生成优化初值中的中间点和每段时长。
    // 这里将 opt_vars.points 直接设为吸引点位置，而 time_stamps 则来自最接近吸引点的引导轨迹采样时间。
    // * 3) Time and waypoint allocation for hot initialization
    // 为每个段间吸引点准备“最近引导轨迹采样点”的匹配缓存。
    // opt_vars.waypoint_attractor.cols(): 段间吸引点个数，通常等于“走廊段数 - 1”。
    // min_dis 记录当前找到的最小距离，min_id 记录该最近采样点在引导轨迹中的下标，
    // time_stamps 则保存“起点时刻 + 各段间吸引点对应时刻 + 终点时刻”。
    VecDf min_dis(opt_vars.waypoint_attractor.cols());
    VecDi min_id(opt_vars.waypoint_attractor.cols());
    VecDf time_stamps(opt_vars.waypoint_attractor.cols() + 2);

    // time_stamps(0) = 0.0 表示起点时刻。
    // time_stamps(last) = opt_vars.guide_t.back() 表示引导轨迹终点时刻。
    // 中间位置 time_stamps(j + 1) 用来存第 j 个吸引点对应的参考时刻。
    time_stamps(0) = 0.0;
    time_stamps(opt_vars.waypoint_attractor.cols() + 1) = opt_vars.guide_t.back();
    min_id.setConstant(0);
    min_dis.setConstant(std::numeric_limits<double>::max());

    // 外层遍历引导轨迹上的每个采样点 guide_path[i]
    // 内层遍历每个段间吸引点 waypoint_attractor.col(j)
    // 计算两者的欧氏距离
    // 如果当前这个引导轨迹点比之前记录的更靠近第 j 个吸引点，就更新匹配结果
    for (int i = 0; i < opt_vars.guide_path.size(); i++) {
        for (int j = 0; j < opt_vars.waypoint_attractor.cols(); j++) {
            const double dis = (opt_vars.guide_path[i] - opt_vars.waypoint_attractor.col(j)).norm();
            if (dis < min_dis[j]) {
                min_dis[j] = dis;
                min_id[j] = i;
                // 表示中间点的位置最终取吸引点本身，不是取 guide_path[i]。
                opt_vars.points.col(j) = opt_vars.waypoint_attractor.col(j);//opt_vars.guide_path[i];
                // 第 j 个中间点采用“与它最近的引导轨迹采样点”的时刻。
                time_stamps(j + 1) = opt_vars.guide_t[i];
            }
        }
    }
    
    // 把前面得到的一串绝对时间戳 time_stamps 转换成每一段轨迹的持续时间 opt_vars.times。
    // 需要注意的是，time_stamps 中存储的是每个段间吸引点对应的绝对时间，而优化变量 opt_vars.times 需要存储每段的持续时间
    for (int i = 1; i < time_stamps.size(); i++) {
        opt_vars.times(i - 1) = time_stamps(i) - time_stamps(i - 1);
        opt_vars.times(i - 1) = std::max(0.01, opt_vars.times(i - 1));
    }

    // 最后一段走廊没有后继交叠区，因此这里只补它自身的顶点缓存。
    if (!geometry_utils::enumerateVs(opt_vars.hPolytopes.back(), curIV)) {
        return false;
    }
    nv = curIV.cols();
    curIOB.resize(3, nv);
    curIOB.col(0) = curIV.col(0);
    curIOB.rightCols(nv - 1) = curIV.rightCols(nv - 1).colwise() - curIV.col(0);
    opt_vars.vPolytopes.push_back(curIOB);
    return true;
}

/**
 * 功能:
 * 在默认初始化分支中，根据初始折线路径各段的几何长度，给优化器生成一组基础的分段时间和中间控制点初值。
 * 这是一种不依赖引导轨迹时间戳的冷启动方式：时间由路径长度和最大速度粗略换算得到，
 * 空间中间点直接取各段之间的吸引点。
 * 输入:
 * 无显式形参。函数依赖成员状态中的 opt_vars.init_path（由起点、段间吸引点、终点拼成的初始折线）、
 * opt_vars.piece_num（轨迹分段数）、cfg_.max_vel（默认时间分配使用的最大速度）以及
 * opt_vars.waypoint_attractor（段间吸引点集合）。
 * 输出:
 * 无显式返回值。函数会直接更新 opt_vars.times（每一段的初始持续时间）和
 * opt_vars.points（优化使用的中间控制点初值），供后续优化变量构造阶段使用。
 */
void ExpTrajOpt::defaultInitialization() {
    // 先计算初始折线路径相邻关键点之间的直线距离，
    // 再按“距离 / 最大速度”的方式为每一段给出一个粗略时长。
    const VecDf dis = (opt_vars.init_path.leftCols(opt_vars.piece_num) -
                       opt_vars.init_path.rightCols(opt_vars.piece_num)).colwise().norm();
    const double speed = cfg_.max_vel;
    opt_vars.times = dis / speed;

    // 默认初始化不再单独搜索新的中间点，直接把段间吸引点作为空间初值。
    opt_vars.points = opt_vars.waypoint_attractor;
}

/**
 * 功能:
 * 在真正进入 L-BFGS 轨迹优化之前，统一完成问题规模计算、走廊预处理结果接入、初始路径与时间构造、
 * 变量维度统计以及 MINCO 边界条件和梯度缓存的分配。
 * 这是整个后端优化的“建模入口”，负责把前端给出的安全走廊、起终点状态和热启动信息整理成可优化的问题。

 * 输入:
 * 无显式形参。函数依赖成员状态中的 opt_vars.hPolytopes（安全走廊半空间约束）、
 * opt_vars.headPVAJ（轨迹起点状态）、opt_vars.tailPVAJ（轨迹终点状态）、
 * opt_vars.default_init（决定走默认初始化还是引导轨迹热启动）以及前面已经写入的走廊缓存和引导轨迹信息。

 * 输出:
 * 返回 bool，表示优化问题是否已成功建立。成功时会更新 opt_vars.piece_num（分段数）、
 * opt_vars.points 和 opt_vars.times（优化初值）、opt_vars.init_path（用于默认初始化和分段统计的折线路径）、
 * opt_vars.vPolyIdx / opt_vars.hPolyIdx（空间变量与走廊约束的索引映射）、opt_vars.temporalDim / opt_vars.spatialDim
 * （优化变量维度）以及 MINCO 与梯度缓存；失败时返回 false，调用方应停止后续优化。
 */
bool ExpTrajOpt::setupProblemAndCheck() {
    // 根据走廊数量确定轨迹分段数，并为时间初值和中间控制点初值预留空间。
    // 分段数等于安全走廊多面体个数，中间控制点个数则比它少 1。
    // init internal variables size;
    opt_vars.piece_num = static_cast<int>(opt_vars.hPolytopes.size());
    opt_vars.times.resize(opt_vars.piece_num);
    opt_vars.points.resize(3, opt_vars.piece_num - 1);


    // 选择走廊预处理与初始化路径构造方式。
    // 当前版本的默认初始化分支会直接抛异常，因此实际可走的是“带引导轨迹的热启动”这一路径；
    // 若预处理失败，则后续优化变量和约束索引都无法建立。
    // Check corridor and init points
    if (opt_vars.default_init) {
        throw std::runtime_error("Not support default init in this version.");
        if (!processCorridor()) {
            return false;
        }
    } else {
        if (!processCorridorWithGuideTraj()) {
            return false;
        }
    }

    // 把“起点 + 各段吸引点 + 终点”拼成一条折线路径。
    // 这条折线一方面可供默认初始化根据几何长度分配时间，
    // 另一方面也作为后续统计每段对应变量数量和走廊索引时的基础几何参考。
    opt_vars.init_path.resize(3, opt_vars.piece_num + 1);
    for (long i = 0; i < opt_vars.piece_num - 1; i++) {
        opt_vars.init_path.col(i + 1) = opt_vars.waypoint_attractor.col(i);
    }
    opt_vars.init_path.col(0) = opt_vars.headPVAJ.col(0);
    opt_vars.init_path.rightCols(1) = opt_vars.tailPVAJ.col(0);

    // 根据初始化模式生成时间初值。
    // 默认模式会重新按折线路径长度估算时长；热启动模式已经在 processCorridorWithGuideTraj 中生成了时长，
    // 这里只做一次统一缩放，让初始时间略微紧凑，便于后续优化继续调整。
    if (opt_vars.default_init) {
        defaultInitialization();
    } else {
        opt_vars.times *= 0.8;
    }

    // 若初始化时间已经出现 NaN，则说明前面的热启动或折线路径构造异常，问题不能继续建立。
    if (std::isnan(opt_vars.times.sum())) {
        cout << YELLOW << " -- [ExpOpt] Init times and point failed." << RESET << endl;
        return false;
    }

    // 统计每段折线对应的几何跨度，并据此生成 pieceIdx。
    // 当前实现通过除以 INFINITY 再加 1 的方式，使每段至少分配 1 份时间变量/空间变量映射单元。
    const Mat3Df deltas = opt_vars.init_path.rightCols(opt_vars.piece_num)
                          - opt_vars.init_path.leftCols(opt_vars.piece_num);
    opt_vars.pieceIdx = (deltas.colwise().norm() / INFINITY).cast<int>().transpose();
    opt_vars.pieceIdx.array() += 1;

    // 初始化优化变量维度和“轨迹段 <-> 走廊缓存”索引表。
    // temporalDim 对应时间变量个数，spatialDim 对应空间变量个数；
    // vPolyIdx / hPolyIdx 用于在后续前向映射和约束评估时找到每一段所关联的顶点缓存和 H 表示走廊。
    opt_vars.temporalDim = opt_vars.piece_num;
    opt_vars.spatialDim = 0;
    opt_vars.vPolyIdx.resize(opt_vars.piece_num - 1);
    opt_vars.hPolyIdx.resize(opt_vars.piece_num);

    // 根据位置约束参数化方式建立空间变量维度。
    // 若使用 WAYPOINT 模式，空间变量就是所有中间控制点坐标，维度固定为 3 * (段数 - 1)；
    // 否则使用走廊内部参数化，空间变量维度取决于每个关联多面体顶点表示的列数之和。
    switch (cfg_.pos_constraint_type) {
        case 1: {
            for (int i = 0, j = 0, k; i < opt_vars.piece_num; i++) {
                k = opt_vars.pieceIdx(i);
                for (int l = 0; l < k; l++, j++) {
                    if (l < k - 1) {
                        opt_vars.vPolyIdx(j) = 2 * i;
                    } else if (i < opt_vars.piece_num - 1) {
                        opt_vars.vPolyIdx(j) = 2 * i + 1;
                    }
                    opt_vars.hPolyIdx(j) = i;
                }
            }
            opt_vars.spatialDim = 3 * (opt_vars.piece_num - 1);
            break;
        }
        default: {
            for (int i = 0, j = 0, k; i < opt_vars.piece_num; i++) {
                k = opt_vars.pieceIdx(i);
                for (int l = 0; l < k; l++, j++) {
                    if (l < k - 1) {
                        opt_vars.vPolyIdx(j) = 2 * i;
                        opt_vars.spatialDim += static_cast<int>(opt_vars.vPolytopes[2 * i].cols());
                    } else if (i < opt_vars.piece_num - 1) {
                        opt_vars.vPolyIdx(j) = 2 * i + 1;
                        opt_vars.spatialDim += static_cast<int>(opt_vars.vPolytopes[2 * i + 1].cols());
                    }
                    opt_vars.hPolyIdx(j) = i;
                }
            }
        }
    }

    // 最后设置 MINCO 的起终点边界条件，并分配后续代价函数评估和梯度传播所需的缓存矩阵。
    opt_vars.minco.setConditions(opt_vars.headPVAJ, opt_vars.tailPVAJ, opt_vars.piece_num);
    opt_vars.gradByPoints.resize(3, opt_vars.piece_num - 1);
    opt_vars.gradByTimes.resize(opt_vars.piece_num);
    opt_vars.partialGradByCoeffs.resize(8 * opt_vars.piece_num, 3);
    opt_vars.partialGradByTimes.resize(opt_vars.piece_num);
    return true;
}

/**
 * 功能:
 * 在优化问题尺寸已经建立之后，用外部提供的中间点和分段时间直接覆盖当前优化器内部的初值。
 * 这相当于手动指定一次热启动结果，使后续 L-BFGS 不再依赖默认初始化或 guide trajectory 自动生成的初值。

 * 输入:
 * init_ps 表示外部给定的中间控制点序列；这些点会依次写入 opt_vars.points，作为空间变量的初始位置。
 * init_ts 表示外部给定的分段时间序列；这些时间会写入 opt_vars.times，作为时间变量的初始值。
 * 两者的尺寸必须与当前问题已分配好的内部缓存一致，否则无法安全覆盖。

 * 输出:
 * 返回 bool，表示外部初值是否成功写入。返回 true 时会关闭默认初始化分支，并更新 opt_vars.points 与
 * opt_vars.times；返回 false 时说明输入尺寸与当前问题规模不匹配，内部初值保持不变。
 */
bool ExpTrajOpt::setInitPsAndTs(const vec_Vec3f &init_ps, const vector<double> &init_ts) {
    // 一旦显式给定外部初值，就不再走默认初始化路径。
    opt_vars.default_init = false;

    // 先检查外部给定的时间数量是否与当前分段数一致。
    // 这里的 opt_vars.times.size() 代表当前优化问题已经分配好的分段时间个数。
    if (opt_vars.times.size() != init_ts.size()) {
        return false;
    }

    // 再检查外部给定的中间点数量是否与当前问题需要的中间控制点个数一致。
    // 这里的 opt_vars.points.cols() 代表“起点和终点之间”需要优化的内部连接点数量。
    if (opt_vars.points.cols() != init_ps.size()) {
        return false;
    }

    // 逐个写入除最后一段之外的时间和对应中间点。
    // points 的每一列对应一个内部连接点，times 的前若干项对应前面的轨迹分段时长。
    for (long i = 0; i < opt_vars.points.cols(); i++) {
        opt_vars.times[i] = init_ts[i];
        opt_vars.points.col(i) = init_ps[i];
    }

    // 最后一段时间没有对应的内部中间点，因此在循环外单独补写。
    opt_vars.times[opt_vars.times.size() - 1] = init_ts.back();
    return true;
}

double ExpTrajOpt::optimize(Trajectory &traj, const double &relCostTol) {
    /* 1) allocate vector for optimization varibles */
    VecDf x(opt_vars.temporalDim + opt_vars.spatialDim);
    /*    creat map for the opt_var vector */
    Eigen::Map<VecDf> tau(x.data(), opt_vars.temporalDim);
    Eigen::Map<VecDf> xi(x.data() + opt_vars.temporalDim, opt_vars.spatialDim);

    opt_vars.penalty_log.resize(8);
    opt_vars.penalty_log.setZero();

    /* 2) check the initial value of the optimization varibles */
    if (opt_vars.times.minCoeff() < 1e-3) {
        cout << YELLOW << " -- [TrajOpt] Error, the init times have zero, force return." << RESET << endl;
        cout << " -- Head PVAJ: " << endl;
        cout << opt_vars.headPVAJ << endl;
        cout << " -- Head PVAJ: " << endl;
        cout << opt_vars.tailPVAJ << endl;
        cout << " -- Times: " << endl;
        cout << opt_vars.times.transpose() << endl;
        return INFINITY;
    }

    if (opt_vars.given_init_ts_and_ps) {
        opt_vars.times = opt_vars.init_ts;
        for (int i = 0; i < opt_vars.init_ps.size(); i++) {
            opt_vars.points.col(i) = opt_vars.init_ps[i];
        }
    }

    /* 3)  construct the initial guess of the optimization varibles*/
    gcopter::backwardMapTToTau(opt_vars.times, tau);
    switch (opt_vars.pos_constraint_type) {
        case 1: {
            MatDf p_e = opt_vars.points;
            xi = Eigen::Map<const VecDf>(p_e.data(), p_e.size());
            break;
        }
        default: {
            gcopter::backwardP(opt_vars.points, opt_vars.vPolyIdx, opt_vars.vPolytopes, xi);
            break;
        }
    }

    /* 4) setup the optimizer's parameters*/
    opt_vars.iter_num = 0;
    double minCostFunctional{0};
    lbfgs::lbfgs_parameter_t lbfgs_params;
    lbfgs_params.mem_size = 256;
    lbfgs_params.past = 3;
    lbfgs_params.min_step = 1.0e-32;
    lbfgs_params.g_epsilon = 0.0;
    lbfgs_params.delta = relCostTol;
    VecDf times_init = opt_vars.times;

    opt_vars.init_ts = opt_vars.times;
    opt_vars.init_ps.clear();
    for (int col = 0; col < opt_vars.points.cols(); col++) {
        opt_vars.init_ps.emplace_back(opt_vars.points.col(col));
    }

    // keep fixed accuracy for
    for (int i = 0; i < opt_vars.waypoint_attractor_dead_d.size(); i++) {
        truncateToSixDecimals(opt_vars.waypoint_attractor_dead_d(i));
        truncateToSixDecimals(opt_vars.waypoint_attractor(0, i));
        truncateToSixDecimals(opt_vars.waypoint_attractor(1, i));
        truncateToSixDecimals(opt_vars.waypoint_attractor(2, i));
    }

    cout << std::fixed << std::setprecision(15);
    auto x0 = x;
    // only for debug
//    cout << " -- [ExpOpt] Start optimization." << x.transpose() << endl;
//    cout << " -- [ExpOpt] minCostFunctional: " << minCostFunctional << endl;
//    cout << " -- [ExpOpt] relCostTol: " << relCostTol << endl;
//    cout << " -- [ExpOpt] weightAtt: " << opt_vars.penaltyWeights(4) << endl;
//    cout << " -- [ExpOpt] waypoint_attractor: " << opt_vars.waypoint_attractor << endl;
//    cout << " -- [ExpOpt] waypoint_attractor_dead_d: " << opt_vars.waypoint_attractor_dead_d.transpose() << endl;
    // TimeConsuming ttt(" -- [ExpTrajOpt]", false);
    opt_vars.iter_num = 0;
    int ret = lbfgs::lbfgs_optimize(x,
                                    minCostFunctional,
                                    &ExpTrajOpt::costFunctional,
                                    nullptr,
                                    nullptr,
                                    &this->opt_vars,
                                    lbfgs_params);
    // double dt = ttt.stop();
    gcopter::forwardMapTauToT(tau, opt_vars.times);
    if (cfg_.print_optimizer_log) {
        cout << " -- [ExpOpt] Opt finish, with iter num: " << opt_vars.iter_num << "\n";
        cout << "\tEnergy: " << opt_vars.penalty_log(0) << endl;
        cout << "\tPos: " << opt_vars.penalty_log(1) << endl;
        cout << "\tVel: " << opt_vars.penalty_log(2) << endl;
        cout << "\tAcc: " << opt_vars.penalty_log(3) << endl;
        cout << "\tJerk: " << opt_vars.penalty_log(4) << endl;
        cout << "\tAttract: " << opt_vars.penalty_log(5) << endl;
        cout << "\tOmg: " << opt_vars.penalty_log(6) << endl;
        cout << "\tThr: " << opt_vars.penalty_log(7) << endl;
        cout << "\tOptimized Time: " << opt_vars.times.transpose() << endl;
    }

    if ((cfg_.penna_pos > 0 && opt_vars.penalty_log(1) > 0.2) ||
        // (cfg_.penna_vel > 0 && opt_vars.penalty_log(2) > cfg_.max_vel * cfg_.penna_margin) ||
        (cfg_.penna_acc > 0 && opt_vars.penalty_log(3) > cfg_.max_acc * cfg_.penna_margin) ||
        (cfg_.penna_omg > 0 && opt_vars.penalty_log(6) > cfg_.max_omg * cfg_.penna_margin) ||
        (cfg_.penna_thr > 0 && opt_vars.penalty_log(7) > cfg_.max_acc * cfg_.penna_margin)) {
        if (cfg_.print_optimizer_log) {
            cout << " -- [ExpOpt] Opt finish, with iter num: " << opt_vars.iter_num << "\n";
            cout << "\tEnergy: " << opt_vars.penalty_log(0) << endl;
            cout << "\tPos: " << opt_vars.penalty_log(1) << endl;
            cout << "\tVel: " << opt_vars.penalty_log(2) << endl;
            cout << "\tAcc: " << opt_vars.penalty_log(3) << endl;
            cout << "\tJerk: " << opt_vars.penalty_log(4) << endl;
            cout << "\tAttract: " << opt_vars.penalty_log(5) << endl;
            cout << "\tOmg: " << opt_vars.penalty_log(6) << endl;
            cout << "\tThr: " << opt_vars.penalty_log(7) << endl;
            cout << "\tOptimized Time: " << opt_vars.times.transpose() << endl;
        }
        ros_ptr_->warn(" -- [ExpOpt] Opt failed, Omg or thr or Pos violation.");
        ret = -1;
    }

    if (ret >= 0) {
        gcopter::forwardMapTauToT(tau, opt_vars.times);
        switch (opt_vars.pos_constraint_type) {
            case 1: {
                VecDf xi_e = xi;
                opt_vars.points = Eigen::Map<Eigen::Matrix<double, 3, Eigen::Dynamic>>(xi_e.data(), 3, xi_e.size() / 3);
                break;
            }
            default: {
                gcopter::forwardP(xi, opt_vars.vPolyIdx,
                                  opt_vars.vPolytopes, opt_vars.points);
                break;
            }
        }
//        opt_vars.minco.setConditions(opt_vars.headPVAJ, opt_vars.tailPVAJ, opt_vars.temporalDim);
        opt_vars.minco.setParameters(opt_vars.points, opt_vars.times);
        opt_vars.minco.getTrajectory(traj);
    } else {
        traj.clear();
        minCostFunctional = INFINITY;
        cout << YELLOW << " -- [MINCO] TrajOpt failed, " << lbfgs::lbfgs_strerror(ret) << RESET << endl;
//        cout << "Init times: " << times_init.transpose() << endl;
    }
    return minCostFunctional + ret;
}

ExpTrajOpt::ExpTrajOpt(const traj_opt::Config &cfg, const ros_interface::RosInterface::Ptr &ros_ptr) :
        cfg_(cfg),
        ros_ptr_(ros_ptr) {
    /// Use time as log file name
    //    auto now = std::chrono::system_clock::now();
    //    std::time_t t = std::chrono::system_clock::to_time_t(now);
    //    std::tm tm = *std::localtime(&t);
    //    std::stringstream ss;
    //    ss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    //    std::string filename = ss.str() + "_exp_opt_log.csv";
    if(cfg_.save_log_en){
        std::string filename = "exp_opt_log.csv";
        failed_traj_log.open(DEBUG_FILE_DIR(filename), std::ios::out | std::ios::trunc);
        penalty_log.open(DEBUG_FILE_DIR("exp_opt_penna.csv"), std::ios::out | std::ios::trunc);
    }

    opt_vars.magnitudeBounds.resize(6);
    opt_vars.penaltyWeights.resize(7);
    opt_vars.magnitudeBounds << cfg_.max_vel, cfg_.max_acc, cfg_.max_jerk,
            cfg_.max_omg, cfg_.min_acc_thr * cfg_.mass, cfg_.max_acc_thr * cfg_.mass;
    opt_vars.penaltyWeights << cfg_.penna_pos, cfg_.penna_vel,
            cfg_.penna_acc, cfg_.penna_jerk,
            cfg_.penna_attract, cfg_.penna_omg,
            cfg_.penna_thr;
    opt_vars.rho = cfg_.penna_t;
    opt_vars.pos_constraint_type = cfg_.pos_constraint_type;
    opt_vars.block_energy_cost = cfg_.block_energy_cost;
    opt_vars.smooth_eps = cfg_.smooth_eps;
    opt_vars.integral_res = cfg_.integral_reso;
    opt_vars.quadrotor_flatness = cfg_.quadrotot_flatness;
}

ExpTrajOpt::~ExpTrajOpt() {
    failed_traj_log.close();
    penalty_log.close();
}


//bool ExpTrajOpt::optimize(const StatePVAJ &headPVAJ, const StatePVAJ &tailPVAJ,
//                          PolytopeVec &sfcs,
//                          Trajectory &out_traj) {
//    /// Check if SFC is valid
//    if (sfcs.empty()) {
//        cout << YELLOW << " -- [TrajOpt] Error, the SFC is empty." << RESET << endl;
//        return false;
//    }
//
//    if (!SimplifySFC(headPVAJ.col(0), tailPVAJ.col(0), sfcs)) {
//        cout << YELLOW << " -- [TrajOpt] Cannot simplify sfcs." << RESET << endl;
//        //        VisualUtils::VisualizePoint(mkr_pub_, headPVAJ.col(0),Color::Pink(),"ill_start",0.5,1);
//        //        VisualUtils::VisualizePoint(mkr_pub_, tailPVAJ.col(0),Color::Pink(),"ill_end",0.5,2);
//        //        cout << "headPVAJ: " << headPVAJ.col(0).transpose() << endl;
//        //        cout << "tailPVAJ: " << tailPVAJ.col(0).transpose() << endl;
//        //        cout << YELLOW << "Killing the node." << RESET << endl;
//        //        exit(-1);
//        return false;
//    }
//
//    for (const auto &poly: sfcs) {
//        if (std::isnan(poly.GetPlanes().sum())) {
//            cout << YELLOW << " -- [TrajOpt] Error, the SFC containes NaN." << RESET << endl;
//            return false;
//        }
//    }
//
//    bool success{true};
//
//    /// Setup optimization problems
//    opt_vars.default_init = true;
//    opt_vars.given_init_ts_and_ps = false;
//    opt_vars.headPVAJ = headPVAJ;
//    opt_vars.tailPVAJ = tailPVAJ;
//    opt_vars.guide_path.clear();
//    opt_vars.guide_t.clear();
//    opt_vars.hPolytopes.resize(sfcs.size());
//    for (long i = 0; i < sfcs.size(); i++) {
//        opt_vars.hPolytopes[i] = sfcs[i].GetPlanes();
//    }
//
//    if (!setupProblemAndCheck()) {
//        cout << YELLOW << " -- [SUPER] Minco corridor preprocess error." << RESET << endl;
//        success = false;
//    }
//
//    if (success && std::isinf(optimize(out_traj, cfg_.opt_accuracy))) {
//        std::cout << YELLOW << " -- [SUPER] in [ExpTrajOpt::optimize]: Optimization failed." << RESET << std::endl;
//        success = false;
//    }
//
//    if(success){
//        out_traj.start_WT = ros_ptr_->getSimTime();
//    }
//
//    if (!success && cfg_.save_log_en) {
//        failed_traj_log << 990419 << endl;
//        failed_traj_log << headPVAJ << endl;
//        failed_traj_log << tailPVAJ << endl;
//        for (long i = 0; i < sfcs.size(); i++) {
//            failed_traj_log << i << endl;
//            failed_traj_log << sfcs[i].GetPlanes() << endl;
//        }
//    }
//
//    return success;
//}

/**
 * 功能:
 * 作为当前主运行路径使用的后端优化入口，接收前端规划器给出的起终点状态、引导路径、引导路径时间戳和安全走廊，
 * 完成输入检查、走廊约束归一化、优化问题初始化，并调用内部核心求解器生成最终轨迹。
 * 与只负责执行 L-BFGS 的 optimize(Trajectory&, relCostTol) 不同，这个重载负责把前端数据整理为 opt_vars 中的完整问题描述。

 * 输入:
 * headPVAJ 表示轨迹起点的位姿、速度、加速度和 jerk 边界状态；
 * tailPVAJ 表示轨迹终点的位姿、速度、加速度和 jerk 边界状态；
 * guide_path 表示前端构造出的离散引导路径，用于热启动初始化；
 * guide_t 表示 guide_path 各采样点对应的相对时间戳，用于热启动时间分配；
 * sfcs 表示沿引导路径生成的安全飞行走廊多面体集合；
 * out_traj 用于输出优化完成后的连续轨迹结果。

 * 输出:
 * 返回 bool，表示整次后端优化是否成功。成功时会把求得的轨迹写入 out_traj，并设置轨迹起始世界时间；
 * 失败时返回 false，out_traj 会保持清空状态，同时在启用日志时写入失败上下文，便于复现问题。
 */
bool ExpTrajOpt::optimize(const StatePVAJ &headPVAJ, const StatePVAJ &tailPVAJ,
                          const vec_E<Vec3f> &guide_path, const vector<double> &guide_t,
                          PolytopeVec &sfcs,
                          Trajectory &out_traj) {
    // 先检查热启动输入是否自洽：
    // 引导路径中的每个离散点都必须有一个对应的相对时间戳，否则无法据此构造时间初值。
    /// Check if hot init is valid
    if (guide_path.size() != guide_t.size()) {
        cout << YELLOW << " -- [TrajOpt] Error, the guide trajectory has wrong path and time stamp." << RESET
             << endl;
        return false;
    }

    // 再检查安全走廊是否非空，并尝试做一次走廊简化。
    // 若简化失败，说明后续既无法构造热启动走廊缓存，也无法建立有效的位置约束。
    /// Check if SFC is valid
    if (sfcs.empty()) {
        cout << YELLOW << " -- [TrajOpt] Error, the SFC is empty." << RESET << endl;
        return false;
    }

    // 先对安全走廊序列做一次简化：
    // 1) 裁掉起点之前和终点之后那些不再需要覆盖本次轨迹的头尾走廊；
    // 2) 在剩余走廊中按“与当前候选走廊是否仍有非空交叠”做贪心远跳，
    //    删除中间可被直接跨越的冗余走廊，只保留保证连续衔接所需的关键走廊段。
    if (!SimplifySFC(headPVAJ.col(0), tailPVAJ.col(0), sfcs)) {
        cout << YELLOW << " -- [TrajOpt] Cannot simplify sfcs." << RESET << endl;
        return false;
    }

    bool success{true};

    // 把本次优化任务的输入统一写入 opt_vars。
    // 这里明确关闭默认初始化，表示后续将走“guide_path + guide_t 驱动的热启动”路径。
    /// Setup optimization problems
    opt_vars.default_init = false;
    opt_vars.given_init_ts_and_ps = false;
    opt_vars.headPVAJ = headPVAJ;
    opt_vars.tailPVAJ = tailPVAJ;
    opt_vars.guide_path = guide_path;
    opt_vars.guide_t = guide_t;
    opt_vars.hPolytopes.resize(sfcs.size());
    // 同时把每个走廊平面的法向量做归一化，便于后面统一计算内部距离、位置罚项和顶点枚举。
    for (long i = 0; i < sfcs.size(); i++) {
        opt_vars.hPolytopes[i] = sfcs[i].GetPlanes();
        const Eigen::ArrayXd norms = opt_vars.hPolytopes[i].leftCols<3>().rowwise().norm();
        opt_vars.hPolytopes[i].array().colwise() /= norms;
    }

    // 建立走廊缓存、热启动初值、变量维度和 MINCO 边界条件。
    // 如果这一步失败，说明优化问题尚未成形，后续不能进入真正的数值求解。
    if (!setupProblemAndCheck()) {
        cout << YELLOW << " -- [SUPER] Minco corridor preprocess error." << RESET << endl;
        success = false;
    }

    // 先清空输出轨迹，避免沿用上一次调用留下的结果。
    out_traj.clear();


    // 真正进入内部核心求解器。
    // 这里的 optimize(out_traj, cfg_.opt_accuracy) 会把当前 opt_vars 中已经准备好的初值和约束送入 L-BFGS。
    if (success && std::isinf(optimize(out_traj, cfg_.opt_accuracy))) {
        cout << YELLOW << " -- [SUPER] Minco exp_traj opt failed." << RESET << endl;
        success = false;
    }

    // 无论成功还是失败，都记录本次约束违反量日志，便于后续分析优化行为。
    penalty_log << opt_vars.penalty_log.transpose() << endl;

    // 成功时补写轨迹在世界时间中的起始时刻，使外部模块可以把这条轨迹放回实时系统中执行。
    if (success) {
        out_traj.start_WT = ros_ptr_->getSimTime();
    }

    // 若优化失败且开启了日志，则把起终点、引导路径、时间戳和走廊平面完整写出，
    // 便于离线复现这次失败输入。
    if (!success && cfg_.save_log_en) {
        failed_traj_log << 990419 << endl;
        failed_traj_log << headPVAJ << endl;
        failed_traj_log << tailPVAJ << endl;
        for (double i: guide_t) {
            failed_traj_log << i << " ";
        }
        failed_traj_log << endl;
        for (const auto &i: guide_path) {
            failed_traj_log << i.transpose() << " ";
        }
        failed_traj_log << endl;
        for (long i = 0; i < sfcs.size(); i++) {
            failed_traj_log << i << endl;
            failed_traj_log << sfcs[i].GetPlanes() << endl;
        }
    }
    return success;
}

/**
 * 功能:
 * 为“外部已经准备好一组中间点和分段时间初值”的场景提供后端优化入口。
 * 这个重载的服务对象不是当前主运行路径中的前端 guide_path/guide_t 直接输入，
 * 而是那些已经在外部构造好了轨迹初值、希望把这组初值强制作为热启动结果送入优化器的调用情形。
 * 为了复用现有的走廊预处理与日志链路，这个函数会先把 init_ps（外部给定的中间控制点）和
 * init_ts（外部给定的每段时间）重新拼成一条离散 guide_path 和对应 guide_t，
 * 同时把 opt_vars.given_init_ts_and_ps 置为 true，使内部核心优化器在真正开始 L-BFGS 前
 * 优先使用这组外部初值覆盖自动热启动结果。
 * 输入:
 * headPVAJ 表示轨迹起点状态；tailPVAJ 表示轨迹终点状态；
 * sfcs 表示安全走廊多面体序列；
 * init_ps 表示外部指定的内部中间点；
 * init_ts 表示外部指定的各段持续时间；
 * out_traj 用于输出优化后的连续轨迹。
 * 输出:
 * 返回 bool，表示基于外部初值的整次后端优化是否成功。成功时会把最终轨迹写入 out_traj；
 * 失败时返回 false，并在启用日志时记录这次输入对应的失败上下文。
 */
bool ExpTrajOpt::optimize(const StatePVAJ &headPVAJ, const StatePVAJ &tailPVAJ,
                          PolytopeVec &sfcs,
                          const vec_Vec3f &init_ps,
                          const VecDf &init_ts,
                          Trajectory &out_traj) {
    vec_Vec3f guide_path;
    guide_path.emplace_back(headPVAJ.col(0));
    for (const auto &i: init_ps) {
        guide_path.emplace_back(i);
    }
    guide_path.emplace_back(tailPVAJ.col(0));
    vector<double> guide_t;
    guide_t.emplace_back(0);
    double accumulate_t = 0;
    for (int i = 0; i < init_ts.size(); i++) {
        accumulate_t += init_ts[i];
        guide_t.emplace_back(accumulate_t);
    }
    /// Check if hot init is valid
    if (guide_path.size() != guide_t.size()) {
        cout << YELLOW << " -- [TrajOpt] Error, the guide trajectory has wrong path and time stamp." << RESET
             << endl;
        return false;
    }
    /// Check if SFC is valid
    if (sfcs.empty()) {
        cout << YELLOW << " -- [TrajOpt] Error, the SFC is empty." << RESET << endl;
        return false;
    }

    if (!SimplifySFC(headPVAJ.col(0), tailPVAJ.col(0), sfcs)) {
        cout << YELLOW << " -- [TrajOpt] Cannot simplify sfcs." << RESET << endl;
        return false;
    }

    bool success{true};

    /// Setup optimization problems
    opt_vars.default_init = false;
    opt_vars.given_init_ts_and_ps = true;
    opt_vars.init_ts = init_ts;
    opt_vars.init_ps = init_ps;
    opt_vars.headPVAJ = headPVAJ;
    opt_vars.tailPVAJ = tailPVAJ;
    opt_vars.guide_path = guide_path;
    opt_vars.guide_t = guide_t;
    opt_vars.hPolytopes.resize(sfcs.size());

    for (long i = 0; i < sfcs.size(); i++) {
        opt_vars.hPolytopes[i] = sfcs[i].GetPlanes();
        const Eigen::ArrayXd norms = opt_vars.hPolytopes[i].leftCols<3>().rowwise().norm();
        opt_vars.hPolytopes[i].array().colwise() /= norms;
    }

    if (!setupProblemAndCheck()) {
        cout << YELLOW << " -- [SUPER] Minco corridor preprocess error." << RESET << endl;
        success = false;
    }

    out_traj.clear();

    if (success && std::isinf(optimize(out_traj, cfg_.opt_accuracy))) {
        cout << YELLOW << " -- [SUPER] Minco exp_traj opt failed." << RESET << endl;
        success = false;
    }
    penalty_log << opt_vars.penalty_log.transpose() << endl;

    if (success) {
        out_traj.start_WT = ros_ptr_->getSimTime();
    }


    if (!success && cfg_.save_log_en) {
        failed_traj_log << 990419 << endl;
        failed_traj_log << headPVAJ << endl;
        failed_traj_log << tailPVAJ << endl;
        for (double i: guide_t) {
            failed_traj_log << i << " ";
        }
        failed_traj_log << endl;
        for (const auto &i: guide_path) {
            failed_traj_log << i.transpose() << " ";
        }
        failed_traj_log << endl;
        for (long i = 0; i < sfcs.size(); i++) {
            failed_traj_log << i << endl;
            failed_traj_log << sfcs[i].GetPlanes() << endl;
        }
    }
    return success;
}
