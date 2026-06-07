# super_planner/src/utils 代码解读

本文面向当前仓库的 ROS2 主流程，解释 `src/SUPER/super_planner/src/utils` 下每个 `.cpp` 文件：

- 它的主要功能是什么
- 它通常被哪个代码文件调用
- 它在运行时大致何时参与
- 它完成的具体任务是什么
- 当前是否建议删除

## 1. 先说结论

`src/utils` 里的大多数文件都不是“没用的杂项文件”，而是：

- 轨迹优化器的底层数学支撑
- 数据结构的实现文件
- corridor / polytope / ellipsoid 的几何支撑
- 轨迹评估与约束检查工具

从 CMake 看，`super_planner` 会把 `src/*.cpp` 全部编进 `super` 库：

```cmake
file(GLOB_RECURSE srcs CONFIGURE_DEPENDS src/*.cpp include/*.h include/*.hpp)
add_library(super ${srcs})
```

所以“文件存在于库中”不等于“每次运行都会显式走到它”，但也绝不意味着它就是无用代码。

## 2. 当前 ROS2 主流程中的调用链

典型主链路如下：

1. `Apps/fsm_node_ros2.cpp`
2. `include/ros_interface/ros2/fsm_ros2.hpp`
3. `src/super_core/fsm.cpp`
4. `src/super_core/super_planner.cpp`
5. `src/super_core/astar.cpp`
6. `src/super_core/corridor_generator.cpp`
7. `src/super_core/ciri.cpp`
8. `src/traj_opt/exp_traj_optimizer_s4.cpp`
9. `src/traj_opt/backup_traj_optimizer_s4.cpp`
10. `src/traj_opt/yaw_traj_opt.cpp`
11. 上述文件再继续调用 `src/utils/*.cpp`

也就是说，`utils` 更像“底层工具层”，一般不会从入口直接被调用，而是通过：

- `super_planner.cpp`
- `corridor_generator.cpp`
- `ciri.cpp`
- `exp_traj_optimizer_s4.cpp`
- `backup_traj_optimizer_s4.cpp`
- `yaw_traj_opt.cpp`

间接参与主流程。

## 3. 总体分类

### 3.1 主链路高频必需

- `geometry_utils.cpp`
- `trajectory.cpp`
- `piece.cpp`
- `polytope.cpp`
- `ellipsoid.cpp`
- `minco.cpp`
- `lbfgs.cpp`
- `banded_system.cpp`
- `mvie.cpp`
- `sdlp.cpp`

### 3.2 主链路条件参与

- `quickhull.cpp`
- `root_finder.cpp`
- `polynomial_interpolation.cpp`
- `optimization_utils.cpp`

### 3.3 当前 ROS2 主流程中不明显直接使用

- `waypoint_trajectory_optimizer.cpp`

注意：即便是“当前主流程不明显直接使用”的文件，也不代表可以随手删除。删除前必须先改 CMake、完整编译并验证运行链路。

---

## 4. 按文件逐个解读

### 4.1 `banded_system.cpp`

**功能**

实现带状矩阵线性系统：

- 创建带状矩阵存储
- LU 分解
- 正向/逆向求解
- 求伴随系统

这是一个纯数值线性代数工具。

**主要调用者**

- `include/traj_opt/minco.h` 对应的 `src/utils/minco.cpp`
- `include/utils/optimization/polynomial_interpolation.h`

**何时被调用**

当系统需要从轨迹边界条件、连续性约束和中间点约束中解出多项式系数时。

**完成什么任务**

为 MINCO 和多项式插值提供高效的带状线性系统求解能力。

**在主流程中的位置**

主轨迹优化、backup 轨迹优化、yaw 轨迹插值时都属于底层依赖。

**是否建议删除**

不建议。它是 `minco.cpp` 和插值模块的基础设施。

---

### 4.2 `ellipsoid.cpp`

**功能**

实现 `Ellipsoid` 数据结构：

- 椭球从 `(C, d)` 或 `(R, r, d)` 构造
- 点到椭球距离
- 最近点、最近障碍点查询
- 坐标系变换相关支撑

**主要调用者**

- `src/super_core/ciri.cpp`
- `src/utils/mvie.cpp`
- `src/utils/polytope.cpp`
- `include/ros_interface/ros2/ros2_interface.hpp` 可视化

**何时被调用**

在 corridor 生成阶段，CIRI 会围绕 seed line / seed point 维护椭球；MVIE 也会求多面体的最大体积内接椭球。

**完成什么任务**

为 corridor 凸分解提供“椭球”这类中间几何对象。

**在主流程中的位置**

不是入口层代码，但 corridor 生成时很核心。

**是否建议删除**

不建议。`ciri.cpp` 和 `mvie.cpp` 都依赖它。

---

### 4.3 `geometry_utils.cpp`

**功能**

这是 `utils` 里最核心的综合几何工具文件之一，内容包括：

- 点到椭圆/椭球距离
- yaw / quaternion / rotation 互转
- 微分平坦性输出转姿态与角速度
- 线段与盒体相交
- path 长度计算
- polytope 内点、重叠、交集深度
- half-space 到顶点枚举
- quickhull 相关支撑
- 时间分配工具，如 `simplePMTimeAllocator`

**主要调用者**

- `src/super_core/super_planner.cpp`
- `src/super_core/corridor_generator.cpp`
- `src/traj_opt/exp_traj_optimizer_s4.cpp`
- `src/traj_opt/backup_traj_optimizer_s4.cpp`
- `src/traj_opt/yaw_traj_opt.cpp`
- `include/ros_interface/ros2/fsm_ros2.hpp`
- `src/utils/polytope.cpp`
- `src/utils/ellipsoid.cpp`

**何时被调用**

几乎整个规划流程都会用到：

- 前端 guide path 时间分配
- corridor 重叠深度判断
- polytope 顶点枚举
- yaw 归一化
- 发布控制命令时的 flatness 转换

**完成什么任务**

它相当于整个项目的几何工具箱。

**在主流程中的位置**

高频核心依赖。

**是否建议删除**

绝对不建议。

---

### 4.4 `lbfgs.cpp`

**功能**

实现 L-BFGS 优化器本体，包括：

- 参数结构
- 主优化循环
- 线搜索
- 返回码解释

**主要调用者**

- `src/traj_opt/exp_traj_optimizer_s4.cpp`
- `src/traj_opt/backup_traj_optimizer_s4.cpp`
- `src/utils/mvie.cpp`
- `src/utils/optimization_utils.cpp`
- `src/utils/waypoint_trajectory_optimizer.cpp`

**何时被调用**

只要开始做连续优化求解，就会调用。

**完成什么任务**

负责求解：

- 主轨迹优化问题
- backup 轨迹优化问题
- MVIE 内部优化问题
- 其他 gcopter 相关优化问题

**在主流程中的位置**

虽然不在顶层出现，但所有优化器都会落到这里。

**是否建议删除**

绝对不建议。

---

### 4.5 `minco.cpp`

**功能**

实现 `MINCO_S2NU / S3NU / S4NU`：

- 根据起终状态与中间点设置问题
- 求多项式系数
- 生成 `Trajectory`
- 计算能量与梯度传播

**主要调用者**

- `src/traj_opt/exp_traj_optimizer_s4.cpp`
- `src/traj_opt/backup_traj_optimizer_s4.cpp`
- `src/utils/waypoint_trajectory_optimizer.cpp`

**何时被调用**

在位置轨迹优化中，优化变量确定之后，需要把结果映射成真正的多项式轨迹。

**完成什么任务**

是“优化变量 -> 多项式轨迹”的核心桥梁。

**在主流程中的位置**

主轨迹和 backup 轨迹优化都会使用。

**是否建议删除**

绝对不建议。

---

### 4.6 `mvie.cpp`

**功能**

实现最大体积内接椭球（MVIE）求解：

- 目标函数
- 平滑 L1 约束
- 优化求解
- 从 polytope 中求最大内接椭球

**主要调用者**

- `src/super_core/ciri.cpp`

**何时被调用**

在 CIRI 的迭代凸分解里，每得到一组约束平面后，会更新其中的内接椭球。

**完成什么任务**

为 corridor 生成中的 ellipsoid / polytope 交替优化提供核心步骤。

**在主流程中的位置**

属于 corridor 生成的内部核心步骤。

**是否建议删除**

不建议。删掉后 `ciri.cpp` 基本无法工作。

---

### 4.7 `optimization_utils.cpp`

**功能**

提供通用优化工具，主要是 `Gcopter` 模板相关：

- 时间变量与优化变量映射
- 梯度传播
- corridor 中点的 forward/backward 映射
- 约束层处理

**主要调用者**

- `exp_traj_optimizer_s4.h/.cpp`
- `backup_traj_optimizer_s4.h/.cpp`
- `yaw_traj_opt.h`
- `waypoint_trajectory_optimizer.h/.cpp`

**何时被调用**

在优化器内部做变量映射和梯度处理时。

**完成什么任务**

相当于这些优化器共用的一套“优化变量操作库”。

**在主流程中的位置**

属于优化器的基础设施。

**是否建议删除**

不建议。

---

### 4.8 `piece.cpp`

**功能**

实现单段多项式轨迹 `Piece`：

- 位置/速度/加速度/jerk 求值
- 最大速度、最大加速度检查
- 单段时间长度管理

**主要调用者**

- `src/utils/trajectory.cpp`
- 通过 `Trajectory` 间接被：
  - `super_planner.cpp`
  - `fsm_ros2.hpp`
  - `exp_traj_optimizer_s4.cpp`
  - `backup_traj_optimizer_s4.cpp`
  - `yaw_traj_opt.cpp`

**何时被调用**

凡是需要对轨迹按时间求状态，最终都会落到某个 `Piece` 上。

**完成什么任务**

它是 `Trajectory` 的基本单元。

**在主流程中的位置**

高频底层依赖。

**是否建议删除**

绝对不建议。

---

### 4.9 `polynomial_interpolation.cpp`

**功能**

这个 `.cpp` 文件本身几乎没有有效实现，当前内容更像占位文件。  
真正的模板实现主要在对应头文件 `include/utils/optimization/polynomial_interpolation.h` 里。

**主要调用者**

- `src/traj_opt/yaw_traj_opt.cpp`

**何时被调用**

在 yaw 轨迹优化中做多项式插值时，真正逻辑主要来自头文件模板。

**完成什么任务**

为偏航轨迹提供插值支撑。

**在主流程中的位置**

属于 yaw 优化的工具模块，但当前 `.cpp` 自身很轻。

**是否建议删除**

不建议直接删。它看起来像占位文件，但是否可删应先确认头文件模板实例化方式和链接行为。

---

### 4.10 `polytope.cpp`

**功能**

实现 `Polytope` 数据结构：

- 设置/存储 half-space 平面
- polytope 交集 `CrossWith`
- 重叠判断
- 点内判定
- 体积估计
- 记录 seed line、ellipsoid 等附加信息

**主要调用者**

- `src/super_core/corridor_generator.cpp`
- `src/super_core/ciri.cpp`
- `src/traj_opt/exp_traj_optimizer_s4.cpp`
- `src/traj_opt/backup_traj_optimizer_s4.cpp`
- ROS 可视化适配层

**何时被调用**

一旦开始生成 corridor、做 corridor overlap 检查、或将 corridor 交给优化器处理时。

**完成什么任务**

它是 Safe Flight Corridor 的核心数据结构实现。

**在主流程中的位置**

高频核心依赖。

**是否建议删除**

绝对不建议。

---

### 4.11 `quickhull.cpp`

**功能**

实现 QuickHull 凸包算法。

**主要调用者**

- `src/utils/polytope.cpp`
- `include/ros_interface/ros2/ros2_adapter.hpp`

**何时被调用**

- 计算 polytope 体积
- polytope 可视化时将顶点组织成凸包网格

**完成什么任务**

把枚举出的顶点点集转成凸包。

**在主流程中的位置**

通常不是每个规划周期的核心算法，但做 polytope 体积/可视化时会参与。

**是否建议删除**

不建议。尤其如果你还要看 polytope 体积或可视化输出。

---

### 4.12 `root_finder.cpp`

**功能**

实现多项式根求解工具：

- 多项式求值
- 区间根个数
- 实根隔离
- 三次/四次方程求解
- safe Newton

**主要调用者**

- `src/utils/piece.cpp`

**何时被调用**

在 `Piece` 做最大速度/最大加速度分析和约束检查时。

**完成什么任务**

帮助判断一段多项式轨迹在时间区间内是否超速度/超加速度。

**在主流程中的位置**

属于轨迹约束检查的底层支撑。

**是否建议删除**

不建议。虽然不是每个主循环都显式调用，但 `Piece` 的边界检查依赖它。

---

### 4.13 `sdlp.cpp`

**功能**

实现小规模线性规划（Seidel 风格）工具。

**主要调用者**

- `src/utils/mvie.cpp`
- `src/utils/geometry_utils.cpp`

**何时被调用**

- 求 polytope 内点
- 求交集深度
- MVIE 内部线性规划子问题

**完成什么任务**

为几何模块提供小规模 LP 能力。

**在主流程中的位置**

属于 corridor / polytope / ellipsoid 几何支撑。

**是否建议删除**

不建议。

---

### 4.14 `trajectory.cpp`

**功能**

实现整条轨迹 `Trajectory`：

- 多段 `Piece` 管理
- 总时长、段数、采样
- 裁剪局部轨迹
- 拼接轨迹
- 查询某个时间的状态
- 最大速度/最大加速度检查

**主要调用者**

- `src/super_core/super_planner.cpp`
- `include/ros_interface/ros2/fsm_ros2.hpp`
- `exp_traj_optimizer_s4.cpp`
- `backup_traj_optimizer_s4.cpp`
- `yaw_traj_opt.cpp`
- ROS 可视化层

**何时被调用**

几乎所有轨迹生成、重规划、命令发布、轨迹裁剪、旧前缀拼接时都会调用。

**完成什么任务**

它是上层规划器真正操作的轨迹容器。

**在主流程中的位置**

高频核心依赖。

**是否建议删除**

绝对不建议。

---

### 4.15 `waypoint_trajectory_optimizer.cpp`

**功能**

实现另一套 waypoint 轨迹优化器：

- `GcopterWayptS3`
- `GcopterWayptS4`

它更像通用 waypoint 轨迹优化工具，而不是当前主规划链路的核心实现。

**主要调用者**

从当前 ROS2 主链路的显式搜索结果看，没有发现 `super_planner.cpp` 主流程直接调用它。

**何时被调用**

更可能用于：

- 实验
- 独立工具
- 备用轨迹优化方案
- 非当前主链路的测试/调参场景

**完成什么任务**

提供以 waypoint 为输入的另一类轨迹优化能力。

**在主流程中的位置**

当前 ROS2 主链路里不是必读核心。

**是否建议删除**

仍然不建议直接删除。虽然当前主链路不明显依赖它，但它被编进库，且可能被工具程序、调参与备用流程使用。

---

## 5. 运行时“谁真正会走到”的简化版

如果你跑的是标准 ROS2 主流程，那么 `src/utils` 可以这么理解：

### 5.1 基本一定会走到

- `geometry_utils.cpp`
- `trajectory.cpp`
- `piece.cpp`
- `polytope.cpp`
- `ellipsoid.cpp`
- `minco.cpp`
- `lbfgs.cpp`
- `banded_system.cpp`
- `mvie.cpp`
- `sdlp.cpp`

### 5.2 很可能会走到，但更偏条件/场景

- `quickhull.cpp`
- `root_finder.cpp`
- `optimization_utils.cpp`
- `polynomial_interpolation.cpp`

### 5.3 当前主链路里可以先不重点读

- `waypoint_trajectory_optimizer.cpp`

---

## 6. 我是否可以删除这些“没用”的文件？

### 6.1 不建议直接删除

原因有三：

1. 它们被全量编进 `super` 库，当前看起来“没显式调用”不代表真的无依赖。
2. 很多文件是间接依赖，不会从入口处直接看到。
3. 有些文件属于条件触发、调试工具、可视化、ROS1/ROS2共用或备用实现。

### 6.2 如果你真的想做精简，正确顺序是

1. 先标记目标文件为“疑似未用”
2. 先从 CMake 中移除，而不是物理删除
3. 完整编译
4. 运行至少这些链路：
   - ROS2 benchmark
   - click demo
   - backup 轨迹会触发的场景
   - 相关日志/可视化工具
5. 确认无回归后，再考虑真正删除

### 6.3 当前最像“主流程可暂时不看”的文件

- `waypoint_trajectory_optimizer.cpp`

但“可暂时不看”不等于“可以直接删”。

---

## 7. 阅读建议

如果你的目标是先吃透主链路，而不是立刻穷尽所有数学细节，建议先按以下顺序读：

1. `geometry_utils.cpp`
2. `trajectory.cpp`
3. `piece.cpp`
4. `polytope.cpp`
5. `ellipsoid.cpp`
6. `minco.cpp`
7. `lbfgs.cpp`
8. `mvie.cpp`
9. `sdlp.cpp`

然后再回头看：

10. `optimization_utils.cpp`
11. `root_finder.cpp`
12. `quickhull.cpp`
13. `waypoint_trajectory_optimizer.cpp`

如果你的目标是“理解运行中的每一个主步骤”，那前 1 到 9 已经覆盖了绝大多数关键依赖。
